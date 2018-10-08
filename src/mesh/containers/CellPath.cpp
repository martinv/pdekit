#include <iostream>

#include "mesh/MeshConstants.hpp"
#include "mesh/containers/CellPath.hpp"

namespace pdekit
{

namespace mesh
{

namespace detail
{

// ----------------------------------------------------------------------------
// Class PathSegment
// ----------------------------------------------------------------------------

PathSegment::PathSegment() : m_value(0u)
{
}

// ----------------------------------------------------------------------------

PathSegment::~PathSegment()
{
}

// ----------------------------------------------------------------------------

void PathSegment::set_entry(const Uint pos, const Uint value)
{
  // Set 4 consequent bits (0xF) and shifts them
  // to given position (a 'shift' is always a multiple of 4)
  // This gives 000011110000 ... 0000
  // Clear mask is the 'negative' of this - hence the '~'
  // clear_mask will therefore be 111100001111 ... 1111
  const uint32_t clear_mask = ~(0xF << (4 * pos));

  // Define the value as 4 bits (for example 1011) and
  // move them to given position by shifting by a multiple of 4
  const uint32_t set_mask = (0xF & value) << (4 * pos);

  // (m_value & clear_mask) will preserve all bits except a block of 4,
  // which will be set to 0
  // Bitwise 'or' with set_mask will then set these 4 previously reset bits
  m_value = (m_value & clear_mask) | set_mask;
}

// ----------------------------------------------------------------------------

Uint PathSegment::entry(const Uint pos) const
{
  return (m_value >> (4 * pos) & 0xF);
}

// ----------------------------------------------------------------------------

} // namespace detail

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const detail::PathSegment &segment)
{
  for (Uint i = 0; i + 1 < segment.max_num_entries; ++i)
  {
    os << segment.entry(i) << "-";
  }
  os << segment.entry(segment.max_num_entries - 1);

  return os;
}

// ----------------------------------------------------------------------------
// Class CellPath
// ----------------------------------------------------------------------------

CellPath::CellPath() : m_zero_level_cell(FlatIdx(INVALID_CELL_ID)), m_num_entries(0)
{
}

// ----------------------------------------------------------------------------

CellPath::CellPath(const CellPath &rhs)
{
  if (rhs.m_path)
  {
    const detail::PathSegment &block0 = rhs.m_path[0];
    const Uint num_blocks             = block0.entry(0);
    reallocate_storage_blocks(num_blocks);

    for (Uint b = 0; b < num_blocks; ++b)
    {
      m_path[b] = rhs.m_path[b];
    }
  }

  m_zero_level_cell = rhs.m_zero_level_cell;
  m_num_entries     = rhs.m_num_entries;
}

// ----------------------------------------------------------------------------

CellPath::CellPath(const FlatIdx zero_level_cell_id, const std::vector<Uint> &entries)
{
  build(zero_level_cell_id, entries);
}

// ----------------------------------------------------------------------------

CellPath::~CellPath()
{
}

// ----------------------------------------------------------------------------

CellPath &CellPath::operator=(const CellPath &rhs)
{
  if (rhs.m_path)
  {
    const detail::PathSegment &block0 = rhs.m_path[0];
    const Uint num_blocks             = block0.entry(0);
    reallocate_storage_blocks(num_blocks);

    for (Uint b = 0; b < num_blocks; ++b)
    {
      m_path[b] = rhs.m_path[b];
    }
  }

  m_zero_level_cell = rhs.m_zero_level_cell;
  m_num_entries     = rhs.m_num_entries;

  return *this;
}

// ----------------------------------------------------------------------------

void CellPath::build(const FlatIdx zero_level_cell_id, const std::vector<Uint> &entries)
{
  m_zero_level_cell = zero_level_cell_id;

  if (entries.size() == 0)
  {
    m_num_entries = 0;
    m_path.reset(nullptr);
    return;
  }

  const Uint blk_size = detail::PathSegment::max_num_entries;

  const Uint nb_needed_blocks = ((entries.size() + 1) % blk_size == 0)
                                    ? (entries.size() + 1) / blk_size
                                    : (entries.size() + 1) / blk_size + 1;

  Uint nb_current_blocks = 0;
  if (m_path)
  {
    nb_current_blocks = m_path[0].entry(0);
  }

  if (nb_current_blocks != nb_needed_blocks)
  {
    reallocate_storage_blocks(nb_needed_blocks);
  }
  m_path[0].set_entry(0, nb_needed_blocks);

  for (Uint i = 0; i < entries.size(); ++i)
  {
    const Uint block_idx    = (i + 1) / detail::PathSegment::max_num_entries;
    const Uint pos_in_block = (i + 1) % detail::PathSegment::max_num_entries;

    detail::PathSegment &block = m_path[block_idx];

    block.set_entry(pos_in_block, entries[i]);
  }

  m_num_entries = entries.size();
}

// ----------------------------------------------------------------------------

Uint CellPath::size() const
{
  return m_num_entries;
}

// ----------------------------------------------------------------------------

Uint CellPath::capacity() const
{
  if (!m_path)
  {
    return 0;
  }

  const detail::PathSegment &block0 = m_path[0];
  return block0.entry(0) * detail::PathSegment::max_num_entries - 1;
}

// ----------------------------------------------------------------------------

void CellPath::reserve(const Uint capacity)
{
  const Uint blk_size = detail::PathSegment::max_num_entries;

  const Uint nb_needed_blocks =
      ((capacity + 1) % blk_size == 0) ? (capacity + 1) / blk_size : (capacity + 1) / blk_size + 1;
  reallocate_storage_blocks(nb_needed_blocks);
  m_num_entries = 0;
}

// ----------------------------------------------------------------------------

void CellPath::set_base_cell_id(const FlatIdx cell_id)
{
  m_zero_level_cell = cell_id;
}

// ----------------------------------------------------------------------------

FlatIdx CellPath::base_cell_id() const
{
  return m_zero_level_cell;
}

// ----------------------------------------------------------------------------

void CellPath::push_back(const Uint value)
{
  if (!m_path)
  {
    reallocate_storage_blocks(1);
    m_path[0].set_entry(1, value);
    m_num_entries = 1;
    return;
  }
  else
  {
    if ((m_num_entries + 1) > capacity())
    {
      const Uint new_nb_blocks = m_path[0].entry(0) + 1;
      reallocate_storage_blocks(new_nb_blocks);
    }

    // The i-the entry of the path is stored on position i+1
    // m_path[0].entry(0) stores the number of blocks in the object
    const Uint store_pos = m_num_entries + 1;

    const Uint block_idx    = store_pos / detail::PathSegment::max_num_entries;
    const Uint pos_in_block = store_pos % detail::PathSegment::max_num_entries;

    detail::PathSegment &block = m_path[block_idx];

    block.set_entry(pos_in_block, value);
    m_num_entries++;
  }
}

// ----------------------------------------------------------------------------

void CellPath::pop()
{
  if (m_num_entries == 0)
    return;

  // The condition here is not
  // if ((m_num_entries - 1) % detail::PathSegment::max_num_entries == 0)
  // but
  // if (m_num_entries % detail::PathSegment::max_num_entries == 0)

  // The reason is that the total number of stored data is
  // m_num_entries + 1, (first value in the path is the number of blocks!)
  // This means that after removing one entry, the total number of values
  // held in the object is 1 + m_num_entries - 1 = m_num_entries

  if (m_num_entries % detail::PathSegment::max_num_entries == 0)
  {
    const Uint num_blocks = m_path[0].entry(0);
    reallocate_storage_blocks(num_blocks - 1);
  }
  m_num_entries--;
}

// ----------------------------------------------------------------------------

Uint CellPath::path_entry(const Uint i) const
{
  const Uint block_idx    = (i + 1) / detail::PathSegment::max_num_entries;
  const Uint pos_in_block = (i + 1) % detail::PathSegment::max_num_entries;

  const detail::PathSegment &block = m_path[block_idx];

  return block.entry(pos_in_block);
}

// ----------------------------------------------------------------------------

void CellPath::print() const
{
  std::cout << "Zero level cell: " << m_zero_level_cell.id() << "\n";
  std::cout << " size: " << size() << ", capacity = " << capacity() << "\n";
  if (!m_path)
  {
    std::cout << " no path allocated";
    return;
  }

  if (size() == 0)
  {
    return;
  }

  const Uint num_segments = m_path[0].entry(0);

  Uint segment_id     = 0;
  Uint pos_in_segment = 2;

  detail::PathSegment curr_segment = m_path[segment_id];

  // Output the first entry in the path
  std::cout << " " << curr_segment.entry(1);

  Uint idx = 1;
  while (idx < size())
  {
    if (pos_in_segment != 0)
    {
      std::cout << "-";
    }

    std::cout << curr_segment.entry(pos_in_segment);
    pos_in_segment = (pos_in_segment + 1) % detail::PathSegment::max_num_entries;

    if (pos_in_segment == 0)
    {
      std::cout << " | ";
      segment_id++;
      if (segment_id < num_segments)
      {
        curr_segment = m_path[segment_id];
      }
    }
    idx++;
  } // while

  std::cout << std::endl;
  return;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const CellPath &path)
{
  os << "[" << path.m_zero_level_cell.id() << "] {";
  if (!path.m_path)
  {
    os << " }";
    return os;
  }

  if (path.size() == 0)
  {
    os << " }";
    return os;
  }

  const Uint num_segments = path.m_path[0].entry(0);

  Uint segment_id     = 0;
  Uint pos_in_segment = 2;

  detail::PathSegment curr_segment = path.m_path[segment_id];

  // Output the first entry in the path
  os << " " << curr_segment.entry(1);

  Uint idx = 1;
  while (idx < path.size())
  {
    if (pos_in_segment != 0)
    {
      os << "-";
    }

    os << curr_segment.entry(pos_in_segment);
    pos_in_segment = (pos_in_segment + 1) % detail::PathSegment::max_num_entries;

    if (pos_in_segment == 0)
    {
      os << " | ";
      segment_id++;
      if (segment_id < num_segments)
      {
        curr_segment = path.m_path[segment_id];
      }
    }
    idx++;
  } // while

  os << " }";

  return os;
}

// ----------------------------------------------------------------------------

void CellPath::reallocate_storage_blocks(const Uint num_blocks)
{
  if (num_blocks == 0)
  {
    m_path.reset(nullptr);
    m_num_entries = 0;
    return;
  }

  if (m_path)
  {
    const detail::PathSegment &block0 = m_path[0];

    const Uint curr_nb_blocks = block0.entry(0);

    if (curr_nb_blocks != num_blocks)
    {
      std::unique_ptr<detail::PathSegment[]> new_blocks(new detail::PathSegment[num_blocks]);
      const Uint min_nb_blocks = std::min(curr_nb_blocks, num_blocks);

      for (Uint b = 0; b < min_nb_blocks; ++b)
      {
        new_blocks[b] = m_path[b];
      }
      std::swap(m_path, new_blocks);
      m_path[0].set_entry(0, num_blocks);
    }
  }
  else
  {
    m_path = std::unique_ptr<detail::PathSegment[]>(new detail::PathSegment[num_blocks]);
    m_path[0].set_entry(0, num_blocks);
    m_num_entries = 0;
  }
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
