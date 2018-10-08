#ifndef PDEKIT_Common_Block_Array_hpp
#define PDEKIT_Common_Block_Array_hpp

#include <algorithm>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#include "common/ArrayView.hpp"
#include "common/Constants.hpp"
#include "common/PDEKit.hpp"
#include "common/Range1D.hpp"
#include "common/VectorUtils.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
class BlockArray;

template <typename T, typename SizeType>
std::ostream &operator<<(std::ostream &os, const BlockArray<T, SizeType> &array);

// ----------------------------------------------------------------------------

template <typename T, typename SizeType = Uint>
class BlockArray
{

  public:
  /// TYPEDEFS
  using value_type      = T;
  using reference       = T &;
  using const_reference = const T &;
  using size_type       = SizeType;
  using difference_type = std::ptrdiff_t;

  /// METHODS
  /// Constructor
  BlockArray();

  /// Copy constructor
  BlockArray(const BlockArray &other_array);

  /// Destructor
  ~BlockArray();

  /// Assignment operator
  BlockArray &operator=(const BlockArray &array_rhs);

  /// Build the block array by copying data
  /// This is convenient but expensive. Use for small arrays.
  /// @param values      ... contains contatenated values of all blocks
  /// @param block_sizes ... contains the size of each respective block
  ///                        the size of 'block_sizes' determines the number
  ///                        of blocks
  void build(const std::vector<T> &values, const std::vector<SizeType> &block_sizes);

  /// Build the block array
  /// The creation of 'values' and 'block_sizes' is more involved, but values
  /// can be moved inside the object without additional copying
  /// @param values      ... contains contatenated values of all blocks
  /// @param block_sizes ... contains the size of each respective block
  ///                        the size of 'block_sizes' determines the number
  ///                        of blocks
  void build(std::unique_ptr<std::vector<T>> &&values,
             std::unique_ptr<std::vector<SizeType>> &&block_sizes);

  /// Build the block array
  /// @param values        ... contains contatenated values of all blocks
  /// @param block_offsets ... contains the offsets determining the position
  ///                          of each respective block
  void build_from_offsets(std::unique_ptr<std::vector<T>> &&values,
                          std::unique_ptr<std::vector<SizeType>> &&block_offsets);

  /// Get the size of the array (i.e. number of all entries)
  size_type const size() const;

  /// Get the size of one block
  size_type const nb_blocks() const;

  /// Check if the array is empty
  bool empty() const;

  /// Get the capacity of the array
  size_type const capacity() const;

  /// Resize this array
  void resize(const size_type nb_values, const size_type nb_blocks);

  /// Reserve storage in this array
  void reserve(const size_type nb_values, const size_type nb_blocks);

  /// Resize blocks
  void resize_blocks(const std::vector<size_type> &new_block_sizes);

  /// Remove certain blocks
  void remove_blocks(const std::vector<size_type> &blocks_to_remove);

  /// Fill the array with the same value
  void fill(const T &value);

  /// Swap with another block array
  void swap(BlockArray<T, SizeType> &rhs);

  /// Set a new block after the currently last block
  void create_back_block(const size_type nb_entries);

  /// Insert values of one block
  void insert_block(const SizeType block_idx,
                    const common::ArrayView<const T, _1D, SizeType> &block_values);

  /// Insert values of one block
  void fill_last_block(const common::ArrayView<const T, _1D, SizeType> &block_values);

  /// Insert one value in block
  void insert_value_in_block(const SizeType block_idx, const SizeType pos_in_block, const T value);

  /// Get one block
  common::ArrayView<const T, _1D, SizeType> const_block(const size_type i) const;

  /// Get one block
  common::ArrayView<T, _1D, SizeType> block(const size_type i);

  /// Get references to raw data
  common::VectorPack<const std::vector<T>, const std::vector<SizeType>> raw_data_crefs() const;

  /// Get raw pointers to all data
  /// Needed by some external libraries
  std::tuple<const T *, const SizeType *> raw_data() const;

  /// Get raw pointers to data of one block
  const std::tuple<const T *, const SizeType> raw_data_block(const size_type i) const;

  /// Return the position of one block, i.e. two indices
  /// indicating the first and last linearized entry of given block
  const common::Range1D<size_type> block_limits(const size_type block_id) const;

  /// Clear all data
  void clear();

  /// Return distance of a pointer from first entry of member data
  difference_type distance_from_origin(const T *ptr) const;

  /// Return reference to the first element of the array
  common::ArrayView<const T, _1D, SizeType> front_block();

  /// Return reference to the first element of the array,
  /// const version
  common::ArrayView<const T, _1D, SizeType> front_block() const;

  /// Return reference to the last element of the array
  common::ArrayView<const T, _1D, SizeType> back_block();

  /// Return reference to the last element of the array,
  /// const version
  common::ArrayView<const T, _1D, SizeType> back_block() const;

  /// Sort the blocks
  template <typename Compare = std::less<T>>
  void sort_blocks(Compare comp = Compare());

  /// Return the memory used
  Real mem_used_mb() const;

  private:
  /// DATA

  /// The actual data stored
  std::vector<T> m_values;

  /// Vector of indexes to delimit boundaries of blocks
  std::vector<size_type> m_offsets;
};

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
BlockArray<T, SizeType>::BlockArray()
{
  m_values.resize(0);
  m_offsets.resize(1);
  m_offsets[0] = 0;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
BlockArray<T, SizeType>::BlockArray(const BlockArray &other_array)
{
  m_values.resize(other_array.m_values.size());
  m_values = other_array.m_values;

  m_offsets.resize(other_array.m_offsets.size());
  m_offsets = other_array.m_offsets;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
BlockArray<T, SizeType>::~BlockArray()
{
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
BlockArray<T, SizeType> &BlockArray<T, SizeType>::operator=(const BlockArray &array_rhs)
{
  m_values.resize(array_rhs.m_values.size());
  m_values = array_rhs.m_values;

  m_offsets.resize(array_rhs.m_offsets.size());
  m_offsets = array_rhs.m_offsets;

  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::build(const std::vector<T> &values,
                                    const std::vector<SizeType> &block_sizes)
{
  m_values.resize(values.size());
  std::copy(values.begin(), values.end(), m_values.begin());

  m_offsets.resize(block_sizes.size() + 1);
  m_offsets[0] = 0;

  for (size_t i = 1; i < m_offsets.size(); ++i)
  {
    m_offsets[i] = m_offsets[i - 1] + block_sizes[i - 1];
  }
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::build(std::unique_ptr<std::vector<T>> &&values,
                                    std::unique_ptr<std::vector<SizeType>> &&block_sizes)
{
  m_values.swap(*values);

  // block_sizes have now length corresponding to the number of blocks
  // The member variable m_offsets, however, needs to have length (nb.
  // blocks+1) Therefore we're going to add one entry into block_sizes, and
  // then turn each of its entries into an offset by accumulating the block
  // sizes When the representation of 'block_sizes' is correct, we swap it
  // with 'm_offsets'

  block_sizes->push_back(0);
  for (SizeType i = block_sizes->size(); i > 0; --i)
  {
    (*block_sizes)[i] = (*block_sizes)[i - 1];
  }
  (*block_sizes)[0] = 0;

  for (SizeType i = 1; i < block_sizes->size(); ++i)
  {
    (*block_sizes)[i] += (*block_sizes)[i - 1];
  }
  m_offsets.swap(*block_sizes);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::build_from_offsets(
    std::unique_ptr<std::vector<T>> &&values,
    std::unique_ptr<std::vector<SizeType>> &&block_offsets)
{
  m_values.swap(*values);
  m_offsets.swap(*block_offsets);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
inline typename BlockArray<T, SizeType>::size_type const BlockArray<T, SizeType>::size() const
{
  return m_values.size();
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
inline typename BlockArray<T, SizeType>::size_type const BlockArray<T, SizeType>::nb_blocks() const
{
  // return (m_offsets.empty()) ? 0 : m_offsets.size() - 1;
  return (m_offsets.size() <= 1) ? 0 : m_offsets.size() - 1;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
inline bool BlockArray<T, SizeType>::empty() const
{
  // return m_offsets.empty();
  return (m_offsets.size() <= 1);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
inline typename BlockArray<T, SizeType>::size_type const BlockArray<T, SizeType>::capacity() const
{
  return m_values.capacity();
}

// ----------------------------------------------------------------------------

// This method does not make sense: resize does not define block offsets!
template <typename T, typename SizeType>
void BlockArray<T, SizeType>::resize(const size_type nb_values, const size_type nb_blocks)
{
  m_values.resize(nb_values);
  m_offsets.resize(nb_blocks + 1);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::reserve(const size_type nb_values, const size_type nb_blocks)
{
  m_values.reserve(nb_values);
  m_offsets.reserve(nb_blocks + 1);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::resize_blocks(const std::vector<size_type> &new_block_sizes)
{
  if ((new_block_sizes.size() + 1) != m_offsets.size())
  {
    std::cerr << "BlockArray::resize_blocks: array of new block sizes has "
                 "invalid length ("
              << new_block_sizes.size() << ")" << std::endl;
    return;
  }

  const size_type new_nb_values =
      std::accumulate(new_block_sizes.begin(), new_block_sizes.end(), 0);
  std::vector<T> new_values;
  new_values.reserve(new_nb_values);
  new_values.resize(0);

  for (size_type b = 0; b < new_block_sizes.size(); ++b)
  {
    const size_type old_block_size = m_offsets[b + 1] - m_offsets[b];
    const size_type new_block_size = new_block_sizes[b];

    const common::ArrayView<const T, _1D, size_type> block_data = const_block(b);

    if (old_block_size < new_block_size)
    {
      for (size_type i = 0; i < old_block_size; ++i)
      {
        new_values.push_back(block_data[i]);
      }
      for (size_type i = old_block_size; i < new_block_size; ++i)
      {
        new_values.push_back(T());
      }
    }
    else
    {
      for (size_type i = 0; i < new_block_size; ++i)
      {
        new_values.push_back(block_data[i]);
      }
    }

  } // Loop over blocks

  m_offsets.assign(m_offsets.size(), 0);
  for (size_type b = 0; b < new_block_sizes.size(); ++b)
  {
    m_offsets[b + 1] = m_offsets[b] + new_block_sizes[b];
  }
  m_values.swap(new_values);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::remove_blocks(const std::vector<size_type> &blocks_to_remove)
{
  const Uint nb_old_blocks = (m_offsets.size() < 1) ? 0 : (m_offsets.size() - 1);

  if (nb_old_blocks == 0)
  {
    std::cerr << "BlockArray::remove_blocks: error, no facets to remove" << std::endl;
    return;
  }

  std::vector<bool> remove_block(nb_old_blocks);
  remove_block.assign(nb_old_blocks, false);

  for (auto id : blocks_to_remove)
  {
    remove_block[id] = true;
  }

  // Position where data can be copied
  Uint fill_pos = 0;
  // Position to first entry that should be preserved
  Uint valid_data_pos = 0;

  Uint nb_new_blocks = 0;

  // I) Modify values array
  Uint nb_entries_after_removal = 0;
  for (Uint b = 0; b < nb_old_blocks; ++b)
  {
    const Uint block_data_len = m_offsets[b + 1] - m_offsets[b];
    // If facet data should be preserved and there were
    // some data removed prior to processing the facet,
    // copy the data on the first available position
    if (!remove_block[b])
    {
      nb_entries_after_removal += block_data_len;
      if (fill_pos != m_offsets[b])
      {
        for (Uint i = 0; i < block_data_len; ++i)
        {
          m_values[fill_pos] = m_values[valid_data_pos];
          fill_pos++;
          valid_data_pos++;
        }
      }
      else
      {
        fill_pos += block_data_len;
        valid_data_pos += block_data_len;
      }

      nb_new_blocks++;
    }
    else
    {
      valid_data_pos += block_data_len;
    }
  }

  m_values.resize(nb_entries_after_removal);

  fill_pos = 0;

  // II) Turn temporarily offsets into block sizes
  // Size of block b is stored on position [b+1]
  // m_offsets[0] will remain = 0
  for (Uint b = nb_old_blocks; b > 0; --b)
  {
    m_offsets[b] = m_offsets[b] - m_offsets[b - 1];
  }

  m_offsets[0] = 0;

  for (Uint b = 0; b < nb_old_blocks; ++b)
  {
    if (!remove_block[b])
    {
      if (fill_pos != b)
      {
        m_offsets[fill_pos + 1] = m_offsets[b + 1];
      }
      fill_pos++;
    }
  }

  // Convert block sizes back to offsets
  for (Uint b = 1; b <= nb_new_blocks; ++b)
  {
    m_offsets[b] += m_offsets[b - 1];
  }

  m_offsets.resize(nb_new_blocks + 1);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::fill(const T &value)
{
  m_values.assign(m_values.size(), value);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::swap(BlockArray<T, SizeType> &rhs)
{
  m_values.swap(rhs.m_values);
  m_offsets.swap(rhs.m_offsets);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::create_back_block(const size_type nb_entries)
{
  if (m_offsets.empty())
  {
    m_offsets.push_back(0);
  }

  const size_type last_offset = m_offsets.back();
  m_offsets.push_back(last_offset + nb_entries);

  for (size_type i = 0; i < nb_entries; ++i)
  {
    m_values.push_back(T());
  }
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::insert_block(
    const SizeType block_idx, const common::ArrayView<const T, _1D, SizeType> &block_values)
{
  const SizeType insert_pos = m_offsets[block_idx];
  for (SizeType i = 0; i < block_values.size(); ++i)
  {
    m_values[insert_pos + i] = block_values[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::fill_last_block(
    const common::ArrayView<const T, _1D, SizeType> &block_values)
{
  const SizeType insert_pos = (m_offsets.size() < 1) ? 0 : m_offsets[nb_blocks() - 1];

  for (SizeType i = 0; i < block_values.size(); ++i)
  {
    m_values[insert_pos + i] = block_values[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::insert_value_in_block(const SizeType block_idx,
                                                    const SizeType pos_in_block, const T value)
{
  m_values[m_offsets[block_idx] + pos_in_block] = value;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
common::ArrayView<const T, _1D, SizeType> BlockArray<T, SizeType>::const_block(
    const size_type i) const
{
  common::ArrayView<const T, _1D, SizeType> b(m_values.data() + m_offsets[i],
                                              m_offsets[i + 1] - m_offsets[i]);
  return b;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
common::ArrayView<T, _1D, SizeType> BlockArray<T, SizeType>::block(const size_type i)
{
  common::ArrayView<T, _1D, SizeType> b(m_values.data() + m_offsets[i],
                                        m_offsets[i + 1] - m_offsets[i]);
  return b;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
common::VectorPack<const std::vector<T>, const std::vector<SizeType>> BlockArray<
    T, SizeType>::raw_data_crefs() const
{
  return common::pack_vectors(m_values, m_offsets);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
std::tuple<const T *, const SizeType *> BlockArray<T, SizeType>::raw_data() const
{
  std::tuple<const T *, const SizeType *> result(m_values.data(), m_offsets.data());
  return result;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
const std::tuple<const T *, const SizeType> BlockArray<T, SizeType>::raw_data_block(
    const size_type i) const
{
  std::tuple<const T *, const SizeType> result(m_values.data() + m_offsets[i],
                                               m_offsets[i + 1] - m_offsets[i]);
  return result;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
inline const common::Range1D<SizeType> BlockArray<T, SizeType>::block_limits(
    const size_type block_id) const
{
  return common::Range1D<SizeType>(m_offsets[block_id], m_offsets[block_id + 1] - 1);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
void BlockArray<T, SizeType>::clear()
{
  m_values.clear();
  m_offsets.clear();
  m_offsets.push_back(0);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
typename BlockArray<T, SizeType>::difference_type BlockArray<T, SizeType>::distance_from_origin(
    const T *ptr) const
{
  return std::distance(m_values.data(), ptr);
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
common::ArrayView<const T, _1D, SizeType> BlockArray<T, SizeType>::front_block()
{
  common::ArrayView<const T, _1D, SizeType> b(m_values.data(), m_offsets[1] - m_offsets[0]);
  return b;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
common::ArrayView<const T, _1D, SizeType> BlockArray<T, SizeType>::front_block() const
{
  common::ArrayView<const T, _1D, SizeType> b(m_values.data(), m_offsets[1] - m_offsets[0]);
  return b;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
common::ArrayView<const T, _1D, SizeType> BlockArray<T, SizeType>::back_block()
{
  const size_type offset_one_before_last = m_offsets[m_offsets.size() - 2];
  common::ArrayView<const T, _1D, SizeType> b(m_values.data() + offset_one_before_last,
                                              m_offsets.back() - offset_one_before_last);
  return b;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
common::ArrayView<const T, _1D, SizeType> BlockArray<T, SizeType>::back_block() const
{
  const size_type offset_one_before_last = m_offsets[m_offsets.size() - 2];
  common::ArrayView<const T, _1D, SizeType> b(m_values.data() + offset_one_before_last,
                                              m_offsets.back() - offset_one_before_last);
  return b;
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
template <typename Compare>
void BlockArray<T, SizeType>::sort_blocks(Compare comp)
{
  if (nb_blocks() == 0)
  {
    return; // Nothing to sort
  }

  for (Uint b = 0; (b + 1) < m_offsets.size(); ++b)
  {
    std::sort(m_values.data() + m_offsets[b], m_values.data() + m_offsets[b + 1], comp);
  }
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
Real BlockArray<T, SizeType>::mem_used_mb() const
{
  return 1.e-6 * (m_values.size() * sizeof(T) + m_offsets.size() * sizeof(SizeType));
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
std::ostream &operator<<(std::ostream &os, const BlockArray<T, SizeType> &array)
{
  for (SizeType ib = 0; ib < array.nb_blocks(); ++ib)
  {
    auto block = array.const_block(ib);
    for (SizeType e = 0; e < block.size(); ++e)
    {
      os << block[e] << " ";
    }
    os << std::endl;
  }
  return os;
}

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit

#endif // PDEKIT_Common_BlockArray_hpp
