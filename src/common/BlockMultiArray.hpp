#ifndef PDEKIT_Common_Block_Multi_Array_hpp
#define PDEKIT_Common_Block_Multi_Array_hpp

#include <algorithm>
#include <iostream>
#include <memory>
#include <numeric>

#include "common/ArrayView.hpp"
#include "common/Constants.hpp"
#include "common/Range1D.hpp"
#include "common/TupleMeta.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

namespace internal
{
template <Uint I, Uint N, typename... T>
struct DataAlgorithm;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
class BlockMultiArray;

template <typename... StoredTypes>
std::ostream &operator<<(std::ostream &os, const BlockMultiArray<StoredTypes...> &array);

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
class BlockMultiArray
{

  public:
  enum
  {
    NFields = sizeof...(StoredTypes)
  };

  /// TYPEDEFS
  using value_types = std::tuple<StoredTypes...>;

  template <Uint N>
  using field_type = typename meta::tuple_nth_elem_type<N, StoredTypes...>::type;

  template <Uint N>
  using reference = typename meta::tuple_nth_elem_type<N, StoredTypes...>::type &;

  template <Uint N>
  using const_reference = const typename meta::tuple_nth_elem_type<N, StoredTypes...>::type &;

  using size_type       = Uint;
  using difference_type = std::ptrdiff_t;

  /// METHODS
  /// Constructor
  BlockMultiArray();

  /// Copy constructor
  BlockMultiArray(const BlockMultiArray &other_array);

  /// Destructor
  ~BlockMultiArray();

  /// Assignment operator
  BlockMultiArray &operator=(const BlockMultiArray &array_rhs);

  /// Build the block array
  /// @param values      ... contains contatenated values of all blocks
  /// @param block_sizes ... contains the size of each respective block
  ///                        the size of 'block_sizes' determines the number
  ///                        of blocks
  void build(std::tuple<std::unique_ptr<std::vector<StoredTypes>>...> &&values,
             std::unique_ptr<std::vector<Uint>> &&block_sizes);

  /// Build the block array
  /// @param values        ... contains contatenated values of all blocks
  /// @param block_offsets ... contains the offsets determining the position
  ///                          of each respective block
  void build_from_offsets(std::tuple<std::unique_ptr<std::vector<StoredTypes>>...> &&values,
                          std::unique_ptr<std::vector<Uint>> &&block_offsets);

  /// Get the size of the array (i.e. number of all entries)
  size_type size() const;

  /// Get the size of one block
  size_type nb_blocks() const;

  /// Check if the array is empty
  bool empty() const;

  /// Get the capacity of the array
  template <Uint I>
  size_type capacity() const;

  /// Resize this array
  void resize(const size_type nb_values, const size_type nb_blocks);

  /// Reserve storage in this array
  void reserve(const size_type nb_values, const size_type nb_blocks);

  /// Resize blocks
  void resize_blocks(const std::vector<size_type> &new_block_sizes);

  /// Remove certain blocks
  void remove_blocks(const std::vector<size_type> &blocks_to_remove);

  /// Fill the array with the same value
  template <Uint I>
  void fill(const field_type<I> &value);

  /// Set a new block after the currently last block
  void create_back_block(const size_type block_size);

  /// Insert values of one block
  template <Uint I>
  void insert_block(const size_type block_idx,
                    const common::ArrayView<const field_type<I>, _1D, Uint> &block_values);

  /// Insert values of the last block
  template <Uint I>
  void fill_last_block(const common::ArrayView<const field_type<I>, _1D, Uint> &block_values);

  /// Insert one value in block
  template <Uint I>
  void insert_value_in_block(const size_type block_idx, const size_type pos_in_block,
                             const field_type<I> value);

  /// Get one block, const version
  template <Uint I>
  common::ArrayView<const field_type<I>, _1D, Uint> const_block(const size_type i) const;

  /// Get one block
  template <Uint I>
  common::ArrayView<field_type<I>, _1D, Uint> block(const size_type i);

  /// Get raw pointers to field data and offset array
  /// Needed by some external libraries
  template <Uint I>
  std::tuple<const field_type<I> *, const Uint *> raw_data() const;

  /// Return the position of one block, i.e. two indices
  /// indicating the first and last linearized entry of given block
  const common::Range1D<size_type> block_limits(const size_type block_id) const;

  /// Clear all data
  void clear();

  /// Return distance of a pointer from first entry of member data
  template <Uint I>
  difference_type distance_from_origin(const field_type<I> *ptr) const;

  /// Return reference to the first element of the array
  template <Uint I>
  common::ArrayView<const field_type<I>, _1D, Uint> front_block();

  /// Return reference to the first element of the array,
  /// const version
  template <Uint I>
  common::ArrayView<const field_type<I>, _1D, Uint> front_block() const;

  /// Return reference to the last element of the array
  template <Uint I>
  common::ArrayView<const field_type<I>, _1D, Uint> back_block();

  /// Return reference to the last element of the array,
  /// const version
  template <Uint I>
  common::ArrayView<const field_type<I>, _1D, Uint> back_block() const;

  /// Sort the blocks
  template <Uint I, typename Compare = std::less<field_type<I>>>
  void sort_blocks(Compare comp = Compare());

  /// Return the memory used
  Real mem_used_mb() const;

  /// Output operator
  friend std::ostream &operator<<<StoredTypes...>(std::ostream &os,
                                                  const BlockMultiArray<StoredTypes...> &array);

  private:
  /// TYPEDEFS
  template <Uint I, Uint N, typename... T>
  using DataAlgo = internal::DataAlgorithm<I, N, T...>;

  using storage_type = typename meta::TransformToVector<StoredTypes...>::type;

  /// DATA
  /// The actual data stored
  storage_type m_values;

  /// Vector of indexes to delimit boundaries of blocks
  std::vector<size_type> m_offsets;
};

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
BlockMultiArray<StoredTypes...>::BlockMultiArray()
{
  // internal::DataAlgorithm<0, NFields - 1, StoredTypes...>::resize(m_values,
  // 0);
  DataAlgo<0, NFields - 1, StoredTypes...>::resize(m_values, 0);

  m_offsets.resize(1);
  m_offsets[0] = 0;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
BlockMultiArray<StoredTypes...>::BlockMultiArray(const BlockMultiArray &other_array)
{
  DataAlgo<0, NFields - 1, StoredTypes...>::copy(m_values, other_array.m_values);

  m_offsets.resize(other_array.m_offsets.size());
  m_offsets = other_array.m_offsets;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
BlockMultiArray<StoredTypes...>::~BlockMultiArray()
{
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
BlockMultiArray<StoredTypes...> &BlockMultiArray<StoredTypes...>::operator=(
    const BlockMultiArray &array_rhs)
{
  DataAlgo<0, NFields - 1, StoredTypes...>::copy(m_values, array_rhs.m_values);

  m_offsets.resize(array_rhs.m_offsets.size());
  m_offsets = array_rhs.m_offsets;

  return *this;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
void BlockMultiArray<StoredTypes...>::build(
    std::tuple<std::unique_ptr<std::vector<StoredTypes>>...> &&values,
    std::unique_ptr<std::vector<Uint>> &&block_sizes)
{
  // Check that the number of arrays holding field data is correct
  /*
  using input_arg_type =
  std::tuple<std::unique_ptr<std::vector<StoredTypes>>...>;
  static_assert(NFields == std::tuple_size<input_arg_type>::value,
                "The number of input fields is not correct.");
  */

  // block_sizes have now length corresponding to the number of blocks
  // The member variable m_offsets, however, needs to have length (nb.
  // blocks+1) Therefore we're going to add one entry into block_sizes, and
  // then turn each of its entries into an offset by accumulating the block
  // sizes When the representation of 'block_sizes' is correct, we swap it
  // with 'm_offsets'

  block_sizes->push_back(0);
  for (Uint i = block_sizes->size(); i > 0; --i)
  {
    (*block_sizes)[i] = (*block_sizes)[i - 1];
  }
  (*block_sizes)[0] = 0;

  for (Uint i = 1; i < block_sizes->size(); ++i)
  {
    (*block_sizes)[i] += (*block_sizes)[i - 1];
  }
  m_offsets.swap(*block_sizes);

  DataAlgo<0, NFields - 1, StoredTypes...>::build_values(m_values, std::move(values));
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
void BlockMultiArray<StoredTypes...>::build_from_offsets(
    std::tuple<std::unique_ptr<std::vector<StoredTypes>>...> &&values,
    std::unique_ptr<std::vector<Uint>> &&block_offsets)
{
  DataAlgo<0, NFields - 1, StoredTypes...>::build_values(m_values, std::move(values));
  m_offsets.swap(*block_offsets);
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
typename BlockMultiArray<StoredTypes...>::size_type BlockMultiArray<StoredTypes...>::size() const
{
  return std::get<0>(m_values).size();
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
typename BlockMultiArray<StoredTypes...>::size_type BlockMultiArray<StoredTypes...>::nb_blocks()
    const
{
  return (m_offsets.size() <= 1) ? 0 : m_offsets.size() - 1;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
inline bool BlockMultiArray<StoredTypes...>::empty() const
{
  // return m_offsets.empty();
  return (m_offsets.size() <= 1);
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
typename BlockMultiArray<StoredTypes...>::size_type BlockMultiArray<StoredTypes...>::capacity()
    const
{
  return std::get<I>(m_values).capacity();
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
void BlockMultiArray<StoredTypes...>::resize(const size_type nb_values, const size_type nb_blocks)
{
  DataAlgo<0, NFields - 1, StoredTypes...>::resize(m_values, nb_values);

  m_offsets.resize(nb_blocks + 1);
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
void BlockMultiArray<StoredTypes...>::reserve(const size_type nb_values, const size_type nb_blocks)
{
  DataAlgo<0, NFields - 1, StoredTypes...>::reserve(m_values, nb_values);

  m_offsets.reserve(nb_blocks + 1);
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
void BlockMultiArray<StoredTypes...>::resize_blocks(const std::vector<size_type> &new_block_sizes)
{
  if ((new_block_sizes.size() + 1) != m_offsets.size())
  {
    std::cerr << "BlockMultiArray::resize_blocks: array of new block sizes has "
                 "invalid length ("
              << new_block_sizes.size() << ")" << std::endl;
    return;
  }

  std::vector<size_type> new_offsets(m_offsets.size());
  new_offsets[0] = 0;

  for (Uint i = 0; i < new_block_sizes.size(); ++i)
  {
    new_offsets[i + 1] = new_offsets[i] + new_block_sizes[i];
  }

  DataAlgo<0, NFields - 1, StoredTypes...>::resize_data_storage(m_offsets, new_offsets, m_values);

  m_offsets.swap(new_offsets);
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
void BlockMultiArray<StoredTypes...>::remove_blocks(const std::vector<size_type> &blocks_to_remove)
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

  // I) Modify values array
  DataAlgo<0, NFields - 1, StoredTypes...>::remove_data_blocks_storage(remove_block, m_offsets,
                                                                       m_values);

  // II) Turn temporarily offsets into block sizes
  // Size of block b is stored on position [b+1]
  // m_offsets[0] will remain = 0
  Uint fill_pos = 0;

  // Count the blocks that will remain after removal
  Uint nb_new_blocks = 0;
  for (auto flag : remove_block)
  {
    if (!flag)
    {
      nb_new_blocks++;
    }
  }

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

template <typename... StoredTypes>
template <Uint I>
void BlockMultiArray<StoredTypes...>::fill(const field_type<I> &value)
{
  std::get<I>(m_values).assign(std::get<I>(m_values).size(), value);
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
void BlockMultiArray<StoredTypes...>::create_back_block(const size_type block_size)
{
  if (m_offsets.empty())
  {
    m_offsets.push_back(0);
  }

  const size_type last_offset = m_offsets.back();
  m_offsets.push_back(last_offset + block_size);

  DataAlgo<0, NFields - 1, StoredTypes...>::extend_data_array(block_size, m_values);
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
void BlockMultiArray<StoredTypes...>::insert_block(
    const typename BlockMultiArray<StoredTypes...>::size_type block_idx,
    const common::ArrayView<const field_type<I>, _1D,
                            typename BlockMultiArray<StoredTypes...>::size_type> &block_values)
{
  const size_type insert_pos = m_offsets[block_idx];
  for (size_type i = 0; i < block_values.size(); ++i)
  {
    std::get<I>(m_values)[insert_pos + i] = block_values[i];
  }
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
void BlockMultiArray<StoredTypes...>::fill_last_block(
    const common::ArrayView<const field_type<I>, _1D,
                            typename BlockMultiArray<StoredTypes...>::size_type> &block_values)
{
  const Uint insert_pos = (m_offsets.size() < 1) ? 0 : m_offsets[nb_blocks() - 1];

  for (Uint i = 0; i < block_values.size(); ++i)
  {
    std::get<I>(m_values)[insert_pos + i] = block_values[i];
  }
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
void BlockMultiArray<StoredTypes...>::insert_value_in_block(
    const typename BlockMultiArray<StoredTypes...>::size_type block_idx,
    const typename BlockMultiArray<StoredTypes...>::size_type pos_in_block,
    const field_type<I> value)
{
  std::get<I>(m_values)[m_offsets[block_idx] + pos_in_block] = value;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
common::ArrayView<const typename BlockMultiArray<StoredTypes...>::template field_type<I>, _1D, Uint>
BlockMultiArray<StoredTypes...>::const_block(
    const typename BlockMultiArray<StoredTypes...>::size_type i) const
{
  common::ArrayView<const field_type<I>, _1D, Uint> b(std::get<I>(m_values).data() + m_offsets[i],
                                                      m_offsets[i + 1] - m_offsets[i]);
  return b;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
common::ArrayView<typename BlockMultiArray<StoredTypes...>::template field_type<I>, _1D, Uint>
BlockMultiArray<StoredTypes...>::block(const typename BlockMultiArray<StoredTypes...>::size_type i)
{
  common::ArrayView<field_type<I>, _1D, Uint> b(std::get<I>(m_values).data() + m_offsets[i],
                                                m_offsets[i + 1] - m_offsets[i]);
  return b;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
std::tuple<const typename BlockMultiArray<StoredTypes...>::template field_type<I> *, const Uint *>
BlockMultiArray<StoredTypes...>::raw_data() const
{
  std::tuple<const field_type<I> *, const Uint *> result(std::get<I>(m_values).data(),
                                                         m_offsets.data());
  return result;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
const common::Range1D<typename BlockMultiArray<StoredTypes...>::size_type> BlockMultiArray<
    StoredTypes...>::block_limits(const size_type block_id) const
{
  return common::Range1D<size_type>(m_offsets[block_id], m_offsets[block_id + 1] - 1);
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
void BlockMultiArray<StoredTypes...>::clear()
{
  DataAlgo<0, NFields - 1, StoredTypes...>::clear(m_values);

  m_offsets.clear();
  m_offsets.push_back(0);
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
typename BlockMultiArray<StoredTypes...>::difference_type BlockMultiArray<
    StoredTypes...>::distance_from_origin(const field_type<I> *ptr) const
{
  return std::distance(std::get<I>(m_values).data(), ptr);
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
common::ArrayView<const typename BlockMultiArray<StoredTypes...>::template field_type<I>, _1D, Uint>
BlockMultiArray<StoredTypes...>::front_block()
{
  common::ArrayView<const field_type<I>, _1D, Uint> b(std::get<I>(m_values).data(),
                                                      m_offsets[1] - m_offsets[0]);
  return b;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
common::ArrayView<const typename BlockMultiArray<StoredTypes...>::template field_type<I>, _1D, Uint>
BlockMultiArray<StoredTypes...>::front_block() const
{
  common::ArrayView<const field_type<I>, _1D, Uint> b(std::get<I>(m_values).data(),
                                                      m_offsets[1] - m_offsets[0]);
  return b;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
common::ArrayView<const typename BlockMultiArray<StoredTypes...>::template field_type<I>, _1D, Uint>
BlockMultiArray<StoredTypes...>::back_block()
{
  const size_type offset_one_before_last = m_offsets[m_offsets.size() - 2];
  common::ArrayView<const field_type<I>, _1D, Uint> b(std::get<I>(m_values).data() +
                                                          offset_one_before_last,
                                                      m_offsets.back() - offset_one_before_last);
  return b;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I>
common::ArrayView<const typename BlockMultiArray<StoredTypes...>::template field_type<I>, _1D, Uint>
BlockMultiArray<StoredTypes...>::back_block() const
{
  const size_type offset_one_before_last = m_offsets[m_offsets.size() - 2];
  common::ArrayView<const field_type<I>, _1D, Uint> b(std::get<I>(m_values).data() +
                                                          offset_one_before_last,
                                                      m_offsets.back() - offset_one_before_last);
  return b;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
template <Uint I, typename Compare>
void BlockMultiArray<StoredTypes...>::sort_blocks(Compare comp)
{
  if (nb_blocks() == 0)
  {
    return; // Nothing to sort
  }

  storage_type reorder_cache;

  // Perform index sort on the I-th field
  // Then reorder the entries of each block based on the obtained indices
  std::vector<Uint> sort_indices;

  for (Uint b = 0; (b + 1) < m_offsets.size(); ++b)
  {
    sort_indices.resize(m_offsets[b + 1] - m_offsets[b]);
    Uint n = 0;
    std::generate(std::begin(sort_indices), std::end(sort_indices), [&] { return n++; });

    common::ArrayView<const field_type<I>, _1D, Uint> data_block(
        std::get<I>(m_values).data() + m_offsets[b], m_offsets[b + 1] - m_offsets[b]);

    // Sort the array of indices based on comparator which operates on
    // data_block extracted from the I-th field
    std::sort(std::begin(sort_indices), std::end(sort_indices),
              [&](Uint i1, Uint i2) { return comp(data_block[i1], data_block[i2]); });

    internal::DataAlgorithm<0, BlockMultiArray<StoredTypes...>::NFields - 1,
                            StoredTypes...>::one_field_sort_cache(m_offsets, sort_indices, b,
                                                                  m_values, reorder_cache);

    /*
    std::cout << "Block: " << data_block << std::endl;
    std::cout << "Sorted indices:";
    for (auto val : sort_indices)
    {
      std::cout << " " << val;
    }
    std::cout << "\n";

    std::sort(std::get<I>(m_values).data() + m_offsets[b],
              std::get<I>(m_values).data() + m_offsets[b + 1], comp);
    */
  }
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
Real BlockMultiArray<StoredTypes...>::mem_used_mb() const
{
  Real mem_used = m_offsets.size() * sizeof(Uint);
  DataAlgo<0, NFields - 1, StoredTypes...>::mem_used_mb(m_values, mem_used);

  return 1.e-6 * mem_used;
}

// ----------------------------------------------------------------------------

template <typename... StoredTypes>
std::ostream &operator<<(std::ostream &os, const BlockMultiArray<StoredTypes...> &array)
{

  internal::DataAlgorithm<0, BlockMultiArray<StoredTypes...>::NFields - 1,
                          StoredTypes...>::print_to_stream(array.m_offsets, array.m_values, os);

  return os;
}

// ----------------------------------------------------------------------------
// Internal data structures for BlockMultiArray
// ----------------------------------------------------------------------------

namespace internal
{

struct DataAlgorithmImpl
{
  template <typename SizeType, typename T>
  static void resize_data_storage(const std::vector<SizeType> &old_offsets,
                                  const std::vector<SizeType> &new_offsets, std::vector<T> &values)
  {
    SizeType new_nb_values = 0;

    for (Uint b = 0; (b + 1) < new_offsets.size(); ++b)
    {
      new_nb_values += new_offsets[b + 1] - new_offsets[b];
    }

    std::vector<T> new_values;
    new_values.reserve(new_nb_values);
    new_values.resize(0);

    for (Uint b = 0; (b + 1) < new_offsets.size(); ++b)
    {
      const SizeType old_block_size = old_offsets[b + 1] - old_offsets[b];
      const SizeType new_block_size = new_offsets[b + 1] - new_offsets[b];

      const common::ArrayView<const T, _1D, Uint> block_data(values.data() + old_offsets[b],
                                                             old_block_size);

      if (old_block_size < new_block_size)
      {
        for (SizeType i = 0; i < old_block_size; ++i)
        {
          new_values.push_back(block_data[i]);
        }
        for (SizeType i = old_block_size; i < new_block_size; ++i)
        {
          new_values.push_back(T());
        }
      }
      else
      {
        for (SizeType i = 0; i < new_block_size; ++i)
        {
          new_values.push_back(block_data[i]);
        }
      }

    } // Loop over blocks

    values.swap(new_values);
  }

  template <typename SizeType, typename T>
  static void remove_data_blocks_storage(const std::vector<bool> &remove_block,
                                         const std::vector<SizeType> &old_offsets,
                                         std::vector<T> &values)
  {
    const Uint nb_old_blocks = (old_offsets.size() < 1) ? 0 : (old_offsets.size() - 1);

    if (nb_old_blocks == 0)
    {
      std::cerr << "BlockMultiArray::remove_blocks: error, no blocks to "
                   "remove\n";
      return;
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
      const Uint block_data_len = old_offsets[b + 1] - old_offsets[b];
      // If facet data should be preserved and there were
      // some data removed prior to processing the facet,
      // copy the data on the first available position
      if (!remove_block[b])
      {
        nb_entries_after_removal += block_data_len;
        if (fill_pos != old_offsets[b])
        {
          for (Uint i = 0; i < block_data_len; ++i)
          {
            values[fill_pos] = values[valid_data_pos];
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

    values.resize(nb_entries_after_removal);
  }

  template <typename SizeType, typename T>
  static void sort_block_data_impl(const std::vector<Uint> &sort_indices,
                                   std::vector<T> &sort_cache,
                                   common::ArrayView<T, _1D, Uint> &data_block)
  {
    sort_cache.resize(sort_indices.size());
    for (Uint i = 0; i < sort_indices.size(); ++i)
    {
      sort_cache[i] = data_block[i];
    }

    for (Uint i = 0; i < sort_indices.size(); ++i)
    {
      data_block[i] = sort_cache[sort_indices[i]];
    }
  }

  template <typename SizeType, typename T>
  static void values_to_stream(const std::vector<SizeType> &offsets, const std::vector<T> &values,
                               std::ostream &os)
  {
    for (SizeType ib = 0; (ib + 1) < offsets.size(); ++ib)
    {
      const common::ArrayView<const T, _1D, Uint> block(values.data() + offsets[ib],
                                                        offsets[ib + 1] - offsets[ib]);
      for (SizeType e = 0; e < block.size(); ++e)
      {
        os << block[e] << " ";
      }
      os << "\n";
    }
  }
};

// ----------------------------------------------------------------------------

template <Uint I, Uint N, typename... T>
struct DataAlgorithm
{
  using vector_tuple_type = typename meta::TransformToVector<T...>::type;
  using field_type        = typename meta::tuple_nth_elem_type<I, T...>::type;

  static void resize(vector_tuple_type &data, const Uint size)
  {
    std::get<I>(data).resize(size);
    DataAlgorithm<I + 1, N, T...>::resize(data, size);
  }

  static void reserve(vector_tuple_type &data, const Uint size)
  {
    std::get<I>(data).reserve(size);
    DataAlgorithm<I + 1, N, T...>::reserve(data, size);
  }

  static void copy(vector_tuple_type &data_to, const vector_tuple_type &data_from)
  {
    std::get<I>(data_to).resize(std::get<I>(data_from).size());
    std::get<I>(data_to) = std::get<I>(data_from);

    DataAlgorithm<I + 1, N, T...>::copy(data_to, data_from);
  }

  static void build_values(vector_tuple_type &data_to,
                           std::tuple<std::unique_ptr<std::vector<T>>...> &&values)
  {
    std::get<I>(data_to).swap(*(std::get<I>(values)));
    DataAlgorithm<I + 1, N, T...>::build_values(data_to, std::move(values));
  }

  static void clear(vector_tuple_type &data)
  {
    std::get<I>(data).clear();
    DataAlgorithm<I + 1, N, T...>::clear(data);
  }

  static void extend_data_array(const Uint num_values, vector_tuple_type &data)
  {
    for (Uint i = 0; i < num_values; ++i)
    {
      std::get<I>(data).push_back(field_type());
    }
    DataAlgorithm<I + 1, N, T...>::extend_data_array(num_values, data);
  }

  template <typename SizeType>
  static void resize_data_storage(const std::vector<SizeType> &old_offsets,
                                  const std::vector<SizeType> &new_offsets, vector_tuple_type &data)
  {
    DataAlgorithmImpl::resize_data_storage(old_offsets, new_offsets, std::get<I>(data));
    DataAlgorithm<I + 1, N, T...>::resize_data_storage(old_offsets, new_offsets, data);
  }

  template <typename SizeType>
  static void print_to_stream(const std::vector<SizeType> &offsets, const vector_tuple_type &data,
                              std::ostream &os)
  {
    os << "BlockMultiArray: Field " << I << "\n";
    DataAlgorithmImpl::values_to_stream(offsets, std::get<I>(data), os);
    DataAlgorithm<I + 1, N, T...>::print_to_stream(offsets, data, os);
  }

  template <typename SizeType>
  static void remove_data_blocks_storage(const std::vector<bool> &remove_block,
                                         const std::vector<SizeType> &old_offsets,
                                         vector_tuple_type &data)
  {
    DataAlgorithmImpl::remove_data_blocks_storage(remove_block, old_offsets, std::get<I>(data));
    DataAlgorithm<I + 1, N, T...>::remove_data_blocks_storage(remove_block, old_offsets, data);
  }

  template <typename SizeType>
  static void one_field_sort_cache(const std::vector<SizeType> &offsets,
                                   const std::vector<Uint> &sort_indices, const Uint block_id,
                                   vector_tuple_type &data_to_sort, vector_tuple_type &sort_cache)
  {
    std::vector<field_type> &one_field_data_to_sort = std::get<I>(data_to_sort);
    std::vector<field_type> &one_field_sort_cache   = std::get<I>(sort_cache);

    common::ArrayView<field_type, _1D, Uint> data_block(one_field_data_to_sort.data() +
                                                            offsets[block_id],
                                                        offsets[block_id + 1] - offsets[block_id]);

    DataAlgorithmImpl::sort_block_data_impl<field_type>(sort_indices, one_field_sort_cache,
                                                        data_block);
    DataAlgorithm<I + 1, N, T...>::one_field_sort_cache(offsets, sort_indices, block_id,
                                                        data_to_sort, sort_cache);
  }

  static void mem_used_mb(const vector_tuple_type &data, Real &mem)
  {
    mem += std::get<I>(data).size() * sizeof(field_type);

    DataAlgorithm<I + 1, N, T...>::mem_used_mb(data, mem);
  }
};

// ----------------------------------------------------------------------------

template <Uint N, typename... T>
struct DataAlgorithm<N, N, T...>
{
  using vector_tuple_type = typename meta::TransformToVector<T...>::type;
  using field_type        = typename meta::tuple_nth_elem_type<N, T...>::type;

  static void resize(vector_tuple_type &data, const Uint size)
  {
    std::get<N>(data).resize(size);
  }

  static void reserve(vector_tuple_type &data, const Uint size)
  {
    std::get<N>(data).reserve(size);
  }

  static void copy(vector_tuple_type &data_to, const vector_tuple_type &data_from)
  {
    std::get<N>(data_to).resize(std::get<N>(data_from).size());
    std::get<N>(data_to) = std::get<N>(data_from);
  }

  static void build_values(vector_tuple_type &data_to,
                           std::tuple<std::unique_ptr<std::vector<T>>...> &&values)
  {
    std::get<N>(data_to).swap(*(std::get<N>(values)));
  }

  static void clear(vector_tuple_type &data)
  {
    std::get<N>(data).clear();
  }

  static void extend_data_array(const Uint num_values, vector_tuple_type &data)
  {
    for (Uint i = 0; i < num_values; ++i)
    {
      std::get<N>(data).push_back(field_type());
    }
  }

  template <typename SizeType>
  static void resize_data_storage(const std::vector<SizeType> &old_offsets,
                                  const std::vector<SizeType> &new_offsets, vector_tuple_type &data)
  {
    DataAlgorithmImpl::resize_data_storage(old_offsets, new_offsets, std::get<N>(data));
  }

  template <typename SizeType>
  static void print_to_stream(const std::vector<SizeType> &offsets, const vector_tuple_type &data,
                              std::ostream &os)
  {
    os << "BlockMultiArray: Field " << N << "\n";
    DataAlgorithmImpl::values_to_stream(offsets, std::get<N>(data), os);
  }

  template <typename SizeType>
  static void remove_data_blocks_storage(const std::vector<bool> &remove_block,
                                         const std::vector<SizeType> &old_offsets,
                                         vector_tuple_type &data)
  {
    DataAlgorithmImpl::remove_data_blocks_storage(remove_block, old_offsets, std::get<N>(data));
  }

  template <typename SizeType>
  static void reorder_block_entries(const std::vector<Uint> &sort_indices,
                                    vector_tuple_type &sort_cache,
                                    common::ArrayView<field_type, _1D, Uint> &block_data)
  {
  }

  template <typename SizeType>
  static void one_field_sort_cache(const std::vector<SizeType> &offsets,
                                   const std::vector<Uint> &sort_indices, const Uint block_id,
                                   vector_tuple_type &data_to_sort, vector_tuple_type &sort_cache)
  {
    std::vector<field_type> &one_field_data_to_sort = std::get<N>(data_to_sort);
    std::vector<field_type> &one_field_sort_cache   = std::get<N>(sort_cache);

    common::ArrayView<field_type, _1D, Uint> data_block(one_field_data_to_sort.data() +
                                                            offsets[block_id],
                                                        offsets[block_id + 1] - offsets[block_id]);

    DataAlgorithmImpl::sort_block_data_impl<field_type>(sort_indices, one_field_sort_cache,
                                                        data_block);
  }

  static void mem_used_mb(const vector_tuple_type &data, Real &mem)
  {
    mem += std::get<N>(data).size() * sizeof(field_type);
  }
};

} // namespace internal

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit

#endif // PDEKIT_Common_BlockMultiArray_hpp
