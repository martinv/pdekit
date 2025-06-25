#ifndef PDEKIT_Trilinos_Tpetra_Crs_Matrix_Diag_Blocks_hpp
#define PDEKIT_Trilinos_Tpetra_Crs_Matrix_Diag_Blocks_hpp

#include <vector>

#include "linear_system/TpetraCrsMatrix.hpp"

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------

template <typename Scalar        = TpetraDefaultTraits::Scalar,
          typename LocalOrdinal  = TpetraDefaultTraits::LocalOrdinal,
          typename GlobalOrdinal = TpetraDefaultTraits::GlobalOrdinal,
          typename Node          = TpetraDefaultTraits::Node,
          const bool classic     = TpetraDefaultTraits::classic>
class TpetraCrsMatrixDiagBlocks
{
  public:
  // TYPEDEFS
  using trilinos_matrix_type = TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;

  typedef typename trilinos_matrix_type::scalar_type scalar_type;
  typedef typename trilinos_matrix_type::local_ordinal_type local_ordinal_type;
  typedef typename trilinos_matrix_type::global_ordinal_type global_ordinal_type;
  typedef typename trilinos_matrix_type::node_type node_type;

  /// Default constructor
  TpetraCrsMatrixDiagBlocks();

  /// Copy constructor
  TpetraCrsMatrixDiagBlocks(const TpetraCrsMatrixDiagBlocks &other);

  /// Destructor
  ~TpetraCrsMatrixDiagBlocks();

  /// Assignment operator
  TpetraCrsMatrixDiagBlocks &operator=(const TpetraCrsMatrixDiagBlocks &rhs);

  /// Set block sizes
  void set_block_sizes(const std::vector<Uint> &block_sizes);

  /// Return the number of blocks
  const Uint num_blocks() const;

  /// Return the size of one block
  const Uint block_size(const Uint block_id) const;

  /// Retrieve the data of one block
  /// Values in output vector are stored row-wise
  void get_block_data(const trilinos_matrix_type &mat, const Uint block_id,
                      std::vector<scalar_type> &values);

  private:
  /// Information about block sizes and positions
  std::vector<Uint> m_block_offsets;

  /// Block indices to retrieve
  std::vector<local_ordinal_type> m_block_indices;
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrixDiagBlocks<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                          classic>::TpetraCrsMatrixDiagBlocks()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrixDiagBlocks<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
    TpetraCrsMatrixDiagBlocks(const TpetraCrsMatrixDiagBlocks &other)
{
  m_block_offsets.resize(other.m_block_offsets.size());
  m_block_offsets = other.m_block_offsets;

  m_block_indices.resize(other.m_block_indices.size());
  m_block_indices = other.m_block_indices;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrixDiagBlocks<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                          classic>::~TpetraCrsMatrixDiagBlocks()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrixDiagBlocks<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>
    &TpetraCrsMatrixDiagBlocks<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::operator=(
        const TpetraCrsMatrixDiagBlocks &rhs)
{
  m_block_offsets.resize(rhs.m_block_offsets.size());
  m_block_offsets = rhs.m_block_offsets;

  m_block_indices.resize(rhs.m_block_indices.size());
  m_block_indices = rhs.m_block_indices;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrixDiagBlocks<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::set_block_sizes(
    const std::vector<Uint> &block_sizes)
{
  m_block_offsets.resize(block_sizes.size() + 1);
  m_block_offsets[0] = 0;

  for (Uint i = 0; i < block_sizes.size(); ++i)
  {
    m_block_offsets[i + 1] = m_block_offsets[i] + block_sizes[i];
  }
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
const Uint TpetraCrsMatrixDiagBlocks<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                                     classic>::num_blocks() const
{
  return (m_block_offsets.size() < 1) ? 0 : m_block_offsets.size() - 1;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
const Uint TpetraCrsMatrixDiagBlocks<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                                     classic>::block_size(const Uint block_id) const
{
  return m_block_offsets[block_id + 1] - m_block_offsets[block_id];
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrixDiagBlocks<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::get_block_data(
    const trilinos_matrix_type &mat, const Uint block_id, std::vector<scalar_type> &block_values)
{
  const Uint block_size = m_block_offsets[block_id + 1] - m_block_offsets[block_id];
  m_block_indices.resize(block_size);
  for (Uint i = 0; i < block_size; ++i)
  {
    m_block_indices[i] = m_block_offsets[block_id] + i;
  }

  math::DenseConstVecView<GlobalOrdinal> col_indices;
  math::DenseConstVecView<Scalar> values;

  block_values.resize(block_size * block_size);

  for (Uint i = 0; i < block_size; ++i)
  {
    mat.row(m_block_indices[i], col_indices, values);
    for (Uint j = 0; j < block_size; ++j)
    {
      block_values[i * block_size + j] = values[j];
    }
  }
}

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_Trilinos_Tpetra_Crs_Matrix_Diag_Blocks_hpp
