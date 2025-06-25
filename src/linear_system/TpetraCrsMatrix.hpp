#ifndef PDEKIT_Trilinos_Tpetra_Crs_Matrix_hpp
#define PDEKIT_Trilinos_Tpetra_Crs_Matrix_hpp

#include <iostream>
#include <memory>
#include <mutex>
#include <vector>

#include "PDEKit_Config.hpp"
#include "common/MPI/MPIEnv.hpp"
#include "common/PDEKit.hpp"
#include "linear_system/TpetraComm.hpp"
#include "linear_system/TpetraDofMap.hpp"
#include "linear_system/TpetraFacadeTools.hpp"
#include "linear_system/TpetraMultiVector.hpp"
#include "linear_system/TpetraTraits.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseMatView.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "Teuchos_RCP.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
#endif

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------
//             Forward declarations to enable friend functions
// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
class TpetraMultiVector;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
class TpetraCrsMatrix;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
std::ostream &operator<<(
    std::ostream &os,
    const TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &tril_crs_mat);

// ----------------------------------------------------------------------------

template <typename Scalar        = TpetraDefaultTraits::Scalar,
          typename LocalOrdinal  = TpetraDefaultTraits::LocalOrdinal,
          typename GlobalOrdinal = TpetraDefaultTraits::GlobalOrdinal,
          typename Node          = TpetraDefaultTraits::Node,
          const bool classic     = TpetraDefaultTraits::classic>
class TpetraCrsMatrix

{
  public:
#if PDEKIT_HAVE_TRILINOS
  using trilinos_matrix_type =
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
#else
  using trilinos_matrix_type = TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
#endif

#if PDEKIT_HAVE_TRILINOS
  using scalar_type         = typename trilinos_matrix_type::scalar_type;
  using local_ordinal_type  = typename trilinos_matrix_type::local_ordinal_type;
  using global_ordinal_type = typename trilinos_matrix_type::global_ordinal_type;
  using node_type           = typename trilinos_matrix_type::node_type;
#else
  using scalar_type          = Scalar;
  using local_ordinal_type   = LocalOrdinal;
  using global_ordinal_type  = GlobalOrdinal;
  using node_type            = Node;
#endif

  /// Default constructor
  TpetraCrsMatrix();

  /// Construct from pointers to bare Tpetra objects

#if PDEKIT_HAVE_TRILINOS
  TpetraCrsMatrix(
      const Teuchos::RCP<Teuchos::Comm<GlobalOrdinal> const> &comm,
      const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
      const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>>
          &matrix);
#endif

  /// Constructor with predefined number of rows
  TpetraCrsMatrix(const Uint global_nb_rows, const Uint nb_entries_per_row = 8);

  /// Constructor with predefined number of columns for each row
  TpetraCrsMatrix(const graph::Graph<GlobalOrdinal> &matrix_nonzeros);

  /// Constructor with predefined number of columns for each row
  TpetraCrsMatrix(const math::MatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros);

  /// Constructor with predefined number of columns for each row
  TpetraCrsMatrix(const math::BlockMatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros);

  /// Destructor
  ~TpetraCrsMatrix();

/// Initialize from (trilinos) shared pointers to Tpetra objects
#if PDEKIT_HAVE_TRILINOS
  void init(
      const Teuchos::RCP<Teuchos::Comm<GlobalOrdinal> const> &comm,
      const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
      const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>>
          &matrix);
#endif

  /// Initialize by providing global number of elements
  void init(const Uint global_nb_rows, const Uint nb_entries_per_row = 8);

  /// Initialize by providing nonzero entries in a sparsity graph
  void init(const graph::Graph<GlobalOrdinal> &matrix_nonzeros);

  /// Initialize by providing nonzero entries in a matrix sparsity pattern
  void init(const math::MatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros);

  /// Initialize by providing nonzero entries in a matrix sparsity pattern
  void init(const math::BlockMatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros);

  /// Get the number of global elements in the matrix
  Uint nb_global_elem() const;

  /// Get the number of rows in the matrix across all processes
  Uint rows_global() const;

  /// Get the number of columns in the matrix across all processes
  Uint cols_global() const;

  /// Get the local number of rows (owned by the calling process)
  Uint rows_local() const;

  /// Get the local number of columns (owned by the calling process)
  Uint cols_local() const;

  /// Lock the matrix (both structure and entries)
  /// This means that no value of the matrix can be changed
  void lock();

  /// Unlock the matrix - this enables adding entries to the
  /// matrix or modifying its structure
  void unlock();

  /// Insert values in one row
  void insert_values_in_row(const GlobalOrdinal row_idx, std::vector<Scalar> &values,
                            std::vector<GlobalOrdinal> &indices);

  /// Accumulate values in one row
  void add_values_to_row(const GlobalOrdinal row_idx, std::vector<Scalar> &values,
                         std::vector<GlobalOrdinal> &indices);

  /// Accumulate a series of values
  void add_values(const std::vector<std::tuple<Uint, Uint, Real>> &buffer);

  /// Get one row of the matrix
  void row(const GlobalOrdinal row_idx, math::DenseConstVecView<GlobalOrdinal> &indices,
           math::DenseConstVecView<Scalar> &values) const;

  /// Get a block of data
  /// @param row_range    ... range of row indices to retrieve
  /// @param bottom_right ... range of column indices to retrieve
  /// @param block_data   ... block of data that stores the output values
  void get_block(const common::Range1D<GlobalOrdinal> row_range,
                 const common::Range1D<GlobalOrdinal> col_range,
                 math::DenseMatView<Scalar> &block_data) const;

  /// Insert a block of data
  /// @param row_range    ... range of row indices to retrieve
  /// @param bottom_right ... range of column indices to retrieve
  /// @param block_data   ... block of data that stores the output values
  void insert_block(const common::Range1D<GlobalOrdinal> row_range,
                    const common::Range1D<GlobalOrdinal> col_range,
                    const math::DenseMatView<Scalar> &block_data);

  /// Get the underlying map
  const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> map() const;

  /// Get the domain map of this matrix
  /// @note This should be used to construct the vectors
  ///       for linear systems involving this matrix!
  const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> domain_map() const;

  /// Get the domain map of this matrix
  /// @note This should be used to construct the vectors
  ///       resulting from applying this matrix as linear operator
  const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> range_map() const;

  /// Apply the matrix as linear operator to vector
  /// Computes Y := beta*Y + alpha*Op(A)*X with Op(A) = A or Op(A) =
  /// transpose(A)
  void apply(const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &X,
             TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &Y,
             const bool transp = false, Scalar alpha = 1.0, Scalar beta = 0.0) const;

  /// Assign the same scalar value to all entries in the sparse matrix
  void fill(const Scalar value);

  /// Print the structure of the matrix to file
  void print_structure_to_file(const std::string &filename) const;

  /// Print info
  void print_info() const
  {
#if PDEKIT_HAVE_TRILINOS
    std::cout << "LSS with trilinos support" << std::endl;
#else
    std::cout << "LSS without trilinos support" << std::endl;
#endif

    if (mat_storage_is_dynamic())
    {
      std::cout << "Tpetra CRS matrix storage profile: dynamic" << std::endl;
    }
    else if (mat_storage_is_static())
    {
      std::cout << "Tpetra CRS matrix storage profile: static" << std::endl;
    }
    else
    {
      std::cout << "Tpetra CRS matrix storage profile: unknown!" << std::endl;
    }
  }

  /// Print the values of the matrix
  friend std::ostream &operator<<<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
      std::ostream &os,
      const TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &tril_crs_mat);

  private:
  /// FRIENDS
  friend class internal::TpetraInternalAccess;

  /// Helper functions to check the storage type of Tpetra matrix
  bool mat_storage_is_dynamic() const;
  bool mat_storage_is_static() const;

#if PDEKIT_HAVE_TRILINOS
  /// TYPES
  using trilinos_map_type = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;

  /// Export pointer to the real tpetra object
  /// This is needed for interaction with other Tpetra wrappers
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> get_matrix()
      const;

  /// Trilinos communicator object
  /// Takes MPI communicator in constructor
  TpetraComm<GlobalOrdinal> m_comm;

  /// Map - needed to create the matrix
  TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> m_map;

  /// Sparse distributed matrix from Trilinos
  Teuchos::RCP<trilinos_matrix_type> m_tpetra_matrix;

  /// Mutex to make some methods thread-safe
  mutable std::mutex m_mutex;
#endif
};

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraCrsMatrix()
    : m_comm(), m_map(), m_tpetra_matrix(nullptr)
{
}
#else
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraCrsMatrix()
{
}
#endif

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraCrsMatrix(
    const Teuchos::RCP<Teuchos::Comm<GlobalOrdinal> const> &comm,
    const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
    const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>>
        &matrix)
    : m_comm(comm), m_map(map), m_tpetra_matrix(matrix)
{
}
#endif

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraCrsMatrix(
    const Uint global_nb_rows, const Uint nb_entries_per_row)
{
#if PDEKIT_HAVE_TRILINOS
  this->init(global_nb_rows, nb_entries_per_row);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraCrsMatrix(
    const graph::Graph<GlobalOrdinal> &matrix_nonzeros)
{
  this->init(matrix_nonzeros);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraCrsMatrix(
    const math::MatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros)
{
#if PDEKIT_HAVE_TRILINOS
  this->init(matrix_nonzeros);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraCrsMatrix(
    const math::BlockMatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros)
{
#if PDEKIT_HAVE_TRILINOS
  this->init(matrix_nonzeros);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::~TpetraCrsMatrix()
{
}

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::init(
    const Teuchos::RCP<Teuchos::Comm<GlobalOrdinal> const> &comm,
    const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
    const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>>
        &matrix)
{
  m_comm          = comm;
  m_map           = map;
  m_tpetra_matrix = matrix;
}
#endif

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::init(
    const Uint global_nb_rows, const Uint nb_entries_per_row)
{
#if PDEKIT_HAVE_TRILINOS
  std::cout << "Initializing Tpetra CRS matrix with dynamic structure ... " << std::endl;
  m_comm = TpetraComm<GlobalOrdinal>(common::mpi::MPIEnv::instance().comm());
  m_map = TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>(static_cast<LocalOrdinal>(global_nb_rows),
                                                          0, m_comm);
  m_tpetra_matrix = Teuchos::RCP<trilinos_matrix_type>(new trilinos_matrix_type(
      internal::TpetraInternalAccess::get_dof_map(m_map), nb_entries_per_row));
  std::cout << "finished initializing Tpetra CRS matrix with dynamic structure." << std::endl;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::init(
    const graph::Graph<GlobalOrdinal> &matrix_nonzeros)
{
#if PDEKIT_HAVE_TRILINOS
  std::cout << "Initializing Tpetra CRS matrix from graph ... " << std::endl;

  m_comm = TpetraComm<GlobalOrdinal>(common::mpi::MPIEnv::instance().comm());
  m_map.reset(static_cast<LocalOrdinal>(matrix_nonzeros.nb_vertices()), 0, m_comm);

  using iaccess = internal::TpetraInternalAccess;
  TpetraFacadeTools::create_matrix_sparsity(iaccess::get_dof_map(m_map), matrix_nonzeros,
                                            m_tpetra_matrix);

  std::cout << "finished initializing Tpetra CRS matrix from graph." << std::endl;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::init(
    const math::MatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros)
{
#if PDEKIT_HAVE_TRILINOS
  std::cout << "Constructing Tpetra CRS matrix from matrix sparsity pattern ... " << std::endl;

  m_comm = TpetraComm<GlobalOrdinal>(common::mpi::MPIEnv::instance().comm());
  m_map.reset(static_cast<LocalOrdinal>(matrix_nonzeros.nb_lines()), 0, m_comm);

  using iaccess = internal::TpetraInternalAccess;
  TpetraFacadeTools::create_matrix_sparsity(iaccess::get_dof_map(m_map), matrix_nonzeros,
                                            m_tpetra_matrix);

  std::cout << "finished constructing Tpetra CRS matrix from sparsity pattern." << std::endl;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::init(
    const math::BlockMatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros)
{
#if PDEKIT_HAVE_TRILINOS
  std::cout << "Constructing Tpetra CRS matrix from block matrix sparsity "
               "pattern ... "
            << std::endl;

  m_comm = TpetraComm<GlobalOrdinal>(common::mpi::MPIEnv::instance().comm());

  const common::ArrayShape<_2D, Uint> mat_shape = matrix_nonzeros.size();
  m_map.reset(static_cast<LocalOrdinal>(mat_shape.size(0)), 0, m_comm);

  using iaccess = internal::TpetraInternalAccess;
  TpetraFacadeTools::create_matrix_sparsity(iaccess::get_dof_map(m_map), matrix_nonzeros,
                                            m_tpetra_matrix);

  std::cout << "finished constructing Tpetra CRS matrix from block sparsity "
               "pattern."
            << std::endl;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Uint TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::nb_global_elem() const
{
#if PDEKIT_HAVE_TRILINOS
  return m_map.global_num_elements();
#else
  return 0;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Uint TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::rows_global() const
{
#if PDEKIT_HAVE_TRILINOS
  return static_cast<Uint>(m_tpetra_matrix->getGlobalNumRows());
#else
  return 0;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Uint TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::cols_global() const
{
#if PDEKIT_HAVE_TRILINOS
  return static_cast<Uint>(m_tpetra_matrix->getGlobalNumCols());
#else
  return 0;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Uint TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::rows_local() const
{
#if PDEKIT_HAVE_TRILINOS
  return static_cast<Uint>(m_tpetra_matrix->getNodeNumRows());
#else
  return 0;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Uint TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::cols_local() const
{
#if PDEKIT_HAVE_TRILINOS
  return static_cast<Uint>(m_tpetra_matrix->getNodeNumCols());
#else
  return 0;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::lock()
{
#if PDEKIT_HAVE_TRILINOS
  m_tpetra_matrix->fillComplete();
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::unlock()
{
#if PDEKIT_HAVE_TRILINOS
  m_tpetra_matrix->resumeFill();
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::insert_values_in_row(
    const GlobalOrdinal row_idx, std::vector<Scalar> &values, std::vector<GlobalOrdinal> &indices)
{
#if PDEKIT_HAVE_TRILINOS
  /*
  const GlobalOrdinal global_row_idx = static_cast<GlobalOrdinal>(row_idx);
  m_index_buffer.resize(indices.size());
  for (Uint i = 0; i < indices.size(); ++i)
  {
    m_index_buffer[i] = static_cast<GlobalOrdinal>(indices[i]);
  }
  */
  Teuchos::ArrayView<GlobalOrdinal> indices_view(indices.data(), indices.size());
  Teuchos::ArrayView<Scalar> values_view(values.data(), values.size());

  if (mat_storage_is_dynamic())
  {
    m_tpetra_matrix->insertGlobalValues(row_idx, indices_view, values_view);
  }
  else // Should be static
  {
    m_tpetra_matrix->replaceGlobalValues(row_idx, indices_view, values_view);
  }
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::add_values_to_row(
    const GlobalOrdinal row_idx, std::vector<Scalar> &values, std::vector<GlobalOrdinal> &indices)
{
#if PDEKIT_HAVE_TRILINOS
  /*
  const GlobalOrdinal global_row_idx = static_cast<GlobalOrdinal>(row_idx);
  m_index_buffer.resize(indices.size());
  for (Uint i = 0; i < indices.size(); ++i)
  {
    m_index_buffer[i] = static_cast<GlobalOrdinal>(indices[i]);
  }
  */

  Teuchos::ArrayView<GlobalOrdinal> indices_view(indices.data(), indices.size());
  Teuchos::ArrayView<Scalar> values_view(values.data(), values.size());

  // This should work with both static and dynamic profiles
  m_tpetra_matrix->sumIntoGlobalValues(row_idx, indices_view, values_view);
// m_tpetra_matrix->sumIntoLocalValues(row_idx, indices_view, values_view);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::add_values(
    const std::vector<std::tuple<Uint, Uint, Real>> &buffer)
{
#if PDEKIT_HAVE_TRILINOS
  std::lock_guard<std::mutex> lock(m_mutex);

  GlobalOrdinal row_idx, col_idx;
  Scalar value;

  Teuchos::ArrayView<GlobalOrdinal> col_index_view(&col_idx, 1);
  Teuchos::ArrayView<Scalar> value_view(&value, 1);

  trilinos_matrix_type &matrix = (*m_tpetra_matrix);

  // This should work with both static and dynamic profiles
  for (Uint i = 0; i < buffer.size(); ++i)
  {
    row_idx = static_cast<GlobalOrdinal>(std::get<0>(buffer[i]));
    col_idx = static_cast<GlobalOrdinal>(std::get<1>(buffer[i]));
    value   = static_cast<Scalar>(std::get<2>(buffer[i]));

    matrix.sumIntoGlobalValues(row_idx, col_index_view, value_view);
    // matrix.sumIntoLocalValues(row_idx, col_index_view, value_view);
  }
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::row(
    const GlobalOrdinal row_idx, math::DenseConstVecView<GlobalOrdinal> &indices,
    math::DenseConstVecView<Scalar> &values) const
{
#if PDEKIT_HAVE_TRILINOS
  Teuchos::ArrayView<const GlobalOrdinal> teuchos_indices;
  Teuchos::ArrayView<const Scalar> teuchos_values;

  m_tpetra_matrix->getLocalRowView(row_idx, teuchos_indices, teuchos_values);
  indices = math::DenseConstVecView<GlobalOrdinal>(
      teuchos_indices.getRawPtr(), teuchos_indices.getRawPtr() + teuchos_indices.size() - 1);
  values = math::DenseConstVecView<Real>(teuchos_values.getRawPtr(),
                                         teuchos_values.getRawPtr() + teuchos_values.size() - 1);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::get_block(
    const common::Range1D<GlobalOrdinal> row_range, const common::Range1D<GlobalOrdinal> col_range,
    math::DenseMatView<Scalar> &block_data) const
{
#if PDEKIT_HAVE_TRILINOS
  TpetraFacadeTools::get_block(*m_tpetra_matrix, row_range, col_range, block_data);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::insert_block(
    const common::Range1D<GlobalOrdinal> row_range, const common::Range1D<GlobalOrdinal> col_range,
    const math::DenseMatView<Scalar> &block_data)
{
#if PDEKIT_HAVE_TRILINOS
  TpetraFacadeTools::insert_block(row_range, col_range, block_data, *m_tpetra_matrix);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> TpetraCrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::map() const
{
#if PDEKIT_HAVE_TRILINOS
  return m_map;
#else
  TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> map;
  return map;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> TpetraCrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::domain_map() const
{
#if PDEKIT_HAVE_TRILINOS
  TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> domain_map(m_tpetra_matrix->getDomainMap());
  return domain_map;
#else
  TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> map;
  return map;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> TpetraCrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::range_map() const
{
#if PDEKIT_HAVE_TRILINOS
  TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> domain_map(m_tpetra_matrix->getRangeMap());
  return domain_map;
#else
  TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> map;
  return map;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::apply(
    const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &X,
    TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &Y, const bool transp,
    Scalar alpha, Scalar beta) const
{
#if PDEKIT_HAVE_TRILINOS
  const Teuchos::ETransp mode = (transp == true) ? Teuchos::TRANS : Teuchos::NO_TRANS;

  m_tpetra_matrix->apply(*internal::TpetraInternalAccess::get_vector_data(X),
                         *internal::TpetraInternalAccess::get_vector_data(Y), mode, alpha, beta);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::fill(const Scalar value)
{
#if PDEKIT_HAVE_TRILINOS
  m_tpetra_matrix->setAllToScalar(value);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::print_structure_to_file(
    const std::string &filename) const
{
#if PDEKIT_HAVE_TRILINOS
  std::cout << "TpetraCrsMatrix::print_structure_to_file:" << std::endl
            << "Tpetra and Ifpack2 do not enable printing to postscript yet" << std::endl;
// Ifpack_PrintSparsity(*m_epetra_matrix, filename.c_str());
#endif
}

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> TpetraCrsMatrix<
    Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::get_matrix() const
{
  return m_tpetra_matrix;
}
#endif

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
inline bool TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                            classic>::mat_storage_is_dynamic() const
{
#if PDEKIT_HAVE_TRILINOS
  return (m_tpetra_matrix->getProfileType() == Tpetra::DynamicProfile);
#else
  return false;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
inline bool TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                            classic>::mat_storage_is_static() const
{
#if PDEKIT_HAVE_TRILINOS
  return (m_tpetra_matrix->getProfileType() == Tpetra::StaticProfile);
#else
  return false;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
std::ostream &operator<<(
    std::ostream &os,
    const TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &tril_crs_mat)
{
  os << *tril_crs_mat.m_tpetra_matrix;
  return os;
}

// ----------------------------------------------------------------------------

// Defer instantiation of Real TpetraCrsMatrix to the cpp file
extern template class TpetraCrsMatrix<Real>;

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_Trilinos_Tpetra_Crs_Matrix_hpp
