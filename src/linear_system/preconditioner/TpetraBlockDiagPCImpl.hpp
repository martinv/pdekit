#ifndef PDEKIT_Linear_System_Tpetra_Block_Diag_PC_Impl_hpp
#define PDEKIT_Linear_System_Tpetra_Block_Diag_PC_Impl_hpp

#include <iostream>

#include "linear_system/TpetraInternalAccess.hpp"
#include "linear_system/preconditioner/TrilinosPC.hpp"
#include "math/DenseDMatArray.hpp"

#if PDEKIT_HAVE_TRILINOS

#include "Ifpack2_Factory.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_oblackholestream.hpp"

namespace pdekit
{

namespace ls
{

namespace detail
{

// ----------------------------------------------------------------------------

// Define a class for our user-defined operator.
// In this case, it is the tridiagonal matrix [-1,2,-1].
// You may define it to be whatever you like.
//
// In general, Trilinos does NOT require the user to deal with MPI
// communication explicitly.  If you want to define your own operator
// though, there's no getting around it.  Fortunately, Trilinos makes
// this relatively straightforward with the use of Map and Import
// objects.  All you have to do is define your initial data
// distribution (which is a block row distribution here), and the data
// distribution you need to perform the operations of your
// matrix-vector multiply.  For instance, when performing a
// matrix-vector multiply with a tridiagonal matrix (with a block row
// distribution), each process needs to know the last element owned by
// the previous process and the first element owned by the next
// process.
//
// If you are only interested in running the code sequentially, you
// may safely ignore everything here regarding Map and Import objects.

template <typename Scalar        = Tpetra::Details::DefaultTypes::scalar_type,
          typename LocalOrdinal  = Tpetra::Details::DefaultTypes::local_ordinal_type,
          typename GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
          typename Node          = Tpetra::Details::DefaultTypes::node_type>

class TpetraBlockDiagPCImpl : public Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{
  public:
  // Tpetra::Operator subclasses should always define these four typedefs.

  using TpetraOp            = Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using scalar_type         = typename TpetraOp::scalar_type;
  using local_ordinal_type  = typename TpetraOp::local_ordinal_type;
  using global_ordinal_type = typename TpetraOp::global_ordinal_type;
  using node_type           = typename TpetraOp::node_type;

  /// The type of the input and output arguments of apply().
  using MV = Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  /// The Map specialization used by this class.
  using map_type = Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>;

  private:
  /// This is an implementation detail; users don't need to see it.
  using import_type = Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type>;

  public:
  /// Constructor
  /// n: Global number of rows and columns in the operator.
  /// comm: The communicator over which to distribute those rows and columns.
  TpetraBlockDiagPCImpl(const global_ordinal_type n,
                        const Teuchos::RCP<const Teuchos::Comm<int>> comm);

  /// Destructor
  ~TpetraBlockDiagPCImpl() override;

  void initialize();

  /// Compute the preconditioner
  void compute_impl(
      const std::shared_ptr<const TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A);

  /// Set preconditioner parameter
  void setParameters(const Teuchos::ParameterList &params);

  // These functions are required since we inherit from Tpetra::Operator

  /// Get the domain Map of this Operator subclass.
  Teuchos::RCP<const map_type> getDomainMap() const override;

  /// Get the range Map of this Operator subclass.
  Teuchos::RCP<const map_type> getRangeMap() const override;

  /// Compute Y := alpha Op X + beta Y.
  ///
  /// We ignore the cases alpha != 1 and beta != 0 for simplicity.
  void apply(const MV &X, MV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
             scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
             scalar_type beta  = Teuchos::ScalarTraits<scalar_type>::zero()) const override;

  /// Build the preconditioner data
  void initialize_structure(
      const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A,
      std::unique_ptr<std::vector<common::Range1D<LocalOrdinal>>> &&blocks);

  /// Print the contents of the preconditioner
  void print_impl() const;

  private:
  Teuchos::RCP<const map_type> m_opMap, m_redistMap;
  Teuchos::RCP<const import_type> m_importer;

  /// Ranges delimiting blocks
  std::vector<common::Range1D<LocalOrdinal>> m_block_ranges;

  /// The actual diagonal blocks
  math::DenseDMatArray<Scalar> m_blocks;
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraBlockDiagPCImpl(
    const global_ordinal_type n, const Teuchos::RCP<const Teuchos::Comm<int>> comm)
{
  TEUCHOS_TEST_FOR_EXCEPTION(comm.is_null(), std::invalid_argument,
                             "MyOp constructor: The input Comm object must be nonnull.");

  // Get the rank of this process and the number of processes
  // We're going to have to do something special with the first and last
  // processes
  const int my_rank   = comm->getRank();
  const int num_procs = comm->getSize();

  /*
  if (myRank == 0)
  {
    std::cout << "TpetraBlockDiagPreconditionerImpl constructor" << std::endl;
  }
  */

  // Construct a map for our block row distribution
  const global_ordinal_type index_base = 0;
  m_opMap                              = Teuchos::rcp(new map_type(n, index_base, comm));

  // Get the local number of rows
  local_ordinal_type nlocal = m_opMap->getNodeNumElements();

  //
  // Define the distribution that you need for the matvec.  When you
  // define this for your own operator, it is helpful to draw
  // pictures on a sheet of paper to keep track of who needs to
  // receive which entries of the source vector.  Here, each process
  // needs to receive one entry from each of its neighbors.
  //

  // All processes but the first will receive one element from the
  // previous process.
  if (my_rank > 0)
  {
    ++nlocal;
  }
  // All processes but the last will receive one element from the
  // next process.
  if (my_rank < num_procs - 1)
  {
    ++nlocal;
  }
  // Construct a list of columns where this process has nonzero
  // elements.  For our tridiagonal matrix, this is
  // firstRowItOwns-1:lastRowItOwns+1.
  std::vector<global_ordinal_type> indices;
  indices.reserve(nlocal);
  // The first process is a special case...
  if (my_rank > 0)
  {
    indices.push_back(m_opMap->getMinGlobalIndex() - 1);
  }
  for (global_ordinal_type i = m_opMap->getMinGlobalIndex(); i <= m_opMap->getMaxGlobalIndex(); ++i)
  {
    indices.push_back(i);
  }
  // So is the last process...
  if (my_rank < num_procs - 1)
  {
    indices.push_back(m_opMap->getMaxGlobalIndex() + 1);
  }

  // Wrap our vector in an array view, which is like a pointer
  Teuchos::ArrayView<const global_ordinal_type> element_list(indices);

  // Make a column Map for handling the redistribution.
  //
  // There will be some redundancies (i.e., some of the entries will
  // be owned by multiple processes).  Those redundancies will help
  // express the communication pattern for the sparse mat-vec.
  const global_ordinal_type num_global_elements = n + 2 * (num_procs - 1);
  m_redistMap = Teuchos::rcp(new map_type(num_global_elements, element_list, index_base, comm));

  // Make an Import object that describes how data will be
  // redistributed.  It takes a Map describing who owns what
  // originally, and a Map that describes who you WANT to own what.
  m_importer = Teuchos::rcp(new import_type(m_opMap, m_redistMap));
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~TpetraBlockDiagPCImpl()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::initialize()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::compute_impl(
    const std::shared_ptr<const TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A)
{
  std::vector<Scalar> wsp_block;
  math::LapackMatrixInverse<Scalar> inverter;

  for (Uint b = 0; b < m_block_ranges.size(); ++b)
  {
    const SUint block_size = static_cast<SUint>(m_block_ranges[b].size());
    wsp_block.resize(block_size * block_size);
    math::DenseMatView<Scalar> block(wsp_block.data(), block_size, block_size, block_size);

    (*A).get_block(m_block_ranges[b], m_block_ranges[b], block);
    math::DenseMatView<Scalar> block_inverse = m_blocks.mat_view(b);

    inverter.compute(wsp_block.data(), &block_inverse(0, 0), block_size);

  } // Loop over blocks
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setParameters(
    const Teuchos::ParameterList &params)
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<
    const typename TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type>
TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap() const
{
  return m_opMap;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<
    const typename TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type>
TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap() const
{
  return m_opMap;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MV &X, MV &Y,
                                                                             Teuchos::ETransp mode,
                                                                             scalar_type alpha,
                                                                             scalar_type beta) const
{
#if 0
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  Teuchos::RCP<const Teuchos::Comm<int>> comm = m_opMap->getComm();
  const int my_rank = comm->getRank();
  const int num_procs = comm->getSize();

  if (my_rank == 0)
  {
    std::cout << "TpetraBlockDiagPreconditioner::apply" << std::endl;
  }

  // We're writing the Operator subclass, so we are responsible for
  // error handling.  You can decide how much error checking you
  // want to do.  Just remember that checking things like Map
  // sameness or compatibility are expensive.
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::invalid_argument,
                             "X and Y do not have the same numbers of vectors (columns).");

  // Let's make sure alpha is 1 and beta is 0...
  // This will throw an exception if that is not the case.
  TEUCHOS_TEST_FOR_EXCEPTION(alpha != STS::one() || beta != STS::zero(), std::logic_error,
                             "MyOp::apply was given alpha != 1 or beta != 0. "
                             "These cases are not implemented.");

  // Get the number of vectors (columns) in X (and Y).
  const size_t num_vecs = X.getNumVectors();

  // Make a temporary multivector for holding the redistributed
  // data.  You could also create this in the constructor and reuse
  // it across different apply() calls, but you would need to be
  // careful to reallocate if it has a different number of vectors
  // than X.  The number of vectors in X can vary across different
  // apply() calls.
  Teuchos::RCP<MV> redist_data = Teuchos::rcp(new MV(m_redistMap, num_vecs));

  // Redistribute the data.
  // This will do all the necessary communication for you.
  // All processes now own enough data to do the matvec.
  redist_data->doImport(X, *m_importer, Tpetra::INSERT);

  // Get the number of local rows in X, on the calling process.
  const local_ordinal_type n_loc_rows = static_cast<local_ordinal_type>(X.getLocalLength());

  // Perform the matvec with the data we now locally own.
  //
  // For each column...
  for (size_t c = 0; c < num_vecs; ++c)
  {
    // Get a view of the desired column
    Teuchos::ArrayRCP<scalar_type> col_view = redist_data->getDataNonConst(c);

    local_ordinal_type offset;
    // Y[0,c] = -colView[0] + 2*colView[1] - colView[2] (using local indices)
    if (my_rank > 0)
    {
      Y.replaceLocalValue(0, c, -col_view[0] + 2 * col_view[1] - col_view[2]);
      offset = 0;
    }
    // Y[0,c] = 2*colView[1] - colView[2] (using local indices)
    else
    {
      Y.replaceLocalValue(0, c, 2 * col_view[0] - col_view[1]);
      offset = 1;
    }

    // Y[r,c] = -colView[r-offset] + 2*colView[r+1-offset] - colView[r+2-offset]
    for (local_ordinal_type r = 1; r < n_loc_rows - 1; ++r)
    {
      const scalar_type new_val =
          -col_view[r - offset] + 2 * col_view[r + 1 - offset] - col_view[r + 2 - offset];
      Y.replaceLocalValue(r, c, new_val);
    }
    // Y[nlocRows-1,c] = -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset]
    //                   - colView[nlocRows+1-offset]
    if (my_rank < num_procs - 1)
    {
      const scalar_type new_val = -col_view[n_loc_rows - 1 - offset] +
                                  2 * col_view[n_loc_rows - offset] -
                                  col_view[n_loc_rows + 1 - offset];
      Y.replaceLocalValue(n_loc_rows - 1, c, new_val);
    }
    // Y[nlocRows-1,c] = -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset]
    else
    {
      const scalar_type new_val =
          -col_view[n_loc_rows - 1 - offset] + 2 * col_view[n_loc_rows - offset];
      Y.replaceLocalValue(n_loc_rows - 1, c, new_val);
    }
  }
#endif

  // Get the number of vectors (columns) in X (and Y).
  const Uint num_vecs = X.getNumVectors();

  for (Uint b = 0; b < m_blocks.size(); ++b)
  {
    const math::DenseConstMatView<scalar_type> block_op = m_blocks.const_mat_view(b);
    const local_ordinal_type row_offset                 = m_block_ranges[b].lbound();

    for (Uint c = 0; c < num_vecs; ++c)
    {
      Teuchos::ArrayRCP<const scalar_type> col_view_X = X.getData(c);
      Teuchos::ArrayRCP<scalar_type> col_view_Y       = Y.getDataNonConst(c);

      // Compute the action of operator on X
      for (Uint r = 0; r < block_op.rows(); ++r)
      {
        scalar_type op_X = scalar_type();

        for (Uint c = 0; c < block_op.cols(); ++c)
        {
          op_X += block_op(r, c) * col_view_X[row_offset + c];
        }
        /// Compute Y := alpha Op X + beta Y.
        col_view_Y[row_offset + r] = alpha * op_X + beta * col_view_Y[row_offset + r];
      }
    }
  }
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::initialize_structure(
    const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A,
    std::unique_ptr<std::vector<common::Range1D<LocalOrdinal>>> &&blocks)
{
  m_block_ranges.swap(*blocks);

  using mat_shape_2D = common::ArrayShape<_2D, SUint>;

  std::unique_ptr<std::vector<mat_shape_2D>> mat_shapes(new std::vector<mat_shape_2D>());
  mat_shapes->resize(m_block_ranges.size());

  for (Uint b = 0; b < m_block_ranges.size(); ++b)
  {
    const SUint block_size = static_cast<SUint>(m_block_ranges[b].size());
    (*mat_shapes)[b]       = mat_shape_2D(block_size, block_size);
  }

  m_blocks.allocate(std::move(mat_shapes));
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_impl() const
{
  for (Uint b = 0; b < m_blocks.size(); ++b)
  {
    const math::DenseConstMatView<Scalar> block = m_blocks.const_mat_view(b);
    std::cout << block << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_HAVE_TRILINOS

#endif
