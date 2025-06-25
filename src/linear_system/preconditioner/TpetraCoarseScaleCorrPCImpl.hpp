#ifndef PDEKIT_Linear_System_Tpetra_Coarse_Scale_Corr_PC_Impl_hpp
#define PDEKIT_Linear_System_Tpetra_Coarse_Scale_Corr_PC_Impl_hpp

#include <iostream>

#include "linear_system/CoarseScaleCorrection.hpp"
#include "linear_system/TpetraDirectSolver.hpp"
#include "linear_system/preconditioner/TrilinosPC.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "Amesos2.hpp"
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

template <typename Scalar        = TpetraDefaultTraits::Scalar,
          typename LocalOrdinal  = TpetraDefaultTraits::LocalOrdinal,
          typename GlobalOrdinal = TpetraDefaultTraits::GlobalOrdinal,
          typename Node          = TpetraDefaultTraits::Node>
class TpetraCoarseScaleCorrPCImpl
    : public Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>
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
  using matrix_type =
      Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using import_type = Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type>;

  public:
  /// Constructor
  /// n: Global number of rows and columns in the operator.
  /// comm: The communicator over which to distribute those rows and columns.
  TpetraCoarseScaleCorrPCImpl(const global_ordinal_type n,
                              const Teuchos::RCP<const Teuchos::Comm<GlobalOrdinal>> comm);

  /// Destructor
  virtual ~TpetraCoarseScaleCorrPCImpl();

  void initialize();

  /// Build the preconditioner data
  void initialize_structure(
      const std::shared_ptr<TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A,
      const common::BlockArray<Uint, Uint> &mesh_dual_graph_crs,
      const std::vector<Uint> &coarse_block_sizes, const std::vector<Uint> &fine_block_sizes,
      const std::vector<Uint> &cell_reordering,
      std::unique_ptr<ls::LocalTransferOps<Real>> &&restriction_ops,
      std::unique_ptr<ls::LocalTransferOps<Real>> &&prolongation_ops);

  /// Compute the preconditioner
  void compute_impl(
      const std::shared_ptr<TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A);

  /// Set preconditioner parameter
  void setParameters(const Teuchos::ParameterList &params);

  // These functions are required since we inherit from Tpetra::Operator

  /// Get the domain Map of this Operator subclass.
  Teuchos::RCP<const map_type> getDomainMap() const;

  /// Get the range Map of this Operator subclass.
  Teuchos::RCP<const map_type> getRangeMap() const;

  /// Compute Y := alpha Op X + beta Y.
  ///
  /// We ignore the cases alpha != 1 and beta != 0 for simplicity.
  void apply(const MV &X, MV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
             scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
             scalar_type beta  = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  /// Print the contents of the preconditioner
  void print_impl() const;

  private:
  using dof_map_wrapper_type = TpetraDofMap<local_ordinal_type, global_ordinal_type, node_type>;
  using crs_mat_wrapper_type =
      TpetraCrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using multivec_wrapper_type =
      TpetraMultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using coarse_scale_solver_type =
      internal::TpetraDirectSolver<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using fine_smoother_type =
      IfpackPC<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;

  Teuchos::RCP<const map_type> m_opMap, m_redistMap;
  Teuchos::RCP<const import_type> m_importer;

  /// Coarse scale correction implementation
  CoarseScaleCorrection m_correction;

  Teuchos::RCP<const Teuchos::Comm<GlobalOrdinal>> m_comm;

  /*
  Teuchos::RCP<const map_type> m_map_coarse;

  /// Coarse-level matrix and vector
  Teuchos::RCP<matrix_type> m_A_coarse;
  Teuchos::RCP<MV> m_b_coarse;
  Teuchos::RCP<MV> m_u_coarse;

  /// Coarse-scale solver
  Teuchos::RCP<Amesos2::Solver<matrix_type, MV>> m_solver_coarse;
  */

  std::shared_ptr<dof_map_wrapper_type> m_map_coarse;

  /// Coarse-level matrix and vector
  std::shared_ptr<crs_mat_wrapper_type> m_A_coarse;
  std::shared_ptr<multivec_wrapper_type> m_b_coarse;
  std::shared_ptr<multivec_wrapper_type> m_u_coarse;

  /// Fine level-residual
  std::shared_ptr<multivec_wrapper_type> m_res_fine;
  /// Fine-level solution
  std::shared_ptr<multivec_wrapper_type> m_u_fine;

  /// Coarse-scale solver
  coarse_scale_solver_type m_solver_coarse;

  /// Fine-level smoother
  std::shared_ptr<crs_mat_wrapper_type> m_A_fine;
  std::shared_ptr<fine_smoother_type> m_smoother_fine;
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraCoarseScaleCorrPCImpl(
    const global_ordinal_type n, const Teuchos::RCP<const Teuchos::Comm<GlobalOrdinal>> comm)
    : m_comm(comm)
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
    std::cout << "TpetraCoarseScaleCorrPreconditionerImpl constructor" <<
  std::endl;
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
TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal,
                            Node>::~TpetraCoarseScaleCorrPCImpl()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::initialize()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::initialize_structure(
    const std::shared_ptr<TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A,
    const common::BlockArray<Uint, Uint> &mesh_dual_graph_crs,
    const std::vector<Uint> &coarse_block_sizes, const std::vector<Uint> &fine_block_sizes,
    const std::vector<Uint> &cell_reordering,
    std::unique_ptr<ls::LocalTransferOps<Real>> &&restriction_ops,
    std::unique_ptr<ls::LocalTransferOps<Real>> &&prolongation_ops)
{
  m_correction.build_block_matrix_sparsity_pattern(
      mesh_dual_graph_crs, coarse_block_sizes, fine_block_sizes, cell_reordering,
      std::move(restriction_ops), std::move(prolongation_ops));

  const math::BlockMatrixSparsityPattern<Int> &coarse_sparsity =
      m_correction.coarse_block_sparsity();

  const common::ArrayShape<_2D, Uint> mat_shape = coarse_sparsity.size();

  /*
  // Create distributed map and matrix of of coarse-scale system
  m_map_coarse = Teuchos::rcp(new map_type(mat_shape.size(0), 0, m_comm));

  // Allocate coarse system matrix and generate coarse matrix sparsity pattern
  TpetraFacadeTools::create_matrix_sparsity(m_map_coarse, coarse_sparsity,
  m_A_coarse);

  // Allocate coarse vectors for left-hand side and solution
  const size_t num_vectors = 1;
  m_b_coarse = Teuchos::rcp(new MV(m_map_coarse, num_vectors));
  m_u_coarse = Teuchos::rcp(new MV(m_map_coarse, num_vectors));

  // Prepare system
  const std::string solver_name = "KLU2";
  // const std::string solver_name = "SuperLU";
  // const std::string solver_name = "Basker";
  m_solver_coarse.reset();
  m_solver_coarse =
      Amesos2::create<matrix_type, MV>(solver_name, m_A_coarse, m_u_coarse,
  m_b_coarse);

  std::cout << "Number of rows in coarse matrix = " <<
  m_A_coarse->getGlobalNumRows() << std::endl;
  */

  // ------------------------------------------------------
  // Create bare Tpetra objects, fill them with initial
  // data and then pass them to constructors
  // of wrapper classes
  // ------------------------------------------------------

  // Create distributed map and matrix of of coarse-scale system
  Teuchos::RCP<map_type const> raw_map_coarse =
      Teuchos::rcp(new map_type(mat_shape.size(0), 0, m_comm));
  m_map_coarse = std::make_shared<dof_map_wrapper_type>(raw_map_coarse);

  // Allocate coarse system matrix and generate coarse matrix sparsity pattern
  // First create just a Teuchos RCP pointer to Tpetra matrix
  Teuchos::RCP<matrix_type> raw_A_coarse;
  TpetraFacadeTools::create_matrix_sparsity(raw_map_coarse, coarse_sparsity, raw_A_coarse);

  // Allocate coarse system matrix and generate coarse matrix sparsity pattern
  m_A_coarse = std::make_shared<crs_mat_wrapper_type>(m_comm, raw_map_coarse, raw_A_coarse);

  // ------------------------
  // Allocate vectors
  // ------------------------
  // Allocate coarse vectors for left-hand side and solution
  const size_t num_vectors = 1;
  m_b_coarse = std::make_shared<multivec_wrapper_type>(m_A_coarse->map(), num_vectors);
  m_u_coarse = std::make_shared<multivec_wrapper_type>(m_A_coarse->map(), num_vectors);

  m_res_fine = std::make_shared<multivec_wrapper_type>(A->map(), num_vectors);
  m_u_fine   = std::make_shared<multivec_wrapper_type>(A->map(), num_vectors);

  // Prepare system
  const bool verbose = false;
  m_solver_coarse.initialize_solver(m_A_coarse, m_u_coarse, m_b_coarse, verbose);
  m_solver_coarse.update_after_mat_graph_change(m_A_coarse, m_u_coarse, m_b_coarse, verbose);

  // Prepare fine-level smoother
  m_smoother_fine = std::make_shared<fine_smoother_type>();
  m_smoother_fine->create("ILUT", A);
  m_smoother_fine->initialize(A);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::compute_impl(
    const std::shared_ptr<TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A)
{
  const bool verbose = false;
  m_correction.fill_coarse_level_matrix(*A, *m_A_coarse);
  m_solver_coarse.update_after_mat_values_change(m_A_coarse, m_u_coarse, m_b_coarse, verbose);

  m_smoother_fine->compute(A);
  m_A_fine = A;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setParameters(
    const Teuchos::ParameterList &params)
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<
    const typename TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type>
TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap() const
{
  return m_opMap;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<
    const typename TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type>
TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap() const
{
  return m_opMap;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(
    const MV &X, MV &Y, Teuchos::ETransp mode, scalar_type alpha, scalar_type beta) const
{
  using iaccess = internal::TpetraInternalAccess;

  const bool verbose = false;

  // Restrict RHS vector
  m_correction.restrict_vec(X, *iaccess::get_vector_data(*m_b_coarse));

  /*
  Teuchos::ArrayRCP<const Scalar> X_fine_view = X.getData(0);
  Teuchos::ArrayRCP<const Scalar> b_coarse_view =
  iaccess::get_vector_data(*m_b_coarse)->getData(0);

  std::cout << "Size of fine view = " << X_fine_view.size() << std::endl;
  std::cout << "Size of coarse view = " << b_coarse_view.size() << std::endl;

  std::cout << "\nVector before restriction: " << std::endl;
  for (Uint i = 0; i < X_fine_view.size(); ++i)
  {
    std::cout << X_fine_view[i] << std::endl;
  }

  std::cout << "\nVector after restriction: " << std::endl;
  for (Uint i = 0; i < b_coarse_view.size(); ++i)
  {
    std::cout << b_coarse_view[i] << std::endl;
  }
  */

  // Solve coarse system
  // m_solver_coarse.solve(100, 1.e-9, verbose);

  // The 'solve' method of coarse solver is not 'const', but the 'apply'
  // method here is It's therefore not possible to call m_solver_coarse.solve(
  // ... ) here and we have to use a workaround
  coarse_scale_solver_type &coarse_solver = const_cast<coarse_scale_solver_type &>(m_solver_coarse);
  coarse_solver.solve(100, 1.e-9, verbose);

  /*
  Teuchos::ArrayRCP<const Scalar> x_coarse_view =
  iaccess::get_vector_data(*m_u_coarse)->getData(0); std::cout << "\nCoarse
  system solution: " << std::endl; for (Uint i = 0; i < x_coarse_view.size();
  ++i)
  {
    std::cout << x_coarse_view[i] << std::endl;
  }
  */

  // Prolongate the solution
  m_u_fine->fill(0.0);
  m_correction.prolongate_vec(*iaccess::get_vector_data(*m_u_coarse),
                              *iaccess::get_vector_data(*m_u_fine));

  /*
  Teuchos::ArrayRCP<const Scalar> x_fine_view =
  iaccess::get_vector_data(*m_u_fine)->getData(0); std::cout << "\nFine system
  solution: " << std::endl; for (Uint i = 0; i < x_fine_view.size(); ++i)
  {
    std::cout << x_fine_view[i] << std::endl;
  }
  */

  // iaccess::get_vector_data(*m_res_fine)->assign(X);
  iaccess::get_vector_data(*m_res_fine)->assign(*iaccess::get_vector_data(*m_u_fine));

  // Compute fine-level residual
  // void apply(X, Y, transp, alpha, beta)
  // apply: Y := beta*Y + alpha*Op(A)*X with Op(A) = A or Op(A) = transpose(A)
  iaccess::get_matrix_data(*m_A_fine)->apply(*iaccess::get_vector_data(*m_u_fine),
                                             *iaccess::get_vector_data(*m_res_fine),
                                             Teuchos::NO_TRANS, -1.0, 1.0);

  /*
  Teuchos::ArrayRCP<const Scalar> res_fine_view =
  iaccess::get_vector_data(*m_res_fine)->getData(0); std::cout << "\nFine
  level residual: " << std::endl; for (Uint i = 0; i < res_fine_view.size();
  ++i)
  {
    std::cout << res_fine_view[i] << std::endl;
  }
  */

  // Update Y with smoothed fine-level residual
  // void apply(X, Y, transp, alpha, beta)
  // apply: Y := beta*Y + alpha*Op(A)*X with Op(A) = A or Op(A) = transpose(A)

  Y.assign(*iaccess::get_vector_data(*m_u_fine));

  /*
  iaccess::get_preconditioner(*m_smoother_fine)
      ->apply(*iaccess::get_vector_data(*m_res_fine), Y, Teuchos::NO_TRANS,
  0.5, 1.0);
  */

  // iaccess::get_preconditioner(*m_smoother_fine)->apply(X, Y,
  // Teuchos::NO_TRANS, 1.0, 0.0);

  // Y.assign(X);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraCoarseScaleCorrPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_impl() const
{
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_HAVE_TRILINOS

#endif
