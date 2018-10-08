#ifndef PDEKIT_Linear_System_Tpetra_Iterative_Solver_hpp
#define PDEKIT_Linear_System_Tpetra_Iterative_Solver_hpp

#include <iostream>

#include "PDEKit_Config.hpp"
#include "common/PDEKit.hpp"
#include "linear_system/TpetraLinSolver.hpp"
#include "linear_system/preconditioner/AMGPC.hpp"
#include "linear_system/preconditioner/IfpackPC.hpp"
#include "linear_system/preconditioner/TpetraBlockDiagPC.hpp"
#include "linear_system/preconditioner/TpetraCoarseScaleCorrPC.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Ifpack2_Factory.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Tpetra_CrsMatrix.hpp"
#endif

namespace pdekit
{

namespace ls
{

namespace internal
{

// ----------------------------------------------------------------------------

template <typename Scalar        = TpetraDefaultTraits::Scalar,
          typename LocalOrdinal  = TpetraDefaultTraits::LocalOrdinal,
          typename GlobalOrdinal = TpetraDefaultTraits::GlobalOrdinal,
          typename Node          = TpetraDefaultTraits::Node,
          const bool classic     = TpetraDefaultTraits::classic>
class TpetraIterativeSolver
    : public TpetraLinSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>
{
  public:
  using base_type           = TpetraLinSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
  using matrix_type         = typename base_type::matrix_type;
  using vector_type         = typename base_type::vector_type;
  using preconditioner_type = typename base_type::preconditioner_type;

  /// Constructor
  TpetraIterativeSolver();

  /// Destructor
  ~TpetraIterativeSolver() override;

#if 0
  /// Configure the solver: set lhs and rhs
  void configure(const Teuchos::RCP<const matrix_type> &mat_lhs,
                 const Teuchos::RCP<vector_type> &vec_x,
                 const Teuchos::RCP<const vector_type> &vec_rhs,
                 const bool recompute_preconditioner, const bool print_output = false);
#endif

  /// Build the solver and its preconditioner
  void initialize_solver(const std::shared_ptr<matrix_type> &mat_lhs,
                         const std::shared_ptr<vector_type> &vec_x,
                         const std::shared_ptr<vector_type> &vec_rhs,
                         const bool print_output = false) override;

  /// Update the preconditioner and other data when the matrix graph changes
  /// (we added or removed elements from the matrix, thus changing the number
  /// of its nonzeros
  void update_after_mat_graph_change(const std::shared_ptr<matrix_type> &mat_lhs,
                                     const std::shared_ptr<vector_type> &vec_x,
                                     const std::shared_ptr<vector_type> &vec_rhs,
                                     const bool print_output = false) override;

  /// Update the preconditioner and other data when the matrix or rhs change
  /// (i.e the matrix graph is constant, but values in the matrix or rhs (or
  /// both) were updated
  void update_after_mat_values_change(const std::shared_ptr<matrix_type> &mat_lhs,
                                      const std::shared_ptr<vector_type> &vec_x,
                                      const std::shared_ptr<vector_type> &vec_rhs,
                                      const bool print_output = false) override;

  /// Run the actual solver
  void solve(const Uint nb_iter, const Real tol, const bool verbose) override;

  /// Connect preconditioner
  void connect_preconditioner(std::shared_ptr<preconditioner_type> &preconditioner) override;

  /// Release preconditioner.
  void release_preconditioner() override;

  private:
#if PDEKIT_HAVE_TRILINOS
  using tpetra_op_type = Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Belos_lin_problem_type =
      typename Belos::LinearProblem<Scalar, typename vector_type::trilinos_multivector_type,
                                    tpetra_op_type>;

  /// METHODS

  void set_parameters_gmres();
  void set_parameters_cg();

  /// --------------------------------
  /// DATA FOR SPARSE ITERATIVE SOLVER
  /// --------------------------------

  /// Linear system
  Teuchos::RCP<Belos_lin_problem_type> m_lin_problem;

  /// Iterative solver - Belos package
  Teuchos::RCP<
      Belos::SolverManager<Scalar, typename vector_type::trilinos_multivector_type, tpetra_op_type>>
      m_iterative_solver;

  /// Parameter list for the iterative linear solver
  Teuchos::RCP<Teuchos::ParameterList> m_iterative_param_list;

  /// Trilinos preconditioner
  std::shared_ptr<TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> m_tpetra_prec;

  /// Type of preconditioner
  LeftRightPrec m_prec_type;
#endif
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraIterativeSolver()
{
#if PDEKIT_HAVE_TRILINOS
  set_parameters_gmres();
  // set_parameters_cg();
  m_prec_type = LeftRightPrec::ePrecRight;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::~TpetraIterativeSolver()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::initialize_solver(
    const std::shared_ptr<matrix_type> &mat_lhs, const std::shared_ptr<vector_type> &vec_x,
    const std::shared_ptr<vector_type> &vec_rhs, const bool print_output)
{
#if PDEKIT_HAVE_TRILINOS
  set_parameters_gmres();
  // set_parameters_cg();

  using iaccess = internal::TpetraInternalAccess;

  Belos::SolverFactory<Scalar, typename vector_type::trilinos_multivector_type, tpetra_op_type>
      belos_factory;
  m_iterative_solver.reset();
  // "GMRES", "FLEXIBLE GMRES", "BLOCK GMRES"
  m_iterative_solver = belos_factory.create("GMRES", m_iterative_param_list);
  // m_iterative_solver = belos_factory.create("CG", m_iterative_param_list);

  m_lin_problem.reset();
  m_lin_problem = Teuchos::RCP<Belos_lin_problem_type>(new Belos_lin_problem_type(
      iaccess::get_matrix_data(*mat_lhs), iaccess::get_vector_data(*vec_x),
      iaccess::get_vector_data(*vec_rhs)));

  // m_tpetra_prec =
  //     std::make_shared<IfpackPreconditioner<Scalar, LocalOrdinal,
  //     GlobalOrdinal, Node, classic>>();

  /*
  m_tpetra_prec = std::make_shared<AMGPreconditioner<Scalar, LocalOrdinal,
  GlobalOrdinal, Node, classic>>();
  */

  // m_tpetra_prec->create("ILUT", mat_lhs);

  // m_lin_problem->setRightPrec(iaccess::get_preconditioner(*m_tpetra_prec));
  // m_lin_problem->setLeftPrec(iaccess::get_preconditioner(*m_tpetra_prec));

  // Tell the solver what problem you want to solve.
  m_iterative_solver->setProblem(m_lin_problem);
#endif
}

// ----------------------------------------------------------------------------
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
    update_after_mat_graph_change(const std::shared_ptr<matrix_type> &mat_lhs,
                                  const std::shared_ptr<vector_type> &vec_x,
                                  const std::shared_ptr<vector_type> &vec_rhs,
                                  const bool print_output)
{
#if PDEKIT_HAVE_TRILINOS
  using iaccess = internal::TpetraInternalAccess;

  if (m_tpetra_prec)
  {
    Teuchos::oblackholestream blackHole;
    std::ostream &out = (print_output) ? std::cout : blackHole;

    Teuchos::RCP<Teuchos::Time> initialize_timer =
        Teuchos::TimeMonitor::getNewCounter("Preconditioner::initialize");

    // Initialize the preconditioner.  If the sparse matrix's _structure_
    // changes, you have to call initialize() and compute() again, in that
    // sequence, before you may use ("apply") the preconditioner.
    out << "-- Initializing preconditioner" << std::endl;
    {
      Teuchos::TimeMonitor mon(*initialize_timer);
      m_tpetra_prec->initialize(mat_lhs);

      // NEEDED TO RESET THE POINTER TO PRECONDITIONER!
      if (m_prec_type == LeftRightPrec::ePrecLeft)
      {
        m_lin_problem->setLeftPrec(iaccess::get_preconditioner(*m_tpetra_prec));
      }
      else if (m_prec_type == LeftRightPrec::ePrecRight)
      {
        m_lin_problem->setRightPrec(iaccess::get_preconditioner(*m_tpetra_prec));
      }
    }
  }
#if 0
  else
  {
    std::cerr << "TpetraIterativeSolver::ERRROR 1: NO PRECONDITIONER!" << std::endl;
  }
#endif
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
    update_after_mat_values_change(const std::shared_ptr<matrix_type> &mat_lhs,
                                   const std::shared_ptr<vector_type> &vec_x,
                                   const std::shared_ptr<vector_type> &vec_rhs,
                                   const bool print_output)
{
#if PDEKIT_HAVE_TRILINOS
  using iaccess = internal::TpetraInternalAccess;

  if (m_tpetra_prec)
  {
    Teuchos::oblackholestream blackHole;
    std::ostream &out = (print_output) ? std::cout : blackHole;

    Teuchos::RCP<Teuchos::Time> compute_timer =
        Teuchos::TimeMonitor::getNewCounter("Preconditioner::compute");

    // Initialize the preconditioner.  If the sparse matrix's _structure_
    // changes, you have to call initialize() and compute() again, in that
    // sequence, before you may use ("apply") the preconditioner.
    out << "-- Computing preconditioner" << std::endl;
    {
      Teuchos::TimeMonitor mon(*compute_timer);
      m_tpetra_prec->compute(mat_lhs);

      // NEEDED TO RESET THE POINTER TO PRECONDITIONER!

      if (m_prec_type == LeftRightPrec::ePrecLeft)
      {
        m_lin_problem->setLeftPrec(iaccess::get_preconditioner(*m_tpetra_prec));
      }
      else if (m_prec_type == LeftRightPrec::ePrecRight)
      {
        m_lin_problem->setRightPrec(iaccess::get_preconditioner(*m_tpetra_prec));
      }
    }
  }
#if 0
  else
  {
    std::cerr << "TpetraIterativeSolver::ERRROR 2: NO PRECONDITIONER!" << std::endl;
  }
#endif
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::solve(
    const Uint nb_iter, const Real tol, const bool verbose)
{
#if PDEKIT_HAVE_TRILINOS
  m_iterative_param_list->set("Maximum Iterations", static_cast<int>(nb_iter));
  m_iterative_param_list->set("Convergence Tolerance", static_cast<double>(tol));
  // m_iterative_param_list->set("Implicit Residual Scaling", "Norm of Initial
  // Residual");

  m_iterative_solver->setParameters(m_iterative_param_list);

  // Tell the linear solver to prepare the problem to solve.
  // This computes, for example, the initial residual vector(s).
  m_iterative_solver->reset(Belos::Problem);

  // Perform solve.
  Belos::ReturnType ret = m_iterative_solver->solve();

  // Did we converge?
  // if (myRank == 0) {

  // Teuchos::oblackholestream blackHole;
  // std::ostream &out = (print_output) ? std::cout : blackHole;

  if (verbose)
  {
    if (ret == Belos::Converged)
    {
      // out << "Belos converged." << std::endl;
      std::cout << "Belos converged." << std::endl;
    }
    else
    {
      // out << "Belos did not converge." << std::endl;
      std::cout << "Belos did not converge." << std::endl;
    }
  }
//}
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
    connect_preconditioner(std::shared_ptr<preconditioner_type> &preconditioner)
{
#if PDEKIT_HAVE_TRILINOS
  using iaccess = internal::TpetraInternalAccess;

  m_tpetra_prec = preconditioner;

  if (m_prec_type == LeftRightPrec::ePrecLeft)
  {
    m_lin_problem->setLeftPrec(iaccess::get_preconditioner(*m_tpetra_prec));
  }
  else if (m_prec_type == LeftRightPrec::ePrecRight)
  {
    m_lin_problem->setRightPrec(iaccess::get_preconditioner(*m_tpetra_prec));
  }
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                           classic>::release_preconditioner()
{
#if PDEKIT_HAVE_TRILINOS
  m_lin_problem->setRightPrec(Teuchos::RCP<tpetra_op_type>(Teuchos::null));
  m_lin_problem->setLeftPrec(Teuchos::RCP<tpetra_op_type>(Teuchos::null));
  m_tpetra_prec.reset();
#endif
}

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                           classic>::set_parameters_gmres()
{
  m_iterative_param_list = Teuchos::rcp(new Teuchos::ParameterList("IterativeSolverGMRESParams"));
  // Blocksize to be used by iterative solver. Trilinos default is 300
  m_iterative_param_list->set("Num Blocks", 300);
  // Maximum number of iterations allowed. Trilinos default is 1000
  m_iterative_param_list->set("Maximum Iterations", 200);
  // Maximum number of restarts. Trilinos default is 20
  m_iterative_param_list->set("Maximum Restarts", 20);
  // Relative convergence tolerance requested
  m_iterative_param_list->set("Convergence Tolerance", 1e-9);

  m_iterative_param_list->set("Orthogonalization", "ICGS");
  m_iterative_param_list->set("Implicit Residual Scaling", "Norm of Initial Residual");

  // int verb = Belos::Errors + Belos::Warnings + Belos::TimingDetails +
  // Belos::FinalSummary; int verb = Belos::Errors + Belos::Warnings +
  // Belos::FinalSummary;
  int verb = Belos::Errors + Belos::Warnings;
  m_iterative_param_list->set("Verbosity", verb);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::set_parameters_cg()
{
  // Permissible parameters and default values:
  /*
    "Convergence Tolerance" : double = 1e-08
    "Maximum Iterations" : int = 1000
    "Assert Positive Definiteness" : bool = 1
    "Verbosity" : int = 0
    "Output Style" : int = 0
    "Output Frequency" : int = -1
    "Deflation Quorum" : int = 1
    "Output Stream" : Teuchos::RCP<std::ostream> =
    Teuchos::RCP<std::ostream>{ptr=0x7fd7f2c6f460,node=0x557a74b8b420,strong_count=4,weak_count=0}
    "Show Maximum Residual Norm Only" : bool = 0
    "Implicit Residual Scaling" : string = Norm of Initial Residual
    "Estimate Condition Number" : bool = 0
    "Residual Scaling" : string = Norm of Initial Residual
    "Timer Label" : string = Belos
  */

  m_iterative_param_list = Teuchos::rcp(new Teuchos::ParameterList("IterativeSolverCGParams"));

  m_iterative_param_list->set("Convergence Tolerance", 1e-9);
  m_iterative_param_list->set("Maximum Iterations", 200);
  m_iterative_param_list->set("Implicit Residual Scaling", "Norm of Initial Residual");
  m_iterative_param_list->set("Residual Scaling", "Norm of Initial Residual");

  // int verb = Belos::Errors + Belos::Warnings + Belos::TimingDetails +
  // Belos::FinalSummary; int verb = Belos::Errors + Belos::Warnings +
  // Belos::FinalSummary;
  int verb = Belos::Errors + Belos::Warnings;
  m_iterative_param_list->set("Verbosity", verb);
}
#endif

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace ls

} // namespace pdekit

#endif
