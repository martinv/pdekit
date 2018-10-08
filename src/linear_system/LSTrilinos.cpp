#include "linear_system/LSTrilinos.hpp"

#if PDEKIT_HAVE_TRILINOS

#include "BelosSolverFactory.hpp"
#include "Ifpack_AdditiveSchwarz.h"

namespace pdekit
{

namespace ls
{

// ============================================================================

LSTrilinos::LSTrilinos()
{
  m_param_list = Teuchos::rcp(new Teuchos::ParameterList());
  m_param_list->set("Block Size",
                    1); // Blocksize to be used by iterative solver
  m_param_list->set("Maximum Iterations",
                    100); // Maximum number of iterations allowed
  m_param_list->set("Convergence Tolerance",
                    1e-9); // Relative convergence tolerance requested
  m_param_list->set("Orthogonalization", "ICGS");
  m_param_list->set("Implicit Residual Scaling", "Norm of Initial Residual");
  int verb = Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary;
  m_param_list->set("Verbosity", verb);
}

// ============================================================================

LSTrilinos::~LSTrilinos()
{
}

// ============================================================================

void LSTrilinos::configure(std::shared_ptr<TrilinosCrsMatrix> &matrix,
                           std::shared_ptr<TrilinosMultiVector> &rhs,
                           std::shared_ptr<TrilinosMultiVector> &x)
{
  // --------------------------------------------
  // Configure operator, left and right-hand side
  // --------------------------------------------

  m_matrix = matrix;
  m_rhs    = rhs;
  m_x      = x;

  // -----------------------------
  // Preconditioner initialization
  // -----------------------------
  Teuchos::ParameterList ifpack_param_list;

  // Allocate an IFPACK factory.  The object contains no data, only
  // the Create() method for creating preconditioners.
  Ifpack factory;

  // Create the preconditioner.  For the list of PrecType values that
  // Create() accepts, please check the IFPACK documentation.
  // Accepted values of preconditioner type:
  // "point relaxation" : returns an instance of
  // Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation> (no Additive Schwarz in
  // serial) "point relaxation stand-alone" : returns an instance of
  // Ifpack_PointRelaxation (value of overlap is ignored). "block relaxation"
  // : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation>
  // (no Additive Schwarz in serial) "block relaxation stand-alone)" : returns
  // an instance of Ifpack_BlockRelaxation. "Amesos" : returns an instance of
  // Ifpack_AdditiveSchwarz<Ifpack_Amesos> (no Additive Schwarz in serial)
  // "Amesos stand-alone" : returns an instance of Ifpack_Amesos.
  // "IC" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_IC> (no
  // Additive Schwarz in serial) "IC stand-alone" : returns an instance of
  // Ifpack_IC. "ICT" : returns an instance of
  // Ifpack_AdditiveSchwarz<Ifpack_ICT> (no Additive Schwarz in serial) "ICT
  // stand-alone" : returns an instance of Ifpack_ICT. "ILU" : returns an
  // instance of Ifpack_AdditiveSchwarz<Ifpack_ILU> (no Additive Schwarz in
  // serial)
  // "ILU stand-alone" : returns an instance of Ifpack_ILU.
  // "ILUT" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_ILUT> (no
  // Additive Schwarz in serial) "ILUT stand-alone" : returns an instance of
  // Ifpack_ILUT.

  const std::string prec_type = "ILU"; // incomplete LU
  int overlap_level           = 1;     // must be >= 0. If Comm.NumProc() == 1,
                                       // it is ignored.

  m_ifpack_prec =
      Teuchos::rcp(factory.Create(prec_type, &*(m_matrix->m_epetra_matrix), overlap_level));
  TEUCHOS_TEST_FOR_EXCEPTION(m_ifpack_prec == Teuchos::null, std::runtime_error,
                             "IFPACK failed to create a preconditioner of type \""
                                 << prec_type << "\" with overlap level " << overlap_level << ".");

  // Specify parameters for ILU.  ILU is local to each MPI process.
  ifpack_param_list.set("fact: drop tolerance", 1e-9);
  ifpack_param_list.set("fact: level-of-fill", 1);

  // IFPACK uses overlapping Schwarz domain decomposition over all
  // participating processes to combine the results of ILU on each
  // process.  IFPACK's Schwarz method can use any of the following
  // combine modes to combine overlapping results:
  //
  // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
  //
  // The Epetra_CombineMode.h header file defines their meaning.
  ifpack_param_list.set("schwarz: combine mode", "Add");
  // Set the parameters.
  // IFPACK_CHK_ERR(m_ifpack_prec->SetParameters(ifpack_param_list));
  m_ifpack_prec->SetParameters(ifpack_param_list);

  // Initialize the preconditioner. At this point the matrix must have
  // been FillComplete()'d, but actual values are ignored.
  // IFPACK_CHK_ERR(m_ifpack_prec->Initialize());
  m_ifpack_prec->Initialize();

  // Build the preconditioner, by looking at the values of the matrix.
  // IFPACK_CHK_ERR(ifpack_prec->Compute());
  m_ifpack_prec->Compute();

  // Create the Belos preconditioned operator from the Ifpack preconditioner.
  // NOTE:  This is necessary because Belos expects an operator to apply the
  //        preconditioner with Apply() NOT ApplyInverse().
  m_belos_prec = Teuchos::rcp(new Belos::EpetraPrecOp(m_ifpack_prec));

  // -----------------------------
  // Set the linear problem
  // -----------------------------

  // m_lin_problem = Teuchos::rcp(new Belos::LinearProblem<double, MV, OP>(
  //    m_matrix->m_epetra_matrix, m_x->m_vector, m_rhs->m_vector));
  m_lin_problem = Teuchos::rcp(new Belos::LinearProblem<double, MV, OP>());

  m_lin_problem->setOperator(m_matrix->m_epetra_matrix);
  bool set = m_lin_problem->setProblem(m_x->m_vector, m_rhs->m_vector);

  TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
                             "*** Belos::LinearProblem failed to set up correctly! ***");

  // // m_lin_problem->setLHS(m_x->m_vector);
  // // m_lin_problem->setRHS(m_rhs->m_vector);

  // Set the IFPACK preconditioner.
  //
  // We're using it as a right preconditioner.  It's better to use a
  // right preconditioner than a left preconditioner in GMRES, because
  // otherwise the projected problem will have a different residual
  // (in exact arithmetic) than the original problem.  This makes it
  // harder for GMRES to tell when it has converged.
  m_lin_problem->setRightPrec(m_belos_prec);
  // bool set = m_lin_problem->setProblem();
  // TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
  //                            "*** Belos::LinearProblem failed to set up
  //                            correctly! ***");

  // Belos::SolverFactory<double, MV, OP> belos_factory;
  // m_solver = belos_factory.create("GMRES", m_param_list);
  m_solver.setProblem(m_lin_problem);
  m_solver.setParameters(m_param_list);
  // std::cout << "Is problem set = " << m_lin_problem->isProblemSet() <<
  // std::endl;
}

// ============================================================================

std::shared_ptr<TrilinosCrsMatrix> LSTrilinos::matrix()
{
  return m_matrix;
}

// ============================================================================

const std::shared_ptr<TrilinosCrsMatrix> LSTrilinos::matrix() const
{
  return m_matrix;
}

// ============================================================================

std::shared_ptr<TrilinosMultiVector> LSTrilinos::rhs()
{
  return m_rhs;
}

// ============================================================================

const std::shared_ptr<TrilinosMultiVector> LSTrilinos::rhs() const
{
  return m_rhs;
}

// ============================================================================

void LSTrilinos::solve(const Uint nb_iter, const Real tol)
{
  m_param_list->set("Maximum Iterations", static_cast<int>(nb_iter));
  m_param_list->set("Convergence Tolerance", static_cast<double>(tol));

  m_solver.setParameters(m_param_list);

  // Perform solve.
  Belos::ReturnType ret = m_solver.solve();

  // Did we converge?
  // if (myRank == 0) {
  if (ret == Belos::Converged)
  {
    std::cout << "Belos converged." << std::endl;
  }
  else
  {
    std::cout << "Belos did not converge." << std::endl;
  }
  //}
}

// ============================================================================

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_HAVE_TRILINOS
