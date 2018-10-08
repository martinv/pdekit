#ifndef PDEKIT_Linear_System_Trilinos_hpp
#define PDEKIT_Linear_System_Trilinos_hpp

#include <iostream>

#include "PDEKit_Config.hpp"
#include "linear_system/TrilinosCrsMatrix.hpp"
#include "linear_system/TrilinosMultiVector.hpp"

#if PDEKIT_HAVE_TRILINOS

#include "BelosBlockGmresSolMgr.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosLinearProblem.hpp"
#include "Epetra_LinearProblem.h"
#include "Ifpack.h"
#include "Teuchos_ParameterList.hpp"

namespace pdekit
{

namespace ls
{

class LSTrilinos
{
  public:
  /// Default constructor
  LSTrilinos();

  /// Destructor
  ~LSTrilinos();

  /// Set the matrix of the LSS, vector of RHS and vector of unknowns
  void configure(std::shared_ptr<TrilinosCrsMatrix> &matrix,
                 std::shared_ptr<TrilinosMultiVector> &rhs,
                 std::shared_ptr<TrilinosMultiVector> &x);

  /// Get the system matrix
  std::shared_ptr<TrilinosCrsMatrix> matrix();

  /// Get the system matrix, const version
  const std::shared_ptr<TrilinosCrsMatrix> matrix() const;

  /// Get the rhs vector
  std::shared_ptr<TrilinosMultiVector> rhs();

  /// Get the rhs vector, const version
  const std::shared_ptr<TrilinosMultiVector> rhs() const;

  void solve(const Uint nb_iter = 100, const Real tol = 1.e-9);

  /// Print info
  void print_info() const
  {
#if PDEKIT_HAVE_TRILINOS
    std::cout << "LSS with Trilinos support" << std::endl;
#else
    std::cout << "LSS without Trilinos support" << std::endl;
#endif
  }

  private:
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;

  /// Pointer to system matrix
  std::shared_ptr<TrilinosCrsMatrix> m_matrix;

  /// Pointer to RHS vector
  std::shared_ptr<TrilinosMultiVector> m_rhs;

  /// Pointer to vector of unknowns
  std::shared_ptr<TrilinosMultiVector> m_x;

  /// Linear system
  Teuchos::RCP<Belos::LinearProblem<double, MV, OP>> m_lin_problem;

  /// Belos solver
  Belos::BlockGmresSolMgr<double, MV, OP> m_solver;

  /// Parameter list for the linear solver
  Teuchos::RCP<Teuchos::ParameterList> m_param_list;

  /// Ifpack preconditioner
  Teuchos::RCP<Ifpack_Preconditioner> m_ifpack_prec;

  /// Belos preconditioner
  Teuchos::RCP<Belos::EpetraPrecOp> m_belos_prec;
};

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_HAVE_TRILINOS

#endif
