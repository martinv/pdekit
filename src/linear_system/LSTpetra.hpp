#ifndef PDEKIT_Linear_System_Tpetra_hpp
#define PDEKIT_Linear_System_Tpetra_hpp

#include <iostream>

#include "PDEKit_Config.hpp"
#include "linear_system/TpetraDirectSolver.hpp"
#include "linear_system/TpetraInternalAccess.hpp"
#include "linear_system/TpetraIterativeSolver.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Ifpack2_Factory.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_oblackholestream.hpp"

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
#endif

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
class LSTpetra
{
  public:
  using matrix_type         = TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
  using multivector_type    = TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
  using preconditioner_type = TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;

  /// Default constructor
  LSTpetra(const SparseSolverType solver_type = SparseSolverType::eIterative);

  /// Destructor
  ~LSTpetra();

#if 0
  /// Set the matrix of the LSS, vector of RHS and vector of unknowns
  void configure(std::shared_ptr<matrix_type> const &matrix,
                 std::shared_ptr<multivector_type> const &rhs,
                 std::shared_ptr<multivector_type> const &x, const bool recompute_preconditioner,
                 const bool print_output);
#endif

  /// Build the solver and its preconditioner
  void initialize_solver(std::shared_ptr<matrix_type> const &matrix,
                         std::shared_ptr<multivector_type> const &rhs,
                         std::shared_ptr<multivector_type> const &x, const bool print_output);

  /// Connect preconditioner
  void connect_preconditioner(std::shared_ptr<preconditioner_type> &preconditioner);

  /// Disconnect preconditioner
  void disconnect_preconditioner();

  /// Update the preconditioner and other data when the matrix graph changes
  /// (we added or removed elements from the matrix, thus changing the number
  /// of its nonzeros
  void update_after_mat_graph_change(std::shared_ptr<matrix_type> const &matrix,
                                     std::shared_ptr<multivector_type> const &rhs,
                                     std::shared_ptr<multivector_type> const &x,
                                     const bool print_output);

  /// Update the preconditioner and other data when the matrix or rhs change
  /// (i.e the matrix graph is constant, but values in the matrix or rhs (or
  /// both) were updated
  void update_after_mat_values_change(std::shared_ptr<matrix_type> const &matrix,
                                      std::shared_ptr<multivector_type> const &rhs,
                                      std::shared_ptr<multivector_type> const &x,
                                      const bool print_output);

  /// Get the system matrix
  std::shared_ptr<matrix_type> matrix();

  /// Get the system matrix, const version
  const std::shared_ptr<matrix_type> matrix() const;

  /// Get the rhs vector
  std::shared_ptr<multivector_type> rhs();

  /// Get the rhs vector, const version
  const std::shared_ptr<multivector_type> rhs() const;

  /// Get the rhs vector
  std::shared_ptr<multivector_type> solution_vector();

  /// Get the rhs vector, const version
  const std::shared_ptr<multivector_type> solution_vector() const;

  /// Get preconditioner
  std::shared_ptr<preconditioner_type> preconditioner();

  /// Get preconditioner, const version
  const std::shared_ptr<preconditioner_type> preconditioner() const;

  /// Set preconditioner
  void set_preconditioner(std::shared_ptr<preconditioner_type> &preconditioner);

  /// Evaluate the system residual, i.e. compute the expression res = b - A *
  /// x
  void residual(multivector_type &res) const;

  void solve(const Uint nb_iter = 100, const Real tol = 1.e-9, const bool verbose = false);

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
  using trilinos_lin_solver_type =
      typename internal::TpetraLinSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
  using trilinos_iterative_solver_type =
      typename internal::TpetraIterativeSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
  using trilinos_direct_solver_type =
      typename internal::TpetraDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;

  /// Type of solver: sparse or direct?
  SparseSolverType m_solver_type;

  /// Pointer to system matrix
  std::shared_ptr<matrix_type> m_matrix;

  /// Pointer to RHS vector
  std::shared_ptr<multivector_type> m_rhs;

  /// Pointer to vector of unknowns
  std::shared_ptr<multivector_type> m_x;

  /// Tpetra linear solver
  std::unique_ptr<trilinos_lin_solver_type> m_linear_solver;

  /// Tpetra preconditioner
  std::shared_ptr<preconditioner_type> m_preconditioner;
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::LSTpetra(
    const SparseSolverType solver_type)
    : m_solver_type(solver_type)
{
  if (solver_type == eIterative)
  {
    m_linear_solver =
        std::unique_ptr<trilinos_lin_solver_type>(new trilinos_iterative_solver_type());
  }
  else if (solver_type == eDirect)
  {
    m_linear_solver = std::unique_ptr<trilinos_lin_solver_type>(new trilinos_direct_solver_type());
  }
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::~LSTpetra()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::initialize_solver(
    std::shared_ptr<matrix_type> const &matrix, std::shared_ptr<multivector_type> const &rhs,
    std::shared_ptr<multivector_type> const &x, const bool print_output)
{
  // --------------------------------------------
  // Configure operator, left and right-hand side
  // --------------------------------------------

  m_matrix = matrix;
  m_rhs    = rhs;
  m_x      = x;

  // --------------------------------------------
  // Configure iterative solver
  // --------------------------------------------
  m_linear_solver->initialize_solver(m_matrix, m_x, m_rhs, print_output);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::connect_preconditioner(
    std::shared_ptr<preconditioner_type> &preconditioner)
{
  m_preconditioner = preconditioner;
  m_linear_solver->connect_preconditioner(m_preconditioner);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::disconnect_preconditioner()
{
  m_linear_solver->release_preconditioner();
  m_preconditioner.reset();
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::update_after_mat_graph_change(
    std::shared_ptr<matrix_type> const &matrix, std::shared_ptr<multivector_type> const &rhs,
    std::shared_ptr<multivector_type> const &x, const bool print_output)
{
  m_linear_solver->update_after_mat_graph_change(m_matrix, m_x, m_rhs, print_output);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::update_after_mat_values_change(
    std::shared_ptr<matrix_type> const &matrix, std::shared_ptr<multivector_type> const &rhs,
    std::shared_ptr<multivector_type> const &x, const bool print_output)
{
  m_linear_solver->update_after_mat_values_change(m_matrix, m_x, m_rhs, print_output);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
std::shared_ptr<TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> LSTpetra<
    Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::matrix()
{
  return m_matrix;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
const std::shared_ptr<TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> LSTpetra<
    Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::matrix() const
{
  return m_matrix;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
std::shared_ptr<TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> LSTpetra<
    Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::rhs()
{
  return m_rhs;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
const std::shared_ptr<TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> LSTpetra<
    Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::rhs() const
{
  return m_rhs;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
std::shared_ptr<
    typename LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::multivector_type>
LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::solution_vector()
{
  return m_x;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
const std::shared_ptr<
    typename LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::multivector_type>
LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::solution_vector() const
{
  return m_x;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
std::shared_ptr<
    typename LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::preconditioner_type>
LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::preconditioner()
{
  return m_preconditioner;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
const std::shared_ptr<
    typename LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::preconditioner_type>
LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::preconditioner() const
{
  return m_preconditioner;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::set_preconditioner(
    std::shared_ptr<preconditioner_type> &preconditioner)
{
  m_preconditioner = preconditioner;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::residual(
    multivector_type &res) const
{
  res.assign(*m_rhs);

  // void apply(X, Y, transp, alpha, beta)
  // apply: Y := beta*Y + alpha*Op(A)*X with Op(A) = A or Op(A) = transpose(A)
  m_matrix->apply(*m_x, res, false, -1.0, 1.0);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void LSTpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::solve(const Uint nb_iter,
                                                                         const Real tol,
                                                                         const bool verbose)
{
  m_linear_solver->solve(nb_iter, tol, verbose);
}

// ----------------------------------------------------------------------------

//  Defer instantiation of Real LSTpetra to the cpp file
extern template class LSTpetra<Real>;

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif
