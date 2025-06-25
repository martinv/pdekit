#ifndef PDEKIT_Linear_System_Tpetra_Direct_Solver_hpp
#define PDEKIT_Linear_System_Tpetra_Direct_Solver_hpp

#include <iostream>

#include "PDEKit_Config.hpp"
#include "common/PDEKit.hpp"
#include "linear_system/TpetraLinSolver.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
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
class TpetraDirectSolver
    : public TpetraLinSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>
{
  public:
  /// TYPEDEFS
  using base_type           = TpetraLinSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
  using matrix_type         = typename base_type::matrix_type;
  using vector_type         = typename base_type::vector_type;
  using preconditioner_type = typename base_type::preconditioner_type;

  /// Constructor
  TpetraDirectSolver();

  /// Destructor
  ~TpetraDirectSolver() override;

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

  void solve(const Uint nb_iter, const Real tol, const bool verbose) override;

  /// Set preconditioner. Does nothing because we don't use preconditioners
  /// with direct solvers
  void connect_preconditioner(std::shared_ptr<preconditioner_type> &preconditioner) override;

  /// Release preconditioner. Does nothing because we don't use
  /// preconditioners with direct solvers
  void release_preconditioner() override;

  private:
#if PDEKIT_HAVE_TRILINOS
  /// -----------------------------
  /// DATA FOR SPARSE DIRECT SOLVER
  /// -----------------------------

  Teuchos::RCP<Amesos2::Solver<typename matrix_type::trilinos_matrix_type,
                               typename vector_type::trilinos_multivector_type>>
      m_solver;
#endif
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraDirectSolver()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::~TpetraDirectSolver()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::initialize_solver(
    const std::shared_ptr<matrix_type> &mat_lhs, const std::shared_ptr<vector_type> &vec_x,
    const std::shared_ptr<vector_type> &vec_rhs, const bool print_output)
{
#if PDEKIT_HAVE_TRILINOS
  using iaccess = internal::TpetraInternalAccess;

  const std::string solver_name = "KLU2";
  // const std::string solver_name = "SuperLU";
  // const std::string solver_name = "Basker";
  m_solver.reset();
  m_solver = Amesos2::create<typename matrix_type::trilinos_matrix_type,
                             typename vector_type::trilinos_multivector_type>(
      solver_name, iaccess::get_matrix_data(*mat_lhs), iaccess::get_vector_data(*vec_x),
      iaccess::get_vector_data(*vec_rhs));
  m_solver->symbolicFactorization();
// m_solver->numericFactorization();
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
    update_after_mat_graph_change(const std::shared_ptr<matrix_type> &mat_lhs,
                                  const std::shared_ptr<vector_type> &vec_x,
                                  const std::shared_ptr<vector_type> &vec_rhs,
                                  const bool print_output)
{
#if PDEKIT_HAVE_TRILINOS
  m_solver->symbolicFactorization();
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
    update_after_mat_values_change(const std::shared_ptr<matrix_type> &mat_lhs,
                                   const std::shared_ptr<vector_type> &vec_x,
                                   const std::shared_ptr<vector_type> &vec_rhs,
                                   const bool print_output)
{
#if PDEKIT_HAVE_TRILINOS
  m_solver->numericFactorization();
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::solve(
    const Uint nb_iter, const Real tol, const bool verbose)
{
#if PDEKIT_HAVE_TRILINOS
  m_solver->solve();
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::connect_preconditioner(
    std::shared_ptr<preconditioner_type> &preconditioner)
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                        classic>::release_preconditioner()
{
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace ls

} // namespace pdekit

#endif
