#ifndef PDEKIT_Linear_System_Tpetra_Lin_Solver_hpp
#define PDEKIT_Linear_System_Tpetra_Lin_Solver_hpp

#include <iostream>

#include "linear_system/TpetraCrsMatrix.hpp"
#include "linear_system/TpetraMultiVector.hpp"
#include "linear_system/preconditioner/TrilinosPC.hpp"

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------

enum LeftRightPrec
{
  ePrecLeft  = 0,
  ePrecRight = 1
};

enum SparseSolverType
{
  eIterative = 0,
  eDirect    = 1
};

// ----------------------------------------------------------------------------

namespace internal
{

// ----------------------------------------------------------------------------

template <typename Scalar        = TpetraDefaultTraits::Scalar,
          typename LocalOrdinal  = TpetraDefaultTraits::LocalOrdinal,
          typename GlobalOrdinal = TpetraDefaultTraits::GlobalOrdinal,
          typename Node          = TpetraDefaultTraits::Node,
          const bool classic     = TpetraDefaultTraits::classic>
class TpetraLinSolver
{
  public:
  /// TYPEDEFS

  using matrix_type         = TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
  using vector_type         = TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
  using preconditioner_type = TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;

  /// Constructor
  TpetraLinSolver();

  /// Destructor
  virtual ~TpetraLinSolver();

  /// Build the solver and its preconditioner
  virtual void initialize_solver(const std::shared_ptr<matrix_type> &mat_lhs,
                                 const std::shared_ptr<vector_type> &vec_x,
                                 const std::shared_ptr<vector_type> &vec_rhs,
                                 const bool print_output = false) = 0;

  /// Update the preconditioner and other data when the matrix graph changes
  /// (we added or removed elements from the matrix, thus changing the number
  /// of its nonzeros
  virtual void update_after_mat_graph_change(const std::shared_ptr<matrix_type> &mat_lhs,
                                             const std::shared_ptr<vector_type> &vec_x,
                                             const std::shared_ptr<vector_type> &vec_rhs,
                                             const bool print_output = false) = 0;

  /// Update the preconditioner and other data when the matrix or rhs change
  /// (i.e the matrix graph is constant, but values in the matrix or rhs (or
  /// both) were updated
  virtual void update_after_mat_values_change(const std::shared_ptr<matrix_type> &mat_lhs,
                                              const std::shared_ptr<vector_type> &vec_x,
                                              const std::shared_ptr<vector_type> &vec_rhs,
                                              const bool print_output = false) = 0;

  /// Run the actual solver
  virtual void solve(const Uint nb_iter, const Real tol, const bool verbose) = 0;

  /// Connect preconditioner
  virtual void connect_preconditioner(std::shared_ptr<preconditioner_type> &preconditioner) = 0;

  /// Release preconditioner
  virtual void release_preconditioner() = 0;

  private:
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraLinSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraLinSolver()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraLinSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::~TpetraLinSolver()
{
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace ls

} // namespace pdekit

#endif
