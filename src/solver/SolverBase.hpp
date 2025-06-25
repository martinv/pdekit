#ifndef PDEKIT_Solver_Solver_Base_hpp
#define PDEKIT_Solver_Solver_Base_hpp

#include <array>
#include <vector>

#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"

namespace pdekit
{

namespace solver
{

// ----------------------------------------------------------------------------

enum SolverOption : unsigned short
{
  RecomputePreconditioner = 0
};

// ----------------------------------------------------------------------------

enum SolverScalFn
{
  time_step = 0
};

enum SolverVecFn
{
  solution  = 0,
  residuals = 1,
  sources   = 2
};

// ----------------------------------------------------------------------------

class SolverBase
{
  public:
  /// Constructor
  SolverBase();

  /// Destructor
  virtual ~SolverBase();

  /// Flag saying whether solver is continuous or not
  virtual bool is_continuous() const = 0;

  /// Flag saying whether solver is explicit
  virtual bool is_explicit() const = 0;

  /// Configure a vector function field
  virtual void set_vec_function(const SolverVecFn vec_fn_id,
                                const interpolation::VectorMeshFunction<Real>::ptr &vec_fn);

  /// Get a scalar mesh function
  interpolation::ScalarMeshFunction<Real>::ptr scal_function(const SolverScalFn scal_fn) const;

  /// Get a vector mesh function
  interpolation::VectorMeshFunction<Real>::ptr vec_function(const SolverVecFn vec_fn) const;

  /// Assemble both lhs and rhs simultaneously
  virtual void assemble_lhs_and_rhs(const Real CFL) = 0;

  /// Assemble only rhs (assuming that left-hand side) is already
  /// assembled
  virtual void assemble_rhs() = 0;

  /// Solve the underlying system
  virtual void solve(std::vector<SolverOption> const &solver_options = {}) = 0;

  /// Set the time step
  void set_cfl(const Real CFL);

  protected:
  /// Time step
  Real m_cfl;

  private:
  std::array<interpolation::ScalarMeshFunction<Real>::ptr, 1> m_scal_functions;
  std::array<interpolation::VectorMeshFunction<Real>::ptr, 3> m_vec_functions;
};

} // namespace solver

} // namespace pdekit

#endif
