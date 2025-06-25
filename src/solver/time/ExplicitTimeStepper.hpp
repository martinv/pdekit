#ifndef PDEKIT_Solver_Time_Explicit_Time_Stepper_hpp
#define PDEKIT_Solver_Time_Explicit_Time_Stepper_hpp

#include <memory>

#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "solver/ExplicitSolver.hpp"
#include "solver/SolverBase.hpp"
#include "solver/rdm/RDTimeUpdate.hpp"

namespace pdekit
{

namespace solver
{

namespace time
{

// ----------------------------------------------------------------------------

class ExplicitTimeStepperBase
{
  public:
  /// Default constructor
  ExplicitTimeStepperBase();

  /// Default destructor
  virtual ~ExplicitTimeStepperBase();

  /// Resize the data
  virtual void setup(const Uint nb_fields, const Uint nb_entries);

  /// Perform one time iteration
  virtual void advance_in_time(SolverBase &solver, const Real CFL) = 0;

  protected:
};

// ----------------------------------------------------------------------------
// Forward Euler time discretization
// ----------------------------------------------------------------------------

class ExplicitEuler : public ExplicitTimeStepperBase
{
  public:
  /// Default constructor
  ExplicitEuler();

  /// Default destructor
  ~ExplicitEuler() override;

  /// Make one time step
  void advance_in_time(SolverBase &scheme, const Real CFL) override;

  private:
};

// ----------------------------------------------------------------------------
// Runge-Kutta time discretization
// ----------------------------------------------------------------------------

class RK3TVD : public ExplicitTimeStepperBase
{
  public:
  /// Default constructor
  RK3TVD();

  /// Default destructor
  ~RK3TVD() override;

  /// Resize the data
  void setup(const Uint nb_fields, const Uint nb_entries) override;

  /// Make one time step
  void advance_in_time(SolverBase &scheme, const Real CFL) override;

  private:
  std::shared_ptr<interpolation::VectorMeshFunction<Real>> res_tmp;
  std::shared_ptr<interpolation::VectorMeshFunction<Real>> u_1;
  std::shared_ptr<interpolation::VectorMeshFunction<Real>> u_2;
};

// ----------------------------------------------------------------------------

} // namespace time

} // namespace solver

} // namespace pdekit

#endif
