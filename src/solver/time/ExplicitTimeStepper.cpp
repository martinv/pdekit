#include "solver/time/ExplicitTimeStepper.hpp"

namespace pdekit
{

namespace solver
{

namespace time
{

// ----------------------------------------------------------------------------

ExplicitTimeStepperBase::ExplicitTimeStepperBase()
{
}

ExplicitTimeStepperBase::~ExplicitTimeStepperBase()
{
}

void ExplicitTimeStepperBase::setup(const Uint nb_fields, const Uint nb_entries)
{
}

// ----------------------------------------------------------------------------

ExplicitEuler::ExplicitEuler()
{
}

ExplicitEuler::~ExplicitEuler()
{
}

void ExplicitEuler::advance_in_time(SolverBase &scheme, const Real CFL)
{
  scheme.assemble_lhs_and_rhs(CFL);
  scheme.solve({});
}

// ----------------------------------------------------------------------------
// 3rd order TVD Runge-Kutta
// ----------------------------------------------------------------------------

RK3TVD::RK3TVD()
    : res_tmp(std::make_shared<interpolation::VectorMeshFunction<Real>>("", "res_tmp")),
      u_1(std::make_shared<interpolation::VectorMeshFunction<Real>>("", "u_1")),
      u_2(std::make_shared<interpolation::VectorMeshFunction<Real>>("", "u_2"))
{
}

RK3TVD::~RK3TVD()
{
}

void RK3TVD::setup(const Uint nb_fields, const Uint nb_entries)
{
  ExplicitTimeStepperBase::setup(nb_fields, nb_entries);

  res_tmp->resize(nb_fields, nb_entries);
  u_1->resize(nb_fields, nb_entries);
  u_2->resize(nb_fields, nb_entries);
}

void RK3TVD::advance_in_time(SolverBase &scheme, const Real CFL)
{
  std::shared_ptr<interpolation::VectorMeshFunction<Real>> u_n =
      scheme.vec_function(SolverVecFn::solution);
  std::shared_ptr<interpolation::VectorMeshFunction<Real>> res_n =
      scheme.vec_function(SolverVecFn::residuals);

  // Assemble lhs, rhs and evaluate the time step
  scheme.assemble_lhs_and_rhs(CFL);

  interpolation::ScalarMeshFunction<Real> const &dt =
      *(scheme.scal_function(SolverScalFn::time_step));

  // -------------------------------
  // Step 1: u_1 = u_n + dt * L(u_n)
  // -------------------------------
  (*u_1) = (*u_n) - dt * (*res_n);

  // -------------------------------------------------------
  // Step 2: u_2 = 3/4 * u_n + 1/4 * u_1 + 1/4 * dt * L(u_1)
  // -------------------------------------------------------
  scheme.set_vec_function(SolverVecFn::solution, u_1);
  scheme.set_vec_function(SolverVecFn::residuals, res_tmp);
  res_tmp->fill(0.0);

  scheme.assemble_rhs();

  (*u_2) = 3. / 4. * (*u_n) + 1. / 4. * (*u_1) - 1. / 4. * (dt * (*res_tmp));

  // -----------------------------------------------------------
  // Step 3: u_(n+1) = 1/3 * u_n + 2/3 * u_2 + 2/3 * dt * L(u_2)
  // -----------------------------------------------------------
  scheme.set_vec_function(SolverVecFn::solution, u_2);
  res_tmp->fill(0.0);

  scheme.assemble_rhs();

  // Use temporarily the value of u_1 to store the result:
  (*u_1) = 1. / 3. * (*u_n) + 2. / 3 * (*u_2) - 2. / 3. * (dt * (*res_tmp));

  // Now put the result in u_n again:
  (*u_n) = (*u_1);

  // Cleanup: let the scheme use correct solution and residuals again:
  scheme.set_vec_function(SolverVecFn::solution, u_n);
  scheme.set_vec_function(SolverVecFn::residuals, res_n);
}

// ----------------------------------------------------------------------------

} // namespace time

} // namespace solver

} // namespace pdekit
