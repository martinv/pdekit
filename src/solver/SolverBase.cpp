#include "solver/SolverBase.hpp"

namespace pdekit
{

namespace solver
{

// ----------------------------------------------------------------------------

SolverBase::SolverBase() : m_cfl(1.0)
{
  m_scal_functions[SolverScalFn::time_step] =
      std::make_shared<interpolation::ScalarMeshFunction<Real>>("", "time_step");
  m_scal_functions[SolverScalFn::time_step]->resize(0);
}

// ----------------------------------------------------------------------------

SolverBase::~SolverBase()
{
}

// ----------------------------------------------------------------------------

void SolverBase::set_vec_function(const SolverVecFn vec_fn_id,
                                  const interpolation::VectorMeshFunction<Real>::ptr &vec_fn)
{
  m_vec_functions[vec_fn_id] = vec_fn;
}

// ----------------------------------------------------------------------------

interpolation::ScalarMeshFunction<Real>::ptr SolverBase::scal_function(
    const SolverScalFn scal_fn) const
{
  return m_scal_functions[scal_fn];
}

// ----------------------------------------------------------------------------

interpolation::VectorMeshFunction<Real>::ptr SolverBase::vec_function(
    const SolverVecFn vec_fn) const
{
  return m_vec_functions[vec_fn];
}

// ----------------------------------------------------------------------------

void SolverBase::set_cfl(const Real CFL)
{
  m_cfl = CFL;
}

// ----------------------------------------------------------------------------

} // namespace solver

} // namespace pdekit
