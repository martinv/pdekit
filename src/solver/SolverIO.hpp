#ifndef PDEKIT_Solver_Solver_IO_hpp
#define PDEKIT_Solver_Solver_IO_hpp

#include "common/PDEKit.hpp"
#include "math/DenseDVec.hpp"

namespace pdekit
{

namespace solver
{

class SolverIO
{
  public:
  /// Default constructor
  SolverIO() = default;

  /// Default destructor
  ~SolverIO() = default;

  static void print_iter_and_res_norm(const Uint iter, const Real CFL,
                                      const math::DenseDVec<Real> &norm);

  static void print_iter_and_res_norm_w_timing(const Uint iter, const Real CFL,
                                               const math::DenseDVec<Real> &norm,
                                               const double cpu_duration);
};

} // namespace solver

} // namespace pdekit

#endif
