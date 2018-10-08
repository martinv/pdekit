#ifndef PDEKIT_RD_Solver_Implicit_Solver_hpp
#define PDEKIT_RD_Solver_Implicit_Solver_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace solver
{

class ImplicitSolver
{
  public:
  /// Default constructor
  ImplicitSolver();

  /// Default destructor
  virtual ~ImplicitSolver();

  virtual void make_implicit_time_step(const Real dt);
};

} // namespace solver

} // namespace pdekit

#endif
