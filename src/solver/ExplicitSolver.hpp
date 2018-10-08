#ifndef PDEKIT_RD_Solver_Explicit_Solver_hpp
#define PDEKIT_RD_Solver_Explicit_Solver_hpp

#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"

namespace pdekit
{

namespace solver
{

class ExplicitSolver
{
  public:
  /// Default constructor
  ExplicitSolver();

  /// Default destructor
  virtual ~ExplicitSolver();
};

} // namespace solver

} // namespace pdekit

#endif
