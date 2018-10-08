#include "physics/scalar/AdvDiff2DFOSProperties.hpp"

namespace pdekit
{

namespace physics
{

// ----------------------------------------------------------------------------

AdvDiff2DFOSProperties::AdvDiff2DFOSProperties()
{
}

// ----------------------------------------------------------------------------

AdvDiff2DFOSProperties::~AdvDiff2DFOSProperties()
{
}

// ----------------------------------------------------------------------------

void AdvDiff2DFOSProperties::print() const
{
  std::cout << "Coordinates [" << coords[X0] << "," << coords[X1] << "]" << std::endl;
  std::cout << "vars   = " << vars << std::endl;
  std::cout << "grad_u = " << grad_vars << std::endl;
  std::cout << "lambda = " << V << std::endl;
  std::cout << "mu     = " << mu << std::endl;
  std::cout << "Lr     = " << Lr << std::endl;
}

} // namespace physics

} // namespace pdekit
