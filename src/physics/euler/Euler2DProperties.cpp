#include "physics/euler/Euler2DProperties.hpp"

namespace pdekit
{

namespace physics
{

// ----------------------------------------------------------------------------

Euler2DProperties::Euler2DProperties()
{
  gamma         = 1.4;
  R             = 287.05;
  gamma_minus_1 = gamma - 1.;
}

// ----------------------------------------------------------------------------

Euler2DProperties::~Euler2DProperties()
{
}

// ----------------------------------------------------------------------------

void Euler2DProperties::print() const
{
  std::cout << "Coordinates [" << coords[X] << "," << coords[Y] << "]" << std::endl;
  std::cout << "u           = " << vars << std::endl;
  std::cout << "grad(u)     = " << grad_vars << std::endl;
  std::cout << "pressure    = " << P << std::endl;
  std::cout << "temperature = " << T << std::endl;
  std::cout << "Mach number = " << Ma << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace physics

} // namespace pdekit
