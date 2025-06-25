#include "physics/euler/Euler3DProperties.hpp"

namespace pdekit
{

namespace physics
{

// ----------------------------------------------------------------------------

Euler3DProperties::Euler3DProperties()
{
  gamma         = 1.4;
  R             = 287.05;
  gamma_minus_1 = gamma - 1.;
}

// ----------------------------------------------------------------------------

Euler3DProperties::~Euler3DProperties()
{
}

// ----------------------------------------------------------------------------

void Euler3DProperties::print() const
{
  std::cout << "Coordinates [" << coords[X] << "," << coords[Y] << "," << coords[Z] << "]"
            << std::endl;
  std::cout << "u           = " << vars << std::endl;
  std::cout << "grad(u)     = " << grad_vars << std::endl;
  std::cout << "pressure    = " << P << std::endl;
  std::cout << "temperature = " << T << std::endl;
  std::cout << "Mach number = " << Ma << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace physics

} // namespace pdekit
