#include "physics/scalar/Adv2DProperties.hpp"

namespace pdekit
{

namespace physics
{

// ----------------------------------------------------------------------------

Adv2DProperties::Adv2DProperties()
{
}

// ----------------------------------------------------------------------------

Adv2DProperties::~Adv2DProperties()
{
}

// ----------------------------------------------------------------------------

void Adv2DProperties::print() const
{
  std::cout << "Coordinates [" << coords[X] << "," << coords[Y] << "]" << std::endl;
  std::cout << "vars   = " << vars << std::endl;
  std::cout << "grad_u = " << grad_vars << std::endl;
  std::cout << "lambda = " << V << std::endl;
  std::cout << "mu     = " << mu << std::endl;
}

} // namespace physics

} // namespace pdekit
