#include "mesh/shape_function/ModalExpansion.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

// ----------------------------------------------------------------------------

ModalExpansion::ModalExpansion() : P(0)
{
}

// ----------------------------------------------------------------------------

ModalExpansion::ModalExpansion(const Uint poly_order) : P(poly_order)
{
}

// ----------------------------------------------------------------------------

ModalExpansion::~ModalExpansion()
{
}

// ----------------------------------------------------------------------------

void ModalExpansion::set_polynomial_order(const Uint poly_order)
{
  P = poly_order;
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit
