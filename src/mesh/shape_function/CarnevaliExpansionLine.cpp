#include "mesh/shape_function/CarnevaliExpansionLine.hpp"

#include "math/DenseConstVecView.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

// ----------------------------------------------------------------------------

CarnevaliExpansionLine::CarnevaliExpansionLine() : ModalExpansion()
{
}

// ----------------------------------------------------------------------------

CarnevaliExpansionLine::CarnevaliExpansionLine(const Uint poly_order) : ModalExpansion(poly_order)
{
}

// ----------------------------------------------------------------------------

CarnevaliExpansionLine::~CarnevaliExpansionLine()
{
}

// ----------------------------------------------------------------------------

/// Return the number of modes
Uint CarnevaliExpansionLine::nb_modes() const
{
  return (P + 1);
}

// ----------------------------------------------------------------------------

Uint CarnevaliExpansionLine::topo_dim() const
{
  return _1D;
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionLine::evaluate_in_one_point(const math::DenseDVec<Real> &point,
                                                   math::DenseDVec<Real> &values)
{
  values.resize(P + 1); // (P+1) is the total number of modes

  const Real xi0 = point[XI0];

  const Real L0 = 0.5 * (1. - xi0);
  const Real L1 = 0.5 * (1. + xi0);

  values[0] = L0;
  values[1] = L1;

  if (P >= 2)
  {
    values[2] = -2.0 * L0 * L1;
  }

  if (P >= 3)
  {
    values[3] = -2.0 * L0 * L1 * (L1 - L0);
  }
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionLine::Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                                                math::DenseDMat<Real> &values)
{
  const Uint nb_pts = coordinates.rows();

  values.resize(nb_pts, P + 1); // (P+1) is the total number of modes

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {
    const Real xi0 = coordinates(pt, XI0);

    const Real L0 = 0.5 * (1. - xi0);
    const Real L1 = 0.5 * (1. + xi0);

    values(pt, 0) = L0;
    values(pt, 1) = L1;
  }

  if (P >= 2)
  {
    for (Uint pt = 0; pt < nb_pts; ++pt)
    {
      const Real xi0 = coordinates(pt, XI0);

      const Real L0 = 0.5 * (1. - xi0);
      const Real L1 = 0.5 * (1. + xi0);

      values(pt, 2) = -2.0 * L0 * L1;
    }
  }

  if (P >= 3)
  {
    for (Uint pt = 0; pt < nb_pts; ++pt)
    {
      const Real xi0 = coordinates(pt, XI0);

      const Real L0 = 0.5 * (1. - xi0);
      const Real L1 = 0.5 * (1. + xi0);

      values(pt, 3) = -2.0 * L0 * L1 * (L1 - L0);
    }
  }
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionLine::evaluate_derivatives_in_one_point(const math::DenseDVec<Real> &point,
                                                               const math::DenseDVec<Real> &values,
                                                               math::DenseDMat<Real> &derivatives)
{
  derivatives.resize(P + 1, 1); // (P+1) is the total number of modes
  derivatives.fill(0.0);

  const Real xi0 = point[XI0];

  const Real L0 = 0.5 * (1. - xi0);
  const Real L1 = 0.5 * (1. + xi0);

  const Real dL0_dxi0 = -0.5;
  const Real dL1_dxi0 = 0.5;

  derivatives(0, XI0) = dL0_dxi0;
  derivatives(1, XI0) = dL1_dxi0;

  if (P >= 2)
  {
    derivatives(2, XI0) = -2.0 * (dL0_dxi0 * L1 + L0 * dL1_dxi0);
  }

  if (P >= 3)
  {
    derivatives(3, XI0) = -2.0 * (dL0_dxi0 * L1 * (L1 - L0) + L0 * dL1_dxi0 * (L1 - L0) +
                                  L0 * L1 * (dL1_dxi0 - dL0_dxi0));
  }
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionLine::Vandermonde_matrix_derivatives(
    const math::DenseDMat<Real> &coordinates, std::vector<math::DenseDMat<Real>> &derivative_values)
{
  const Uint nb_pts = coordinates.rows();

  derivative_values.resize(_1D);
  derivative_values[XI0].resize(nb_pts, P + 1);
  for (Uint pt = 0; pt < nb_pts; ++pt)
  {
    const Real dL0_dxi0 = -0.5;
    const Real dL1_dxi0 = 0.5;

    derivative_values[XI0](pt, 0) = dL0_dxi0;
    derivative_values[XI0](pt, 1) = dL1_dxi0;
  }

  if (P >= 2)
  {
    for (Uint pt = 0; pt < nb_pts; ++pt)
    {
      const Real xi0 = coordinates(pt, XI0);

      const Real L0 = 0.5 * (1. - xi0);
      const Real L1 = 0.5 * (1. + xi0);

      const Real dL0_dxi0 = -0.5;
      const Real dL1_dxi0 = 0.5;

      derivative_values[XI0](pt, 2) = -2.0 * (dL0_dxi0 * L1 + L0 * dL1_dxi0);
    }
  }

  if (P >= 3)
  {
    for (Uint pt = 0; pt < nb_pts; ++pt)
    {
      const Real xi0 = coordinates(pt, XI0);

      const Real L0 = 0.5 * (1. - xi0);
      const Real L1 = 0.5 * (1. + xi0);

      const Real dL0_dxi0 = -0.5;
      const Real dL1_dxi0 = 0.5;

      derivative_values[XI0](pt, 3) =
          -2.0 *
          (dL0_dxi0 * L1 * (L1 - L0) + L0 * dL1_dxi0 * (L1 - L0) + L0 * L1 * (dL1_dxi0 - dL0_dxi0));
    }
  }
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionLine::is_leading_expansion_term(math::DenseDVec<bool> &is_leading_term)
{
  is_leading_term.resize(P + 1);
  is_leading_term.fill(false);

  if (P == P0)
  {
    is_leading_term[0] = true;
    return;
  }

  if (P == P1)
  {
    // In this case, we have two linear modes
    is_leading_term[0] = true;
    is_leading_term[1] = true;
    return;
  }

  if (P == P2)
  {
    is_leading_term[2] = true;
  }

  if (P == P3)
  {
    is_leading_term[3] = true;
  }
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionLine::mode_poly_deg(math::DenseDVec<Uint> &poly_deg)
{
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit
