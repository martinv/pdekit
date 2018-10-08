#include "mesh/shape_function/CarnevaliExpansionTriag.hpp"

#include "math/DenseConstVecView.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

// ----------------------------------------------------------------------------

CarnevaliExpansionTriag::CarnevaliExpansionTriag() : ModalExpansion()
{
}

// ----------------------------------------------------------------------------

CarnevaliExpansionTriag::CarnevaliExpansionTriag(const Uint poly_order) : ModalExpansion(poly_order)
{
}

// ----------------------------------------------------------------------------

CarnevaliExpansionTriag::~CarnevaliExpansionTriag()
{
}

// ----------------------------------------------------------------------------

/// Return the number of modes
Uint CarnevaliExpansionTriag::nb_modes() const
{
  return (P + 1) * (P + 2) / 2;
}

// ----------------------------------------------------------------------------

Uint CarnevaliExpansionTriag::topo_dim() const
{
  return _2D;
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionTriag::evaluate_in_one_point(const math::DenseDVec<Real> &point,
                                                    math::DenseDVec<Real> &values)
{
  values.resize((P + 1) * (P + 2) / 2); // (P+1)*(P+2)/2 is the total number of modes

  const Real xi0 = point[XI0];
  const Real xi1 = point[XI1];

  const Real L0 = -0.5 * (xi0 + xi1);
  const Real L1 = 0.5 * (1.0 + xi0);
  const Real L2 = 0.5 * (1.0 + xi1);

  values[0] = L0;
  values[1] = L1;
  values[2] = L2;

  if (P >= 2)
  {
    values[3] = -2.0 * L0 * L1;
    values[4] = -2.0 * L1 * L2;
    values[5] = -2.0 * L2 * L0;
  }

  if (P >= 3)
  {
    values[6] = -2.0 * L0 * L1 * (L1 - L0);
    values[7] = -2.0 * L1 * L2 * (L2 - L1);
    values[8] = -2.0 * L2 * L0 * (L0 - L2);
    values[9] = L0 * L1 * L2;
  }
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionTriag::Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                                                 math::DenseDMat<Real> &values)
{
  const Uint nb_pts = coordinates.rows();

  values.resize(nb_pts, (P + 1) * (P + 2) / 2); // (P+1)*(P+2)/2 is the total number of modes

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {
    const Real xi0 = coordinates(pt, XI0);
    const Real xi1 = coordinates(pt, XI1);
    values(pt, 0)  = -0.5 * (xi0 + xi1);
    values(pt, 1)  = 0.5 * (1.0 + xi0);
    values(pt, 2)  = 0.5 * (1.0 + xi1);
  }

  if (P >= 2)
  {
    for (Uint pt = 0; pt < nb_pts; ++pt)
    {
      const Real xi0 = coordinates(pt, XI0);
      const Real xi1 = coordinates(pt, XI1);

      const Real L0 = -0.5 * (xi0 + xi1);
      const Real L1 = 0.5 * (1.0 + xi0);
      const Real L2 = 0.5 * (1.0 + xi1);

      values(pt, 3) = -2.0 * L0 * L1;
      values(pt, 4) = -2.0 * L1 * L2;
      values(pt, 5) = -2.0 * L2 * L0;
    }
  }

  if (P >= 3)
  {
    for (Uint pt = 0; pt < nb_pts; ++pt)
    {
      const Real xi0 = coordinates(pt, XI0);
      const Real xi1 = coordinates(pt, XI1);

      const Real L0 = -0.5 * (xi0 + xi1);
      const Real L1 = 0.5 * (1.0 + xi0);
      const Real L2 = 0.5 * (1.0 + xi1);

      values(pt, 6) = -2.0 * L0 * L1 * (L1 - L0);
      values(pt, 7) = -2.0 * L1 * L2 * (L2 - L1);
      values(pt, 8) = -2.0 * L2 * L0 * (L0 - L2);
      values(pt, 9) = L0 * L1 * L2;
    }
  }
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionTriag::evaluate_derivatives_in_one_point(const math::DenseDVec<Real> &point,
                                                                const math::DenseDVec<Real> &values,
                                                                math::DenseDMat<Real> &derivatives)
{
  derivatives.resize((P + 1) * (P + 2) / 2,
                     2); // (P+1)*(P+2)/2 is the total number of modes
  derivatives.fill(0.0);

  const Real xi0 = point[XI0];
  const Real xi1 = point[XI1];

  const Real L0 = -0.5 * (xi0 + xi1);
  const Real L1 = 0.5 * (1.0 + xi0);
  const Real L2 = 0.5 * (1.0 + xi1);

  const Real dL0_dxi0 = -0.5;
  const Real dL0_dxi1 = -0.5;
  const Real dL1_dxi0 = 0.5;
  const Real dL1_dxi1 = 0.0;
  const Real dL2_dxi0 = 0.0;
  const Real dL2_dxi1 = 0.5;

  derivatives(0, XI0) = dL0_dxi0;
  derivatives(0, XI1) = dL0_dxi1;
  derivatives(1, XI0) = dL1_dxi0;
  derivatives(1, XI1) = dL1_dxi1;
  derivatives(2, XI0) = dL2_dxi0;
  derivatives(2, XI1) = dL2_dxi1;

  if (P >= 2)
  {
    derivatives(3, XI0) = -2.0 * (dL0_dxi0 * L1 + L0 * dL1_dxi0);
    derivatives(3, XI1) = -2.0 * (dL0_dxi1 * L1 + L0 * dL1_dxi1);
    derivatives(4, XI0) = -2.0 * (dL1_dxi0 * L2 + L1 * dL2_dxi0);
    derivatives(4, XI1) = -2.0 * (dL1_dxi1 * L2 + L1 * dL2_dxi1);
    derivatives(5, XI0) = -2.0 * (dL2_dxi0 * L0 + L2 * dL0_dxi0);
    derivatives(5, XI1) = -2.0 * (dL2_dxi1 * L0 + L2 * dL0_dxi1);
  }

  if (P >= 3)
  {
    derivatives(6, XI0) = -2.0 * (dL0_dxi0 * L1 * (L1 - L0) + L0 * dL1_dxi0 * (L1 - L0) +
                                  L0 * L1 * (dL1_dxi0 - dL0_dxi0));
    derivatives(6, XI1) = -2.0 * (dL0_dxi1 * L1 * (L1 - L0) + L0 * dL1_dxi1 * (L1 - L0) +
                                  L0 * L1 * (dL1_dxi1 - dL0_dxi1));

    derivatives(7, XI0) = -2.0 * (dL1_dxi0 * L2 * (L2 - L1) + L1 * dL2_dxi0 * (L2 - L1) +
                                  L1 * L2 * (dL2_dxi0 - dL1_dxi0));
    derivatives(7, XI1) = -2.0 * (dL1_dxi1 * L2 * (L2 - L1) + L1 * dL2_dxi1 * (L2 - L1) +
                                  L1 * L2 * (dL2_dxi1 - dL1_dxi1));

    derivatives(8, XI0) = -2.0 * (dL2_dxi0 * L0 * (L0 - L2) + L2 * dL0_dxi0 * (L0 - L2) +
                                  L2 * L0 * (dL0_dxi0 - dL2_dxi0));
    derivatives(8, XI1) = -2.0 * (dL2_dxi1 * L0 * (L0 - L2) + L2 * dL0_dxi1 * (L0 - L2) +
                                  L2 * L0 * (dL0_dxi1 - dL2_dxi1));

    derivatives(9, XI0) = dL0_dxi0 * L1 * L2 + L0 * dL1_dxi0 * L2 + L0 * L1 * dL2_dxi0;
    derivatives(9, XI1) = dL0_dxi1 * L1 * L2 + L0 * dL1_dxi1 * L2 + L0 * L1 * dL2_dxi1;
  }
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionTriag::Vandermonde_matrix_derivatives(
    const math::DenseDMat<Real> &coordinates, std::vector<math::DenseDMat<Real>> &derivative_values)
{
  const Uint nb_pts = coordinates.rows();

  derivative_values.resize(_2D);
  derivative_values[XI0].resize(nb_pts, (P + 1) * (P + 2) / 2);
  derivative_values[XI1].resize(nb_pts, (P + 1) * (P + 2) / 2);

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {
    const Real dL0_dxi0 = -0.5;
    const Real dL0_dxi1 = -0.5;
    const Real dL1_dxi0 = 0.5;
    const Real dL1_dxi1 = 0.0;
    const Real dL2_dxi0 = 0.0;
    const Real dL2_dxi1 = 0.5;

    derivative_values[XI0](pt, 0) = dL0_dxi0;
    derivative_values[XI1](pt, 0) = dL0_dxi1;
    derivative_values[XI0](pt, 1) = dL1_dxi0;
    derivative_values[XI1](pt, 1) = dL1_dxi1;
    derivative_values[XI0](pt, 2) = dL2_dxi0;
    derivative_values[XI1](pt, 2) = dL2_dxi1;
  }

  if (P >= 2)
  {
    for (Uint pt = 0; pt < nb_pts; ++pt)
    {
      const Real xi0 = coordinates(pt, XI0);
      const Real xi1 = coordinates(pt, XI1);

      const Real L0 = -0.5 * (xi0 + xi1);
      const Real L1 = 0.5 * (1.0 + xi0);
      const Real L2 = 0.5 * (1.0 + xi1);

      const Real dL0_dxi0 = -0.5;
      const Real dL0_dxi1 = -0.5;
      const Real dL1_dxi0 = 0.5;
      const Real dL1_dxi1 = 0.0;
      const Real dL2_dxi0 = 0.0;
      const Real dL2_dxi1 = 0.5;

      derivative_values[XI0](pt, 3) = -2.0 * (dL0_dxi0 * L1 + L0 * dL1_dxi0);
      derivative_values[XI1](pt, 3) = -2.0 * (dL0_dxi1 * L1 + L0 * dL1_dxi1);
      derivative_values[XI0](pt, 4) = -2.0 * (dL1_dxi0 * L2 + L1 * dL2_dxi0);
      derivative_values[XI1](pt, 4) = -2.0 * (dL1_dxi1 * L2 + L1 * dL2_dxi1);
      derivative_values[XI0](pt, 5) = -2.0 * (dL2_dxi0 * L0 + L2 * dL0_dxi0);
      derivative_values[XI1](pt, 5) = -2.0 * (dL2_dxi1 * L0 + L2 * dL0_dxi1);
    }
  }

  if (P >= 3)
  {
    for (Uint pt = 0; pt < nb_pts; ++pt)
    {
      const Real xi0 = coordinates(pt, XI0);
      const Real xi1 = coordinates(pt, XI1);

      const Real L0 = -0.5 * (xi0 + xi1);
      const Real L1 = 0.5 * (1.0 + xi0);
      const Real L2 = 0.5 * (1.0 + xi1);

      const Real dL0_dxi0 = -0.5;
      const Real dL0_dxi1 = -0.5;
      const Real dL1_dxi0 = 0.5;
      const Real dL1_dxi1 = 0.0;
      const Real dL2_dxi0 = 0.0;
      const Real dL2_dxi1 = 0.5;

      derivative_values[XI0](pt, 6) =
          -2.0 *
          (dL0_dxi0 * L1 * (L1 - L0) + L0 * dL1_dxi0 * (L1 - L0) + L0 * L1 * (dL1_dxi0 - dL0_dxi0));
      derivative_values[XI1](pt, 6) =
          -2.0 *
          (dL0_dxi1 * L1 * (L1 - L0) + L0 * dL1_dxi1 * (L1 - L0) + L0 * L1 * (dL1_dxi1 - dL0_dxi1));

      derivative_values[XI0](pt, 7) =
          -2.0 *
          (dL1_dxi0 * L2 * (L2 - L1) + L1 * dL2_dxi0 * (L2 - L1) + L1 * L2 * (dL2_dxi0 - dL1_dxi0));
      derivative_values[XI1](pt, 7) =
          -2.0 *
          (dL1_dxi1 * L2 * (L2 - L1) + L1 * dL2_dxi1 * (L2 - L1) + L1 * L2 * (dL2_dxi1 - dL1_dxi1));

      derivative_values[XI0](pt, 8) =
          -2.0 *
          (dL2_dxi0 * L0 * (L0 - L2) + L2 * dL0_dxi0 * (L0 - L2) + L2 * L0 * (dL0_dxi0 - dL2_dxi0));
      derivative_values[XI1](pt, 8) =
          -2.0 *
          (dL2_dxi1 * L0 * (L0 - L2) + L2 * dL0_dxi1 * (L0 - L2) + L2 * L0 * (dL0_dxi1 - dL2_dxi1));

      derivative_values[XI0](pt, 9) = dL0_dxi0 * L1 * L2 + L0 * dL1_dxi0 * L2 + L0 * L1 * dL2_dxi0;
      derivative_values[XI1](pt, 9) = dL0_dxi1 * L1 * L2 + L0 * dL1_dxi1 * L2 + L0 * L1 * dL2_dxi1;
    }
  }
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionTriag::is_leading_expansion_term(math::DenseDVec<bool> &is_leading_term)
{
  is_leading_term.resize((P + 1) * (P + 2) / 2);
  is_leading_term.fill(false);

  if (P == P0)
  {
    is_leading_term[0] = true;
    return;
  }

  if (P == P1)
  {
    // In this case, we have three linear modes
    is_leading_term[0] = true;
    is_leading_term[1] = true;
    is_leading_term[2] = true;
    return;
  }

  // In what follows, we exploit the fact that in simplex expansions, the
  // highes-order modes are those modes (p,q) for which p + q == P

  // is_leading_term[mode_id(P, 0)] = true;
  // is_leading_term[mode_id(P - 1, 1)] = true;

  for (Uint p = 0; p < P; ++p)
  {
    for (Uint q = 1; q < (P - p); ++q)
    {
      if ((p + q + 1) == P)
      {
        // Terms (p,q+1)
        // is_leading_term[mode_id(p, q + 1)] = true;
      }
    }
  }
}

// ----------------------------------------------------------------------------

void CarnevaliExpansionTriag::mode_poly_deg(math::DenseDVec<Uint> &poly_deg)
{
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit
