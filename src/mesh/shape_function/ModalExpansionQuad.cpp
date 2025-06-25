#include "mesh/shape_function/ModalExpansionQuad.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

// ----------------------------------------------------------------------------

ModalExpansionQuad::ModalExpansionQuad() : ModalExpansion()
{
}

// ----------------------------------------------------------------------------

ModalExpansionQuad::ModalExpansionQuad(const Uint poly_order) : ModalExpansion(poly_order)
{
}

// ----------------------------------------------------------------------------

ModalExpansionQuad::~ModalExpansionQuad()
{
}

// ----------------------------------------------------------------------------

/// Return the number of modes
Uint ModalExpansionQuad::nb_modes() const
{
  return ((P + 1) * (P + 1));
}

// ----------------------------------------------------------------------------

Uint ModalExpansionQuad::topo_dim() const
{
  return _2D;
}

// ----------------------------------------------------------------------------

void ModalExpansionQuad::evaluate_in_one_point(const math::DenseDVec<Real> &point,
                                               math::DenseDVec<Real> &values)
{
  values.resize((P + 1) * (P + 1)); // ((P+1)* (P+1)) is the total number of modes

  // If the polynomial order is equal to 0, we have one (constant) mode
  if (P == 0)
  {
    values[0] = 0.5; // (1/sqrt(2) * 1/sqrt(2) )
    return;
  }

  // Else we have at least four linear modes (corner modes)

  Uint mode_id = 0;

  for (Uint m = 0; m < (P + 1); ++m)
  {
    const Real gamma_m  = 2. / (2. * m + 1.);
    const Real mode_xi0 = m_JP(m, 0, 0, point[XI0]);

    for (Uint n = 0; n < (P + 1); ++n)
    {
      const Real gamma_n  = 2. / (2. * n + 1.);
      const Real mode_xi1 = m_JP(n, 0, 0, point[XI1]);
      values[mode_id++]   = mode_xi0 * mode_xi1 / (std::sqrt(gamma_m * gamma_n));
    }
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionQuad::Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                                            math::DenseDMat<Real> &values)
{
  const Uint nb_pts = coordinates.rows();

  values.resize(nb_pts,
                (P + 1) * (P + 1)); // (P+1)*(P+1) is the total number of modes

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {

    // If the polynomial order is equal to 0, we have one (constant) mode
    if (P == 0)
    {
      values(pt, 0) = 0.5; // (1/sqrt(2) * 1/sqrt(2) )
    }

    else
    {
      // Else we have at least four linear modes (corner modes)

      const Real xi0 = coordinates(pt, XI0);
      const Real xi1 = coordinates(pt, XI1);

      // Interior modes

      Uint mode_id = 0;

      for (Uint m = 0; m < (P + 1); ++m)
      {
        const Real gamma_m  = 2. / (2. * m + 1.);
        const Real mode_xi0 = m_JP(m, 0, 0, xi0);

        for (Uint n = 0; n < (P + 1); ++n)
        {
          const Real gamma_n    = 2. / (2. * n + 1.);
          const Real mode_xi1   = m_JP(n, 0, 0, xi1);
          values(pt, mode_id++) = (mode_xi0 * mode_xi1) / (std::sqrt(gamma_m * gamma_n));
        }
      }
    } // else
  }   // Loop over all coordinate points
}

// ----------------------------------------------------------------------------

void ModalExpansionQuad::evaluate_derivatives_in_one_point(const math::DenseDVec<Real> &point,
                                                           const math::DenseDVec<Real> &values,
                                                           math::DenseDMat<Real> &derivatives)
{
  derivatives.resize((P + 1) * (P + 1),
                     2); // ((P+1)*(P+1)) is the total number of modes
  derivatives.fill(0.0);

  // If the polynomial order is equal to 0, we have one (constant) mode

  if (P == 0)
  {
    derivatives(0, XI0) = 0.0;
    derivatives(0, XI1) = 0.0;
    return;
  }

  // Else we have at least four linear modes (corner modes)

  // The normalization factor for each mode is 1/sqrt(gamma(n)) with gamma(n)
  // = 2/(2n+1) Hence gamma(0) = 2, gamma(1) = 2/3 etc. So for psi_0(xi) *
  // psi_0(XI1)  we have factor = 1/sqrt(gamma(0)*gamma(1)) = 1/2 For psi0 *
  // psi_1 the factor is 1/( sqrt(gamma(0)*gamma(1) ) = 1/sqrt(2*2/3)) =
  // sqrt(3)/2

  Uint mode_id = 0;

  for (Uint m = 0; m < (P + 1); ++m)
  {
    const Real gamma_m      = 2. / (2. * m + 1.);
    const Real mode_xi0     = m_JP(m, 0, 0, point[XI0]);
    const Real der_mode_xi0 = m_JP.dx(m, 0, 0, point[XI0]);

    for (Uint n = 0; n < (P + 1); ++n)
    {
      const Real gamma_n      = 2. / (2. * n + 1.);
      const Real mode_xi1     = m_JP(n, 0, 0, point[XI1]);
      const Real der_mode_xi1 = m_JP.dx(n, 0, 0, point[XI1]);

      const Real factor = 1.0 / (std::sqrt(gamma_m * gamma_n));

      derivatives(mode_id, XI0) = factor * der_mode_xi0 * mode_xi1;
      derivatives(mode_id, XI1) = factor * mode_xi0 * der_mode_xi1;
      mode_id++;
    }
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionQuad::Vandermonde_matrix_derivatives(
    const math::DenseDMat<Real> &coordinates, std::vector<math::DenseDMat<Real>> &derivative_values)
{
  const Uint nb_pts = coordinates.rows();

  derivative_values.resize(2);
  derivative_values[XI0].resize(nb_pts, (P + 1) * (P + 1)); // (P+1) is the total number of modes
  derivative_values[XI1].resize(nb_pts, (P + 1) * (P + 1));

  math::DenseDMat<Real> &dVdXi0 = derivative_values[XI0];
  math::DenseDMat<Real> &dVdXi1 = derivative_values[XI1];

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {
    // If the polynomial order is equal to 0, we have one (constant) mode

    if (P == 0)
    {
      dVdXi0(pt, 0) = 0.0;
      dVdXi1(pt, 0) = 0.0;
    }

    else
    {
      // Else we have at least four linear modes (corner modes)

      const Real xi0 = coordinates(pt, XI0);
      const Real xi1 = coordinates(pt, XI1);

      // The normalization factor for each mode is 1/sqrt(gamma(n)) with
      // gamma(n) = 2/(2n+1)
      // Hence gamma(0) = 2, gamma(1) = 2/3 etc.
      // So for psi_0(xi) * psi_0(XI1)  we have factor =
      // 1/sqrt(gamma(0)*gamma(1)) = 1/2
      // For psi0 * psi_1 the factor is 1/( sqrt(gamma(0)*gamma(1) ) =
      // 1/sqrt(2*2/3)) = sqrt(3)/2

      Uint mode_id = 0;

      for (Uint m = 0; m < (P + 1); ++m)
      {
        const Real gamma_m      = 2. / (2. * m + 1.);
        const Real mode_xi0     = m_JP(m, 0, 0, xi0);
        const Real der_mode_xi0 = m_JP.dx(m, 0, 0, xi0);

        for (Uint n = 0; n < (P + 1); ++n)
        {
          const Real gamma_n      = 2. / (2. * n + 1.);
          const Real mode_xi1     = m_JP(n, 0, 0, xi1);
          const Real der_mode_xi1 = m_JP.dx(n, 0, 0, xi1);

          const Real factor = 1.0 / (std::sqrt(gamma_m * gamma_n));

          dVdXi0(pt, mode_id) = factor * der_mode_xi0 * mode_xi1;
          dVdXi1(pt, mode_id) = factor * mode_xi0 * der_mode_xi1;
          mode_id++;
        }
      }

    } // else
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionQuad::is_leading_expansion_term(math::DenseDVec<bool> &is_leading_term)
{
  is_leading_term.resize((P + 1) * (P + 1));
  is_leading_term.fill(false);

  Uint mode_id = 0;

  for (Uint m = 0; m < (P + 1); ++m)
  {
    for (Uint n = 0; n < (P + 1); ++n)
    {
      if ((m == P) || (n == P))
      {
        is_leading_term[mode_id++] = true;
      }
    }
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionQuad::mode_poly_deg(math::DenseDVec<Uint> &poly_deg)
{
  poly_deg.resize((P + 1) * (P + 1));

  Uint mode_id = 0;

  for (Uint m = 0; m < (P + 1); ++m)
  {
    for (Uint n = 0; n < (P + 1); ++n)
    {
      poly_deg[mode_id++] = m + n;
    }
  }
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit
