#include "mesh/shape_function/ModalExpansionHexa.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

// ----------------------------------------------------------------------------

ModalExpansionHexa::ModalExpansionHexa() : ModalExpansion()
{
}

// ----------------------------------------------------------------------------

ModalExpansionHexa::ModalExpansionHexa(const Uint poly_order) : ModalExpansion(poly_order)
{
}

// ----------------------------------------------------------------------------

ModalExpansionHexa::~ModalExpansionHexa()
{
}

// ----------------------------------------------------------------------------

/// Return the number of modes
Uint ModalExpansionHexa::nb_modes() const
{
  return ((P + 1) * (P + 1) * (P + 1));
}

// ----------------------------------------------------------------------------

Uint ModalExpansionHexa::topo_dim() const
{
  return _3D;
}

// ----------------------------------------------------------------------------

void ModalExpansionHexa::evaluate_in_one_point(const math::DenseDVec<Real> &point,
                                               math::DenseDVec<Real> &values)
{
  // ((P+1) * (P+1) * (P+1)) is the total number of modes
  values.resize((P + 1) * (P + 1) * (P + 1));

  // If the polynomial order is equal to 0, we have one (constant) mode
  if (P == 0)
  {
    values[0] = 0.25 * std::sqrt(2.); // (1/sqrt(2) * 1/sqrt(2) * 1/sqrt(2) )
    return;
  }

  // Else we have at least eight linear modes (corner modes)

  Uint mode_id = 0;

  for (Uint r = 0; r < (P + 1); ++r)
  {
    const Real gamma_r  = 2. / (2. * r + 1.);
    const Real mode_xi0 = m_JP(r, 0, 0, point[XI0]);

    for (Uint s = 0; s < (P + 1); ++s)
    {
      const Real gamma_s  = 2. / (2. * s + 1.);
      const Real mode_xi1 = m_JP(s, 0, 0, point[XI1]);

      for (Uint t = 0; t < (P + 1); ++t)
      {
        const Real gamma_t  = 2. / (2. * t + 1.);
        const Real mode_xi2 = m_JP(t, 0, 0, point[XI2]);
        values[mode_id++] =
            mode_xi0 * mode_xi1 * mode_xi2 / (std::sqrt(gamma_r * gamma_s * gamma_t));
      }
    }
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionHexa::Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                                            math::DenseDMat<Real> &values)
{
  const Uint nb_pts = coordinates.rows();

  // ((P+1) * (P+1) * (P+1)) is the total number of modes
  values.resize(nb_pts, (P + 1) * (P + 1) * (P + 1));

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {

    // If the polynomial order is equal to 0, we have one (constant) mode
    if (P == 0)
    {
      values(pt, 0) = 0.25 * std::sqrt(2.); // (1/sqrt(2) * 1/sqrt(2) * 1/sqrt(2) )
    }

    else
    {
      // Else we have at least eight linear modes (corner modes)

      const Real xi0 = coordinates(pt, XI0);
      const Real xi1 = coordinates(pt, XI1);
      const Real xi2 = coordinates(pt, XI2);

      // Interior modes

      Uint mode_id = 0;

      for (Uint r = 0; r < (P + 1); ++r)
      {
        const Real gamma_r  = 2. / (2. * r + 1.);
        const Real mode_xi0 = m_JP(r, 0, 0, xi0);

        for (Uint s = 0; s < (P + 1); ++s)
        {
          const Real gamma_s  = 2. / (2. * s + 1.);
          const Real mode_xi1 = m_JP(s, 0, 0, xi1);

          for (Uint t = 0; t < (P + 1); ++t)
          {
            const Real gamma_t  = 2. / (2. * t + 1.);
            const Real mode_xi2 = m_JP(t, 0, 0, xi2);
            values(pt, mode_id++) =
                (mode_xi0 * mode_xi1 * mode_xi2) / (std::sqrt(gamma_r * gamma_s * gamma_t));
          }
        }
      }
    } // else
  }   // Loop over all coordinate points
}

// ----------------------------------------------------------------------------

void ModalExpansionHexa::evaluate_derivatives_in_one_point(const math::DenseDVec<Real> &point,
                                                           const math::DenseDVec<Real> &values,
                                                           math::DenseDMat<Real> &derivatives)
{
  // ((P+1)*(P+1)*(P+1)) is the total number of modes
  derivatives.resize((P + 1) * (P + 1) * (P + 1), 3);
  derivatives.fill(0.0);

  // If the polynomial order is equal to 0, we have one (constant) mode

  if (P == 0)
  {
    derivatives(0, XI0) = 0.0;
    derivatives(0, XI1) = 0.0;
    derivatives(0, XI2) = 0.0;
    return;
  }

  // Else we have at least eight linear modes (corner modes)

  // The normalization factor for each mode is 1/sqrt(gamma(n)) with gamma(n)
  // = 2/(2n+1) Hence gamma(0) = 2, gamma(1) = 2/3 etc. So for psi_0(xi) *
  // psi_0(eta) * psi_0(zeta)  we have factor = 1/sqrt(gamma(0) * gamma(1) *
  // gamma(2)) =  1/(2 * sqrt(2)) = 0.25 * sqrt(2) For psi0 * psi_1 * psi_2
  // the factor is 1/( sqrt(gamma(0)*gamma(1)*gamma(2) ) = 1/sqrt(2*2/3 * 2))
  // = sqrt(3/2)/2

  Uint mode_id = 0;

  for (Uint r = 0; r < (P + 1); ++r)
  {
    const Real gamma_r      = 2. / (2. * r + 1.);
    const Real mode_xi0     = m_JP(r, 0, 0, point[XI0]);
    const Real der_mode_xi0 = m_JP.dx(r, 0, 0, point[XI0]);

    for (Uint s = 0; s < (P + 1); ++s)
    {
      const Real gamma_s      = 2. / (2. * s + 1.);
      const Real mode_xi1     = m_JP(s, 0, 0, point[XI1]);
      const Real der_mode_xi1 = m_JP.dx(s, 0, 0, point[XI1]);

      for (Uint t = 0; t < (P + 1); ++t)
      {
        const Real gamma_t      = 2. / (2. * t + 1.);
        const Real mode_xi2     = m_JP(t, 0, 0, point[XI2]);
        const Real der_mode_xi2 = m_JP.dx(t, 0, 0, point[XI2]);

        const Real factor = 1.0 / (std::sqrt(gamma_r * gamma_s * gamma_t));

        derivatives(mode_id, XI0) = factor * der_mode_xi0 * mode_xi1 * mode_xi2;
        derivatives(mode_id, XI1) = factor * mode_xi0 * der_mode_xi1 * mode_xi2;
        derivatives(mode_id, XI2) = factor * mode_xi0 * mode_xi1 * der_mode_xi2;
        mode_id++;
      } // Loop over order t
    }   // Loop over order s
  }     // Loop over order r
}

// ----------------------------------------------------------------------------

void ModalExpansionHexa::Vandermonde_matrix_derivatives(
    const math::DenseDMat<Real> &coordinates, std::vector<math::DenseDMat<Real>> &derivative_values)
{
  const Uint nb_pts = coordinates.rows();

  // (P+1)*(P+1)*(P+1) is the total number of modes
  derivative_values.resize(3);
  derivative_values[XI0].resize(nb_pts, (P + 1) * (P + 1) * (P + 1));
  derivative_values[XI1].resize(nb_pts, (P + 1) * (P + 1) * (P + 1));
  derivative_values[XI2].resize(nb_pts, (P + 1) * (P + 1) * (P + 1));

  math::DenseDMat<Real> &dVdXi0 = derivative_values[XI0];
  math::DenseDMat<Real> &dVdXi1 = derivative_values[XI1];
  math::DenseDMat<Real> &dVdXi2 = derivative_values[XI2];

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {
    // If the polynomial order is equal to 0, we have one (constant) mode

    if (P == 0)
    {
      dVdXi0(pt, 0) = 0.0;
      dVdXi1(pt, 0) = 0.0;
      dVdXi2(pt, 0) = 0.0;
    }

    else
    {
      // Else we have at least eight linear modes (corner modes)

      const Real xi0 = coordinates(pt, XI0);
      const Real xi1 = coordinates(pt, XI1);
      const Real xi2 = coordinates(pt, XI2);

      // The normalization factor for each mode is 1/sqrt(gamma(n)) with
      // gamma(n) = 2/(2n+1) Hence gamma(0) = 2, gamma(1) = 2/3 etc.

      Uint mode_id = 0;

      for (Uint r = 0; r < (P + 1); ++r)
      {
        const Real gamma_r      = 2. / (2. * r + 1.);
        const Real mode_xi0     = m_JP(r, 0, 0, xi0);
        const Real der_mode_xi0 = m_JP.dx(r, 0, 0, xi0);

        for (Uint s = 0; s < (P + 1); ++s)
        {
          const Real gamma_s      = 2. / (2. * s + 1.);
          const Real mode_xi1     = m_JP(s, 0, 0, xi1);
          const Real der_mode_xi1 = m_JP.dx(s, 0, 0, xi1);

          for (Uint t = 0; t < (P + 1); ++t)
          {
            const Real gamma_t      = 2. / (2. * t + 1);
            const Real mode_xi2     = m_JP(t, 0, 0, xi2);
            const Real der_mode_xi2 = m_JP.dx(t, 0, 0, xi2);

            const Real factor = 1.0 / (std::sqrt(gamma_r * gamma_s * gamma_t));

            dVdXi0(pt, mode_id) = factor * der_mode_xi0 * mode_xi1 * mode_xi2;
            dVdXi1(pt, mode_id) = factor * mode_xi0 * der_mode_xi1 * mode_xi2;
            dVdXi2(pt, mode_id) = factor * mode_xi0 * mode_xi1 * der_mode_xi2;
            mode_id++;

          } // Loop over order t
        }   // Loop over order s
      }     // Loop over order r

    } // else
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionHexa::is_leading_expansion_term(math::DenseDVec<bool> &is_leading_term)
{
  is_leading_term.resize((P + 1) * (P + 1) * (P + 1));
  is_leading_term.fill(false);

  Uint mode_id = 0;

  for (Uint r = 0; r < (P + 1); ++r)
  {
    for (Uint s = 0; s < (P + 1); ++s)
    {
      for (Uint t = 0; t < (P + 1); ++t)
      {
        if ((r == P) || (s == P) || (t == P))
        {
          is_leading_term[mode_id++] = true;
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionHexa::mode_poly_deg(math::DenseDVec<Uint> &poly_deg)
{
  poly_deg.resize((P + 1) * (P + 1) * (P + 1));

  Uint mode_id = 0;

  for (Uint r = 0; r < (P + 1); ++r)
  {
    for (Uint s = 0; s < (P + 1); ++s)
    {
      for (Uint t = 0; t < (P + 1); ++t)
      {
        poly_deg[mode_id++] = r + s + t;
      }
    }
  }
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit
