#include "mesh/shape_function/ModalExpansionLine.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

// ----------------------------------------------------------------------------

ModalExpansionLine::ModalExpansionLine() : ModalExpansion()
{
}

// ----------------------------------------------------------------------------

ModalExpansionLine::ModalExpansionLine(const Uint poly_order) : ModalExpansion(poly_order)
{
}

// ----------------------------------------------------------------------------

ModalExpansionLine::~ModalExpansionLine()
{
}

// ----------------------------------------------------------------------------

/// Return the number of modes
Uint ModalExpansionLine::nb_modes() const
{
  return P + 1;
}

// ----------------------------------------------------------------------------

Uint ModalExpansionLine::topo_dim() const
{
  return _1D;
}

// ----------------------------------------------------------------------------

void ModalExpansionLine::evaluate_in_one_point(const math::DenseDVec<Real> &point,
                                               math::DenseDVec<Real> &values)
{
  values.resize(P + 1); // (P+1) is the total number of modes

  for (Uint n = 0; n < (P + 1); ++n)
  {
    const Real gamma = 2. / (2. * n + 1.);
    const Real jp    = m_JP(n, 0, 0, point[KSI]);
    values[n]        = jp / std::sqrt(gamma);
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionLine::Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                                            math::DenseDMat<Real> &values)
{
  const Uint nb_pts = coordinates.rows();

  values.resize(nb_pts, P + 1); // (P+1) is the total number of modes

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {
    const Real xi = coordinates(pt, KSI);

    for (Uint n = 0; n < (P + 1); ++n)
    {
      const Real gamma = 2. / (2. * n + 1.);
      const Real jp    = m_JP(n, 0, 0, xi);
      values(pt, n)    = jp / std::sqrt(gamma);
    }
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionLine::evaluate_derivatives_in_one_point(const math::DenseDVec<Real> &point,
                                                           const math::DenseDVec<Real> &values,
                                                           math::DenseDMat<Real> &derivatives)
{

  derivatives.resize(P + 1, 1); // (P+1) is the total number of modes
  derivatives.fill(0.0);

  for (Uint n = 0; n < (P + 1); ++n)
  {
    const Real gamma    = 2. / (2. * n + 1.);
    const Real jp_deriv = m_JP.dx(n, 0, 0, point[KSI]);
    derivatives(n, 0)   = jp_deriv / std::sqrt(gamma);
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionLine::Vandermonde_matrix_derivatives(
    const math::DenseDMat<Real> &coordinates, std::vector<math::DenseDMat<Real>> &derivative_values)
{
  const Uint nb_pts = coordinates.rows();

  derivative_values.resize(1);
  derivative_values[KSI].resize(nb_pts,
                                P + 1); // (P+1) is the total number of modes

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {
    const Real xi = coordinates(pt, KSI);

    for (Uint n = 0; n < (P + 1); ++n)
    {
      const Real gamma              = 2. / (2. * n + 1.);
      const Real jp_deriv           = m_JP.dx(n, 0, 0, xi);
      derivative_values[KSI](pt, n) = jp_deriv / std::sqrt(gamma);
    }
  }
}

// ----------------------------------------------------------------------------

void ModalExpansionLine::is_leading_expansion_term(math::DenseDVec<bool> &is_leading_term)
{
  is_leading_term.resize(P + 1);
  is_leading_term.fill(false);
  is_leading_term[P] = true;
}

// ----------------------------------------------------------------------------

void ModalExpansionLine::mode_poly_deg(math::DenseDVec<Uint> &poly_deg)
{
  poly_deg.resize(P + 1);
  for (Uint i = 0; i < (P + 1); ++i)
  {
    poly_deg[i] = i;
  }
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit
