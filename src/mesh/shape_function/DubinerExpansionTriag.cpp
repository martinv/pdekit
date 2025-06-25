#include "mesh/shape_function/DubinerExpansionTriag.hpp"

#include "math/DenseConstVecView.hpp"
#include "math/polynomials/JacobiPolynomial.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

// ----------------------------------------------------------------------------

DubinerExpansionTriag::DubinerExpansionTriag() : ModalExpansion()
{
}

// ----------------------------------------------------------------------------

DubinerExpansionTriag::DubinerExpansionTriag(const Uint poly_order) : ModalExpansion(poly_order)
{
}

// ----------------------------------------------------------------------------

DubinerExpansionTriag::~DubinerExpansionTriag()
{
}

// ----------------------------------------------------------------------------

/// Return the number of modes
Uint DubinerExpansionTriag::nb_modes() const
{
  return (P + 1) * (P + 2) / 2;
}

// ----------------------------------------------------------------------------

Uint DubinerExpansionTriag::topo_dim() const
{
  return _2D;
}

// ----------------------------------------------------------------------------

void DubinerExpansionTriag::evaluate_in_one_point(const math::DenseDVec<Real> &point,
                                                  math::DenseDVec<Real> &values)
{
#if USE_SINGULARITY_FREE_EVALUATION_ON_TRIAG

  values.resize((P + 1) * (P + 2) / 2); // (P+1)*(P+2)/2 is the total number of modes

  const Real xi0 = point[KSI];
  const Real xi1 = point[ETA];

  values[0]             = 1.;
  values[mode_id(1, 0)] = 0.5 * (1 + 2. * xi0 + xi1);
  for (Uint p = 1; p < P; ++p)
  {
    values[mode_id(p + 1, 0)] =
        0.5 * (2. * p + 1.) / (p + 1.) * (1. + 2. * xi0 + xi1) * values[mode_id(p, 0)] -
        0.25 * p / (p + 1.) * (1. - xi1) * (1. - xi1) * values[mode_id(p - 1, 0)];
  }

  for (Uint p = 0; p < P; ++p)
  {
    values[mode_id(p, 1)] = 0.5 * values[mode_id(p, 0)] * (1. + 2. * p + (3. + 2. * p) * xi1);
  }

  for (Uint p = 0; p < P; ++p)
  {
    for (Uint q = 1; q < (P - p); ++q)
    {
      values[mode_id(p, q + 1)] =
          (a(q, 2. * p + 1, 0) * xi1 + b(q, 2. * p + 1, 0)) * values[mode_id(p, q)] -
          c(q, 2. * p + 1, 0) * values[mode_id(p, q - 1)];
    }
  }

  // This is here for better conditioning of the basis when evaluated
  // in Vandermonde matrix?
  //  for(Uint p = 0; p < (P+1); ++p)
  //  {
  //    for(Uint q = 0; q < (P-p+1); ++q)
  //    {
  //      values[mode_id(p,q)] *= std::sqrt((p+0.5)*(p+q+1.0));
  //    }
  //  }

#else

  math::JacobiPolynomial jp;

  values.resize((P + 1) * (P + 2) / 2); // (P+1)*(P+2)/2 is the total number of modes

  // Variables transformed from triangle to square
  Real eta0 = 0.0;
  Real eta1 = 0.0;

  // const Real sqrt_2 = std::sqrt(2.0);

  for (Uint p = 0; p < (P + 1); ++p)
  {
    for (Uint q = 0; q < (P + 1 - p); ++q)
    {
      const Real xi0 = point[KSI];
      const Real xi1 = point[ETA];

      if (((1.0 - 1e-6) < xi1) && (xi1 < 1.0 + 1e-6))
      {
        eta0 = -1.0;
      }
      else
      {
        eta0 = 2. * (1. + xi0) / (1. - xi1) - 1.;
      }
      eta1 = xi1;

      // values[sf_idx] = sqrt_2 *
      // jp(i,0,0,eta0)*jp(j,2*i+1,0,eta1)*std::pow(1.-eta1,i);
      values[mode_id(p, q)] =
          jp(p, 0, 0, eta0) * jp(q, 2 * p + 1, 0, eta1) * std::pow(0.5 * (1. - eta1), p);
    }
  }

#endif
}

// ----------------------------------------------------------------------------

void DubinerExpansionTriag::Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                                               math::DenseDMat<Real> &values)
{
  const Uint nb_pts = coordinates.rows();

  values.resize(nb_pts, (P + 1) * (P + 2) / 2); // (P+1)*(P+2)/2 is the total number of modes

#if USE_SINGULARITY_FREE_EVALUATION_ON_TRIAG

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {
    const Real xi0 = coordinates(pt, KSI);
    const Real xi1 = coordinates(pt, ETA);

    values(pt, 0)             = 1.;
    values(pt, mode_id(1, 0)) = 0.5 * (1 + 2. * xi0 + xi1);

    for (Uint p = 1; p < P; ++p)
    {
      values(pt, mode_id(p + 1, 0)) =
          0.5 * (2. * p + 1.) / (p + 1.) * (1. + 2. * xi0 + xi1) * values(pt, mode_id(p, 0)) -
          0.25 * p / (p + 1.) * (1. - xi1) * (1. - xi1) * values(pt, mode_id(p - 1, 0));
    }

    for (Uint p = 0; p < P; ++p)
    {
      values(pt, mode_id(p, 1)) =
          0.5 * values(pt, mode_id(p, 0)) * (1. + 2. * p + (3. + 2. * p) * xi1);
    }

    for (Uint p = 0; p < P; ++p)
    {
      for (Uint q = 1; q < (P - p); ++q)
      {
        values(pt, mode_id(p, q + 1)) =
            (a(q, 2. * p + 1, 0) * xi1 + b(q, 2 * p + 1, 0)) * values(pt, mode_id(p, q)) -
            c(q, 2. * p + 1, 0) * values(pt, mode_id(p, q - 1));
      }
    }

    // This is here for better conditioning of the matrix?
    //    for(Uint p = 0; p < (P+1); ++p)
    //    {
    //      for(Uint q = 0; q < (P-p+1); ++q)
    //      {
    //        values(pt,mode_id(p,q)) *= std::sqrt((p+0.5)*(p+q+1.0));
    //      }
    //    }
  }

#else

  math::JacobiPolynomial jp;

  // Variables transformed from triangle to square
  Real eta0 = 0.0;
  Real eta1 = 0.0;

  // const Real sqrt_2 = std::sqrt(2.0);

  for (Uint p = 0; p < (P + 1); ++p)
  {
    for (Uint q = 0; q < (P + 1 - p); ++q)
    {

      for (Uint pt = 0; pt < nb_pts; ++pt)
      {
        const Real xi0 = coordinates(pt, KSI);
        const Real xi1 = coordinates(pt, ETA);

        if (((1.0 - 1e-6) < xi1) && (xi1 < 1.0 + 1e-6))
        {
          eta0 = -1.0;
        }
        else
        {
          eta0 = 2. * (1. + xi0) / (1. - xi1) - 1.;
        }
        eta1 = xi1;

        // values(pt,sf_idx) = sqrt_2 *
        // jp(i,0,0,eta0)*jp(j,2*i+1,0,eta1)*std::pow(1.-eta1,i);
        values(pt, mode_id(p, q)) =
            jp(p, 0, 0, eta0) * jp(q, 2 * p + 1, 0, eta1) * std::pow(0.5 * (1. - eta1), p);

      } // Loop over all pts
    }
  }

#endif
}

// ----------------------------------------------------------------------------

void DubinerExpansionTriag::evaluate_derivatives_in_one_point(const math::DenseDVec<Real> &point,
                                                              const math::DenseDVec<Real> &values,
                                                              math::DenseDMat<Real> &derivatives)
{
#if USE_SINGULARITY_FREE_EVALUATION_ON_TRIAG

  derivatives.resize((P + 1) * (P + 2) / 2,
                     2); // (P+1)*(P+2)/2 is the total number of modes
  derivatives.fill(0.0);

  const Real xi0 = point[KSI];
  const Real xi1 = point[ETA];

  derivatives(0, KSI) = 0.;
  derivatives(0, ETA) = 0.;

  derivatives(mode_id(1, 0), KSI) = 1.;
  derivatives(mode_id(1, 0), ETA) = 0.5;

  for (Uint p = 1; p < P; ++p)
  {

    derivatives(mode_id(p + 1, 0), KSI) =
        (2. * p + 1.) / (p + 1.) *
            (values[mode_id(p, 0)] +
             0.5 * (1. + 2. * xi0 + xi1) * derivatives(mode_id(p, 0), KSI)) -
        0.25 * p / (p + 1.) * (1. - xi1) * (1. - xi1) * derivatives(mode_id(p - 1, 0), KSI);

    derivatives(mode_id(p + 1, 0), ETA) =
        (2. * p + 1.) / (p + 1.) *
            (0.5 * values[mode_id(p, 0)] +
             0.5 * (1. + 2. * xi0 + xi1) * derivatives(mode_id(p, 0), ETA)) -
        p / (p + 1.) *
            (-0.5 * (1. - xi1) * values[mode_id(p - 1, 0)] +
             0.25 * (1. - xi1) * (1. - xi1) * derivatives(mode_id(p - 1, 0), ETA));
  }

  for (Uint p = 0; p < P; ++p)
  {
    derivatives(mode_id(p, 1), KSI) =
        0.5 * derivatives(mode_id(p, 0), KSI) * (1. + 2. * p + (3. + 2. * p) * xi1);

    derivatives(mode_id(p, 1), ETA) =
        0.5 * derivatives(mode_id(p, 0), ETA) * (1. + 2. * p + (3. + 2. * p) * xi1) +
        0.5 * (3. + 2. * p) * values[mode_id(p, 0)];
  }

  for (Uint p = 0; p < P; ++p)
  {
    for (Uint q = 1; q < (P - p); ++q)
    {
      derivatives(mode_id(p, q + 1), KSI) =
          (a(q, 2. * p + 1, 0) * xi1 + b(q, 2. * p + 1, 0)) * derivatives(mode_id(p, q), KSI) -
          c(q, 2. * p + 1, 0) * derivatives(mode_id(p, q - 1), KSI);

      derivatives(mode_id(p, q + 1), ETA) =
          a(q, 2. * p + 1, 0) * values[mode_id(p, q)] +
          (a(q, 2. * p + 1, 0) * xi1 + b(q, 2. * p + 1, 0)) * derivatives(mode_id(p, q), ETA) -
          c(q, 2. * p + 1, 0) * derivatives(mode_id(p, q - 1), ETA);
    }
  }

  // This is here for better conditioning of the basis when evaluated
  // in Vandermonde matrix?
  //  for(Uint p = 0; p < (P+1); ++p)
  //  {
  //    for(Uint q = 0; q < (P-p+1); ++q)
  //    {
  //      derivatives(mode_id(p,q),KSI) *= std::sqrt((p+0.5)*(p+q+1.0));
  //      derivatives(mode_id(p,q),ETA) *= std::sqrt((p+0.5)*(p+q+1.0));
  //    }
  //  }

#else

  math::JacobiPolynomial jp;

  derivatives.resize((P + 1) * (P + 2) / 2,
                     2); // (P+1)*(P+2)/2 is the total number of modes

  derivatives.fill(0.0);

  const Real xi0 = point[KSI];
  const Real xi1 = point[ETA];

  // Variables transformed from triangle to square
  Real eta0 = 0.0;
  Real eta1 = 0.0;

  if (((1.0 - 1e-6) < xi1) && (xi1 < 1.0 + 1e-6))
  {
    eta0 = -1.0;
  }
  else
  {
    eta0 = 2. * (1. + xi0) / (1. - xi1) - 1.;
  }
  eta1 = xi1;

  for (Uint p = 0; p < (P + 1); ++p)
  {
    for (Uint q = 0; q < (P + 1 - p); ++q)
    {
      const Real P_p  = jp(p, 0, 0, eta0);
      const Real dP_p = jp.dx(p, 0, 0, eta0);
      const Real P_q  = jp(q, 2 * p + 1, 0, eta1);
      const Real dP_q = jp.dx(q, 2 * p + 1, 0, eta1);

      derivatives(mode_id(p, q), KSI) = dP_p * P_q;

      if (p != 0)
      {
        derivatives(mode_id(p, q), KSI) *= std::pow(0.5 * (1. - eta1), p - 1);
        derivatives(mode_id(p, q), ETA) =
            dP_p * 2. * (1. + xi0) / ((1. - xi1) * (1. - xi1)) * std::pow(0.5 * (1. - eta1), p) *
                P_q +
            P_p * (-0.5 * p * std::pow(0.5 * (1. - eta1), p - 1) * P_q +
                   std::pow(0.5 * (1. - eta1), p) * dP_q);
      }
      else
      {
        derivatives(mode_id(p, q), KSI) *= 2. / (1. - eta1);
        derivatives(mode_id(p, q), ETA) =
            dP_p * 2. * (1. + xi0) / ((1. - xi1) * (1. - xi1)) * P_q + P_p * dP_q;
      }
    }
  }

#endif
}

// ----------------------------------------------------------------------------

void DubinerExpansionTriag::Vandermonde_matrix_derivatives(
    const math::DenseDMat<Real> &coordinates, std::vector<math::DenseDMat<Real>> &derivative_values)
{
  math::DenseDVec<Real> point(_2D);
  math::DenseDMat<Real> derivatives_at_pt((P + 1) * (P + 2) / 2,
                                          _2D); // (P+1)*(P+2)/2 is the total number of modes
  math::DenseDVec<Real> values;

  const Uint nb_pts = coordinates.rows();

  derivative_values.resize(_2D);
  derivative_values[KSI].resize(nb_pts, (P + 1) * (P + 2) / 2);
  derivative_values[ETA].resize(nb_pts, (P + 1) * (P + 2) / 2);

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {

    const Real xi0 = coordinates(pt, KSI);
    const Real xi1 = coordinates(pt, ETA);

    point[KSI] = xi0;
    point[ETA] = xi1;

    derivatives_at_pt.fill(0.0);

    evaluate_in_one_point(point, values);

    derivatives_at_pt(0, KSI) = 0.;
    derivatives_at_pt(0, ETA) = 0.;

    derivatives_at_pt(mode_id(1, 0), KSI) = 1.;
    derivatives_at_pt(mode_id(1, 0), ETA) = 0.5;

    for (Uint p = 1; p < P; ++p)
    {

      derivatives_at_pt(mode_id(p + 1, 0), KSI) =
          (2. * p + 1.) / (p + 1.) *
              (values[mode_id(p, 0)] +
               0.5 * (1. + 2. * xi0 + xi1) * derivatives_at_pt(mode_id(p, 0), KSI)) -
          0.25 * p / (p + 1.) * (1. - xi1) * (1. - xi1) * derivatives_at_pt(mode_id(p - 1, 0), KSI);

      derivatives_at_pt(mode_id(p + 1, 0), ETA) =
          (2. * p + 1.) / (p + 1.) *
              (0.5 * values[mode_id(p, 0)] +
               0.5 * (1. + 2. * xi0 + xi1) * derivatives_at_pt(mode_id(p, 0), ETA)) -
          p / (p + 1.) *
              (-0.5 * (1. - xi1) * values[mode_id(p - 1, 0)] +
               0.25 * (1. - xi1) * (1. - xi1) * derivatives_at_pt(mode_id(p - 1, 0), ETA));
    }

    for (Uint p = 0; p < P; ++p)
    {
      derivatives_at_pt(mode_id(p, 1), KSI) =
          0.5 * derivatives_at_pt(mode_id(p, 0), KSI) * (1. + 2. * p + (3. + 2. * p) * xi1);

      derivatives_at_pt(mode_id(p, 1), ETA) =
          0.5 * derivatives_at_pt(mode_id(p, 0), ETA) * (1. + 2. * p + (3. + 2. * p) * xi1) +
          0.5 * (3. + 2. * p) * values[mode_id(p, 0)];
    }

    for (Uint p = 0; p < P; ++p)
    {
      for (Uint q = 1; q < (P - p); ++q)
      {
        derivatives_at_pt(mode_id(p, q + 1), KSI) =
            (a(q, 2. * p + 1, 0) * xi1 + b(q, 2. * p + 1, 0)) *
                derivatives_at_pt(mode_id(p, q), KSI) -
            c(q, 2. * p + 1, 0) * derivatives_at_pt(mode_id(p, q - 1), KSI);

        derivatives_at_pt(mode_id(p, q + 1), ETA) =
            a(q, 2. * p + 1, 0) * values[mode_id(p, q)] +
            (a(q, 2. * p + 1, 0) * xi1 + b(q, 2. * p + 1, 0)) *
                derivatives_at_pt(mode_id(p, q), ETA) -
            c(q, 2. * p + 1, 0) * derivatives_at_pt(mode_id(p, q - 1), ETA);
      }
    }

    // Insert the values in derivative matrices
    derivative_values[KSI].insert_row(pt, derivatives_at_pt.const_col(KSI));
    derivative_values[ETA].insert_row(pt, derivatives_at_pt.const_col(ETA));

  } // Loop over all points
}

// ----------------------------------------------------------------------------

void DubinerExpansionTriag::is_leading_expansion_term(math::DenseDVec<bool> &is_leading_term)
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

  is_leading_term[mode_id(P, 0)]     = true;
  is_leading_term[mode_id(P - 1, 1)] = true;

  for (Uint p = 0; p < P; ++p)
  {
    for (Uint q = 1; q < (P - p); ++q)
    {
      if ((p + q + 1) == P)
      {
        // Terms (p,q+1)
        is_leading_term[mode_id(p, q + 1)] = true;
      }
    }
  }
}

// ----------------------------------------------------------------------------

void DubinerExpansionTriag::mode_poly_deg(math::DenseDVec<Uint> &poly_deg)
{
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit
