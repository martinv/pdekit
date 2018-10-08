#include "mesh/shape_function/DubinerExpansionTetra.hpp"

#include "math/DenseConstVecView.hpp"
#include "math/polynomials/JacobiPolynomial.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

// ----------------------------------------------------------------------------

DubinerExpansionTetra::DubinerExpansionTetra() : ModalExpansion()
{
}

// ----------------------------------------------------------------------------

DubinerExpansionTetra::DubinerExpansionTetra(const Uint poly_order) : ModalExpansion(poly_order)
{
}

// ----------------------------------------------------------------------------

DubinerExpansionTetra::~DubinerExpansionTetra()
{
}

// ----------------------------------------------------------------------------

Uint DubinerExpansionTetra::nb_modes() const
{
  return (P + 1) * (P + 2) * (P + 3) / 6;
}

// ----------------------------------------------------------------------------

Uint DubinerExpansionTetra::topo_dim() const
{
  return _3D;
}

// ----------------------------------------------------------------------------

void DubinerExpansionTetra::evaluate_in_one_point(const math::DenseDVec<Real> &point,
                                                  math::DenseDVec<Real> &values)
{
#if USE_SINGULARITY_FREE_EVALUATION_ON_TETRA

  values.resize((P + 1) * (P + 2) * (P + 3) / 6); // (P+1)*(P+2)*(P+3)/6 is
                                                  // the total number of modes

  const Real xi0 = point[XI0];
  const Real xi1 = point[XI1];
  const Real xi2 = point[XI2];

  const Real F1 = 0.5 * (2. + 2. * xi0 + xi1 + xi2);
  const Real F2 = 0.25 * (xi1 + xi2) * (xi1 + xi2);
  const Real F3 = 0.5 * (2. + 3. * xi1 + xi2);
  const Real F4 = 0.5 * (1. + 2. * xi1 + xi2);
  const Real F5 = 0.5 * (1. - xi2);

  values[0]                = 1.;
  values[mode_id(1, 0, 0)] = F1;

  for (Uint p = 1; p < P; ++p)
  {
    values[mode_id(p + 1, 0, 0)] = (2. * p + 1.) / (p + 1.) * F1 * values[mode_id(p, 0, 0)] -
                                   p / (p + 1.) * F2 * values[mode_id(p - 1, 0, 0)];
  }

  for (Uint p = 0; p < P; ++p)
  {
    values[mode_id(p, 1, 0)] = (p * (1. + xi1) + F3) * values[mode_id(p, 0, 0)];
  }

  for (Uint p = 0; p < (P - 1); ++p)
  {
    for (Uint q = 1; q < (P - p); ++q)
    {
      values[mode_id(p, q + 1, 0)] =
          (a(q, 2. * p + 1., 0) * F4 + b(q, 2. * p + 1., 0) * F5) * values[mode_id(p, q, 0)] -
          c(q, 2. * p + 1., 0) * F5 * F5 * values[mode_id(p, q - 1, 0)];
    }
  }

  for (Uint p = 0; p < P; ++p)
  {
    for (Uint q = 0; q < (P - p); ++q)
    {
      values[mode_id(p, q, 1)] = (1. + p + q + (2. + q + p) * xi2) * values[mode_id(p, q, 0)];
    }
  }

  for (Uint p = 0; p < (P - 1); ++p)
  {
    for (Uint q = 0; q < (P - p - 1); ++q)
    {
      for (Uint r = 1; r < (P - p - q); ++r)
      {
        values[mode_id(p, q, r + 1)] =
            (a(r, 2 * p + 2 * q + 2, 0) * xi2 + b(r, 2 * p + 2 * q + 2, 0)) *
                values[mode_id(p, q, r)] -
            c(r, 2 * p + 2 * q + 2, 0) * values[mode_id(p, q, r - 1)];
      }
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

  values.resize((P + 1) * (P + 2) * (P + 3) / 6); // (P+1)*(P+2)*(P+3)/6 is
                                                  // the total number of modes

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

void DubinerExpansionTetra::Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                                               math::DenseDMat<Real> &values)
{
  const Uint nb_pts = coordinates.rows();

  values.resize(nb_pts,
                (P + 1) * (P + 2) * (P + 3) / 6); // (P+1)*(P+2)*(P+3)/6
                                                  // is the total number
                                                  // of modes

#if USE_SINGULARITY_FREE_EVALUATION_ON_TETRA

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {
    const Real xi0 = coordinates(pt, XI0);
    const Real xi1 = coordinates(pt, XI1);
    const Real xi2 = coordinates(pt, XI2);

    const Real F1 = 0.5 * (2. + 2. * xi0 + xi1 + xi2);
    const Real F2 = 0.25 * (xi1 + xi2) * (xi1 + xi2);
    const Real F3 = 0.5 * (2. + 3. * xi1 + xi2);
    const Real F4 = 0.5 * (1. + 2. * xi1 + xi2);
    const Real F5 = 0.5 * (1. - xi2);

    values(pt, 0)                = 1.;
    values(pt, mode_id(1, 0, 0)) = F1;

    for (Uint p = 1; p < P; ++p)
    {
      values(pt, mode_id(p + 1, 0, 0)) =
          (2. * p + 1.) / (p + 1.) * F1 * values(pt, mode_id(p, 0, 0)) -
          p / (p + 1.) * F2 * values(pt, mode_id(p - 1, 0, 0));
    }

    for (Uint p = 0; p < P; ++p)
    {
      values(pt, mode_id(p, 1, 0)) = (p * (1. + xi1) + F3) * values(pt, mode_id(p, 0, 0));
    }

    for (Uint p = 0; p < (P - 1); ++p)
    {
      for (Uint q = 1; q < (P - p); ++q)
      {
        values(pt, mode_id(p, q + 1, 0)) =
            (a(q, 2. * p + 1., 0) * F4 + b(q, 2. * p + 1., 0) * F5) * values(pt, mode_id(p, q, 0)) -
            c(q, 2. * p + 1., 0) * F5 * F5 * values(pt, mode_id(p, q - 1, 0));
      }
    }

    for (Uint p = 0; p < P; ++p)
    {
      for (Uint q = 0; q < (P - p); ++q)
      {
        values(pt, mode_id(p, q, 1)) =
            (1. + p + q + (2. + q + p) * xi2) * values(pt, mode_id(p, q, 0));
      }
    }

    for (Uint p = 0; p < (P - 1); ++p)
    {
      for (Uint q = 0; q < (P - p - 1); ++q)
      {
        for (Uint r = 1; r < (P - p - q); ++r)
        {
          values(pt, mode_id(p, q, r + 1)) =
              (a(r, 2 * p + 2 * q + 2, 0) * xi2 + b(r, 2 * p + 2 * q + 2, 0)) *
                  values(pt, mode_id(p, q, r)) -
              c(r, 2 * p + 2 * q + 2, 0) * values(pt, mode_id(p, q, r - 1));
        }
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

void DubinerExpansionTetra::evaluate_derivatives_in_one_point(const math::DenseDVec<Real> &point,
                                                              const math::DenseDVec<Real> &values,
                                                              math::DenseDMat<Real> &derivatives)
{
#if USE_SINGULARITY_FREE_EVALUATION_ON_TETRA

  // (P+1)*(P+2)*(P+3)/6 is the total number of modes
  derivatives.resize((P + 1) * (P + 2) * (P + 3) / 6, 3);
  derivatives.fill(0.0);

  const Real xi0 = point[XI0];
  const Real xi1 = point[XI1];
  const Real xi2 = point[XI2];

  std::array<Real, 5> F;
  std::array<Real, 5> dFdxi0;
  std::array<Real, 5> dFdxi1;
  std::array<Real, 5> dFdxi2;

  F[0] = 0.5 * (2. + 2. * xi0 + xi1 + xi2);
  F[1] = 0.25 * (xi1 + xi2) * (xi1 + xi2);
  F[2] = 0.5 * (2. + 3. * xi1 + xi2);
  F[3] = 0.5 * (1. + 2. * xi1 + xi2);
  F[4] = 0.5 * (1. - xi2);

  dFdxi0[0] = 1.0;
  dFdxi0[1] = 0.0;
  dFdxi0[2] = 0.0;
  dFdxi0[3] = 0.0;
  dFdxi0[4] = 0.0;

  dFdxi1[0] = 0.5;
  dFdxi1[1] = 0.5 * (xi1 + xi2);
  dFdxi1[2] = 1.5;
  dFdxi1[3] = 1.0;
  dFdxi1[4] = 0.0;

  dFdxi2[0] = 0.5;
  dFdxi2[1] = 0.5 * (xi1 + xi2);
  dFdxi2[2] = 0.5;
  dFdxi2[3] = 0.5;
  dFdxi2[4] = -0.5;

  derivatives(0, XI0) = 0.;
  derivatives(0, XI1) = 0.;
  derivatives(0, XI2) = 0.;

  derivatives(mode_id(1, 0, 0), XI0) = dFdxi0[0];
  derivatives(mode_id(1, 0, 0), XI1) = dFdxi1[0];
  derivatives(mode_id(1, 0, 0), XI2) = dFdxi2[0];

  for (Uint p = 1; p < P; ++p)
  {

    derivatives(mode_id(p + 1, 0, 0), XI0) =
        (2. * p + 1) / (p + 1.) *
            (dFdxi0[0] * values[mode_id(p, 0, 0)] + F[0] * derivatives(mode_id(p, 0, 0), XI0)) -
        p / (p + 1.) *
            (dFdxi0[1] * values[mode_id(p - 1, 0, 0)] +
             F[1] * derivatives(mode_id(p - 1, 0, 0), XI0));

    derivatives(mode_id(p + 1, 0, 0), XI1) =
        (2. * p + 1) / (p + 1.) *
            (dFdxi1[0] * values[mode_id(p, 0, 0)] + F[0] * derivatives(mode_id(p, 0, 0), XI1)) -
        p / (p + 1.) *
            (dFdxi1[1] * values[mode_id(p - 1, 0, 0)] +
             F[1] * derivatives(mode_id(p - 1, 0, 0), XI1));

    derivatives(mode_id(p + 1, 0, 0), XI2) =
        (2. * p + 1) / (p + 1.) *
            (dFdxi2[0] * values[mode_id(p, 0, 0)] + F[0] * derivatives(mode_id(p, 0, 0), XI2)) -
        p / (p + 1.) *
            (dFdxi2[1] * values[mode_id(p - 1, 0, 0)] +
             F[1] * derivatives(mode_id(p - 1, 0, 0), XI2));
  }

  for (Uint p = 0; p < P; ++p)
  {
    derivatives(mode_id(p, 1, 0), XI0) =
        (p * (1. + xi1) + F[2]) * derivatives(mode_id(p, 0, 0), XI0);

    derivatives(mode_id(p, 1, 0), XI1) =
        (p + dFdxi1[2]) * values[mode_id(p, 0, 0)] +
        (p * (1. + xi1) + F[2]) * derivatives(mode_id(p, 0, 0), XI1);

    derivatives(mode_id(p, 1, 0), XI2) =
        dFdxi2[2] * values[mode_id(p, 0, 0)] +
        (p * (1. + xi1) + F[2]) * derivatives(mode_id(p, 0, 0), XI2);
  }

  for (Uint p = 0; p < (P - 1); ++p)
  {
    for (Uint q = 1; q < (P - p); ++q)
    {
      derivatives(mode_id(p, q + 1, 0), XI0) =
          (a(q, 2 * p + 1, 0) * F[3] + b(q, 2 * p + 1, 0) * F[4]) *
              derivatives(mode_id(p, q, 0), XI0) -
          c(q, 2 * p + 1, 0) * F[4] * F[4] * derivatives(mode_id(p, q - 1, 0), XI0);

      derivatives(mode_id(p, q + 1, 0), XI1) =
          a(q, 2 * p + 1, 0) * dFdxi1[3] * values[mode_id(p, q, 0)] +
          (a(q, 2 * p + 1, 0) * F[3] + b(q, 2 * p + 1, 0) * F[4]) *
              derivatives(mode_id(p, q, 0), XI1) -
          c(q, 2 * p + 1, 0) * F[4] * F[4] * derivatives(mode_id(p, q - 1, 0), XI1);

      derivatives(mode_id(p, q + 1, 0), XI2) =
          (a(q, 2 * p + 1, 0) * dFdxi2[3] + b(q, 2 * p + 1, 0) * dFdxi2[4]) *
              values[mode_id(p, q, 0)] +
          (a(q, 2 * p + 1, 0) * F[3] + b(q, 2 * p + 1, 0) * F[4]) *
              derivatives(mode_id(p, q, 0), XI2) -
          c(q, 2 * p + 1, 0) * (2. * F[4] * dFdxi2[4] * values[mode_id(p, q - 1, 0)] +
                                F[4] * F[4] * derivatives(mode_id(p, q - 1, 0), XI2));
    }
  }

  for (Uint p = 0; p < P; ++p)
  {
    for (Uint q = 0; q < (P - p); ++q)
    {
      derivatives(mode_id(p, q, 1), XI0) =
          (1. + p + q + (2. + q + p) * xi2) * derivatives(mode_id(p, q, 0), XI0);

      derivatives(mode_id(p, q, 1), XI1) =
          (1. + p + q + (2. + q + p) * xi2) * derivatives(mode_id(p, q, 0), XI1);

      derivatives(mode_id(p, q, 1), XI2) =
          (1. + p + q + (2. + q + p) * xi2) * derivatives(mode_id(p, q, 0), XI2) +
          (2. + q + p) * values[mode_id(p, q, 0)];
    }
  }

  for (Uint p = 0; p < (P - 1); ++p)
  {
    for (Uint q = 0; q < (P - p - 1); ++q)
    {
      for (Uint r = 1; r < (P - p - q); ++r)
      {
        derivatives(mode_id(p, q, r + 1), XI0) =
            (a(r, 2 * p + 2 * q + 2, 0) * xi2 + b(r, 2 * p + 2 * q + 2, 0)) *
                derivatives(mode_id(p, q, r), XI0) -
            c(r, 2 * p + 2 * q + 2, 0) * derivatives(mode_id(p, q, r - 1), XI0);

        derivatives(mode_id(p, q, r + 1), XI1) =
            (a(r, 2 * p + 2 * q + 2, 0) * xi2 + b(r, 2 * p + 2 * q + 2, 0)) *
                derivatives(mode_id(p, q, r), XI1) -
            c(r, 2 * p + 2 * q + 2, 0) * derivatives(mode_id(p, q, r - 1), XI1);

        derivatives(mode_id(p, q, r + 1), XI2) =
            (a(r, 2 * p + 2 * q + 2, 0) * xi2 + b(r, 2 * p + 2 * q + 2, 0)) *
                derivatives(mode_id(p, q, r), XI2) +
            a(r, 2 * p + 2 * q + 2, 0) * values[mode_id(p, q, r)] -
            c(r, 2 * p + 2 * q + 2, 0) * derivatives(mode_id(p, q, r - 1), XI2);
      }
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

  derivatives.resize((P + 1) * (P + 2) * (P + 3) / 6,
                     2); // (P+1)*(P+2)*(P+3)/6 is the total number of modes

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

void DubinerExpansionTetra::Vandermonde_matrix_derivatives(
    const math::DenseDMat<Real> &coordinates, std::vector<math::DenseDMat<Real>> &derivative_values)
{
  math::DenseDVec<Real> point(_3D);
  math::DenseDMat<Real> derivatives_at_pt((P + 1) * (P + 2) * (P + 3) / 6,
                                          _3D); // (P+1)*(P+2)*(P+3)/6 is the total number of modes
  math::DenseDVec<Real> values;

  std::array<Real, 5> F;
  std::array<Real, 5> dFdxi0;
  std::array<Real, 5> dFdxi1;
  std::array<Real, 5> dFdxi2;

  const Uint nb_pts = coordinates.rows();

  derivative_values.resize(_3D);
  derivative_values[KSI].resize(nb_pts, (P + 1) * (P + 2) * (P + 3) / 6);
  derivative_values[ETA].resize(nb_pts, (P + 1) * (P + 2) * (P + 3) / 6);
  derivative_values[ZTA].resize(nb_pts, (P + 1) * (P + 2) * (P + 3) / 6);

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {

    const Real xi0 = coordinates(pt, KSI);
    const Real xi1 = coordinates(pt, ETA);
    const Real xi2 = coordinates(pt, ZTA);

    point[KSI] = xi0;
    point[ETA] = xi1;
    point[ZTA] = xi2;

    derivatives_at_pt.fill(0.0);

    evaluate_in_one_point(point, values);

    F[0] = 0.5 * (2. + 2. * xi0 + xi1 + xi2);
    F[1] = 0.25 * (xi1 + xi2) * (xi1 + xi2);
    F[2] = 0.5 * (2. + 3. * xi1 + xi2);
    F[3] = 0.5 * (1. + 2. * xi1 + xi2);
    F[4] = 0.5 * (1. - xi2);

    dFdxi0[0] = 1.0;
    dFdxi0[1] = 0.0;
    dFdxi0[2] = 0.0;
    dFdxi0[3] = 0.0;
    dFdxi0[4] = 0.0;

    dFdxi1[0] = 0.5;
    dFdxi1[1] = 0.5 * (xi1 + xi2);
    dFdxi1[2] = 1.5;
    dFdxi1[3] = 1.0;
    dFdxi1[4] = 0.0;

    dFdxi2[0] = 0.5;
    dFdxi2[1] = 0.5 * (xi1 + xi2);
    dFdxi2[2] = 0.5;
    dFdxi2[3] = 0.5;
    dFdxi2[4] = -0.5;

    derivatives_at_pt(0, XI0) = 0.;
    derivatives_at_pt(0, XI1) = 0.;
    derivatives_at_pt(0, XI2) = 0.;

    derivatives_at_pt(mode_id(1, 0, 0), XI0) = dFdxi0[0];
    derivatives_at_pt(mode_id(1, 0, 0), XI1) = dFdxi1[0];
    derivatives_at_pt(mode_id(1, 0, 0), XI2) = dFdxi2[0];

    for (Uint p = 1; p < P; ++p)
    {

      derivatives_at_pt(mode_id(p + 1, 0, 0), XI0) =
          (2. * p + 1) / (p + 1.) *
              (dFdxi0[0] * values[mode_id(p, 0, 0)] +
               F[0] * derivatives_at_pt(mode_id(p, 0, 0), XI0)) -
          p / (p + 1.) *
              (dFdxi0[1] * values[mode_id(p - 1, 0, 0)] +
               F[1] * derivatives_at_pt(mode_id(p - 1, 0, 0), XI0));

      derivatives_at_pt(mode_id(p + 1, 0, 0), XI1) =
          (2. * p + 1) / (p + 1.) *
              (dFdxi1[0] * values[mode_id(p, 0, 0)] +
               F[0] * derivatives_at_pt(mode_id(p, 0, 0), XI1)) -
          p / (p + 1.) *
              (dFdxi1[1] * values[mode_id(p - 1, 0, 0)] +
               F[1] * derivatives_at_pt(mode_id(p - 1, 0, 0), XI1));

      derivatives_at_pt(mode_id(p + 1, 0, 0), XI2) =
          (2. * p + 1) / (p + 1.) *
              (dFdxi2[0] * values[mode_id(p, 0, 0)] +
               F[0] * derivatives_at_pt(mode_id(p, 0, 0), XI2)) -
          p / (p + 1.) *
              (dFdxi2[1] * values[mode_id(p - 1, 0, 0)] +
               F[1] * derivatives_at_pt(mode_id(p - 1, 0, 0), XI2));
    }

    for (Uint p = 0; p < P; ++p)
    {
      derivatives_at_pt(mode_id(p, 1, 0), XI0) =
          (p * (1. + xi1) + F[2]) * derivatives_at_pt(mode_id(p, 0, 0), XI0);

      derivatives_at_pt(mode_id(p, 1, 0), XI1) =
          (p + dFdxi1[2]) * values[mode_id(p, 0, 0)] +
          (p * (1. + xi1) + F[2]) * derivatives_at_pt(mode_id(p, 0, 0), XI1);

      derivatives_at_pt(mode_id(p, 1, 0), XI2) =
          dFdxi2[2] * values[mode_id(p, 0, 0)] +
          (p * (1. + xi1) + F[2]) * derivatives_at_pt(mode_id(p, 0, 0), XI2);
    }

    for (Uint p = 0; p < (P - 1); ++p)
    {
      for (Uint q = 1; q < (P - p); ++q)
      {
        derivatives_at_pt(mode_id(p, q + 1, 0), XI0) =
            (a(q, 2 * p + 1, 0) * F[3] + b(q, 2 * p + 1, 0) * F[4]) *
                derivatives_at_pt(mode_id(p, q, 0), XI0) -
            c(q, 2 * p + 1, 0) * F[4] * F[4] * derivatives_at_pt(mode_id(p, q - 1, 0), XI0);

        derivatives_at_pt(mode_id(p, q + 1, 0), XI1) =
            a(q, 2 * p + 1, 0) * dFdxi1[3] * values[mode_id(p, q, 0)] +
            (a(q, 2 * p + 1, 0) * F[3] + b(q, 2 * p + 1, 0) * F[4]) *
                derivatives_at_pt(mode_id(p, q, 0), XI1) -
            c(q, 2 * p + 1, 0) * F[4] * F[4] * derivatives_at_pt(mode_id(p, q - 1, 0), XI1);

        derivatives_at_pt(mode_id(p, q + 1, 0), XI2) =
            (a(q, 2 * p + 1, 0) * dFdxi2[3] + b(q, 2 * p + 1, 0) * dFdxi2[4]) *
                values[mode_id(p, q, 0)] +
            (a(q, 2 * p + 1, 0) * F[3] + b(q, 2 * p + 1, 0) * F[4]) *
                derivatives_at_pt(mode_id(p, q, 0), XI2) -
            c(q, 2 * p + 1, 0) * (2. * F[4] * dFdxi2[4] * values[mode_id(p, q - 1, 0)] +
                                  F[4] * F[4] * derivatives_at_pt(mode_id(p, q - 1, 0), XI2));
      }
    }

    for (Uint p = 0; p < P; ++p)
    {
      for (Uint q = 0; q < (P - p); ++q)
      {
        derivatives_at_pt(mode_id(p, q, 1), XI0) =
            (1. + p + q + (2. + q + p) * xi2) * derivatives_at_pt(mode_id(p, q, 0), XI0);

        derivatives_at_pt(mode_id(p, q, 1), XI1) =
            (1. + p + q + (2. + q + p) * xi2) * derivatives_at_pt(mode_id(p, q, 0), XI1);

        derivatives_at_pt(mode_id(p, q, 1), XI2) =
            (1. + p + q + (2. + q + p) * xi2) * derivatives_at_pt(mode_id(p, q, 0), XI2) +
            (2. + q + p) * values[mode_id(p, q, 0)];
      }
    }

    for (Uint p = 0; p < (P - 1); ++p)
    {
      for (Uint q = 0; q < (P - p - 1); ++q)
      {
        for (Uint r = 1; r < (P - p - q); ++r)
        {
          derivatives_at_pt(mode_id(p, q, r + 1), XI0) =
              (a(r, 2 * p + 2 * q + 2, 0) * xi2 + b(r, 2 * p + 2 * q + 2, 0)) *
                  derivatives_at_pt(mode_id(p, q, r), XI0) -
              c(r, 2 * p + 2 * q + 2, 0) * derivatives_at_pt(mode_id(p, q, r - 1), XI0);

          derivatives_at_pt(mode_id(p, q, r + 1), XI1) =
              (a(r, 2 * p + 2 * q + 2, 0) * xi2 + b(r, 2 * p + 2 * q + 2, 0)) *
                  derivatives_at_pt(mode_id(p, q, r), XI1) -
              c(r, 2 * p + 2 * q + 2, 0) * derivatives_at_pt(mode_id(p, q, r - 1), XI1);

          derivatives_at_pt(mode_id(p, q, r + 1), XI2) =
              (a(r, 2 * p + 2 * q + 2, 0) * xi2 + b(r, 2 * p + 2 * q + 2, 0)) *
                  derivatives_at_pt(mode_id(p, q, r), XI2) +
              a(r, 2 * p + 2 * q + 2, 0) * values[mode_id(p, q, r)] -
              c(r, 2 * p + 2 * q + 2, 0) * derivatives_at_pt(mode_id(p, q, r - 1), XI2);
        }
      }
    }

    // Insert the values in derivative matrices

    derivative_values[KSI].insert_row(pt, derivatives_at_pt.const_col(KSI));
    derivative_values[ETA].insert_row(pt, derivatives_at_pt.const_col(ETA));
    derivative_values[ZTA].insert_row(pt, derivatives_at_pt.const_col(ZTA));

  } // Loop over all points
}

// ----------------------------------------------------------------------------

void DubinerExpansionTetra::is_leading_expansion_term(math::DenseDVec<bool> &is_leading_term)
{
  is_leading_term.resize((P + 1) * (P + 2) * (P + 3) / 6);
  is_leading_term.fill(false);

  if (P == P0)
  {
    is_leading_term[0] = true;
    return;
  }

  if (P == P1)
  {
    // In this case, we have four linear modes
    is_leading_term[0] = true;
    is_leading_term[1] = true;
    is_leading_term[2] = true;
    is_leading_term[3] = true;
    return;
  }

  // In what follows, we exploit the fact that in simplex expansions, the
  // highes-order modes are those modes (p,q,r) for which p + q + r == P

  is_leading_term[mode_id(P, 0, 0)]     = true;
  is_leading_term[mode_id(P - 1, 1, 0)] = true;

  for (Uint p = 0; p < (P - 1); ++p)
  {
    for (Uint q = 1; q < (P - p); ++q)
    {
      if ((p + q + 1) == P)
      {
        // Terms (p,q+1,0)
        is_leading_term[mode_id(p, q + 1, 0)] = true;
      }
    }
  }

  for (Uint p = 0; p < P; ++p)
  {
    for (Uint q = 0; q < (P - p); ++q)
    {
      if ((p + q + 1) == P)
      {
        // Terms (p,q,1)
        is_leading_term[mode_id(p, q, 1)] = true;
      }
    }
  }

  for (Uint p = 0; p < (P - 1); ++p)
  {
    for (Uint q = 0; q < (P - p - 1); ++q)
    {
      for (Uint r = 1; r < (P - p - q); ++r)
      {
        if ((p + q + r + 1) == P)
        {
          // Terms (p,q,r+1)
          is_leading_term[mode_id(p, q, r + 1)] = true;
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------

void DubinerExpansionTetra::mode_poly_deg(math::DenseDVec<Uint> &poly_deg)
{
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit
