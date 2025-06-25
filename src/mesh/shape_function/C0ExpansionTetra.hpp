#ifndef PDEKIT_Mesh_Shape_Function_C0Expansion_Tetra_hpp
#define PDEKIT_Mesh_Shape_Function_C0Expansion_Tetra_hpp

#include <limits>

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/vector.hpp>

#include "common/Meta.hpp"
#include "common/PDEKit.hpp"
#include "math/polynomials/JacobiPolynomial.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

// ----------------------------------------------------------------------------
// C0 continuous Dubiner orthogonal basis on tetrahedra
// ----------------------------------------------------------------------------

template <Uint Order>
class C0ExpansionTetra
{
  public:
  enum
  {
    P = Order
  };

  enum
  {
    nmodes = (P + 1) * (P + 2) * (P + 3) / 6
  };

  /// Default constructor
  C0ExpansionTetra();

  /// Destructor
  ~C0ExpansionTetra();

  /// Compute the values of prime basis modes at one point
  void evaluate_in_one_point(const math::DenseDVec<Real> &point, math::DenseDVec<Real> &values);

  /// Evaluate the prime basis
  void Vandermonde_matrix(const math::DenseDMat<Real> &coordinates, math::DenseDMat<Real> &values);

  /// Compute the derivatives of the prime basis
  void Vandermonde_matrix_derivatives(const math::DenseDMat<Real> &coordinates,
                                      math::DenseDMat<Real> &values);

  private:
};

// ----------------------------------------------------------------------------

template <Uint Order>
C0ExpansionTetra<Order>::C0ExpansionTetra()
{
}

// ----------------------------------------------------------------------------

template <Uint Order>
C0ExpansionTetra<Order>::~C0ExpansionTetra()
{
}

// ----------------------------------------------------------------------------

template <Uint Order>
void C0ExpansionTetra<Order>::evaluate_in_one_point(const math::DenseDVec<Real> &point,
                                                    math::DenseDVec<Real> &values)
{
  math::JacobiPolynomial jp;

  values.resize(nmodes);

  // Variables transformed from tetrahedra to cube
  Real eta0 = 0.0;
  Real eta1 = 0.0;
  Real eta2 = 0.0;

  const Real ksi0 = point[KSI];
  const Real ksi1 = point[ETA];
  const Real ksi2 = point[ZTA];

  if ((-1.e-6 < (ksi1 + ksi2)) && ((ksi1 + ksi2) < 1.e-6))
  {
    eta0 = -1.0;
  }
  else
  {
    eta0 = -2. * (1. + ksi0) / (ksi1 + ksi2) - 1.;
  }

  if (((1.0 - 1e-6) < ksi2) && (ksi2 < 1.0 + 1e-6))
  {
    eta1 = -1.0;
  }
  else
  {
    eta1 = 2. * (1. + ksi1) / (1. - ksi2) - 1.;
  }
  eta2 = ksi2;

  /// First assign vertex modes

  values[0] = 0.125 * (1. - eta0) * (1. - eta1) * (1. - eta2);
  values[1] = 0.125 * (1. + eta0) * (1. - eta1) * (1. - eta2);
  values[2] = 0.125 * (1. - eta0) * (1. + eta1) * (1. - eta2);
  values[3] = 0.5 * (1. + eta2);

  Uint mode_id = 4;

  /// Modes on edge 0-1 (A-B)
  for (Uint p = 1; p < P; ++p)
  {
    values[mode_id] = 0.25 * (1. - eta0) * (1. + eta0) * jp(p - 1, 1, 1, eta0) *
                      std::pow(0.5 * (1. - eta1), p + 1) * std::pow(0.5 * (1. - eta2), p + 1);
    mode_id++;
  }

  /// Modes on edge 2-0 (A-C)
  for (Uint q = 1; q < P; ++q)
  {
    values[mode_id] = 0.125 * (1. - eta0) * (1. - eta1) * (1. + eta1) * jp(q - 1, 1, 1, eta1) *
                      std::pow(0.5 * (1. - eta2), q + 1);
    mode_id++;
  }

  /// Modes on edge 1-2 (B-C)
  for (Uint q = 1; q < P; ++q)
  {
    values[mode_id] = 0.125 * (1. + eta0) * (1. - eta1) * (1. + eta1) * jp(q - 1, 1, 1, eta1) *
                      std::pow(0.5 * (1. - eta2), q + 1);
    mode_id++;
  }

  /// Modes on edge 0-3 (A-D)
  for (Uint r = 1; r < P; ++r)
  {
    values[mode_id] =
        0.0625 * (1. - eta0) * (1. - eta1) * (1. - eta2) * (1. + eta2) * jp(r - 1, 1, 1, eta2);
    mode_id++;
  }

  /// Modes on edge 1-3 (B-D)
  for (Uint r = 1; r < P; ++r)
  {
    values[mode_id] =
        0.0625 * (1. + eta0) * (1. - eta1) * (1. - eta2) * (1. + eta2) * jp(r - 1, 1, 1, eta2);
    mode_id++;
  }

  /// Modes on edge 2-3 (C-D)
  for (Uint r = 1; r < P; ++r)
  {
    values[mode_id] = 0.125 * (1. + eta1) * (1. - eta2) * (1. + eta2) * jp(r - 1, 1, 1, eta2);
    mode_id++;
  }

  /// Face modes

  /// Face A-B-C (0-2-1)
  for (Uint p = 1; p < P; ++p)
  {
    for (Uint q = 1; q < (P - p); ++q)
    {
      values[mode_id] = 0.125 * (1. - eta0) * (1. + eta0) * jp(p - 1, 1, 1, eta0) * (1. + eta1) *
                        jp(q - 1, 2 * p + 1, 1, eta1) * std::pow(0.5 * (1. - eta2), p + q + 1);
      mode_id++;
    }
  }

  /// Face A-B-E (0-1-3)
  for (Uint p = 1; p < P; ++p)
  {
    for (Uint r = 1; r < (P - p); ++r)
    {
      values[mode_id] = 0.125 * (1. - eta0) * (1. + eta0) * jp(p - 1, 1, 1, eta0) *
                        std::pow(0.5 * (1. - eta1), p + 1) * std::pow(0.5 * (1. - eta2), p + 1) *
                        (1. + eta2) * jp(r - 1, 2 * p + 1, 1, eta2);
      mode_id++;
    }
  }

  /// Face A-C-E (0-3-2)
  for (Uint q = 1; q < P; ++q)
  {
    for (Uint r = 1; r < (P - q); ++r)
    {
      values[mode_id] = 0.0625 * (1. - eta0) * (1. - eta1) * (1. + eta1) * jp(q - 1, 1, 1, eta1) *
                        std::pow(0.5 * (1. - eta2), q + 1) * (1. + eta2) *
                        jp(r - 1, 2 * q + 1, 1, eta2);
      mode_id++;
    }
  }

  /// Face B-C-E (3-1-2)
  for (Uint q = 1; q < P; ++q)
  {
    for (Uint r = 1; r < (P - q); ++r)
    {
      values[mode_id] = 0.0625 * (1. + eta0) * (1. - eta1) * (1. + eta1) * jp(q - 1, 1, 1, eta1) *
                        std::pow(0.5 * (1. - eta2), q + 1) * (1. + eta2) *
                        jp(r - 1, 2 * q + 1, 1, eta2);
      mode_id++;
    }
  }

  /// Interior modes

  for (Uint p = 1; p < P; ++p)
  {
    for (Uint q = 1; q < (P - p); ++q)
    {
      for (Uint r = 1; r < (P - (p + q)); ++r)
      {
        values[mode_id] = 0.0625 * (1. - eta0) * (1. + eta0) * jp(p - 1, 1, 1, eta0) *
                          std::pow(0.5 * (1. - eta1), p + 1) * (1. + eta1) *
                          jp(q - 1, 2 * p + 1, 1, eta1) * std::pow(0.5 * (1. - eta2), p + q + 1) *
                          (1. + eta2) * jp(r - 1, 2 * p + 2 * q + 1, 1, eta2);
        mode_id++;
      }
    }
  }
}

// ----------------------------------------------------------------------------

template <Uint Order>
void C0ExpansionTetra<Order>::Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                                                 math::DenseDMat<Real> &values)

{
  math::JacobiPolynomial jp;

  const Uint nb_pts = coordinates.rows();

  values.resize(nb_pts, nmodes);

  // Variables transformed from tetra ref. element to hexa
  Real eta0 = 0.0;
  Real eta1 = 0.0;
  Real eta2 = 0.0;

  Uint mode_id = 0;

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {

    const Real ksi0 = coordinates(pt, KSI);
    const Real ksi1 = coordinates(pt, ETA);
    const Real ksi2 = coordinates(pt, ZTA);

    if ((-1.e-6 < (ksi1 + ksi2)) && ((ksi1 + ksi2) < 1.e-6))
    {
      eta0 = -1.0;
    }
    else
    {
      eta0 = -2. * (1. + ksi0) / (ksi1 + ksi2) - 1.;
    }

    if (((1.0 - 1e-6) < ksi2) && (ksi2 < 1.0 + 1e-6))
    {
      eta1 = -1.0;
    }
    else
    {
      eta1 = 2. * (1. + ksi1) / (1. - ksi2) - 1.;
    }
    eta2 = ksi2;

    /// First assign vertex modes

    values(pt, 0) = 0.125 * (1. - eta0) * (1. - eta1) * (1. - eta2);
    values(pt, 1) = 0.125 * (1. + eta0) * (1. - eta1) * (1. - eta2);
    values(pt, 2) = 0.125 * (1. - eta0) * (1. + eta1) * (1. - eta2);
    values(pt, 3) = 0.5 * (1. + eta2);

    mode_id = 4;

    /// Modes on edge 0-1 (A-B)
    for (Uint p = 1; p < P; ++p)
    {
      values(pt, mode_id) = 0.25 * (1. - eta0) * (1. + eta0) * jp(p - 1, 1, 1, eta0) *
                            std::pow(0.5 * (1. - eta1), p + 1) * std::pow(0.5 * (1. - eta2), p + 1);
      mode_id++;
    }

    /// Modes on edge 2-0 (A-C)
    for (Uint q = 1; q < P; ++q)
    {
      values(pt, mode_id) = 0.125 * (1. - eta0) * (1. - eta1) * (1. + eta1) *
                            jp(q - 1, 1, 1, eta1) * std::pow(0.5 * (1. - eta2), q + 1);
      mode_id++;
    }

    /// Modes on edge 1-2 (B-C)
    for (Uint q = 1; q < P; ++q)
    {
      values(pt, mode_id) = 0.125 * (1. + eta0) * (1. - eta1) * (1. + eta1) *
                            jp(q - 1, 1, 1, eta1) * std::pow(0.5 * (1. - eta2), q + 1);
      mode_id++;
    }

    /// Modes on edge 0-3 (A-D)
    for (Uint r = 1; r < P; ++r)
    {
      values(pt, mode_id) =
          0.0625 * (1. - eta0) * (1. - eta1) * (1. - eta2) * (1. + eta2) * jp(r - 1, 1, 1, eta2);
      mode_id++;
    }

    /// Modes on edge 1-3 (B-D)
    for (Uint r = 1; r < P; ++r)
    {
      values(pt, mode_id) =
          0.0625 * (1. + eta0) * (1. - eta1) * (1. - eta2) * (1. + eta2) * jp(r - 1, 1, 1, eta2);
      mode_id++;
    }

    /// Modes on edge 2-3 (C-D)
    for (Uint r = 1; r < P; ++r)
    {
      values(pt, mode_id) = 0.125 * (1. + eta1) * (1. - eta2) * (1. + eta2) * jp(r - 1, 1, 1, eta2);
      mode_id++;
    }

    /// Face modes

    /// Face A-B-C (0-2-1)
    for (Uint p = 1; p < P; ++p)
    {
      for (Uint q = 1; q < (P - p); ++q)
      {
        values(pt, mode_id) = 0.125 * (1. - eta0) * (1. + eta0) * jp(p - 1, 1, 1, eta0) *
                              (1. + eta1) * jp(q - 1, 2 * p + 1, 1, eta1) *
                              std::pow(0.5 * (1. - eta2), p + q + 1);
        mode_id++;
      }
    }

    /// Face A-B-E (0-1-3)
    for (Uint p = 1; p < P; ++p)
    {
      for (Uint r = 1; r < (P - p); ++r)
      {
        values(pt, mode_id) = 0.125 * (1. - eta0) * (1. + eta0) * jp(p - 1, 1, 1, eta0) *
                              std::pow(0.5 * (1. - eta1), p + 1) *
                              std::pow(0.5 * (1. - eta2), p + 1) * (1. + eta2) *
                              jp(r - 1, 2 * p + 1, 1, eta2);
        mode_id++;
      }
    }

    /// Face A-C-E (0-3-2)
    for (Uint q = 1; q < P; ++q)
    {
      for (Uint r = 1; r < (P - q); ++r)
      {
        values(pt, mode_id) = 0.0625 * (1. - eta0) * (1. - eta1) * (1. + eta1) *
                              jp(q - 1, 1, 1, eta1) * std::pow(0.5 * (1. - eta2), q + 1) *
                              (1. + eta2) * jp(r - 1, 2 * q + 1, 1, eta2);
        mode_id++;
      }
    }

    /// Face B-C-E (3-1-2)
    for (Uint q = 1; q < P; ++q)
    {
      for (Uint r = 1; r < (P - q); ++r)
      {
        values(pt, mode_id) = 0.0625 * (1. + eta0) * (1. - eta1) * (1. + eta1) *
                              jp(q - 1, 1, 1, eta1) * std::pow(0.5 * (1. - eta2), q + 1) *
                              (1. + eta2) * jp(r - 1, 2 * q + 1, 1, eta2);
        mode_id++;
      }
    }

    /// Interior modes

    for (Uint p = 1; p < P; ++p)
    {
      for (Uint q = 1; q < (P - p); ++q)
      {
        for (Uint r = 1; r < (P - p - q); ++r)
        {
          values(pt, mode_id) = 0.0625 * (1. - eta0) * (1. + eta0) * jp(p - 1, 1, 1, eta0) *
                                std::pow(0.5 * (1. - eta1), p + 1) * (1. + eta1) *
                                jp(q - 1, 2 * p + 1, 1, eta1) *
                                std::pow(0.5 * (1. - eta2), p + q + 1) * (1. + eta2) *
                                jp(r - 1, 2 * p + 2 * q + 1, 1, eta2);
          mode_id++;
        }
      }
    }

  } // Loop over points 'pt'
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit

#endif
