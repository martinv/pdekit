#ifndef PDEKIT_Mesh_Shape_Function_C0Expansion_Triag_hpp
#define PDEKIT_Mesh_Shape_Function_C0Expansion_Triag_hpp

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
/// C0 continuous Dubiner orthogonal basis on triangles
// ----------------------------------------------------------------------------

template <Uint Order>
class C0ExpansionTriag
{
  public:
  enum
  {
    P = Order
  };

  enum
  {
    nmodes = (P + 1) * (P + 2) / 2
  };

  // enum { nmodes = 3 + 3 * (P-1) + (P-2)*(P-1)/2 };

  /// Default constructor
  C0ExpansionTriag();

  /// Destructor
  ~C0ExpansionTriag();

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
C0ExpansionTriag<Order>::C0ExpansionTriag()
{
}

// ----------------------------------------------------------------------------

template <Uint Order>
C0ExpansionTriag<Order>::~C0ExpansionTriag()
{
}

// ----------------------------------------------------------------------------

template <Uint Order>
void C0ExpansionTriag<Order>::evaluate_in_one_point(const math::DenseDVec<Real> &point,
                                                    math::DenseDVec<Real> &values)
{
  math::JacobiPolynomial jp;

  values.resize(nmodes);

  // Variables transformed from triangle to square
  Real eta0 = 0.0;
  Real eta1 = 0.0;

  const Real ksi0 = point[KSI];
  const Real ksi1 = point[ETA];

  if (((1.0 - 1e-6) < ksi1) && (ksi1 < 1.0 + 1e-6))
  {
    eta0 = -1.0;
  }
  else
  {
    eta0 = 2. * (1. + ksi0) / (1. - ksi1) - 1.;
  }
  eta1 = ksi1;

  /// First assign vertex modes

  values[0] = 0.25 * (1. - eta0) * (1. - eta1);
  values[1] = 0.25 * (1. + eta0) * (1. - eta1);
  values[2] = 0.5 * (1. + eta1);

  Uint mode_id = 3;

  /// Modes on edge 0-1
  for (Uint p = 1; p < P; ++p)
  {
    values[mode_id] = 0.25 * (1. - eta0) * (1. + eta0) * jp(p - 1, 1, 1, eta0) *
                      std::pow(0.5 * (1. - eta1), p + 1);
    mode_id++;
  }

  /// Modes on edge 2-0
  for (Uint q = 1; q < P; ++q)
  {
    values[mode_id] = 0.125 * (1. - eta0) * (1. - eta1) * (1. + eta1) * jp(q - 1, 1, 1, eta1);
    mode_id++;
  }

  /// Modes on edge 1-2
  for (Uint q = 1; q < P; ++q)
  {
    values[mode_id] = 0.125 * (1. + eta0) * (1. - eta1) * (1. + eta1) * jp(q - 1, 1, 1, eta1);
    mode_id++;
  }

  /// Interior modes

  for (Uint p = 1; p < P; ++p)
  {
    for (Uint q = 1; q < (P - p); ++q)
    {
      const Real z    = std::pow(0.5 * (1. - eta1), p + 1);
      values[mode_id] = 0.125 * (1. - eta0) * (1. + eta0) * (1. + eta1) * jp(p - 1, 1, 1, eta0) *
                        jp(q - 1, 2 * p + 1, 1, eta1) * z;
      mode_id++;
    }
  }
}

// ----------------------------------------------------------------------------

template <Uint Order>
void C0ExpansionTriag<Order>::Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                                                 math::DenseDMat<Real> &values)

{
  math::JacobiPolynomial jp;

  const Uint nb_pts = coordinates.rows();

  values.resize(nb_pts, nmodes);

  // Variables transformed from triangle to square
  Real eta0 = 0.0;
  Real eta1 = 0.0;

  Uint mode_id = 0;

  for (Uint pt = 0; pt < nb_pts; ++pt)
  {

    const Real ksi0 = coordinates(pt, KSI);
    const Real ksi1 = coordinates(pt, ETA);

    if (((1.0 - 1e-6) < ksi1) && (ksi1 < 1.0 + 1e-6))
    {
      eta0 = -1.0;
    }
    else
    {
      eta0 = 2. * (1. + ksi0) / (1. - ksi1) - 1.;
    }
    eta1 = ksi1;

    /// First assign vertex modes

    values(pt, 0) = 0.25 * (1. - eta0) * (1. - eta1);
    values(pt, 1) = 0.25 * (1. + eta0) * (1. - eta1);
    values(pt, 2) = 0.5 * (1. + eta1);

    mode_id = 3;

    /// Modes on edge 0-1
    for (Uint p = 1; p < P; ++p)
    {
      values(pt, mode_id) = 0.25 * (1. - eta0) * (1. + eta0) * jp(p - 1, 1, 1, eta0) *
                            std::pow(0.5 * (1. - eta1), p + 1);
      mode_id++;
    }

    /// Modes on edge 2-0
    for (Uint q = 1; q < P; ++q)
    {
      values(pt, mode_id) = 0.125 * (1. - eta0) * (1. - eta1) * (1. + eta1) * jp(q - 1, 1, 1, eta1);
      mode_id++;
    }

    /// Modes on edge 1-2
    for (Uint q = 1; q < P; ++q)
    {
      values(pt, mode_id) = 0.125 * (1. + eta0) * (1. - eta1) * (1. + eta1) * jp(q - 1, 1, 1, eta1);
      mode_id++;
    }

    /// Interior modes

    for (Uint p = 1; p < P; ++p)
    {
      for (Uint q = 1; q < (P - p); ++q)
      {
        const Real z        = std::pow(0.5 * (1. - eta1), p + 1);
        values(pt, mode_id) = 0.125 * (1. - eta0) * (1. + eta0) * (1. + eta1) *
                              jp(p - 1, 1, 1, eta0) * jp(q - 1, 2 * p + 1, 1, eta1) * z;
        std::cout << "mode id= " << mode_id << " ";
        mode_id++;
      }
    }

    std::cout << std::endl;
  } // Loop over points 'pt'
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit

#endif
