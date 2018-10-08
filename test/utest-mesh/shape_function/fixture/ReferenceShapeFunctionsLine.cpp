#include "test/utest-mesh/shape_function/fixture/ReferenceShapeFunctionsLine.hpp"

namespace pdekit
{

namespace mesh
{

namespace utest_fixture
{

// ============================================================================

const Uint p1_line_Lagrange_sf::poly_order = P1;

// ============================================================================

void p1_line_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                               math::DenseDVec<Real> &values)
{
  values.resize(2);
  values[0] = 0.5 * (1.0 - ref_coords[KSI]);
  values[1] = 0.5 * (1.0 + ref_coords[KSI]);
}

// ============================================================================

void p1_line_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                     math::DenseDMat<Real> &values)
{
  values(0, KSI) = -0.5;
  values(1, KSI) = 0.5;
}

// ============================================================================

const Uint p2_line_Lagrange_sf::poly_order = P2;

// ============================================================================

void p2_line_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                               math::DenseDVec<Real> &values)
{
  values.resize(3);
  values[0] = -0.5 * ref_coords[KSI] * (1.0 - ref_coords[KSI]);
  values[1] = 0.5 * ref_coords[KSI] * (1.0 + ref_coords[KSI]);
  values[2] = 1.0 - ref_coords[KSI] * ref_coords[KSI];
}

// ============================================================================

void p2_line_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                     math::DenseDMat<Real> &values)
{
  values(0, KSI) = -0.5 + ref_coords[KSI];
  values(1, KSI) = 0.5 + ref_coords[KSI];
  values(2, KSI) = -2.0 * ref_coords[KSI];
}

// ============================================================================

const Uint p3_line_Lagrange_sf::poly_order = P3;

void p3_line_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                               math::DenseDVec<Real> &values)
{
  values.resize(4);
  values[0] = -9. / 16. * (ref_coords[KSI] + 1. / 3.) * (ref_coords[KSI] - 1. / 3.) *
              (ref_coords[KSI] - 1.);
  values[1] =
      9. / 16. * (ref_coords[KSI] + 1.) * (ref_coords[KSI] + 1. / 3.) * (ref_coords[KSI] - 1. / 3.);
  values[2] =
      27. / 16. * (ref_coords[KSI] + 1.) * (ref_coords[KSI] - 1. / 3.) * (ref_coords[KSI] - 1.);
  values[3] =
      -27. / 16. * (ref_coords[KSI] + 1.) * (ref_coords[KSI] + 1. / 3.) * (ref_coords[KSI] - 1.);
}

// ============================================================================

void p3_line_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                     math::DenseDMat<Real> &values)
{
  values(0, KSI) = -9. / 16. *
                   ((ref_coords[KSI] - 1. / 3.) * (ref_coords[KSI] - 1.0) +
                    (ref_coords[KSI] + 1. / 3.) * (ref_coords[KSI] - 1.0) +
                    (ref_coords[KSI] + 1. / 3.) * (ref_coords[KSI] - 1. / 3.));

  values(1, KSI) = 9. / 16. *
                   ((ref_coords[KSI] + 1. / 3.) * (ref_coords[KSI] - 1. / 3.) +
                    (ref_coords[KSI] + 1.) * (ref_coords[KSI] - 1. / 3.) +
                    (ref_coords[KSI] + 1.) * (ref_coords[KSI] + 1. / 3.));

  values(2, KSI) = 27. / 16. *
                   ((ref_coords[KSI] - 1. / 3.) * (ref_coords[KSI] - 1.0) +
                    (ref_coords[KSI] + 1.) * (ref_coords[KSI] - 1.0) +
                    (ref_coords[KSI] + 1.) * (ref_coords[KSI] - 1. / 3.));

  values(3, KSI) = -27. / 16. *
                   ((ref_coords[KSI] + 1. / 3.) * (ref_coords[KSI] - 1.0) +
                    (ref_coords[KSI] + 1.) * (ref_coords[KSI] - 1.0) +
                    (ref_coords[KSI] + 1.) * (ref_coords[KSI] + 1. / 3.));
}

// ============================================================================

} // namespace utest_fixture

} // namespace mesh

} // namespace pdekit
