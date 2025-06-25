#include "test/utest-mesh/shape_function/fixture/ReferenceShapeFunctionsTriag.hpp"

namespace pdekit
{

namespace mesh
{

namespace utest_fixture
{

// ----------------------------------------------------------------------------

const Uint p1_triag_Lagrange_sf::poly_order = P1;

// ----------------------------------------------------------------------------

void p1_triag_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                                math::DenseDVec<Real> &values)
{
  values.resize(3);
  values[0] = -0.5 * (ref_coords[KSI] + ref_coords[ETA]);
  values[1] = 0.5 * (1.0 + ref_coords[KSI]);
  values[2] = 0.5 * (1.0 + ref_coords[ETA]);
}

// ----------------------------------------------------------------------------

void p1_triag_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                      math::DenseDMat<Real> &values)
{
  values(0, KSI) = -0.5;
  values(0, ETA) = -0.5;
  values(1, KSI) = 0.5;
  values(1, ETA) = 0.0;
  values(2, KSI) = 0.0;
  values(2, ETA) = 0.5;
}

// ----------------------------------------------------------------------------

const Uint p2_triag_Lagrange_sf::poly_order = P2;

// ----------------------------------------------------------------------------

void p2_triag_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                                math::DenseDVec<Real> &values)
{
  values.resize(6);
  values[0] = 0.5 * (ref_coords[KSI] + ref_coords[ETA]) * (1. + ref_coords[KSI] + ref_coords[ETA]);
  values[1] = 0.5 * ref_coords[KSI] * (1. + ref_coords[KSI]);
  values[2] = 0.5 * ref_coords[ETA] * (1. + ref_coords[ETA]);
  values[3] = -(ref_coords[KSI] + 1.) * (ref_coords[KSI] + ref_coords[ETA]);
  values[4] = (1. + ref_coords[KSI]) * (1. + ref_coords[ETA]);
  values[5] = -(ref_coords[ETA] + 1.) * (ref_coords[KSI] + ref_coords[ETA]);
}

// ----------------------------------------------------------------------------

void p2_triag_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                      math::DenseDMat<Real> &values)
{
  values(0, KSI) = 0.5 + ref_coords[KSI] + ref_coords[ETA];
  values(0, ETA) = 0.5 + ref_coords[KSI] + ref_coords[ETA];
  values(1, KSI) = 0.5 + ref_coords[KSI];
  values(1, ETA) = 0.0;
  values(2, KSI) = 0.0;
  values(2, ETA) = 0.5 + ref_coords[ETA];
  values(3, KSI) = -1. - ref_coords[ETA] - 2. * ref_coords[KSI];
  values(3, ETA) = -1. - ref_coords[KSI];
  values(4, KSI) = 1. + ref_coords[ETA];
  values(4, ETA) = 1. + ref_coords[KSI];
  values(5, KSI) = -1. - ref_coords[ETA];
  values(5, ETA) = -1. - ref_coords[KSI] - 2. * ref_coords[ETA];
}

// ----------------------------------------------------------------------------

const Uint p3_triag_Lagrange_sf::poly_order = P3;

void p3_triag_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                                math::DenseDVec<Real> &values)
{
  values.resize(10);
  const Real L0 = -0.5 * (ref_coords[KSI] + ref_coords[ETA]);
  const Real L1 = 0.5 * (1.0 + ref_coords[KSI]);
  const Real L2 = 0.5 * (1.0 + ref_coords[ETA]);

  values[0] = 0.5 * (3. * L0 - 1.) * (3. * L0 - 2.) * L0;
  values[1] = 0.5 * (3. * L1 - 1.) * (3. * L1 - 2.) * L1;
  values[2] = 0.5 * (3. * L2 - 1.) * (3. * L2 - 2.) * L2;

  values[3] = 9. / 2. * L0 * L1 * (3. * L0 - 1.);
  values[4] = 9. / 2. * L0 * L1 * (3. * L1 - 1.);

  values[5] = 9. / 2. * L1 * L2 * (3. * L1 - 1.);
  values[6] = 9. / 2. * L1 * L2 * (3. * L2 - 1.);

  values[7] = 9. / 2. * L2 * L0 * (3. * L2 - 1.);
  values[8] = 9. / 2. * L2 * L0 * (3. * L0 - 1.);

  values[9] = 27. * L0 * L1 * L2;
}

// ----------------------------------------------------------------------------

void p3_triag_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                      math::DenseDMat<Real> &values)
{
  const Real L0 = -0.5 * (ref_coords[KSI] + ref_coords[ETA]);
  const Real L1 = 0.5 * (1.0 + ref_coords[KSI]);
  const Real L2 = 0.5 * (1.0 + ref_coords[ETA]);

  const Real dL0dxi  = -0.5;
  const Real dL0deta = -0.5;
  const Real dL1dxi  = 0.5;
  const Real dL1deta = 0.0;
  const Real dL2dxi  = 0.0;
  const Real dL2deta = 0.5;

  values(0, KSI) = 0.5 * dL0dxi * (27. * L0 * L0 - 18. * L0 + 2.);
  values(1, KSI) = 0.5 * dL1dxi * (27. * L1 * L1 - 18. * L1 + 2.);
  values(2, KSI) = 0.5 * dL2dxi * (27. * L2 * L2 - 18. * L2 + 2.);
  values(3, KSI) = 4.5 * (dL0dxi * (6. * L0 * L1 - L1) + dL1dxi * L0 * (3. * L0 - 1.));
  values(4, KSI) = 4.5 * (dL1dxi * (6. * L0 * L1 - L0) + dL0dxi * L1 * (3. * L1 - 1.));
  values(5, KSI) = 4.5 * (dL1dxi * (6. * L1 * L2 - L2) + dL2dxi * L1 * (3. * L1 - 1.));
  values(6, KSI) = 4.5 * (dL2dxi * (6. * L1 * L2 - L1) + dL1dxi * L2 * (3. * L2 - 1.));
  values(7, KSI) = 4.5 * (dL2dxi * (6. * L0 * L2 - L0) + dL0dxi * L2 * (3. * L2 - 1.));
  values(8, KSI) = 4.5 * (dL0dxi * (6. * L0 * L2 - L2) + dL2dxi * L0 * (3. * L0 - 1.));
  values(9, KSI) = 27.0 * (dL0dxi * L1 * L2 + L0 * dL1dxi * L2 + L0 * L1 * dL2dxi);

  values(0, ETA) = 0.5 * dL0deta * (27. * L0 * L0 - 18. * L0 + 2.);
  values(1, ETA) = 0.5 * dL1deta * (27. * L1 * L1 - 18. * L1 + 2.);
  values(2, ETA) = 0.5 * dL2deta * (27. * L2 * L2 - 18. * L2 + 2.);
  values(3, ETA) = 4.5 * (dL0deta * (6. * L0 * L1 - L1) + dL1deta * L0 * (3. * L0 - 1.));
  values(4, ETA) = 4.5 * (dL1deta * (6. * L0 * L1 - L0) + dL0deta * L1 * (3. * L1 - 1.));
  values(5, ETA) = 4.5 * (dL1deta * (6. * L1 * L2 - L2) + dL2deta * L1 * (3. * L1 - 1.));
  values(6, ETA) = 4.5 * (dL2deta * (6. * L1 * L2 - L1) + dL1deta * L2 * (3. * L2 - 1.));
  values(7, ETA) = 4.5 * (dL2deta * (6. * L0 * L2 - L0) + dL0deta * L2 * (3. * L2 - 1.));
  values(8, ETA) = 4.5 * (dL0deta * (6. * L0 * L2 - L2) + dL2deta * L0 * (3. * L0 - 1.));
  values(9, ETA) = 27.0 * (dL0deta * L1 * L2 + L0 * dL1deta * L2 + L0 * L1 * dL2deta);
}

// ----------------------------------------------------------------------------

} // namespace utest_fixture

} // namespace mesh

} // namespace pdekit
