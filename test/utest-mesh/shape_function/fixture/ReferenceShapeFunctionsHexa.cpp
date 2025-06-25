#include "test/utest-mesh/shape_function/fixture/ReferenceShapeFunctionsHexa.hpp"

namespace pdekit
{

namespace mesh
{

namespace utest_fixture
{

/// ===========================================================================

const Uint p1_hexa_Lagrange_sf::poly_order = P1;

// ============================================================================

void p1_hexa_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                               math::DenseDVec<Real> &values)
{
  values.resize(8);
  values[0] = 0.125 * (1.0 - ref_coords[XI0]) * (1.0 - ref_coords[XI1]) * (1.0 - ref_coords[XI2]);
  values[1] = 0.125 * (1.0 + ref_coords[XI0]) * (1.0 - ref_coords[XI1]) * (1.0 - ref_coords[XI2]);
  values[2] = 0.125 * (1.0 + ref_coords[XI0]) * (1.0 + ref_coords[XI1]) * (1.0 - ref_coords[XI2]);
  values[3] = 0.125 * (1.0 - ref_coords[XI0]) * (1.0 + ref_coords[XI1]) * (1.0 - ref_coords[XI2]);
  values[4] = 0.125 * (1.0 - ref_coords[XI0]) * (1.0 - ref_coords[XI1]) * (1.0 + ref_coords[XI2]);
  values[5] = 0.125 * (1.0 + ref_coords[XI0]) * (1.0 - ref_coords[XI1]) * (1.0 + ref_coords[XI2]);
  values[6] = 0.125 * (1.0 + ref_coords[XI0]) * (1.0 + ref_coords[XI1]) * (1.0 + ref_coords[XI2]);
  values[7] = 0.125 * (1.0 - ref_coords[XI0]) * (1.0 + ref_coords[XI1]) * (1.0 + ref_coords[XI2]);
}

// ============================================================================

void p1_hexa_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                     math::DenseDMat<Real> &values)
{
  values(0, XI0) = -0.125 * (1. - ref_coords[XI1]) * (1.0 - ref_coords[XI2]);
  values(1, XI0) = 0.125 * (1. - ref_coords[XI1]) * (1.0 - ref_coords[XI2]);
  values(2, XI0) = 0.125 * (1. + ref_coords[XI1]) * (1.0 - ref_coords[XI2]);
  values(3, XI0) = -0.125 * (1. + ref_coords[XI1]) * (1.0 - ref_coords[XI2]);
  values(4, XI0) = -0.125 * (1. - ref_coords[XI1]) * (1.0 + ref_coords[XI2]);
  values(5, XI0) = 0.125 * (1. - ref_coords[XI1]) * (1.0 + ref_coords[XI2]);
  values(6, XI0) = 0.125 * (1. + ref_coords[XI1]) * (1.0 + ref_coords[XI2]);
  values(7, XI0) = -0.125 * (1. + ref_coords[XI1]) * (1.0 + ref_coords[XI2]);

  values(0, XI1) = -0.125 * (1. - ref_coords[XI0]) * (1.0 - ref_coords[XI2]);
  values(1, XI1) = -0.125 * (1. + ref_coords[XI0]) * (1.0 - ref_coords[XI2]);
  values(2, XI1) = 0.125 * (1. + ref_coords[XI0]) * (1.0 - ref_coords[XI2]);
  values(3, XI1) = 0.125 * (1. - ref_coords[XI0]) * (1.0 - ref_coords[XI2]);
  values(4, XI1) = -0.125 * (1. - ref_coords[XI0]) * (1.0 + ref_coords[XI2]);
  values(5, XI1) = -0.125 * (1. + ref_coords[XI0]) * (1.0 + ref_coords[XI2]);
  values(6, XI1) = 0.125 * (1. + ref_coords[XI0]) * (1.0 + ref_coords[XI2]);
  values(7, XI1) = 0.125 * (1. - ref_coords[XI0]) * (1.0 + ref_coords[XI2]);

  values(0, XI2) = -0.125 * (1. - ref_coords[XI0]) * (1.0 - ref_coords[XI1]);
  values(1, XI2) = -0.125 * (1. + ref_coords[XI0]) * (1.0 - ref_coords[XI1]);
  values(2, XI2) = -0.125 * (1. + ref_coords[XI0]) * (1.0 + ref_coords[XI1]);
  values(3, XI2) = -0.125 * (1. - ref_coords[XI0]) * (1.0 + ref_coords[XI1]);
  values(4, XI2) = 0.125 * (1. - ref_coords[XI0]) * (1.0 - ref_coords[XI1]);
  values(5, XI2) = 0.125 * (1. + ref_coords[XI0]) * (1.0 - ref_coords[XI1]);
  values(6, XI2) = 0.125 * (1. + ref_coords[XI0]) * (1.0 + ref_coords[XI1]);
  values(7, XI2) = 0.125 * (1. - ref_coords[XI0]) * (1.0 + ref_coords[XI1]);
}

/// ===========================================================================

const Uint p2_hexa_Lagrange_sf::poly_order = P2;

// ============================================================================

void p2_hexa_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                               math::DenseDVec<Real> &values)
{
  values.resize(27);
  const Real ksi = ref_coords[XI0];
  const Real eta = ref_coords[XI1];
  const Real zta = ref_coords[XI2];

  const Real L0_ksi = -0.5 * ksi * (1. - ksi);
  const Real L1_ksi = 0.5 * ksi * (1. + ksi);
  const Real L2_ksi = 1. - ksi * ksi;

  const Real L0_eta = -0.5 * eta * (1. - eta);
  const Real L1_eta = 0.5 * eta * (1. + eta);
  const Real L2_eta = 1. - eta * eta;

  const Real L0_zta = -0.5 * zta * (1. - zta);
  const Real L1_zta = 0.5 * zta * (1. + zta);
  const Real L2_zta = 1. - zta * zta;

  values[0]  = L0_ksi * L0_eta * L0_zta;
  values[1]  = L1_ksi * L0_eta * L0_zta;
  values[2]  = L1_ksi * L1_eta * L0_zta;
  values[3]  = L0_ksi * L1_eta * L0_zta;
  values[4]  = L0_ksi * L0_eta * L1_zta;
  values[5]  = L1_ksi * L0_eta * L1_zta;
  values[6]  = L1_ksi * L1_eta * L1_zta;
  values[7]  = L0_ksi * L1_eta * L1_zta;
  values[8]  = L2_ksi * L0_eta * L0_zta;
  values[9]  = L0_ksi * L2_eta * L0_zta;
  values[10] = L0_ksi * L0_eta * L2_zta;
  values[11] = L1_ksi * L2_eta * L0_zta;
  values[12] = L1_ksi * L0_eta * L2_zta;
  values[13] = L2_ksi * L1_eta * L0_zta;
  values[14] = L1_ksi * L1_eta * L2_zta;
  values[15] = L0_ksi * L1_eta * L2_zta;
  values[16] = L2_ksi * L0_eta * L1_zta;
  values[17] = L0_ksi * L2_eta * L1_zta;
  values[18] = L1_ksi * L2_eta * L1_zta;
  values[19] = L2_ksi * L1_eta * L1_zta;
  values[20] = L2_ksi * L2_eta * L0_zta;
  values[21] = L2_ksi * L0_eta * L2_zta;
  values[22] = L0_ksi * L2_eta * L2_zta;
  values[23] = L1_ksi * L2_eta * L2_zta;
  values[24] = L2_ksi * L1_eta * L2_zta;
  values[25] = L2_ksi * L2_eta * L1_zta;
  values[26] = L2_ksi * L2_eta * L2_zta;
}

// ============================================================================

void p2_hexa_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                     math::DenseDMat<Real> &values)
{
  const Real ksi = ref_coords[XI0];
  const Real eta = ref_coords[XI1];
  const Real zta = ref_coords[XI2];

  const Real L0_ksi = -0.5 * ksi * (1. - ksi);
  const Real L1_ksi = 0.5 * ksi * (1. + ksi);
  const Real L2_ksi = 1. - ksi * ksi;

  const Real dL0_ksi = -0.5 + ksi;
  const Real dL1_ksi = 0.5 + ksi;
  const Real dL2_ksi = -2.0 * ksi;

  const Real L0_eta = -0.5 * eta * (1. - eta);
  const Real L1_eta = 0.5 * eta * (1. + eta);
  const Real L2_eta = 1. - eta * eta;

  const Real dL0_eta = -0.5 + eta;
  const Real dL1_eta = 0.5 + eta;
  const Real dL2_eta = -2.0 * eta;

  const Real L0_zta = -0.5 * zta * (1. - zta);
  const Real L1_zta = 0.5 * zta * (1. + zta);
  const Real L2_zta = 1. - zta * zta;

  const Real dL0_zta = -0.5 + zta;
  const Real dL1_zta = 0.5 + zta;
  const Real dL2_zta = -2.0 * zta;

  values(0, XI0)  = dL0_ksi * L0_eta * L0_zta;
  values(1, XI0)  = dL1_ksi * L0_eta * L0_zta;
  values(2, XI0)  = dL1_ksi * L1_eta * L0_zta;
  values(3, XI0)  = dL0_ksi * L1_eta * L0_zta;
  values(4, XI0)  = dL0_ksi * L0_eta * L1_zta;
  values(5, XI0)  = dL1_ksi * L0_eta * L1_zta;
  values(6, XI0)  = dL1_ksi * L1_eta * L1_zta;
  values(7, XI0)  = dL0_ksi * L1_eta * L1_zta;
  values(8, XI0)  = dL2_ksi * L0_eta * L0_zta;
  values(9, XI0)  = dL0_ksi * L2_eta * L0_zta;
  values(10, XI0) = dL0_ksi * L0_eta * L2_zta;
  values(11, XI0) = dL1_ksi * L2_eta * L0_zta;
  values(12, XI0) = dL1_ksi * L0_eta * L2_zta;
  values(13, XI0) = dL2_ksi * L1_eta * L0_zta;
  values(14, XI0) = dL1_ksi * L1_eta * L2_zta;
  values(15, XI0) = dL0_ksi * L1_eta * L2_zta;
  values(16, XI0) = dL2_ksi * L0_eta * L1_zta;
  values(17, XI0) = dL0_ksi * L2_eta * L1_zta;
  values(18, XI0) = dL1_ksi * L2_eta * L1_zta;
  values(19, XI0) = dL2_ksi * L1_eta * L1_zta;
  values(20, XI0) = dL2_ksi * L2_eta * L0_zta;
  values(21, XI0) = dL2_ksi * L0_eta * L2_zta;
  values(22, XI0) = dL0_ksi * L2_eta * L2_zta;
  values(23, XI0) = dL1_ksi * L2_eta * L2_zta;
  values(24, XI0) = dL2_ksi * L1_eta * L2_zta;
  values(25, XI0) = dL2_ksi * L2_eta * L1_zta;
  values(26, XI0) = dL2_ksi * L2_eta * L2_zta;

  values(0, XI1)  = L0_ksi * dL0_eta * L0_zta;
  values(1, XI1)  = L1_ksi * dL0_eta * L0_zta;
  values(2, XI1)  = L1_ksi * dL1_eta * L0_zta;
  values(3, XI1)  = L0_ksi * dL1_eta * L0_zta;
  values(4, XI1)  = L0_ksi * dL0_eta * L1_zta;
  values(5, XI1)  = L1_ksi * dL0_eta * L1_zta;
  values(6, XI1)  = L1_ksi * dL1_eta * L1_zta;
  values(7, XI1)  = L0_ksi * dL1_eta * L1_zta;
  values(8, XI1)  = L2_ksi * dL0_eta * L0_zta;
  values(9, XI1)  = L0_ksi * dL2_eta * L0_zta;
  values(10, XI1) = L0_ksi * dL0_eta * L2_zta;
  values(11, XI1) = L1_ksi * dL2_eta * L0_zta;
  values(12, XI1) = L1_ksi * dL0_eta * L2_zta;
  values(13, XI1) = L2_ksi * dL1_eta * L0_zta;
  values(14, XI1) = L1_ksi * dL1_eta * L2_zta;
  values(15, XI1) = L0_ksi * dL1_eta * L2_zta;
  values(16, XI1) = L2_ksi * dL0_eta * L1_zta;
  values(17, XI1) = L0_ksi * dL2_eta * L1_zta;
  values(18, XI1) = L1_ksi * dL2_eta * L1_zta;
  values(19, XI1) = L2_ksi * dL1_eta * L1_zta;
  values(20, XI1) = L2_ksi * dL2_eta * L0_zta;
  values(21, XI1) = L2_ksi * dL0_eta * L2_zta;
  values(22, XI1) = L0_ksi * dL2_eta * L2_zta;
  values(23, XI1) = L1_ksi * dL2_eta * L2_zta;
  values(24, XI1) = L2_ksi * dL1_eta * L2_zta;
  values(25, XI1) = L2_ksi * dL2_eta * L1_zta;
  values(26, XI1) = L2_ksi * dL2_eta * L2_zta;

  values(0, XI2)  = L0_ksi * L0_eta * dL0_zta;
  values(1, XI2)  = L1_ksi * L0_eta * dL0_zta;
  values(2, XI2)  = L1_ksi * L1_eta * dL0_zta;
  values(3, XI2)  = L0_ksi * L1_eta * dL0_zta;
  values(4, XI2)  = L0_ksi * L0_eta * dL1_zta;
  values(5, XI2)  = L1_ksi * L0_eta * dL1_zta;
  values(6, XI2)  = L1_ksi * L1_eta * dL1_zta;
  values(7, XI2)  = L0_ksi * L1_eta * dL1_zta;
  values(8, XI2)  = L2_ksi * L0_eta * dL0_zta;
  values(9, XI2)  = L0_ksi * L2_eta * dL0_zta;
  values(10, XI2) = L0_ksi * L0_eta * dL2_zta;
  values(11, XI2) = L1_ksi * L2_eta * dL0_zta;
  values(12, XI2) = L1_ksi * L0_eta * dL2_zta;
  values(13, XI2) = L2_ksi * L1_eta * dL0_zta;
  values(14, XI2) = L1_ksi * L1_eta * dL2_zta;
  values(15, XI2) = L0_ksi * L1_eta * dL2_zta;
  values(16, XI2) = L2_ksi * L0_eta * dL1_zta;
  values(17, XI2) = L0_ksi * L2_eta * dL1_zta;
  values(18, XI2) = L1_ksi * L2_eta * dL1_zta;
  values(19, XI2) = L2_ksi * L1_eta * dL1_zta;
  values(20, XI2) = L2_ksi * L2_eta * dL0_zta;
  values(21, XI2) = L2_ksi * L0_eta * dL2_zta;
  values(22, XI2) = L0_ksi * L2_eta * dL2_zta;
  values(23, XI2) = L1_ksi * L2_eta * dL2_zta;
  values(24, XI2) = L2_ksi * L1_eta * dL2_zta;
  values(25, XI2) = L2_ksi * L2_eta * dL1_zta;
  values(26, XI2) = L2_ksi * L2_eta * dL2_zta;
}

/// ===========================================================================

} // namespace utest_fixture

} // namespace mesh

} // namespace pdekit
