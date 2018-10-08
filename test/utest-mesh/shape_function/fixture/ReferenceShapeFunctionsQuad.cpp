#include "test/utest-mesh/shape_function/fixture/ReferenceShapeFunctionsQuad.hpp"

namespace pdekit
{

namespace mesh
{

namespace utest_fixture
{

/// ===========================================================================

const Uint p1_quad_Lagrange_sf::poly_order = P1;

// ============================================================================

void p1_quad_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                               math::DenseDVec<Real> &values)
{
  values.resize(4);
  values[0] = 0.25 * (1.0 - ref_coords[KSI]) * (1.0 - ref_coords[ETA]);
  values[1] = 0.25 * (1.0 + ref_coords[KSI]) * (1.0 - ref_coords[ETA]);
  values[2] = 0.25 * (1.0 + ref_coords[KSI]) * (1.0 + ref_coords[ETA]);
  values[3] = 0.25 * (1.0 - ref_coords[KSI]) * (1.0 + ref_coords[ETA]);
}

// ============================================================================

void p1_quad_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                     math::DenseDMat<Real> &values)
{
  values(0, KSI) = -0.25 * (1. - ref_coords[ETA]);
  values(1, KSI) = 0.25 * (1. - ref_coords[ETA]);
  values(2, KSI) = 0.25 * (1. + ref_coords[ETA]);
  values(3, KSI) = -0.25 * (1. + ref_coords[ETA]);

  values(0, ETA) = -0.25 * (1. - ref_coords[KSI]);
  values(1, ETA) = -0.25 * (1. + ref_coords[KSI]);
  values(2, ETA) = 0.25 * (1. + ref_coords[KSI]);
  values(3, ETA) = 0.25 * (1. - ref_coords[KSI]);
}

/// ===========================================================================

const Uint p2_quad_Lagrange_sf::poly_order = P2;

// ============================================================================

void p2_quad_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                               math::DenseDVec<Real> &values)
{
  values.resize(9);
  const Real ksi = ref_coords[KSI];
  const Real eta = ref_coords[ETA];

  const Real ksi2    = ksi * ksi;
  const Real eta2    = eta * eta;
  const Real ksi_eta = ksi * eta;

  values[0] = 0.25 * (1.0 - ksi) * (1.0 - eta) * ksi_eta;
  values[1] = -0.25 * (1.0 + ksi) * (1.0 - eta) * ksi_eta;
  values[2] = 0.25 * (1.0 + ksi) * (1.0 + eta) * ksi_eta;
  values[3] = -0.25 * (1.0 - ksi) * (1.0 + eta) * ksi_eta;
  values[4] = -0.5 * (1.0 - ksi2) * (1.0 - eta) * eta;
  values[5] = 0.5 * (1.0 + ksi) * (1.0 - eta2) * ksi;
  values[6] = 0.5 * (1.0 - ksi2) * (1.0 + eta) * eta;
  values[7] = -0.5 * (1.0 - ksi) * (1.0 - eta2) * ksi;
  values[8] = (1.0 - ksi2) * (1.0 - eta2);
}

// ============================================================================

void p2_quad_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                     math::DenseDMat<Real> &values)
{
  const Real ksi      = ref_coords[KSI];
  const Real eta      = ref_coords[ETA];
  const Real ksi2     = ksi * ksi;
  const Real eta2     = eta * eta;
  const Real ksi_eta  = ksi * eta;
  const Real ksi2_eta = ksi2 * eta;
  const Real ksi_eta2 = ksi * eta2;

  values(0, KSI) = 0.25 * (eta - 2. * ksi_eta - eta2 + 2. * ksi_eta2);
  values(1, KSI) = -0.25 * (eta + 2. * ksi_eta - eta2 - 2. * ksi_eta2);
  values(2, KSI) = 0.25 * (eta + 2. * ksi_eta + eta2 + 2. * ksi_eta2);
  values(3, KSI) = -0.25 * (eta - 2. * ksi_eta + eta2 - 2. * ksi_eta2);
  values(4, KSI) = -0.5 * (-2. * ksi_eta + 2. * ksi_eta2);
  values(5, KSI) = 0.5 * (1. - eta2 + 2. * ksi - 2. * ksi_eta2);
  values(6, KSI) = 0.5 * (-2. * ksi_eta - 2. * ksi_eta2);
  values(7, KSI) = -0.5 * (1. - eta2 - 2. * ksi + 2. * ksi_eta2);
  values(8, KSI) = 2. * ksi_eta2 - 2. * ksi;

  values(0, ETA) = 0.25 * (ksi - ksi2 - 2. * ksi_eta + 2. * ksi2_eta);
  values(1, ETA) = -0.25 * (ksi + ksi2 - 2. * ksi_eta - 2. * ksi2_eta);
  values(2, ETA) = 0.25 * (ksi + ksi2 + 2. * ksi_eta + 2. * ksi2_eta);
  values(3, ETA) = -0.25 * (ksi - ksi2 + 2. * ksi_eta - 2. * ksi2_eta);
  values(4, ETA) = -0.5 * (1. - ksi2 - 2. * eta + 2. * ksi2_eta);
  values(5, ETA) = 0.5 * (-2. * ksi_eta - 2. * ksi2_eta);
  values(6, ETA) = 0.5 * (1. - ksi2 + 2. * eta - 2. * ksi2_eta);
  values(7, ETA) = -0.5 * (-2. * ksi_eta + 2. * ksi2_eta);
  values(8, ETA) = 2. * ksi2_eta - 2. * eta;
}

/// ===========================================================================

const Uint p3_quad_Lagrange_sf::poly_order = P3;

void p3_quad_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                               math::DenseDVec<Real> &values)
{
  values.resize(16);
  const Real ksi = ref_coords[KSI];
  const Real eta = ref_coords[ETA];

  const Real onethird = 1.0 / 3.0;

  const Real L0ksi = -9.0 / 16.0 * (ksi + onethird) * (ksi - onethird) * (ksi - 1.0);
  const Real L1ksi = 27.0 / 16.0 * (ksi + 1.0) * (ksi - onethird) * (ksi - 1.0);
  const Real L2ksi = -27.0 / 16.0 * (ksi + 1.0) * (ksi + onethird) * (ksi - 1.0);
  const Real L3ksi = 9.0 / 16.0 * (ksi + 1.0) * (ksi + onethird) * (ksi - onethird);

  const Real L0eta = -9.0 / 16.0 * (eta + onethird) * (eta - onethird) * (eta - 1.0);
  const Real L1eta = 27.0 / 16.0 * (eta + 1.0) * (eta - onethird) * (eta - 1.0);
  const Real L2eta = -27.0 / 16.0 * (eta + 1.0) * (eta + onethird) * (eta - 1.0);
  const Real L3eta = 9.0 / 16.0 * (eta + 1.0) * (eta + onethird) * (eta - onethird);

  values[0]  = L0ksi * L0eta;
  values[1]  = L3ksi * L0eta;
  values[2]  = L3ksi * L3eta;
  values[3]  = L0ksi * L3eta;
  values[4]  = L1ksi * L0eta;
  values[5]  = L2ksi * L0eta;
  values[6]  = L3ksi * L1eta;
  values[7]  = L3ksi * L2eta;
  values[8]  = L2ksi * L3eta;
  values[9]  = L1ksi * L3eta;
  values[10] = L0ksi * L2eta;
  values[11] = L0ksi * L1eta;
  values[12] = L1ksi * L1eta;
  values[13] = L2ksi * L1eta;
  values[14] = L2ksi * L2eta;
  values[15] = L1ksi * L2eta;
}

// ============================================================================

void p3_quad_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                     math::DenseDMat<Real> &values)
{
  const Real ksi = ref_coords[KSI];
  const Real eta = ref_coords[ETA];

  const Real onethird = 1.0 / 3.0;

  const Real L0ksi = -9.0 / 16.0 * (ksi + onethird) * (ksi - onethird) * (ksi - 1.0);
  const Real L1ksi = 27.0 / 16.0 * (ksi + 1.0) * (ksi - onethird) * (ksi - 1.0);
  const Real L2ksi = -27.0 / 16.0 * (ksi + 1.0) * (ksi + onethird) * (ksi - 1.0);
  const Real L3ksi = 9.0 / 16.0 * (ksi + 1.0) * (ksi + onethird) * (ksi - onethird);

  const Real L0eta = -9.0 / 16.0 * (eta + onethird) * (eta - onethird) * (eta - 1.0);
  const Real L1eta = 27.0 / 16.0 * (eta + 1.0) * (eta - onethird) * (eta - 1.0);
  const Real L2eta = -27.0 / 16.0 * (eta + 1.0) * (eta + onethird) * (eta - 1.0);
  const Real L3eta = 9.0 / 16.0 * (eta + 1.0) * (eta + onethird) * (eta - onethird);

  const Real dL0dksi = -9.0 / 16.0 *
                       ((ksi - onethird) * (ksi - 1.0) + (ksi + onethird) * (ksi - 1.0) +
                        (ksi + onethird) * (ksi - onethird));
  const Real dL1dksi =
      27.0 / 16.0 *
      ((ksi - onethird) * (ksi - 1.0) + (ksi + 1.0) * (ksi - 1.0) + (ksi + 1.0) * (ksi - onethird));
  const Real dL2dksi =
      -27.0 / 16.0 *
      ((ksi + onethird) * (ksi - 1.0) + (ksi + 1.0) * (ksi - 1.0) + (ksi + 1.0) * (ksi + onethird));
  const Real dL3dksi = 9.0 / 16.0 *
                       ((ksi + onethird) * (ksi - onethird) + (ksi + 1.0) * (ksi - onethird) +
                        (ksi + 1.0) * (ksi + onethird));

  const Real dL0deta = -9.0 / 16.0 *
                       ((eta - onethird) * (eta - 1.0) + (eta + onethird) * (eta - 1.0) +
                        (eta + onethird) * (eta - onethird));
  const Real dL1deta =
      27.0 / 16.0 *
      ((eta - onethird) * (eta - 1.0) + (eta + 1.0) * (eta - 1.0) + (eta + 1.0) * (eta - onethird));
  const Real dL2deta =
      -27.0 / 16.0 *
      ((eta + onethird) * (eta - 1.0) + (eta + 1.0) * (eta - 1.0) + (eta + 1.0) * (eta + onethird));
  const Real dL3deta = 9.0 / 16.0 *
                       ((eta + onethird) * (eta - onethird) + (eta + 1.0) * (eta - onethird) +
                        (eta + 1.0) * (eta + onethird));

  values(0, KSI)  = dL0dksi * L0eta;
  values(1, KSI)  = dL3dksi * L0eta;
  values(2, KSI)  = dL3dksi * L3eta;
  values(3, KSI)  = dL0dksi * L3eta;
  values(4, KSI)  = dL1dksi * L0eta;
  values(5, KSI)  = dL2dksi * L0eta;
  values(6, KSI)  = dL3dksi * L1eta;
  values(7, KSI)  = dL3dksi * L2eta;
  values(8, KSI)  = dL2dksi * L3eta;
  values(9, KSI)  = dL1dksi * L3eta;
  values(10, KSI) = dL0dksi * L2eta;
  values(11, KSI) = dL0dksi * L1eta;
  values(12, KSI) = dL1dksi * L1eta;
  values(13, KSI) = dL2dksi * L1eta;
  values(14, KSI) = dL2dksi * L2eta;
  values(15, KSI) = dL1dksi * L2eta;

  values(0, ETA)  = L0ksi * dL0deta;
  values(1, ETA)  = L3ksi * dL0deta;
  values(2, ETA)  = L3ksi * dL3deta;
  values(3, ETA)  = L0ksi * dL3deta;
  values(4, ETA)  = L1ksi * dL0deta;
  values(5, ETA)  = L2ksi * dL0deta;
  values(6, ETA)  = L3ksi * dL1deta;
  values(7, ETA)  = L3ksi * dL2deta;
  values(8, ETA)  = L2ksi * dL3deta;
  values(9, ETA)  = L1ksi * dL3deta;
  values(10, ETA) = L0ksi * dL2deta;
  values(11, ETA) = L0ksi * dL1deta;
  values(12, ETA) = L1ksi * dL1deta;
  values(13, ETA) = L2ksi * dL1deta;
  values(14, ETA) = L2ksi * dL2deta;
  values(15, ETA) = L1ksi * dL2deta;
}

/// ===========================================================================

const Uint p4_quad_Lagrange_sf::poly_order = P4;

void p4_quad_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                               math::DenseDVec<Real> &values)
{
  values.resize(25);

  const Real ksi = ref_coords[KSI];
  const Real eta = ref_coords[ETA];

  const Real onesixth = 1.0 / 6.0;

  const Real h0ksi =
      onesixth * (ksi - ksi * ksi - 4. * ksi * ksi * ksi + 4. * ksi * ksi * ksi * ksi);
  const Real h1ksi =
      onesixth * (-8. * ksi + 16. * ksi * ksi + 8. * ksi * ksi * ksi - 16. * ksi * ksi * ksi * ksi);
  const Real h2ksi = onesixth * (6. - 30. * ksi * ksi + 24. * ksi * ksi * ksi * ksi);
  const Real h3ksi =
      onesixth * (8. * ksi + 16. * ksi * ksi - 8. * ksi * ksi * ksi - 16. * ksi * ksi * ksi * ksi);
  const Real h4ksi =
      onesixth * (-ksi - ksi * ksi + 4. * ksi * ksi * ksi + 4. * ksi * ksi * ksi * ksi);

  const Real h0eta =
      onesixth * (eta - eta * eta - 4. * eta * eta * eta + 4. * eta * eta * eta * eta);
  const Real h1eta =
      onesixth * (-8. * eta + 16. * eta * eta + 8. * eta * eta * eta - 16. * eta * eta * eta * eta);
  const Real h2eta = onesixth * (6. - 30. * eta * eta + 24. * eta * eta * eta * eta);
  const Real h3eta =
      onesixth * (8. * eta + 16. * eta * eta - 8. * eta * eta * eta - 16. * eta * eta * eta * eta);
  const Real h4eta =
      onesixth * (-eta - eta * eta + 4. * eta * eta * eta + 4. * eta * eta * eta * eta);

  values[0]  = h0ksi * h0eta;
  values[1]  = h4ksi * h0eta;
  values[2]  = h4ksi * h4eta;
  values[3]  = h0ksi * h4eta;
  values[4]  = h1ksi * h0eta;
  values[5]  = h2ksi * h0eta;
  values[6]  = h3ksi * h0eta;
  values[7]  = h4ksi * h1eta;
  values[8]  = h4ksi * h2eta;
  values[9]  = h4ksi * h3eta;
  values[10] = h3ksi * h4eta;
  values[11] = h2ksi * h4eta;
  values[12] = h1ksi * h4eta;
  values[13] = h0ksi * h3eta;
  values[14] = h0ksi * h2eta;
  values[15] = h0ksi * h1eta;
  values[16] = h1ksi * h1eta;
  values[17] = h3ksi * h1eta;
  values[18] = h3ksi * h3eta;
  values[19] = h1ksi * h3eta;
  values[20] = h2ksi * h1eta;
  values[21] = h3ksi * h2eta;
  values[22] = h2ksi * h3eta;
  values[23] = h1ksi * h2eta;
  values[24] = h2ksi * h2eta;
}

// ============================================================================

void p4_quad_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                     math::DenseDMat<Real> &values)
{
  const Real ksi = ref_coords[KSI];
  const Real eta = ref_coords[ETA];

  const Real onesixth = 1.0 / 6.0;

  const Real h0ksi =
      onesixth * (ksi - ksi * ksi - 4. * ksi * ksi * ksi + 4. * ksi * ksi * ksi * ksi);
  const Real h1ksi =
      onesixth * (-8. * ksi + 16. * ksi * ksi + 8. * ksi * ksi * ksi - 16. * ksi * ksi * ksi * ksi);
  const Real h2ksi = onesixth * (6. - 30. * ksi * ksi + 24. * ksi * ksi * ksi * ksi);
  const Real h3ksi =
      onesixth * (8. * ksi + 16. * ksi * ksi - 8. * ksi * ksi * ksi - 16. * ksi * ksi * ksi * ksi);
  const Real h4ksi =
      onesixth * (-ksi - ksi * ksi + 4. * ksi * ksi * ksi + 4. * ksi * ksi * ksi * ksi);

  const Real dh0dksi = onesixth * (1. - 2. * ksi - 12. * ksi * ksi + 16. * ksi * ksi * ksi);
  const Real dh1dksi = onesixth * (-8. + 32. * ksi + 24. * ksi * ksi - 64. * ksi * ksi * ksi);
  const Real dh2dksi = onesixth * (-60. * ksi + 96. * ksi * ksi * ksi);
  const Real dh3dksi = onesixth * (8. + 32. * ksi - 24. * ksi * ksi - 64. * ksi * ksi * ksi);
  const Real dh4dksi = onesixth * (-1. - 2. * ksi + 12. * ksi * ksi + 16. * ksi * ksi * ksi);

  const Real h0eta =
      onesixth * (eta - eta * eta - 4. * eta * eta * eta + 4. * eta * eta * eta * eta);
  const Real h1eta =
      onesixth * (-8. * eta + 16. * eta * eta + 8. * eta * eta * eta - 16. * eta * eta * eta * eta);
  const Real h2eta = onesixth * (6. - 30. * eta * eta + 24. * eta * eta * eta * eta);
  const Real h3eta =
      onesixth * (8. * eta + 16. * eta * eta - 8. * eta * eta * eta - 16. * eta * eta * eta * eta);
  const Real h4eta =
      onesixth * (-eta - eta * eta + 4. * eta * eta * eta + 4. * eta * eta * eta * eta);

  const Real dh0deta = onesixth * (1. - 2. * eta - 12. * eta * eta + 16. * eta * eta * eta);
  const Real dh1deta = onesixth * (-8. + 32. * eta + 24. * eta * eta - 64. * eta * eta * eta);
  const Real dh2deta = onesixth * (-60. * eta + 96. * eta * eta * eta);
  const Real dh3deta = onesixth * (8. + 32. * eta - 24. * eta * eta - 64. * eta * eta * eta);
  const Real dh4deta = onesixth * (-1. - 2. * eta + 12. * eta * eta + 16. * eta * eta * eta);

  values(0, KSI)  = dh0dksi * h0eta;
  values(1, KSI)  = dh4dksi * h0eta;
  values(2, KSI)  = dh4dksi * h4eta;
  values(3, KSI)  = dh0dksi * h4eta;
  values(4, KSI)  = dh1dksi * h0eta;
  values(5, KSI)  = dh2dksi * h0eta;
  values(6, KSI)  = dh3dksi * h0eta;
  values(7, KSI)  = dh4dksi * h1eta;
  values(8, KSI)  = dh4dksi * h2eta;
  values(9, KSI)  = dh4dksi * h3eta;
  values(10, KSI) = dh3dksi * h4eta;
  values(11, KSI) = dh2dksi * h4eta;
  values(12, KSI) = dh1dksi * h4eta;
  values(13, KSI) = dh0dksi * h3eta;
  values(14, KSI) = dh0dksi * h2eta;
  values(15, KSI) = dh0dksi * h1eta;
  values(16, KSI) = dh1dksi * h1eta;
  values(17, KSI) = dh3dksi * h1eta;
  values(18, KSI) = dh3dksi * h3eta;
  values(19, KSI) = dh1dksi * h3eta;
  values(20, KSI) = dh2dksi * h1eta;
  values(21, KSI) = dh3dksi * h2eta;
  values(22, KSI) = dh2dksi * h3eta;
  values(23, KSI) = dh1dksi * h2eta;
  values(24, KSI) = dh2dksi * h2eta;

  values(0, ETA)  = h0ksi * dh0deta;
  values(1, ETA)  = h4ksi * dh0deta;
  values(2, ETA)  = h4ksi * dh4deta;
  values(3, ETA)  = h0ksi * dh4deta;
  values(4, ETA)  = h1ksi * dh0deta;
  values(5, ETA)  = h2ksi * dh0deta;
  values(6, ETA)  = h3ksi * dh0deta;
  values(7, ETA)  = h4ksi * dh1deta;
  values(8, ETA)  = h4ksi * dh2deta;
  values(9, ETA)  = h4ksi * dh3deta;
  values(10, ETA) = h3ksi * dh4deta;
  values(11, ETA) = h2ksi * dh4deta;
  values(12, ETA) = h1ksi * dh4deta;
  values(13, ETA) = h0ksi * dh3deta;
  values(14, ETA) = h0ksi * dh2deta;
  values(15, ETA) = h0ksi * dh1deta;
  values(16, ETA) = h1ksi * dh1deta;
  values(17, ETA) = h3ksi * dh1deta;
  values(18, ETA) = h3ksi * dh3deta;
  values(19, ETA) = h1ksi * dh3deta;
  values(20, ETA) = h2ksi * dh1deta;
  values(21, ETA) = h3ksi * dh2deta;
  values(22, ETA) = h2ksi * dh3deta;
  values(23, ETA) = h1ksi * dh2deta;
  values(24, ETA) = h2ksi * dh2deta;
}

/// ===========================================================================

} // namespace utest_fixture

} // namespace mesh

} // namespace pdekit
