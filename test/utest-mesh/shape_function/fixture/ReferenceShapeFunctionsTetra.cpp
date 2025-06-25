#include <array>

#include "test/utest-mesh/shape_function/fixture/ReferenceShapeFunctionsTetra.hpp"

namespace pdekit
{

namespace mesh
{

namespace utest_fixture
{

// ----------------------------------------------------------------------------

const Uint p1_tetra_Lagrange_sf::poly_order = P1;

// ----------------------------------------------------------------------------

void p1_tetra_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                                math::DenseDVec<Real> &values)
{
  values.resize(4);
  values[0] = -.5 * (1. + ref_coords[KSI] + ref_coords[ETA] + ref_coords[ZTA]);
  values[1] = 0.5 * (1. + ref_coords[KSI]);
  values[2] = 0.5 * (1. + ref_coords[ETA]);
  values[3] = 0.5 * (1. + ref_coords[ZTA]);
}

// ----------------------------------------------------------------------------

void p1_tetra_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                      math::DenseDMat<Real> &values)
{
  values(0, KSI) = -0.5;
  values(0, ETA) = -0.5;
  values(0, ZTA) = -0.5;
  values(1, KSI) = 0.5;
  values(1, ETA) = 0.0;
  values(1, ZTA) = 0.0;
  values(2, KSI) = 0.0;
  values(2, ETA) = 0.5;
  values(2, ZTA) = 0.0;
  values(3, KSI) = 0.0;
  values(3, ETA) = 0.0;
  values(3, ZTA) = 0.5;
}

// ----------------------------------------------------------------------------

const Uint p2_tetra_Lagrange_sf::poly_order = P2;

// ----------------------------------------------------------------------------

void p2_tetra_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                                math::DenseDVec<Real> &values)
{
  values.resize(10);
  const Real L0 = -.5 * (1. + ref_coords[KSI] + ref_coords[ETA] + ref_coords[ZTA]);
  const Real L1 = 0.5 * (1. + ref_coords[KSI]);
  const Real L2 = 0.5 * (1. + ref_coords[ETA]);
  const Real L3 = 0.5 * (1. + ref_coords[ZTA]);

  values[0] = (2. * L0 - 1.) * L0;
  values[1] = (2. * L1 - 1.) * L1;
  values[2] = (2. * L2 - 1.) * L2;
  values[3] = (2. * L3 - 1.) * L3;
  values[4] = 4. * L0 * L1;
  values[5] = 4. * L1 * L2;
  values[6] = 4. * L0 * L2;
  values[7] = 4. * L0 * L3;
  values[8] = 4. * L2 * L3;
  values[9] = 4. * L1 * L3;
}

// ----------------------------------------------------------------------------

void p2_tetra_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                      math::DenseDMat<Real> &values)
{
  values(0, XI0) = 1.5 + ref_coords[XI0] + ref_coords[XI1] + ref_coords[XI2];
  values(0, XI1) = 1.5 + ref_coords[XI0] + ref_coords[XI1] + ref_coords[XI2];
  values(0, XI2) = 1.5 + ref_coords[XI0] + ref_coords[XI1] + ref_coords[XI2];

  values(1, XI0) = 0.5 + ref_coords[XI0];
  values(1, XI1) = 0.0;
  values(1, XI2) = 0.0;

  values(2, XI0) = 0.0;
  values(2, XI1) = 0.5 + ref_coords[XI1];
  values(2, XI2) = 0.0;

  values(3, XI0) = 0.0;
  values(3, XI1) = 0.0;
  values(3, XI2) = 0.5 + ref_coords[XI2];

  values(4, XI0) = -(2. + 2. * ref_coords[XI0] + ref_coords[XI1] + ref_coords[XI2]);
  values(4, XI1) = -(1. + ref_coords[XI0]);
  values(4, XI2) = -(1. + ref_coords[XI0]);

  values(5, XI0) = 1. + ref_coords[XI1];
  values(5, XI1) = 1. + ref_coords[XI0];
  values(5, XI2) = 0.0;

  values(6, XI0) = -(1. + ref_coords[XI1]);
  values(6, XI1) = -(2. + ref_coords[XI0] + 2. * ref_coords[XI1] + ref_coords[XI2]);
  values(6, XI2) = -(1. + ref_coords[XI1]);

  values(7, XI0) = -(1. + ref_coords[XI2]);
  values(7, XI1) = -(1. + ref_coords[XI2]);
  values(7, XI2) = -(2. + ref_coords[XI0] + ref_coords[XI1] + 2. * ref_coords[XI2]);

  values(8, XI0) = 0.0;
  values(8, XI1) = 1. + ref_coords[XI2];
  values(8, XI2) = 1. + ref_coords[XI1];

  values(9, XI0) = 1. + ref_coords[XI2];
  values(9, XI1) = 0.0;
  values(9, XI2) = 1. + ref_coords[XI0];
}

// ----------------------------------------------------------------------------

const Uint p3_tetra_Lagrange_sf::poly_order = P3;

// ----------------------------------------------------------------------------

void p3_tetra_Lagrange_sf::eval(const math::DenseDVec<Real> &ref_coords,
                                math::DenseDVec<Real> &values)
{
  values.resize(20);
  const Real L0 = -.5 * (1. + ref_coords[KSI] + ref_coords[ETA] + ref_coords[ZTA]);
  const Real L1 = 0.5 * (1. + ref_coords[KSI]);
  const Real L2 = 0.5 * (1. + ref_coords[ETA]);
  const Real L3 = 0.5 * (1. + ref_coords[ZTA]);

  values[0] = 0.5 * (3. * L0 - 1.) * (3. * L0 - 2.) * L0;
  values[1] = 0.5 * (3. * L1 - 1.) * (3. * L1 - 2.) * L1;
  values[2] = 0.5 * (3. * L2 - 1.) * (3. * L2 - 2.) * L2;
  values[3] = 0.5 * (3. * L3 - 1.) * (3. * L3 - 2.) * L3;

  // Edge 0-1
  values[4] = 4.5 * L0 * L1 * (3. * L0 - 1.);
  values[5] = 4.5 * L0 * L1 * (3. * L1 - 1.);

  // Edge 1-2
  values[6] = 4.5 * L1 * L2 * (3. * L1 - 1.);
  values[7] = 4.5 * L1 * L2 * (3. * L2 - 1.);

  // Edge 2-0
  values[8] = 4.5 * L2 * L0 * (3. * L2 - 1.);
  values[9] = 4.5 * L2 * L0 * (3. * L0 - 1.);

  // Edge 3-0
  values[10] = 4.5 * L3 * L0 * (3. * L3 - 1.);
  values[11] = 4.5 * L3 * L0 * (3. * L0 - 1.);

  // Edge 3-2
  values[12] = 4.5 * L3 * L2 * (3. * L3 - 1.);
  values[13] = 4.5 * L3 * L2 * (3. * L2 - 1.);

  // Edge 3-1
  values[14] = 4.5 * L3 * L1 * (3. * L3 - 1.);
  values[15] = 4.5 * L3 * L1 * (3. * L1 - 1.);

  // Face 0-2-1
  values[16] = 27. * L0 * L2 * L1;

  // Face 0-1-3
  values[17] = 27. * L0 * L1 * L3;

  // Face 0-3-2
  values[18] = 27. * L0 * L3 * L2;

  // Face 3-1-2
  values[19] = 27. * L3 * L1 * L2;
}

// ----------------------------------------------------------------------------

void p3_tetra_Lagrange_sf::eval_deriv(const math::DenseDVec<Real> &ref_coords,
                                      math::DenseDMat<Real> &values)
{
  const std::array<Real, 4> L{-.5 * (1. + ref_coords[KSI] + ref_coords[ETA] + ref_coords[ZTA]),
                              0.5 * (1. + ref_coords[KSI]), 0.5 * (1. + ref_coords[ETA]),
                              0.5 * (1. + ref_coords[ZTA])};

  const std::array<Real, 4> dLdXi0{-0.5, 0.5, 0.0, 0.0};
  const std::array<Real, 4> dLdXi1{-0.5, 0.0, 0.5, 0.0};
  const std::array<Real, 4> dLdXi2{-0.5, 0.0, 0.0, 0.5};

  values(0, XI0) =
      0.5 * ((3. * L[0] - 1.) * (3. * L[0] - 2.) * dLdXi0[0] +
             3. * L[0] * (3. * L[0] - 2.) * dLdXi0[0] + 3. * L[0] * (3. * L[0] - 1.) * dLdXi0[0]);

  values(0, XI1) =
      0.5 * ((3. * L[0] - 1.) * (3. * L[0] - 2.) * dLdXi1[0] +
             3. * L[0] * (3. * L[0] - 2.) * dLdXi1[0] + 3. * L[0] * (3. * L[0] - 1.) * dLdXi1[0]);

  values(0, XI2) =
      0.5 * ((3. * L[0] - 1.) * (3. * L[0] - 2.) * dLdXi2[0] +
             3. * L[0] * (3. * L[0] - 2.) * dLdXi2[0] + 3. * L[0] * (3. * L[0] - 1.) * dLdXi2[0]);

  values(1, XI0) =
      0.5 * ((3. * L[1] - 1.) * (3. * L[1] - 2.) * dLdXi0[1] +
             3. * L[1] * (3. * L[1] - 2.) * dLdXi0[1] + 3. * L[1] * (3. * L[1] - 1.) * dLdXi0[1]);

  values(1, XI1) =
      0.5 * ((3. * L[1] - 1.) * (3. * L[1] - 2.) * dLdXi1[1] +
             3. * L[1] * (3. * L[1] - 2.) * dLdXi1[1] + 3. * L[1] * (3. * L[1] - 1.) * dLdXi1[1]);

  values(1, XI2) =
      0.5 * ((3. * L[1] - 1.) * (3. * L[1] - 2.) * dLdXi2[1] +
             3. * L[1] * (3. * L[1] - 2.) * dLdXi2[1] + 3. * L[1] * (3. * L[1] - 1.) * dLdXi2[1]);

  values(2, XI0) =
      0.5 * ((3. * L[2] - 1.) * (3. * L[2] - 2.) * dLdXi0[2] +
             3. * L[2] * (3. * L[2] - 2.) * dLdXi0[2] + 3. * L[2] * (3. * L[2] - 1.) * dLdXi0[2]);

  values(2, XI1) =
      0.5 * ((3. * L[2] - 1.) * (3. * L[2] - 2.) * dLdXi1[2] +
             3. * L[2] * (3. * L[2] - 2.) * dLdXi1[2] + 3. * L[2] * (3. * L[2] - 1.) * dLdXi1[2]);

  values(2, XI2) =
      0.5 * ((3. * L[2] - 1.) * (3. * L[2] - 2.) * dLdXi2[2] +
             3. * L[2] * (3. * L[2] - 2.) * dLdXi2[2] + 3. * L[2] * (3. * L[2] - 1.) * dLdXi2[2]);

  values(3, XI0) =
      0.5 * ((3. * L[3] - 1.) * (3. * L[3] - 2.) * dLdXi0[3] +
             3. * L[3] * (3. * L[3] - 2.) * dLdXi0[3] + 3. * L[3] * (3. * L[3] - 1.) * dLdXi0[3]);

  values(3, XI1) =
      0.5 * ((3. * L[3] - 1.) * (3. * L[3] - 2.) * dLdXi1[3] +
             3. * L[3] * (3. * L[3] - 2.) * dLdXi1[3] + 3. * L[3] * (3. * L[3] - 1.) * dLdXi1[3]);

  values(3, XI2) =
      0.5 * ((3. * L[3] - 1.) * (3. * L[3] - 2.) * dLdXi2[3] +
             3. * L[3] * (3. * L[3] - 2.) * dLdXi2[3] + 3. * L[3] * (3. * L[3] - 1.) * dLdXi2[3]);

  values(4, XI0) = 4.5 * (L[1] * (3. * L[0] - 1.) * dLdXi0[0] + 3. * L[0] * L[1] * dLdXi0[0] +
                          L[0] * (3. * L[0] - 1.) * dLdXi0[1]);

  values(4, XI1) = 4.5 * (L[1] * (3. * L[0] - 1.) * dLdXi1[0] + 3. * L[0] * L[1] * dLdXi1[0] +
                          L[0] * (3. * L[0] - 1.) * dLdXi1[1]);

  values(4, XI2) = 4.5 * (L[1] * (3. * L[0] - 1.) * dLdXi2[0] + 3. * L[0] * L[1] * dLdXi2[0] +
                          L[0] * (3. * L[0] - 1.) * dLdXi2[1]);

  values(5, XI0) = 4.5 * (L[0] * (3. * L[1] - 1.) * dLdXi0[1] + 3. * L[1] * L[0] * dLdXi0[1] +
                          L[1] * (3. * L[1] - 1.) * dLdXi0[0]);

  values(5, XI1) = 4.5 * (L[0] * (3. * L[1] - 1.) * dLdXi1[1] + 3. * L[1] * L[0] * dLdXi1[1] +
                          L[1] * (3. * L[1] - 1.) * dLdXi1[0]);

  values(5, XI2) = 4.5 * (L[0] * (3. * L[1] - 1.) * dLdXi2[1] + 3. * L[1] * L[0] * dLdXi2[1] +
                          L[1] * (3. * L[1] - 1.) * dLdXi2[0]);

  values(6, XI0) = 4.5 * (L[2] * (3. * L[1] - 1.) * dLdXi0[1] + 3. * L[1] * L[2] * dLdXi0[1] +
                          L[1] * (3. * L[1] - 1.) * dLdXi0[2]);

  values(6, XI1) = 4.5 * (L[2] * (3. * L[1] - 1.) * dLdXi1[1] + 3. * L[1] * L[2] * dLdXi1[1] +
                          L[1] * (3. * L[1] - 1.) * dLdXi1[2]);

  values(6, XI2) = 4.5 * (L[2] * (3. * L[1] - 1.) * dLdXi2[1] + 3. * L[1] * L[2] * dLdXi2[1] +
                          L[1] * (3. * L[1] - 1.) * dLdXi2[2]);

  values(7, XI0) = 4.5 * (L[1] * (3. * L[2] - 1.) * dLdXi0[2] + 3. * L[2] * L[1] * dLdXi0[2] +
                          L[2] * (3. * L[2] - 1.) * dLdXi0[1]);

  values(7, XI1) = 4.5 * (L[1] * (3. * L[2] - 1.) * dLdXi1[2] + 3. * L[2] * L[1] * dLdXi1[2] +
                          L[2] * (3. * L[2] - 1.) * dLdXi1[1]);

  values(7, XI2) = 4.5 * (L[1] * (3. * L[2] - 1.) * dLdXi2[2] + 3. * L[2] * L[1] * dLdXi2[2] +
                          L[2] * (3. * L[2] - 1.) * dLdXi2[1]);

  values(8, XI0) = 4.5 * (L[0] * (3. * L[2] - 1.) * dLdXi0[2] + 3. * L[2] * L[0] * dLdXi0[2] +
                          L[2] * (3. * L[2] - 1.) * dLdXi0[0]);

  values(8, XI1) = 4.5 * (L[0] * (3. * L[2] - 1.) * dLdXi1[2] + 3. * L[2] * L[0] * dLdXi1[2] +
                          L[2] * (3. * L[2] - 1.) * dLdXi1[0]);

  values(8, XI2) = 4.5 * (L[0] * (3. * L[2] - 1.) * dLdXi2[2] + 3. * L[2] * L[0] * dLdXi2[2] +
                          L[2] * (3. * L[2] - 1.) * dLdXi2[0]);

  values(9, XI0) = 4.5 * (L[2] * (3. * L[0] - 1.) * dLdXi0[0] + 3. * L[2] * L[0] * dLdXi0[0] +
                          L[0] * (3. * L[0] - 1.) * dLdXi0[2]);

  values(9, XI1) = 4.5 * (L[2] * (3. * L[0] - 1.) * dLdXi1[0] + 3. * L[2] * L[0] * dLdXi1[0] +
                          L[0] * (3. * L[0] - 1.) * dLdXi1[2]);

  values(9, XI2) = 4.5 * (L[2] * (3. * L[0] - 1.) * dLdXi2[0] + 3. * L[2] * L[0] * dLdXi2[0] +
                          L[0] * (3. * L[0] - 1.) * dLdXi2[2]);

  values(10, XI0) = 4.5 * (L[0] * (3. * L[3] - 1.) * dLdXi0[3] + 3. * L[3] * L[0] * dLdXi0[3] +
                           L[3] * (3. * L[3] - 1.) * dLdXi0[0]);

  values(10, XI1) = 4.5 * (L[0] * (3. * L[3] - 1.) * dLdXi1[3] + 3. * L[3] * L[0] * dLdXi1[3] +
                           L[3] * (3. * L[3] - 1.) * dLdXi1[0]);

  values(10, XI2) = 4.5 * (L[0] * (3. * L[3] - 1.) * dLdXi2[3] + 3. * L[3] * L[0] * dLdXi2[3] +
                           L[3] * (3. * L[3] - 1.) * dLdXi2[0]);

  values(11, XI0) = 4.5 * (L[3] * (3. * L[0] - 1.) * dLdXi0[0] + 3. * L[3] * L[0] * dLdXi0[0] +
                           L[0] * (3. * L[0] - 1.) * dLdXi0[3]);

  values(11, XI1) = 4.5 * (L[3] * (3. * L[0] - 1.) * dLdXi1[0] + 3. * L[3] * L[0] * dLdXi1[0] +
                           L[0] * (3. * L[0] - 1.) * dLdXi1[3]);

  values(11, XI2) = 4.5 * (L[3] * (3. * L[0] - 1.) * dLdXi2[0] + 3. * L[3] * L[0] * dLdXi2[0] +
                           L[0] * (3. * L[0] - 1.) * dLdXi2[3]);

  values(12, XI0) = 4.5 * (L[2] * (3. * L[3] - 1.) * dLdXi0[3] + 3. * L[3] * L[2] * dLdXi0[3] +
                           L[3] * (3. * L[3] - 1.) * dLdXi0[2]);

  values(12, XI1) = 4.5 * (L[2] * (3. * L[3] - 1.) * dLdXi1[3] + 3. * L[3] * L[2] * dLdXi1[3] +
                           L[3] * (3. * L[3] - 1.) * dLdXi1[2]);

  values(12, XI2) = 4.5 * (L[2] * (3. * L[3] - 1.) * dLdXi2[3] + 3. * L[3] * L[2] * dLdXi2[3] +
                           L[3] * (3. * L[3] - 1.) * dLdXi2[2]);

  values(13, XI0) = 4.5 * (L[3] * (3. * L[2] - 1.) * dLdXi0[2] + 3. * L[2] * L[3] * dLdXi0[2] +
                           L[2] * (3. * L[2] - 1.) * dLdXi0[3]);

  values(13, XI1) = 4.5 * (L[3] * (3. * L[2] - 1.) * dLdXi1[2] + 3. * L[2] * L[3] * dLdXi1[2] +
                           L[2] * (3. * L[2] - 1.) * dLdXi1[3]);

  values(13, XI2) = 4.5 * (L[3] * (3. * L[2] - 1.) * dLdXi2[2] + 3. * L[2] * L[3] * dLdXi2[2] +
                           L[2] * (3. * L[2] - 1.) * dLdXi2[3]);

  values(14, XI0) = 4.5 * (L[1] * (3. * L[3] - 1.) * dLdXi0[3] + 3. * L[3] * L[1] * dLdXi0[3] +
                           L[3] * (3. * L[3] - 1.) * dLdXi0[1]);

  values(14, XI1) = 4.5 * (L[1] * (3. * L[3] - 1.) * dLdXi1[3] + 3. * L[3] * L[1] * dLdXi1[3] +
                           L[3] * (3. * L[3] - 1.) * dLdXi1[1]);

  values(14, XI2) = 4.5 * (L[1] * (3. * L[3] - 1.) * dLdXi2[3] + 3. * L[3] * L[1] * dLdXi2[3] +
                           L[3] * (3. * L[3] - 1.) * dLdXi2[1]);

  values(15, XI0) = 4.5 * (L[3] * (3. * L[1] - 1.) * dLdXi0[1] + 3. * L[1] * L[3] * dLdXi0[1] +
                           L[1] * (3. * L[1] - 1.) * dLdXi0[3]);

  values(15, XI1) = 4.5 * (L[3] * (3. * L[1] - 1.) * dLdXi1[1] + 3. * L[1] * L[3] * dLdXi1[1] +
                           L[1] * (3. * L[1] - 1.) * dLdXi1[3]);

  values(15, XI2) = 4.5 * (L[3] * (3. * L[1] - 1.) * dLdXi2[1] + 3. * L[1] * L[3] * dLdXi2[1] +
                           L[1] * (3. * L[1] - 1.) * dLdXi2[3]);

  values(16, XI0) =
      27.0 * (L[1] * L[2] * dLdXi0[0] + L[2] * L[0] * dLdXi0[1] + L[0] * L[1] * dLdXi0[2]);

  values(16, XI1) =
      27.0 * (L[1] * L[2] * dLdXi1[0] + L[2] * L[0] * dLdXi1[1] + L[0] * L[1] * dLdXi1[2]);

  values(16, XI2) =
      27.0 * (L[1] * L[2] * dLdXi2[0] + L[2] * L[0] * dLdXi2[1] + L[0] * L[1] * dLdXi2[2]);

  values(17, XI0) =
      27.0 * (L[1] * L[3] * dLdXi0[0] + L[3] * L[0] * dLdXi0[1] + L[0] * L[1] * dLdXi0[3]);

  values(17, XI1) =
      27.0 * (L[1] * L[3] * dLdXi1[0] + L[3] * L[0] * dLdXi1[1] + L[0] * L[1] * dLdXi1[3]);

  values(17, XI2) =
      27.0 * (L[1] * L[3] * dLdXi2[0] + L[3] * L[0] * dLdXi2[1] + L[0] * L[1] * dLdXi2[3]);

  values(18, XI0) =
      27.0 * (L[2] * L[3] * dLdXi0[0] + L[3] * L[0] * dLdXi0[2] + L[0] * L[2] * dLdXi0[3]);

  values(18, XI1) =
      27.0 * (L[2] * L[3] * dLdXi1[0] + L[3] * L[0] * dLdXi1[2] + L[0] * L[2] * dLdXi1[3]);

  values(18, XI2) =
      27.0 * (L[2] * L[3] * dLdXi2[0] + L[3] * L[0] * dLdXi2[2] + L[0] * L[2] * dLdXi2[3]);

  values(19, XI0) =
      27.0 * (L[2] * L[3] * dLdXi0[1] + L[3] * L[1] * dLdXi0[2] + L[1] * L[2] * dLdXi0[3]);

  values(19, XI1) =
      27.0 * (L[2] * L[3] * dLdXi1[1] + L[3] * L[1] * dLdXi1[2] + L[1] * L[2] * dLdXi1[3]);

  values(19, XI2) =
      27.0 * (L[2] * L[3] * dLdXi2[1] + L[3] * L[1] * dLdXi2[2] + L[1] * L[2] * dLdXi2[3]);
}

// ----------------------------------------------------------------------------

} // namespace utest_fixture

} // namespace mesh

} // namespace pdekit
