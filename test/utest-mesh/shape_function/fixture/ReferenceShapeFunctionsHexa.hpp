#ifndef PDEKIT_Mesh_UTest_Fixture_Reference_Shape_Functions_Hexa_hpp
#define PDEKIT_Mesh_UTest_Fixture_Reference_Shape_Functions_Hexa_hpp

#include "common/Constants.hpp"
#include "common/PDEKit.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"

namespace pdekit
{

namespace mesh
{

namespace utest_fixture
{

// Hard-coded computation of Lagrange P1 shape functions on triangles
struct p1_hexa_Lagrange_sf
{
  static const Uint poly_order;

  static void eval(const math::DenseDVec<Real> &ref_coords, math::DenseDVec<Real> &values);

  static void eval_deriv(const math::DenseDVec<Real> &ref_coords, math::DenseDMat<Real> &values);
};

/// ===========================================================================

// P2 Lagrange shape functions on triangles
struct p2_hexa_Lagrange_sf
{
  static const Uint poly_order;

  static void eval(const math::DenseDVec<Real> &ref_coords, math::DenseDVec<Real> &values);

  static void eval_deriv(const math::DenseDVec<Real> &ref_coords, math::DenseDMat<Real> &values);
};

/// ===========================================================================

// P3 Lagrange shape functions on triangles
struct p3_hexa_Lagrange_sf
{
  static const Uint poly_order;

  static void eval(const math::DenseDVec<Real> &ref_coords, math::DenseDVec<Real> &values);

  static void eval_deriv(const math::DenseDVec<Real> &ref_coords, math::DenseDMat<Real> &values);
};

/// ===========================================================================

// P4 Lagrange shape functions on triangles
struct p4_hexa_Lagrange_sf
{
  static const Uint poly_order;

  static void eval(const math::DenseDVec<Real> &ref_coords, math::DenseDVec<Real> &values);

  static void eval_deriv(const math::DenseDVec<Real> &ref_coords, math::DenseDMat<Real> &values);
};

/// ===========================================================================

} // namespace utest_fixture

} // namespace mesh

} // namespace pdekit

#endif
