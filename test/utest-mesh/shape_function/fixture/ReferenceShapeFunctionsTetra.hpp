#ifndef PDEKIT_Mesh_UTest_Fixture_Reference_Shape_Functions_Tetra_hpp
#define PDEKIT_Mesh_UTest_Fixture_Reference_Shape_Functions_Tetra_hpp

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

// ----------------------------------------------------------------------------

// P1 Lagrange shape functions on tetrahedra
struct p1_tetra_Lagrange_sf
{
  static const Uint poly_order;

  static void eval(const math::DenseDVec<Real> &ref_coords, math::DenseDVec<Real> &values);

  static void eval_deriv(const math::DenseDVec<Real> &ref_coords, math::DenseDMat<Real> &values);
};

// ----------------------------------------------------------------------------

struct p2_tetra_Lagrange_sf
{
  static const Uint poly_order;

  static void eval(const math::DenseDVec<Real> &ref_coords, math::DenseDVec<Real> &values);

  static void eval_deriv(const math::DenseDVec<Real> &ref_coords, math::DenseDMat<Real> &values);
};

// ----------------------------------------------------------------------------

struct p3_tetra_Lagrange_sf
{
  static const Uint poly_order;

  static void eval(const math::DenseDVec<Real> &ref_coords, math::DenseDVec<Real> &values);

  static void eval_deriv(const math::DenseDVec<Real> &ref_coords, math::DenseDMat<Real> &values);
};

// ----------------------------------------------------------------------------

} // namespace utest_fixture

} // namespace mesh

} // namespace pdekit

#endif
