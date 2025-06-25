#include "mesh/shape_function/TetraP1Lagrange3D.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

/// =======================================================================================================

const Uint ShapeFunctionT<Tetra, P1, Lagrange, _3D>::NbNodes = 4;
const Uint ShapeFunctionT<Tetra, P1, Lagrange, _3D>::Order   = P1;

/// =======================================================================================================

ShapeFunctionT<Tetra, P1, Lagrange, _3D>::ShapeFunctionT() : LagrangeShapeFunction<Tetra, P1, _3D>()
{
}

/// =======================================================================================================

ShapeFunctionT<Tetra, P1, Lagrange, _3D>::~ShapeFunctionT()
{
}

/// =======================================================================================================

const Uint ShapeFunctionT<Tetra, P1, Lagrange, _3D>::nb_dof() const
{
  return 4;
}

/// =======================================================================================================

void ShapeFunctionT<Tetra, P1, Lagrange, _3D>::ref_coords(TCoord &coords) const
{
  coords.resize(NbNodes, TopoDim);

  coords(0, KSI) = -1.0;
  coords(0, ETA) = -1.0;
  coords(0, ZTA) = -1.0;

  coords(1, KSI) = 1.0;
  coords(1, ETA) = -1.0;
  coords(1, ZTA) = -1.0;

  coords(2, KSI) = -1.0;
  coords(2, ETA) = 1.0;
  coords(2, ZTA) = -1.0;

  coords(3, KSI) = -1.0;
  coords(3, ETA) = -1.0;
  coords(3, ZTA) = 1.0;
}

/// =======================================================================================================

void ShapeFunctionT<Tetra, P1, Lagrange, _3D>::compute_ref_values(
    const math::VectorBlock<Real> &ref_coords, math::DynamicVector<Real> &values)
{
  values[0] = -0.5 * (1. + ref_coords[KSI] + ref_coords[ETA] + ref_coords[ZTA]);
  values[1] = 0.5 * (1.0 + ref_coords[KSI]);
  values[2] = 0.5 * (1.0 + ref_coords[ETA]);
  values[3] = 0.5 * (1.0 + ref_coords[ZTA]);
}

/// =======================================================================================================

void ShapeFunctionT<Tetra, P1, Lagrange, _3D>::compute_ref_derivatives(
    const math::VectorBlock<Real> &ref_coords, math::DynamicMatrix<Real> &values)
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

/// =======================================================================================================

// const ElementType&
// ShapeFunctionT<Tetra::value,P1::value,SFType::Lagrange,_3D>::face_type(const
// Uint face)
//{

//}

/// =======================================================================================================

void ShapeFunctionT<Tetra, P1, Lagrange, _3D>::print() const
{
}

/// =======================================================================================================

} // namespace sf

} // namespace mesh

} // namespace pdekit
