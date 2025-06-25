#ifndef Tetra_P1_Lagrange_3D_hpp
#define Tetra_P1_Lagrange_3D_hpp

#include "mesh/shape_function/LagrangeShapeFunction.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

template <>
class ShapeFunctionT<Tetra, P1, Lagrange, _3D> : public LagrangeShapeFunction<Tetra, P1, _3D>
{
  public:
  /// Default constructor
  ShapeFunctionT();

  /// Destructor
  ~ShapeFunctionT();

  /// Return the name of this shape function
  static std::string type_name()
  {
    return "Triag-P1-Lagrange-2D";
  }

  /// Return the number of nodes (degrees of freedom) for this sf
  const Uint nb_dof() const;

  /// Fill a matrix containing the reference coordinates
  void ref_coords(TCoord &coords) const;

  /// Compute the values of all shape functions at one point:
  void compute_ref_values(const math::VectorBlock<Real> &ref_coords,
                          math::DynamicVector<Real> &values);

  /// Compute the derivatives of all shape functions at one point:
  void compute_ref_derivatives(const math::VectorBlock<Real> &ref_coords,
                               math::DynamicMatrix<Real> &values);

  void print() const;

  private:
  /// TYPEDFES
  typedef LagrangeShapeFunction<P1, Tetra, _3D> BaseShapeFunction;

  /// STATIC VARIABLES
  static const Uint NbNodes;
  static const Uint Order;
};

typedef ShapeFunctionT<Tetra, P1, Lagrange, _3D> TetraP1Lagrange3D;

} // namespace sf

} // namespace mesh

} // namespace pdekit

#endif
