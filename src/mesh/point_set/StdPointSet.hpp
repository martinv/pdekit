#ifndef PDEKIT_Mesh_Std_Point_Set_hpp
#define PDEKIT_Mesh_Std_Point_Set_hpp

#include "common/Flyweight.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "mesh/EntityRealignCode.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

namespace detail
{

class StdPointSetInstance
{
  public:
  typedef mesh::PointSetTag tag_type;

  static const tag_type undefined;

  /// Construct one instance of quadrature based on tag
  static void construct(const mesh::PointSetTag &tag, StdPointSetInstance &qinstance);

  /// Default constructor
  StdPointSetInstance();

  /// Default destructor
  ~StdPointSetInstance();

  /// Get the tag of this quadrature
  const mesh::PointSetTag &tag() const;

  /// Get the dimension of the reference element for which
  /// this quadrature rule can be used
  Uint dim() const;

  /// Get the codimension of the quadrature rule. For example,
  /// if the quadrature rule defines quadrature on the (2D) faces of
  /// reference hexahedron (which is a 3D element), then
  /// codim = 3 - 2 = 1
  Uint codim() const;

  /// Return the number of local entities (faces, edges)
  /// on which this quadrature has points
  Uint nb_local_entities() const;

  /// Get the number of quadrature points
  Uint size(const Uint local_idx = 0) const;

  /// Get the quadrature coordinates
  const math::DenseDMat<Real> &coordinates(const Uint local_idx = 0) const;

  /// Get the quadrature weights
  const math::DenseDVec<Real> &weights(const Uint local_idx = 0) const;

  private:
  /// Tag of this quadrature
  mesh::PointSetTag m_tag;

  /// Dimension of this quadrature rule
  Uint m_dim;

  /// Codimension of this quadrature rule
  Uint m_codim;

  /// Quadrature coordinates
  std::vector<math::DenseDMat<Real>> m_coords;

  /// Quadrature weights
  std::vector<math::DenseDVec<Real>> m_weights;
};

struct StdPointSetFlyweightPolicy
{
  typedef StdPointSetInstance::tag_type key_type;
};

} // namespace detail

typedef common::Flyweight<detail::StdPointSetInstance, detail::StdPointSetFlyweightPolicy>
    StdPointSet;

} // namespace mesh

} // namespace pdekit

#endif
