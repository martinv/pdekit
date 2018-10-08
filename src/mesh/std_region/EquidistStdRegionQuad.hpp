#ifndef PDEKIT_Mesh_Interpolation_Equidist_Std_Region_Quad_hpp
#define PDEKIT_Mesh_Interpolation_Equidist_Std_Region_Quad_hpp

#include <array>

#include "common/BlockArray.hpp"
#include "math/DenseDMat.hpp"
#include "mesh/ElementTopology.hpp"
#include "mesh/EntityRealignCode.hpp"
#include "mesh/MeshConstants.hpp"
#include "mesh/std_region/StdRegionEntity.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// EQUIDISTANT STANDARD REGION ON QUADRILATERALS
// ----------------------------------------------------------------------------

class EquidistStdRegionQuad
{
  public:
  static constexpr PointSetID Type     = PointSetID::Equidist;
  static constexpr ElemShape GeomShape = ElemShape::Quad;
  static constexpr Uint TopoDim        = ElementTopology<ElemShape::Quad>::TopoDim;
  static constexpr Uint NbVert         = ElementTopology<ElemShape::Quad>::NbVerts;
  static const std::array<Uint, 4> TopologyStorage;
  static constexpr Real ref_measure = 4.0;

  EquidistStdRegionQuad() = default;

  ~EquidistStdRegionQuad() = default;

  static Uint nb_dof(const Uint poly_order);

  static void fill_topology_relations(
      std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences);

  static void fill_edge_to_node_connectivity(
      const Uint poly_order,
      common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec);

  static void fill_face_to_node_connectivity(
      const Uint poly_order,
      common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec);

  static void fill_volume_to_node_connectivity(
      const Uint poly_order,
      common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec);

  static void fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord);

  static void fill_facet_normals(math::DenseDMat<Real> &normals);

  static void fill_permutation(const Uint poly_order, const EntityRealignCode &permutation_code,
                               std::vector<Uint> &permutation_vec);

  private:
  friend class WarpblendStdRegionQuad;

  static void fill_rotation_permutation(const Uint poly_order, std::vector<Uint> &permutation_vec);

  static void fill_flip_permutation(const Uint poly_order, std::vector<Uint> &permutation_vec);
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
