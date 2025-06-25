#ifndef PDEKIT_Mesh_Interpolation_Warpblend_Std_Region_Tetra_hpp
#define PDEKIT_Mesh_Interpolation_Warpblend_Std_Region_Tetra_hpp

#include <array>

#include "common/BlockArray.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseSVec.hpp"
#include "mesh/ElementTopology.hpp"
#include "mesh/EntityRealignCode.hpp"
#include "mesh/MeshConstants.hpp"
#include "mesh/std_region/StdRegionEntity.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// WARPBLEND STANDARD REGION ON TETRAHEDRA
// ----------------------------------------------------------------------------

// LOCAL TOPOLOGY:
// The triangle (0,1,2) forms the base of the tetra in the xi-eta (x-y) plane.
// The vertex 3
// is on the vertical (z,zeta) axis
//                                           ORDER OF EDGES:
//              ETA                          1) 0-1
//             /                             2) 1-2
//            /                              3) 2-0
//           2                               4) 3-0
//          /||                              5) 3-2
//         / | |                             6) 3-1
//        /  |  |
//       /   |   |                           ORDER OF FACES:
//      /    |    |                          1) 0-2-1
//     0.....|.....1-----> XI                2) 0-1-3
//      |    |    /                          3) 0-3-2
//       |   |   /                           4) 3-1-2
//        |  |  /
//         | | /
//          ||/
//           3
//            |
//             |
//            ZETA

class WarpblendStdRegionTetra
{
  public:
  static constexpr PointSetID Type     = PointSetID::Warpblend;
  static constexpr ElemShape GeomShape = ElemShape::Tetra;
  static constexpr Uint TopoDim        = ElementTopology<ElemShape::Tetra>::TopoDim;
  static constexpr Uint NbVert         = ElementTopology<ElemShape::Tetra>::NbVerts;
  static const std::array<Uint, 4> TopologyStorage;
  static constexpr Real ref_measure = 4.0 / 3.0;

  WarpblendStdRegionTetra() = default;

  ~WarpblendStdRegionTetra() = default;

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
  /// METHODS
  static void warp_shift_face_3D(const Uint poly_order, const Real alpha,
                                 const math::DenseDVec<Real> &L1, const math::DenseDVec<Real> &L2,
                                 const math::DenseDVec<Real> &L3, math::DenseDVec<Real> &warp1,
                                 math::DenseDVec<Real> &warp2);

  static void evalwarp(const Uint poly_order, const math::DenseDVec<Real> &xnodes,
                       const math::DenseDVec<Real> &xout, math::DenseDVec<Real> &warp);

  /// DATA
  static const math::DenseSVec<Real, 15> alpopt;
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
