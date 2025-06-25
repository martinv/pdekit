#ifndef PDEKIT_Mesh_Interpolation_Equidist_Std_Region_Pyramid_hpp
#define PDEKIT_Mesh_Interpolation_Equidist_Std_Region_Pyramid_hpp

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
// EQUIDISTANT STANDARD REGION ON PYRAMIDS
// ----------------------------------------------------------------------------

// LOCAL TOPOLOGY:
// The quadrilateral (0,1,2,3) forms the base of the pyramid in the xi-eta (x-y)
// plane. The vertex 4 is directly above vertex 0, not above the center of the
// base (0,1,2,3)
//
//
//                    ZETA                          ORDER OF EDGES:
//                      ^                           1) 0-1
//                      |                           2) 0-3
//         4            |                           3) 0-4
//         |\ \         |                           4) 1-2
//         |:\  \       |                           5) 1-4
//         | \ \   \    |                           6) 2-3
//         | :   \    \ |                           7) 2-4
//         |  :   .\    |\                          8) 3-4
//         |  \      \ .|   \.
//         0--`--------\|-----3                     ORDER OF FACES
//         `\  \        |\    `\                    1) 0-1-4
//          `\ `        |  \   `\                   2) 3-0-4
//           `\ \       |----\--`\-------> ETA      3) 1-2-4
//            `\ `      \      \ `\                 4) 2-3-4
//             `\|       \        \\                5) 0-3-2-1
//               1------------------2
//                         \.
//                          \.
//                          XI
//

class EquidistStdRegionPyramid
{
  public:
  static constexpr PointSetID Type     = PointSetID::Equidist;
  static constexpr ElemShape GeomShape = ElemShape::Pyramid;
  static constexpr Uint TopoDim        = ElementTopology<ElemShape::Pyramid>::TopoDim;
  static constexpr Uint NbVert         = ElementTopology<ElemShape::Pyramid>::NbVerts;
  static const std::array<Uint, 4> TopologyStorage;
  static constexpr Real ref_measure = 8.0 / 3.0;

  EquidistStdRegionPyramid() = default;

  ~EquidistStdRegionPyramid() = default;

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
  // -------------------------------------------
  struct Edge01
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 0;
      if (i == 1)
        return 1;
      return 2 + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge03
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 0;
      if (i == 1)
        return 3;
      return 2 + (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge04
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 0;
      if (i == 1)
        return 4;
      return 2 + 2 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge12
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 1;
      if (i == 1)
        return 2;
      return 2 + 3 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge14
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 1;
      if (i == 1)
        return 4;
      return 2 + 4 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge23
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 2;
      if (i == 1)
        return 3;
      return 2 + 5 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge24
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 2;
      if (i == 1)
        return 4;
      return 2 + 6 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge34
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 3;
      if (i == 1)
        return 4;
      return 2 + 7 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
