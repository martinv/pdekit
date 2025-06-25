#ifndef PDEKIT_Mesh_Interpolation_Equidist_Std_Region_Hexa_hpp
#define PDEKIT_Mesh_Interpolation_Equidist_Std_Region_Hexa_hpp

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
// EQUIDISTANT STANDARD REGION ON HEXAHEDRA
// Check MHexahedron.h in Gmsh source code for local numbering
// ----------------------------------------------------------------------------

// LOCAL TOPOLOGY:
//
// ORDER OF EDGES:     ORDER OF FACES:
// 1)  0-1             1) 0-3-2-1
// 2)  0-3             2) 0-1-5-4
// 3)  0-4             3) 0-4-7-3
// 4)  1-2             4) 1-2-6-5
// 5)  1-5             5) 2-3-7-6
// 6)  2-3             6) 4-5-6-7
// 7)  2-6
// 8)  3-7
// 9)  4-5
// 10) 4-7
// 11) 5-6
// 12) 6-7
//
/*           w
//    4----------7            4----17----7
//    |\     ^   |\           |\         |\
//    | \    |   | \          | 16       | 19
//    |  \   |   |  \        10  \       15 \
//    |   5------+---6        |   5----18+---6
//    |   |  +-- |-- | -> v   |   |      |   |
//    0---+---\--3   |        0---+-9----3   |
//     \  |    \  \  |         \  12      \  14
//      \ |     \  \ |          8 |        13|
//       \|      u  \|           \|         \|
//        1----------2            1----11----2
*/

class EquidistStdRegionHexa
{
  public:
  static constexpr PointSetID Type     = PointSetID::Equidist;
  static constexpr ElemShape GeomShape = ElemShape::Hexa;
  static constexpr Uint TopoDim        = ElementTopology<ElemShape::Hexa>::TopoDim;
  static constexpr Uint NbVert         = ElementTopology<ElemShape::Hexa>::NbVerts;
  static const std::array<Uint, 4> TopologyStorage;
  static constexpr Real ref_measure = 8.0;

  EquidistStdRegionHexa() = default;

  ~EquidistStdRegionHexa() = default;

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
      // If the position index i is not 0 or 1, it is at
      // least 2, which indicates a higher-order
      // node on the edge. In that case, the vertex
      // number must be at least 8 (vertices 0-7 are
      // corner vertices of the reference hexahedron).
      // Therefore the number of such vertex will be 6+i
      return 6 + i;
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
      return 6 + (order - 1) + i;
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
      return 6 + 2 * (order - 1) + i;
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
      return 6 + 3 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge15
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 1;
      if (i == 1)
        return 5;
      return 6 + 4 * (order - 1) + i;
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
      return 6 + 5 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge26
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 2;
      if (i == 1)
        return 6;
      return 6 + 6 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge37
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 3;
      if (i == 1)
        return 7;
      return 6 + 7 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge45
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 4;
      if (i == 1)
        return 5;
      return 6 + 8 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge47
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 4;
      if (i == 1)
        return 7;
      return 6 + 9 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge56
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 5;
      if (i == 1)
        return 6;
      return 6 + 10 * (order - 1) + i;
    }

    static Uint nb_nodes(const Uint order)
    {
      return order + 1;
    }
  };
  // -------------------------------------------
  struct Edge67
  {
    static Uint node_id(const Uint order, const Uint i)
    {
      if (i == 0)
        return 6;
      if (i == 1)
        return 7;
      return 6 + 11 * (order - 1) + i;
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
