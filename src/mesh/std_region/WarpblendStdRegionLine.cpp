#include "mesh/std_region/WarpblendStdRegionLine.hpp"

#include "math/DenseDVec.hpp"
#include "math/polynomials/PolyLib.hpp"
#include "mesh/std_region/EquidistStdRegionLine.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

const std::array<Uint, 4> WarpblendStdRegionLine::TopologyStorage = {1, 1, 1, 1};

// ----------------------------------------------------------------------------

Uint WarpblendStdRegionLine::nb_dof(const Uint poly_order)
{
  return poly_order + 1;
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionLine::fill_topology_relations(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences)
{
  EquidistStdRegionLine::fill_topology_relations(incidences);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionLine::fill_edge_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  EquidistStdRegionLine::fill_edge_to_node_connectivity(poly_order, entity_vec);

  StdRegionEntity &current_entity = *entity_vec[0];
  current_entity.change_point_set_type(PointSetID::Warpblend);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionLine::fill_face_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  StdRegionEntity &current_entity = *entity_vec[0];

  const common::ArrayView<const Uint, _1D, Uint> vert_ids_view;
  const common::ArrayView<const Uint, _1D, Uint> p1_vert_flags_view;

  current_entity.construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined),
                           vert_ids_view, p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionLine::fill_volume_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  StdRegionEntity &current_entity = *entity_vec[0];

  const common::ArrayView<const Uint, _1D, Uint> vert_ids_view;
  const common::ArrayView<const Uint, _1D, Uint> p1_vert_flags_view;

  current_entity.construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined),
                           vert_ids_view, p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionLine::fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord)
{
  const Uint nb_nodes = poly_order + 1;
  coord.resize(nb_nodes, TopoDim);

  // First vertex  (index 0)
  coord(0, KSI) = -1.0;
  // Second vertex (index 1)
  coord(1, KSI) = 1.0;

  if (poly_order == P2)
  {
    coord(2, KSI) = 0.0;
  }

  if (poly_order > P2)
  {
    const Real alpha = 0.0;
    const Real beta  = 0.0;
    math::DenseDVec<Real> z(nb_nodes - 2);
    math::jac_zeros(nb_nodes - 2, z, alpha + 1, beta + 1);

    // Node coordinates of the reference interval [-1,1] - internal nodes
    Uint inode = 2;
    for (Uint i = 0; i < (nb_nodes - 2); ++i)
    {
      coord(inode, KSI) = z[i]; // x-coordinate
      inode++;
    }

  } // If polynomial order > 1
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionLine::fill_facet_normals(math::DenseDMat<Real> &normals)
{
  EquidistStdRegionLine::fill_facet_normals(normals);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionLine::fill_permutation(const Uint poly_order,
                                              const EntityRealignCode &permutation_code,
                                              std::vector<Uint> &permutation_vec)
{
  EquidistStdRegionLine::fill_permutation(poly_order, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionLine::fill_rotation_permutation(const Uint poly_order,
                                                       std::vector<Uint> &permutation_vec)
{
  EquidistStdRegionLine::fill_rotation_permutation(poly_order, permutation_vec);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionLine::fill_flip_permutation(const Uint poly_order,
                                                   std::vector<Uint> &permutation_vec)
{
  EquidistStdRegionLine::fill_flip_permutation(poly_order, permutation_vec);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
