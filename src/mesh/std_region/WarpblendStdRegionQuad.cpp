#include "mesh/std_region/WarpblendStdRegionQuad.hpp"

#include "math/polynomials/PolyLib.hpp"
#include "mesh/std_region/EquidistStdRegionQuad.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

/// A quadrilateral had 4 edges, 1 face and  0 volumes (but we store one default
/// 'invalid' entity as its volume)
const std::array<Uint, 4> WarpblendStdRegionQuad::TopologyStorage = {1, 4, 1, 1};

// ----------------------------------------------------------------------------

Uint WarpblendStdRegionQuad::nb_dof(const Uint poly_order)
{
  return (poly_order + 1) * (poly_order + 1);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionQuad::fill_topology_relations(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences)
{
  EquidistStdRegionQuad::fill_topology_relations(incidences);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionQuad::fill_edge_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  EquidistStdRegionQuad::fill_edge_to_node_connectivity(poly_order, entity_vec);

  // Override the type of the entities
  for (Uint e = 0; e < 4; ++e)
  {
    entity_vec[e]->change_point_set_type(PointSetID::Warpblend);
  }
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionQuad::fill_face_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  EquidistStdRegionQuad::fill_face_to_node_connectivity(poly_order, entity_vec);

  // Override the topology type of the face
  StdRegionEntity &face = *entity_vec[0];
  face.change_point_set_type(PointSetID::Warpblend);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionQuad::fill_volume_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  StdRegionEntity &current_entity = *entity_vec[0];
  // This entity will not be resized - it has length 0 by default

  const common::ArrayView<const Uint, _1D, Uint> vert_ids_view;
  const common::ArrayView<const Uint, _1D, Uint> p1_vert_flags_view;

  current_entity.construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined),
                           vert_ids_view, p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionQuad::fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord)
{
  if (poly_order < P3)
  {
    EquidistStdRegionQuad::fill_coordinates(poly_order, coord);
    return;
  }

  const Uint nb_pts_on_interval = poly_order + 1;
  math::DenseDVec<Real> tmp(nb_pts_on_interval - 2);
  math::DenseDVec<Real> pt_coords_1d(nb_pts_on_interval);

  math::jac_zeros(nb_pts_on_interval - 2, tmp, 1, 1);

  pt_coords_1d[0]                      = -1.0;
  pt_coords_1d[nb_pts_on_interval - 1] = 1.0;
  for (Uint i = 0; i < nb_pts_on_interval - 2; ++i)
  {
    pt_coords_1d[i + 1] = tmp[i];
  }

  for (Uint i = 0; i < pt_coords_1d.size(); ++i)
  {
    if (std::abs(pt_coords_1d[i]) < 1.e-16)
    {
      pt_coords_1d[i] = 0.0;
    }
  }

  // std::cout << "Points on biunit inteval: " << std::endl << pt_coords_1d <<
  // std::endl;

  /// The coordinates have to be written in such a manner that first come
  /// vertex nodes, then edge nodes. For the remaining interior nodes, we
  /// repeat the procedure recursively. Imagine we remove the edges of the
  /// quadrilateral. We will be left with a smaller interior quad. We apply
  /// the same numbering to it again: first vertices, then internal edge
  /// nodes. After, we peel its edges off and so on until no nodes are left
  const Uint nb_nodes = (poly_order + 1) * (poly_order + 1);

  coord.resize(nb_nodes, _2D);

  // First vertex  (index 0)
  coord(0, KSI) = -1.0;
  coord(0, ETA) = -1.0;
  // Second vertex (index 1)
  coord(1, KSI) = 1.0;
  coord(1, ETA) = -1.0;
  // Third vertex  (index 2)
  coord(2, KSI) = 1.0;
  coord(2, ETA) = 1.0;
  // Fourth vertex (index 3)
  coord(3, KSI) = -1.0;
  coord(3, ETA) = 1.0;

  if (poly_order == P2)
  {
    coord(4, KSI) = 0.0;
    coord(4, ETA) = -1.0;
    coord(5, KSI) = 1.0;
    coord(5, ETA) = 0.0;
    coord(6, KSI) = 0.0;
    coord(6, ETA) = 1.0;
    coord(7, KSI) = -1.0;
    coord(7, ETA) = 0.0;
  }

  if (poly_order > P2)
  {
    Uint inode = 4;

    // First finish the outermost edges - points on the contour of the
    // reference element

    for (Uint i = 1; i < poly_order; ++i)
    {
      coord(inode, KSI) = pt_coords_1d[i]; // x-coordinate
      coord(inode, ETA) = -1.0;            // y-coordinate
      inode++;
    }

    // Node coordinates of the second edge
    for (Uint i = 1; i < poly_order; ++i)
    {
      coord(inode, KSI) = 1.;
      coord(inode, ETA) = pt_coords_1d[i];
      inode++;
    }

    // Node coordinates of the third edge
    for (Uint i = 1; i < poly_order; ++i)
    {
      coord(inode, KSI) = pt_coords_1d[poly_order - i];
      coord(inode, ETA) = 1.;
      inode++;
    }

    // Node coordinates of the fourth edge
    for (Uint i = 1; i < poly_order; ++i)
    {
      coord(inode, KSI) = -1.;
      coord(inode, ETA) = pt_coords_1d[poly_order - i];
      inode++;
    }

    // Offset denotes how many layers of nodes from the boundary of the main
    // quad are we (how many contours have we 'peeled off' so far)
    Uint offset = 1;
    // Polynomial order to which the remaining interior nodes correspond
    Int current_poly_order = poly_order - 2;

    while (current_poly_order >= 0)
    {
      // First corner node, this can als0 be center node (e.g. in P3 quad)
      // Therefore we might have to stop the filling after this node.
      coord(inode, KSI) = pt_coords_1d[offset];
      coord(inode, ETA) = pt_coords_1d[offset];
      inode++;
      if (current_poly_order == _0D)
        break;

      // Second corner
      coord(inode, KSI) = pt_coords_1d[nb_pts_on_interval - offset - 1];
      coord(inode, ETA) = pt_coords_1d[offset];
      inode++;

      // Third corner
      coord(inode, KSI) = pt_coords_1d[nb_pts_on_interval - offset - 1];
      coord(inode, ETA) = pt_coords_1d[nb_pts_on_interval - offset - 1];
      inode++;

      // Fourth corner
      coord(inode, KSI) = pt_coords_1d[offset];
      coord(inode, ETA) = pt_coords_1d[nb_pts_on_interval - offset - 1];
      inode++;

      if (current_poly_order == _1D)
        break;

      // Print high-order edge nodes:

      // Node coordinates of the first sub-quad edge
      for (Int i = 1; i < current_poly_order; ++i)
      {
        coord(inode, KSI) = pt_coords_1d[offset + i]; // x-coordinate
        coord(inode, ETA) = pt_coords_1d[offset];     // y-coordinate
        inode++;
      }

      // Node coordinates of the second edge
      for (Int i = 1; i < current_poly_order; ++i)
      {
        coord(inode, KSI) = pt_coords_1d[nb_pts_on_interval - offset - 1];
        coord(inode, ETA) = pt_coords_1d[offset + i];
        inode++;
      }

      // Node coordinates of the third edge
      for (Int i = 1; i < current_poly_order; ++i)
      {
        coord(inode, KSI) = pt_coords_1d[nb_pts_on_interval - offset - i - 1];
        coord(inode, ETA) = pt_coords_1d[nb_pts_on_interval - offset - 1];
        inode++;
      }

      // Node coordinates of the fourth edge
      for (Int i = 1; i < current_poly_order; ++i)
      {
        coord(inode, KSI) = pt_coords_1d[offset];
        coord(inode, ETA) = pt_coords_1d[nb_pts_on_interval - offset - i - 1];
        inode++;
      }

      current_poly_order -= 2;
      offset++;
    }

  } // If polynomial order > 2
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionQuad::fill_facet_normals(math::DenseDMat<Real> &normals)
{
  EquidistStdRegionQuad::fill_facet_normals(normals);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionQuad::fill_permutation(const Uint poly_order,
                                              const EntityRealignCode &permutation_code,
                                              std::vector<Uint> &permutation_vec)
{
  EquidistStdRegionQuad::fill_permutation(poly_order, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionQuad::fill_rotation_permutation(const Uint poly_order,
                                                       std::vector<Uint> &permutation_vec)
{
  EquidistStdRegionQuad::fill_rotation_permutation(poly_order, permutation_vec);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionQuad::fill_flip_permutation(const Uint poly_order,
                                                   std::vector<Uint> &permutation_vec)
{
  EquidistStdRegionQuad::fill_flip_permutation(poly_order, permutation_vec);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
