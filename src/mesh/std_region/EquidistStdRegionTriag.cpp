#include "mesh/std_region/EquidistStdRegionTriag.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

/// A triangle 3 edges, 1 face and  0 volumes (but we store one default
/// 'invalid' entity as its volume)
const std::array<Uint, 4> EquidistStdRegionTriag::TopologyStorage = {1, 3, 1, 1};

// ----------------------------------------------------------------------------

Uint EquidistStdRegionTriag::nb_dof(const Uint poly_order)
{
  return (poly_order + 1) * (poly_order + 2) / 2;
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTriag::fill_topology_relations(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences)
{
  // Node to edge incidences - position (0,1)

  // edge to node incidences - position (1,0)
  common::BlockArray<SUint, SUint> &edge_to_node = incidences[(_3D + 1) * _1D + _0D];

  // Edge 0 is incident to nodes 0 and 1
  // Edge 1 is incident to nodes 1 and 2
  // Edge 2 is incident to nodes 2 and 0
  edge_to_node.build({0, 1, 1, 2, 2, 0}, {2, 2, 2});

  // edge to edge incidences - position (1,1)
  common::BlockArray<SUint, SUint> &edge_to_edge = incidences[(_3D + 1) * _1D + _1D];
  edge_to_edge.build({2, 1, 0, 2, 1, 0}, {2, 2, 2});

  // face-to-node incidence should be stored at position (2,0)
  common::BlockArray<SUint, SUint> &face_to_node = incidences[(_3D + 1) * _2D + _0D];
  // One entity - the triangle face
  face_to_node.build({0}, {1});

  // face-to-edge incidence should be stored at position (2,1)
  common::BlockArray<SUint, SUint> &face_to_edge = incidences[(_3D + 1) * _2D + _1D];
  // Face 0 (the only face present) is incident to edge 0
  // Face 0 (the only face present) is incident to edge 1
  // Face 0 (the only face present) is incident to edge 2
  face_to_edge.build({0, 1, 2}, {3});

  // face-to-face incidence should be stored at position (2,2)
  common::BlockArray<SUint, SUint> &face_to_face = incidences[(_3D + 1) * _2D + _2D];
  // Face 0 (the only face present) is incident only to itself
  // This incidence has to be filled, otherwise 'local_transform'
  // in MeshEntity will fail when the entity is transformed to
  // itself!
  face_to_face.build({0}, {1});

  // volume-to-node incidence should be stored at position (3,0)
  common::BlockArray<SUint, SUint> &volume_to_node = incidences[(_3D + 1) * _3D + _0D];
  // One entity - default (invalid) volume
  volume_to_node.build({0}, {1});
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTriag::fill_edge_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{

  /*
  // Example: numbering of a P4 equidistant point set

  //  2
  //  |\
  //  | \
  //  9  8
  //  |   \
  //  |    \
  // 10 14  7
  //  |      \
  //  |       \
  // 11 12 13  6
  //  |         \
  //  |          \
  //  0--3--4--5--1
  */

  // The triangle has 3 edges ...
  // entity_vec.resize(3);

  std::vector<Uint> vert_ids;
  std::vector<Uint> p1_vert_flags;

  // Set correctly the order and dimension of the entities

  // First edge: [ 0 - 3 - 4 - 5 - ....  - 1 ] - 0 and 1 are the end vertices,
  // the number 2
  // is assigned to the third vertex of the triangle, so the numbers 3, 4, ...
  // are the
  // internal nodes of the first edge

  StdRegionEntity &edge0 = *entity_vec[0];

  vert_ids.resize(poly_order + 1);
  p1_vert_flags.resize(poly_order + 1);
  p1_vert_flags.assign(poly_order + 1, 0);

  vert_ids[0] = 0;
  vert_ids[1] = 1;

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  // Where to place it in the reference entity
  Uint index = 2;
  // Number of the node (its 'value')
  Uint node_nr = 3;

  for (Uint i = 0; i < poly_order - 1; ++i)
  {
    vert_ids[index++] = node_nr++;
  }

  auto vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  auto p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge0.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 0);

  // Second edge: 1 and 2 are the end vertices

  StdRegionEntity &edge1 = *entity_vec[1];

  vert_ids.resize(poly_order + 1);
  p1_vert_flags.resize(poly_order + 1);
  p1_vert_flags.assign(poly_order + 1, 0);

  vert_ids[0] = 1;
  vert_ids[1] = 2;

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  index   = 2;
  node_nr = 2 + (poly_order - 1) + 1;

  for (Uint i = 0; i < poly_order - 1; ++i)
  {
    vert_ids[index++] = node_nr++;
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge1.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 1);

  // Third edge: 2 and 0 are the end nodes

  StdRegionEntity &edge2 = *entity_vec[2];
  // edge2.resize(Order+1);
  vert_ids.resize(poly_order + 1);
  p1_vert_flags.resize(poly_order + 1);
  p1_vert_flags.assign(poly_order + 1, 0);

  vert_ids[0] = 2;
  vert_ids[1] = 0;

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  index   = 2;
  node_nr = 2 + 2 * (poly_order - 1) + 1;

  for (Uint i = 0; i < poly_order - 1; ++i)
  {
    vert_ids[index++] = node_nr++;
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge2.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 2);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTriag::fill_face_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // entity_vec.resize(1);
  StdRegionEntity &face = *entity_vec[0];

  const Uint nb_nodes = (poly_order + 1) * (poly_order + 2) / 2;

  std::vector<Uint> vert_ids(nb_nodes);
  std::iota(vert_ids.begin(), vert_ids.end(), 0);

  std::vector<Uint> p1_vert_flags(nb_nodes);
  p1_vert_flags.assign(nb_nodes, 0);
  for (Uint n = 0; n < 3; ++n)
  {
    p1_vert_flags[n] = 1;
  }

  const auto vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  const auto p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face.construct(PointSetTag(ElemShape::Triag, poly_order, PointSetID::Equidist), vert_ids_view,
                 p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTriag::fill_volume_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // entity_vec.resize(1);
  StdRegionEntity &volume = *entity_vec[0];

  // This entity will not be resized - it has length 0 by default

  const common::ArrayView<const Uint, _1D, Uint> vert_ids_view;
  const common::ArrayView<const Uint, _1D, Uint> p1_vert_flags_view;

  volume.construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined), vert_ids_view,
                   p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTriag::fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord)
{
  /// The coordinates have to be written in such a manner that first come
  /// vertex nodes, then edge nodes. For the remaining interior nodes, we
  /// repeat the procedure recursively. Imagine we remove the edges of the
  /// triangle. We will be left with a smaller interior triangle. We apply the
  /// same numbering to it again: first vertices, then internal edge nodes.
  /// After, we peel its edges off and so on until no nodes are left

  const Uint nb_nodes = (poly_order + 1) * (poly_order + 2) / 2;
  coord.resize(nb_nodes, _2D);

  // First vertex  (index 0)
  coord(0, KSI) = -1.0;
  coord(0, ETA) = -1.0;
  // Second vertex (index 1)
  coord(1, KSI) = 1.0;
  coord(1, ETA) = -1.0;
  // Third vertex  (index 2)
  coord(2, KSI) = -1.0;
  coord(2, ETA) = 1.0;

  if (poly_order > P1)
  {

    const Real dx = 2. / poly_order;
    const Real dy = 2. / poly_order;

    // Node coordinates of the first edge
    Uint inode = 3;
    for (Uint i = 1; i < poly_order; ++i)
    {
      coord(inode, KSI) = -1. + i * dx; // x-coordinate
      coord(inode, ETA) = -1.0;         // y-coordinate
      inode++;
    }

    // Node coordinates of the second edge
    for (Uint i = 1; i < poly_order; ++i)
    {
      coord(inode, KSI) = 1. - i * dx;
      coord(inode, ETA) = -1. + i * dy;
      inode++;
    }

    // Node coordinates of the third edge
    for (Uint i = 1; i < poly_order; ++i)
    {
      coord(inode, KSI) = -1.;
      coord(inode, ETA) = 1. - i * dy;
      inode++;
    }

    if (poly_order > P2)
    {
      // Offset denotes how many layers of nodes from the boundary of the
      // main triangle are we (how many contours have we 'peeled off' so
      // far)
      Uint offset = 1;
      // Polynomial order to which the remaining interior nodes correspond
      Int current_poly_order = poly_order - 3;

      while (current_poly_order >= 0)
      {
        // First corner node, this can alos be center node (e.g. in P3
        // triangle) Therefore we might have to stop the filling after
        // this node.
        coord(inode, KSI) = -1. + offset * dx;
        coord(inode, ETA) = -1. + offset * dy;
        inode++;
        if (current_poly_order == _0D)
          break;

        // Second corner
        coord(inode, KSI) = 1. - (2 * offset) * dx;
        coord(inode, ETA) = -1. + offset * dy;
        inode++;

        // Third corner
        coord(inode, KSI) = -1. + offset * dx;
        coord(inode, ETA) = 1. - (2 * offset) * dy;
        inode++;

        if (current_poly_order == _1D)
          break;

        // Print high-order edge nodes:

        // Node coordinates of the first sub-triag edge
        for (Int i = 1; i < current_poly_order; ++i)
        {
          coord(inode, KSI) = -1. + (offset + i) * dx; // x-coordinate
          coord(inode, ETA) = -1. + offset * dy;       // y-coordinate
          inode++;
        }

        // Node coordinates of the second edge
        for (Int i = 1; i < current_poly_order; ++i)
        {
          coord(inode, KSI) = 1. - (2 * offset + i) * dx;
          coord(inode, ETA) = -1. + (offset + i) * dy;
          inode++;
        }

        // Node coordinates of the third edge
        for (Int i = 1; i < current_poly_order; ++i)
        {
          coord(inode, KSI) = -1. + offset * dx;
          coord(inode, ETA) = 1. - (2 * offset + i) * dy;
          inode++;
        }

        current_poly_order -= 3;
        offset++;
      }

    } // If polynomial order > 2

  } // If polynomial order > 1
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTriag::fill_facet_normals(math::DenseDMat<Real> &normals)
{
  normals.resize(3, 2);

  // Face 1
  normals(0, XI0) = 0.0;
  normals(0, XI1) = -1.0;

  // Face 2
  normals(1, XI0) = 1.0 / std::sqrt(2.0);
  normals(1, XI1) = 1.0 / std::sqrt(2.0);

  // Face 3
  normals(2, XI0) = -1.0;
  normals(2, XI1) = 0.0;
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTriag::fill_permutation(const Uint poly_order,
                                              const EntityRealignCode &permutation_code,
                                              std::vector<Uint> &permutation_vec)
{
  const Uint nb_nodes = (poly_order + 1) * (poly_order + 2) / 2;
  // First set permutation_vec as identity permutation: p[i] = i;
  permutation_vec.resize(nb_nodes);
  for (Uint i = 0; i < permutation_vec.size(); ++i)
  {
    permutation_vec[i] = i;
  }

  // If this permutation is identity, we are done
  if (permutation_code.is_identity(ElemShape::Triag))
  {
    return;
  }

  // Temporary vector
  std::vector<Uint> old_permutation(nb_nodes);
  std::vector<Uint> new_permutation(nb_nodes);

  // Apply all flips first
  for (Uint i = 0; i < permutation_code.nb_flips(); ++i)
  {
    old_permutation.swap(permutation_vec);

    fill_flip_permutation(poly_order, new_permutation);

    for (Uint j = 0; j < new_permutation.size(); ++j)
    {
      permutation_vec[j] = old_permutation[new_permutation[j]];
    }
  }

  // Then apply all rotations
  for (Uint i = 0; i < permutation_code.nb_rotations(); ++i)
  {
    old_permutation.swap(permutation_vec);

    fill_rotation_permutation(poly_order, new_permutation);

    for (Uint j = 0; j < new_permutation.size(); ++j)
    {
      permutation_vec[j] = old_permutation[new_permutation[j]];
    }
  }
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTriag::fill_rotation_permutation(const Uint poly_order,
                                                       std::vector<Uint> &permutation_vec)
{
  const Uint nb_nodes = (poly_order + 1) * (poly_order + 2) / 2;
  permutation_vec.resize(nb_nodes);

  // Offset denotes what is the position of the first node in the
  // permutation_vec
  // for the triangle that we are currently filling.
  Uint offset = 0;
  // Polynomial order to which the remaining interior nodes correspond
  Int current_poly_order = poly_order;

  while (current_poly_order >= 0)
  {
    // If the current triangle is of order 0, it has one center node
    // In that case, insert the center node and exit the while loop
    if (current_poly_order == _0D)
    {
      permutation_vec[offset] = offset;
    }

    // This should be executed when the current polynomial order is >= 1
    if (current_poly_order >= _1D)
    {
      permutation_vec[offset]     = offset + 2;
      permutation_vec[offset + 1] = offset + 0u;
      permutation_vec[offset + 2] = offset + 1;

      // If we have some edge-internal nodes, let's insert them in
      // permutation_vec
      if (current_poly_order >= _2D)
      {
        const Uint nb_internal_pts_on_edge = current_poly_order - 1;

        // Edge 2-0 will become first edge
        for (Uint i = 0; i < nb_internal_pts_on_edge; ++i)
        {
          permutation_vec[offset + 3 + i] = offset + 3 + 2 * (nb_internal_pts_on_edge) + i;
        }

        // Edge 0-1 will become second edge
        for (Uint i = 0; i < nb_internal_pts_on_edge; ++i)
        {
          permutation_vec[offset + 3 + nb_internal_pts_on_edge + i] = offset + 3 + i;
        }

        // Edge 1-2 will become third edge
        for (Uint i = 0; i < nb_internal_pts_on_edge; ++i)
        {
          permutation_vec[offset + 3 + 2 * nb_internal_pts_on_edge + i] =
              offset + 3 + nb_internal_pts_on_edge + i;
        }
      } // Polynomial order >= 2

    } // Polynomial order >= 1

    // For order P, we have P+1 nodes on edge and 3*P nodes on the contour
    // of the triangle
    offset += 3 * current_poly_order;
    current_poly_order -= 3;

  } // While polynomial order >= 0

} // End of method

// ----------------------------------------------------------------------------

void EquidistStdRegionTriag::fill_flip_permutation(const Uint poly_order,
                                                   std::vector<Uint> &permutation_vec)
{
  const Uint nb_nodes = (poly_order + 1) * (poly_order + 2) / 2;
  permutation_vec.resize(nb_nodes);

  // Offset denotes what is the position of the first node in the
  // permutation_vec
  // for the triangle that we are currently filling.
  Uint offset = 0;
  // Polynomial order to which the remaining interior nodes correspond
  Int current_poly_order = poly_order;

  while (current_poly_order >= 0)
  {
    // If the current triangle is of order 0, it has one center node
    // In that case, insert the center node and exit the while loop
    if (current_poly_order == _0D)
    {
      permutation_vec[offset] = offset;
    }

    // This should be executed when the current polynomial order is >= 1
    if (current_poly_order >= _1D)
    {
      permutation_vec[offset]     = offset + 1;
      permutation_vec[offset + 1] = offset + 0u;
      permutation_vec[offset + 2] = offset + 2;

      // If we have some edge-internal nodes, let's insert them in
      // permutation_vec
      if (current_poly_order >= _2D)
      {
        const Uint nb_internal_pts_on_edge = current_poly_order - 1;

        // Edge 0-1 will become edge 1-0
        for (Uint i = 0; i < nb_internal_pts_on_edge; ++i)
        {
          permutation_vec[offset + 3 + i] = offset + 2 + nb_internal_pts_on_edge - i;
        }

        // Edge 0-2
        for (Uint i = 0; i < nb_internal_pts_on_edge; ++i)
        {
          permutation_vec[offset + 3 + nb_internal_pts_on_edge + i] =
              offset + 2 + 3 * nb_internal_pts_on_edge - i;
        }

        // Edge 1-2 will become third edge: 2-1
        for (Uint i = 0; i < nb_internal_pts_on_edge; ++i)
        {
          permutation_vec[offset + 3 + 2 * nb_internal_pts_on_edge + i] =
              offset + 2 + 2 * nb_internal_pts_on_edge - i;
        }
      } // Polynomial order >= 2

    } // Polynomial order >= 1

    // For order P, we have P+1 nodes on edge and 3*P nodes on the contour
    // of the triangle
    offset += 3 * current_poly_order;
    current_poly_order -= 3;

  } // While polynomial order >= 0

} // End of method

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
