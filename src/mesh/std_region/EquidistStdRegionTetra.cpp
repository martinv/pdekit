#include <cassert>

#include "mesh/std_region/EquidistStdRegionTetra.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

/// A tetra has 1 volume, 4 triangular faces and 6 edges
const std::array<Uint, 4> EquidistStdRegionTetra::TopologyStorage = {1, 6, 4, 1};

// ----------------------------------------------------------------------------

/// A TETRAHEDRON OF POLYNOMIAL ORDER P HAS (P+1)*(P+2)*(P+3)/6 NODES
Uint EquidistStdRegionTetra::nb_dof(const Uint poly_order)
{
  return ((poly_order + 1) * (poly_order + 2) * (poly_order + 3) / 6);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTetra::fill_topology_relations(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences)
{
  // edge to node incidences - position (1,0)
  common::BlockArray<SUint, SUint> &edge_to_node = incidences[(_3D + 1) * _1D + _0D];
  // Edge 0-1           Six entities - the edges of the tetra
  // Edge 1-2           The values 0,1,2, ... on the rhs are POSITIONS
  // Edge 2-0           of reference entities as they will be stored
  // Edge 3-0           in entity vector for entities of dim 1
  // Edge 3-2
  // Edge 3-1
  edge_to_node.build({0, 1, 2, 3, 4, 5}, {6});

  // edge to edge incidences - position (1,1)
  common::BlockArray<SUint, SUint> &edge_to_edge = incidences[(_3D + 1) * _1D + _1D];
  edge_to_edge.build({1, 2, 3, 5, 0, 2, 4, 5, 0, 1, 3, 4, 0, 2, 4, 5, 1, 2, 3, 5, 0, 1, 3, 4},
                     {4, 4, 4, 4, 4, 4});

  // edge to face incidences - position (1,2)
  // In a tetrahedron, each edge is incident to all faces of the element
  common::BlockArray<SUint, SUint> &edge_to_face = incidences[(_3D + 1) * _1D + _2D];
  edge_to_face.build({0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3},
                     {4, 4, 4, 4, 4, 4});

  // face-to-node incidence should be stored at position (2,0)
  common::BlockArray<SUint, SUint> &face_to_node = incidences[(_3D + 1) * _2D + _0D];
  // Four triangular faces
  face_to_node.build({0, 1, 2, 3}, {4});

  // face-to-edge incidence - should be stored at position (2,1)
  common::BlockArray<SUint, SUint> &face_to_edge = incidences[(_3D + 1) * _2D + _1D];
  face_to_edge.build({2, 1, 0, 0, 5, 3, 3, 4, 2, 5, 1, 4}, {3, 3, 3, 3});

  // face-to-face incidence - position (2,2)
  // Each face is incident to all remaining faces
  common::BlockArray<SUint, SUint> &face_to_face = incidences[(_3D + 1) * _2D + _2D];
  face_to_face.build({1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2}, {3, 3, 3, 3});

  // volume-to-node incidence should be stored at position (3,0)
  common::BlockArray<SUint, SUint> &volume_to_node = incidences[(_3D + 1) * _3D + _0D];
  // One entity - volume of the reference element
  volume_to_node.build({0}, {1});

  // volume-to-edge incidence: (3,1)
  common::BlockArray<SUint, SUint> &volume_to_edge = incidences[(_3D + 1) * _3D + _1D];
  // Only one row (we have 1 volume) with 6 columns (for 6 edges)
  volume_to_edge.build({0, 1, 2, 3, 4, 5}, {6});

  // volume-to-face incidence is stored at position (3,2)
  common::BlockArray<SUint, SUint> &volume_to_face = incidences[(_3D + 1) * _3D + _2D];
  volume_to_face.build({0, 1, 2, 3}, {4});

  // volume-to-volume incidence should be stored at position (3,3)
  common::BlockArray<SUint, SUint> &volume_to_volume = incidences[(_3D + 1) * _3D + _3D];
  // Volume 0 (the only volume present) is incident only to itself
  // This incidence has to be filled, otherwise 'local_transform'
  // in MeshEntity will fail when the entity is transformed to
  // itself!
  volume_to_volume.build({0}, {1});
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTetra::fill_edge_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // A tetrahedron has 6 edges ...

  std::vector<Uint> vert_ids;
  std::vector<Uint> p1_vert_flags;

  // First edge: [ 0 - 4 - 5 - 6 - ....  - 1 ] - 0 and 1 are the end vertices,
  // the number 4 is assigned to the third vertex of the edge,
  // so the numbers 4, 5, 6, ... are the internal nodes of the first edge

  StdRegionEntity &edge0 = *entity_vec[0];

  vert_ids.resize(Edge01::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge01::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge01::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge01::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  auto vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  auto p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge0.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 0);

  // ---

  StdRegionEntity &edge1 = *entity_vec[1];

  vert_ids.resize(Edge12::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge12::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge12::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge12::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge1.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 1);

  // ---

  StdRegionEntity &edge2 = *entity_vec[2];

  vert_ids.resize(Edge20::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge20::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge20::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge20::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge2.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 2);

  // ---

  StdRegionEntity &edge3 = *entity_vec[3];

  vert_ids.resize(Edge30::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge30::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge30::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge30::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge3.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 3);

  // ---

  StdRegionEntity &edge4 = *entity_vec[4];

  vert_ids.resize(Edge32::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge32::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge32::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge32::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge4.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 4);

  // ---

  StdRegionEntity &edge5 = *entity_vec[5];

  vert_ids.resize(Edge31::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge31::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge31::nb_nodes(poly_order); ++i)
  {
    // const Uint reverse_idx_on_edge = Edge13::nb_nodes() - i -1;
    // edge5.set_vertex( i, Edge13::node_id(reverse_idx_on_edge) );
    vert_ids[i]      = Edge31::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge5.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 5);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTetra::fill_face_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // A tetrahedron has 4 faces ...
  // entity_vec.resize(4);
  // Set correct ids, dimension and element type of all face entities

  assert(entity_vec.size() == 4);

  std::vector<Uint> vert_ids;
  std::vector<Uint> p1_vert_flags;

  /// **********************************
  /// FACE 0-2-1
  /// **********************************
  StdRegionEntity &face0 = *entity_vec[0];

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.resize((poly_order + 1) * (poly_order + 2) / 2);

  vert_ids[0] = 0;
  vert_ids[1] = 2;
  vert_ids[2] = 1;

  p1_vert_flags.assign(p1_vert_flags.size(), 0);
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;

  Uint entry_pos = 3;
  Uint face_node_id;

  if (poly_order > _1D)
  {
    // Fill interior nodes of edge 0-2
    for (Uint i = 2; i < Edge20::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge20::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge20::node_id(poly_order, reverse_idx_on_edge);
    }

    // Fill interior nodes of edge 2-1
    for (Uint i = 2; i < Edge12::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge12::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge12::node_id(poly_order, reverse_idx_on_edge);
    }

    // Fill interior nodes of edge 1-0
    for (Uint i = 2; i < Edge01::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge01::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge01::node_id(poly_order, reverse_idx_on_edge);
    }

    if (poly_order > _2D)
    {
      face_node_id                   = 3 + 6 * (poly_order - 1) + 1;
      const Uint nb_edge_to_node_pts = poly_order + 1;
      for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 3) / 2); ++i)
      {
        vert_ids[entry_pos++] = face_node_id++;
      }
    }
  }

  auto vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  auto p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face0.construct(PointSetTag(ElemShape::Triag, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 0);

  /// **********************************
  /// FACE 0-1-3
  /// **********************************
  StdRegionEntity &face1 = *entity_vec[1];

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.resize((poly_order + 1) * (poly_order + 2) / 2);

  vert_ids[0] = 0;
  vert_ids[1] = 1;
  vert_ids[2] = 3;

  p1_vert_flags.assign(p1_vert_flags.size(), 0);
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;

  entry_pos = 3;

  if (poly_order > _1D)
  {
    // Fill interior nodes of edge 0-1
    for (Uint i = 2; i < Edge01::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge01::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 1-3
    for (Uint i = 2; i < Edge31::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge31::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge31::node_id(poly_order, reverse_idx_on_edge);
    }

    // Fill interior nodes of edge 3-0
    for (Uint i = 2; i < Edge30::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge30::node_id(poly_order, i);
    }

    if (poly_order > _2D)
    {
      //           vertices           edges         1 previous face
      face_node_id = 3 + 6 * (poly_order - 1) + (poly_order - 2) * (poly_order - 1) / 2 + 1;
      const Uint nb_edge_to_node_pts = poly_order + 1;
      for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 3) / 2); ++i)
      {
        vert_ids[entry_pos++] = face_node_id++;
      }
    }
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face1.construct(PointSetTag(ElemShape::Triag, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 1);

  /// **********************************
  /// FACE 0-3-2
  /// **********************************
  StdRegionEntity &face2 = *entity_vec[2];

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.resize((poly_order + 1) * (poly_order + 2) / 2);

  vert_ids[0] = 0;
  vert_ids[1] = 3;
  vert_ids[2] = 2;

  p1_vert_flags.assign(p1_vert_flags.size(), 0);
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;

  entry_pos = 3;

  if (poly_order > _1D)
  {
    // Fill interior nodes of edge 0-3
    for (Uint i = 2; i < Edge30::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge30::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge30::node_id(poly_order, reverse_idx_on_edge);
    }

    // Fill interior nodes of edge 3-2
    for (Uint i = 2; i < Edge32::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge32::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 2-0
    for (Uint i = 2; i < Edge20::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge20::node_id(poly_order, i);
    }

    if (poly_order > _2D)
    {
      //           vertices           edges         2 previous faces
      face_node_id = 3 + 6 * (poly_order - 1) + 2 * (poly_order - 2) * (poly_order - 1) / 2 + 1;
      const Uint nb_edge_to_node_pts = poly_order + 1;
      for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 3) / 2); ++i)
      {
        vert_ids[entry_pos++] = face_node_id++;
      }
    }
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face2.construct(PointSetTag(ElemShape::Triag, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 2);

  /// **********************************
  /// FACE 3-1-2
  /// **********************************
  StdRegionEntity &face3 = *entity_vec[3];
  // face2.resize( (poly_order+1)*(poly_order+2)/2 );

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.resize((poly_order + 1) * (poly_order + 2) / 2);

  vert_ids[0] = 3;
  vert_ids[1] = 1;
  vert_ids[2] = 2;

  p1_vert_flags.assign(p1_vert_flags.size(), 0);
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;

  entry_pos = 3;

  if (poly_order > _1D)
  {
    // Fill interior nodes of edge 3-1
    for (Uint i = 2; i < Edge31::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge31::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 1-2
    for (Uint i = 2; i < Edge12::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge12::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 2-3
    for (Uint i = 2; i < Edge32::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge32::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge32::node_id(poly_order, reverse_idx_on_edge);
    }

    if (poly_order > _2D)
    {
      //           vertices           edges         3 previous faces
      face_node_id = 3 + 6 * (poly_order - 1) + 3 * (poly_order - 2) * (poly_order - 1) / 2 + 1;
      const Uint nb_edge_to_node_pts = poly_order + 1;
      for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 3) / 2); ++i)
      {
        vert_ids[entry_pos++] = face_node_id++;
      }
    }
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face3.construct(PointSetTag(ElemShape::Triag, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 3);

} // end of method

// ----------------------------------------------------------------------------

void EquidistStdRegionTetra::fill_volume_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // entity_vec.resize(1);
  StdRegionEntity &volume = *entity_vec[0];

  const Uint nb_nodes = ((poly_order + 1) * (poly_order + 2) * (poly_order + 3) / 6);
  std::vector<Uint> vert_ids(nb_nodes);
  std::iota(vert_ids.begin(), vert_ids.end(), 0);

  std::vector<Uint> p1_vert_flags(nb_nodes);
  p1_vert_flags.assign(nb_nodes, 0);
  for (Uint n = 0; n < 4; ++n)
  {
    p1_vert_flags[n] = 1;
  }

  const auto vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  const auto p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  volume.construct(PointSetTag(ElemShape::Tetra, poly_order, PointSetID::Equidist), vert_ids_view,
                   p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTetra::fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord)
{
  const Uint nb_nodes = ((poly_order + 1) * (poly_order + 2) * (poly_order + 3) / 6);
  coord.resize(nb_nodes, _3D);

  const Real dx = 2. / poly_order;
  const Real dy = 2. / poly_order;
  const Real dz = 2. / poly_order;

  // Coordinates of local vertex zero
  Real ksi_0 = -1.;
  Real eta_0 = -1.;
  Real zta_0 = -1.;

  Real ksi_1 = 1.;
  Real eta_1 = -1.;
  Real zta_1 = -1.;

  Real ksi_2 = -1.;
  Real eta_2 = 1.;
  Real zta_2 = -1.;

  Real ksi_3 = -1.;
  Real eta_3 = -1.;
  Real zta_3 = 1.;

  // Polynomial order to which the remaining interior nodes correspond
  Int current_poly_order = poly_order;

  Uint inode = 0;
  while (current_poly_order >= 0)
  {
    if (current_poly_order == P0)
    {
      coord(inode, KSI) = ksi_0;
      coord(inode, ETA) = eta_0;
      coord(inode, ZTA) = zta_0;
      inode++;
    }

    if (current_poly_order > P0)
    {
      // First corner
      coord(inode, KSI) = ksi_0;
      coord(inode, ETA) = eta_0;
      coord(inode, ZTA) = zta_0;
      inode++;

      // Second corner
      coord(inode, KSI) = ksi_1;
      coord(inode, ETA) = eta_1;
      coord(inode, ZTA) = zta_1;
      inode++;

      // Third corner
      coord(inode, KSI) = ksi_2;
      coord(inode, ETA) = eta_2;
      coord(inode, ZTA) = zta_2;
      inode++;

      // Fourth corner
      coord(inode, KSI) = ksi_3;
      coord(inode, ETA) = eta_3;
      coord(inode, ZTA) = zta_3;
      inode++;
    }

    // If there are internal nodes on edges, let's fill them
    if (current_poly_order > P1)
    {
      const Uint nb_interior_pts_on_edge = current_poly_order - 1;

      // Edge 0-1
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_0 + (i + 1) * dx;
        coord(inode, ETA) = eta_0;
        coord(inode, ZTA) = zta_0;
        inode++;
      }

      // Edge 1-2
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_1 - (i + 1) * dx;
        coord(inode, ETA) = eta_1 + (i + 1) * dy;
        coord(inode, ZTA) = zta_1;
        inode++;
      }

      // Edge 2-0
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_2;
        coord(inode, ETA) = eta_2 - (i + 1) * dy;
        coord(inode, ZTA) = zta_2;
        inode++;
      }

      // Edge 3-0
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_3;
        coord(inode, ETA) = eta_3;
        coord(inode, ZTA) = zta_3 - (i + 1) * dz;
        inode++;
      }

      // Edge 3-2
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_3;
        coord(inode, ETA) = eta_3 + (i + 1) * dy;
        coord(inode, ZTA) = zta_3 - (i + 1) * dz;
        inode++;
      }

      // Edge 3-1
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_3 + (i + 1) * dx;
        coord(inode, ETA) = eta_3;
        coord(inode, ZTA) = zta_3 - (i + 1) * dz;
        inode++;
      }

    } // If TETRA polynomial order > 1

    if (current_poly_order > P2)
    {
      Real tri_ksi_0, tri_eta_0, tri_zta_0, tri_ksi_1, tri_eta_1, tri_zta_1, tri_ksi_2, tri_eta_2,
          tri_zta_2;

      Int face_poly_order = current_poly_order - 3;

      // ------------------------------
      // Node coordinates of face 0-2-1
      // ------------------------------
      tri_ksi_0 = ksi_0 + dx;
      tri_eta_0 = eta_0 + dy;
      tri_zta_0 = zta_0;

      tri_ksi_1 = ksi_2 + dx;
      tri_eta_1 = eta_2 - 2. * dy;
      tri_zta_1 = zta_2;

      tri_ksi_2 = ksi_1 - 2. * dx;
      tri_eta_2 = eta_1 + dy;
      tri_zta_2 = zta_1;

      while (face_poly_order >= 0)
      {

        if (face_poly_order == P0)
        {
          coord(inode, KSI) = tri_ksi_0;
          coord(inode, ETA) = tri_eta_0;
          coord(inode, ZTA) = tri_zta_0;
          inode++;
        }

        if (face_poly_order > P0)
        {
          coord(inode, KSI) = tri_ksi_0;
          coord(inode, ETA) = tri_eta_0;
          coord(inode, ZTA) = tri_zta_0;
          inode++;

          coord(inode, KSI) = tri_ksi_1;
          coord(inode, ETA) = tri_eta_1;
          coord(inode, ZTA) = tri_zta_1;
          inode++;

          coord(inode, KSI) = tri_ksi_2;
          coord(inode, ETA) = tri_eta_2;
          coord(inode, ZTA) = tri_zta_2;
          inode++;
        }

        if (face_poly_order > P1)
        {
          // Edge 0-1 of the local triangle embedded in tetra face
          // 0-2-1
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_0;
            coord(inode, ETA) = tri_eta_0 + i * dy;
            coord(inode, ZTA) = tri_zta_0;
            inode++;
          }
          // Edge 1-2 of the local triangle embedded in tetra face
          // 0-2-1
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_1 + i * dx;
            coord(inode, ETA) = tri_eta_1 - i * dy;
            coord(inode, ZTA) = tri_zta_1;
            inode++;
          }

          // Edge 2-0 of the local triangle embedded in tetra face
          // 0-2-1
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_2 - i * dx;
            coord(inode, ETA) = tri_eta_2;
            coord(inode, ZTA) = tri_zta_2;
            inode++;
          }
        }

        tri_ksi_0 += dx;
        tri_eta_0 += dy;
        // tri_zta_0 remains unchanged

        tri_ksi_1 += dx;      //-= 2. * dx;
        tri_eta_1 -= 2. * dy; //+= dy;
        // tri_zta_1 remains unchanged

        tri_ksi_2 -= 2. * dx; // += dx;
        tri_eta_2 += dy;      // -= 2. * dy;
        // tri_zta_2 remains unchanged

        face_poly_order -= 3;
      }

      // ------------------------------
      // Node coordinates of face 0-1-3
      // ------------------------------
      face_poly_order = current_poly_order - 3;

      tri_ksi_0 = ksi_0 + dx;
      tri_eta_0 = eta_0;
      tri_zta_0 = zta_0 + dz;

      tri_ksi_1 = ksi_1 - 2. * dx;
      tri_eta_1 = eta_1;
      tri_zta_1 = zta_1 + dz;

      tri_ksi_2 = ksi_3 + dx;
      tri_eta_2 = eta_3;
      tri_zta_2 = zta_3 - 2. * dz;

      while (face_poly_order >= 0)
      {

        if (face_poly_order == P0)
        {
          coord(inode, KSI) = tri_ksi_0;
          coord(inode, ETA) = tri_eta_0;
          coord(inode, ZTA) = tri_zta_0;
          inode++;
        }

        if (face_poly_order > P0)
        {
          coord(inode, KSI) = tri_ksi_0;
          coord(inode, ETA) = tri_eta_0;
          coord(inode, ZTA) = tri_zta_0;
          inode++;

          coord(inode, KSI) = tri_ksi_1;
          coord(inode, ETA) = tri_eta_1;
          coord(inode, ZTA) = tri_zta_1;
          inode++;

          coord(inode, KSI) = tri_ksi_2;
          coord(inode, ETA) = tri_eta_2;
          coord(inode, ZTA) = tri_zta_2;
          inode++;
        }

        if (face_poly_order > P1)
        {
          // Edge 0-1 of the local triangle embedded in tetra face
          // 0-1-3
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_0 + i * dx;
            coord(inode, ETA) = tri_eta_0;
            coord(inode, ZTA) = tri_zta_0;
            inode++;
          }
          // Edge 1-2 of the local triangle embedded in tetra face
          // 0-1-3
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_1 - i * dx;
            coord(inode, ETA) = tri_eta_1;
            coord(inode, ZTA) = tri_zta_1 + i * dz;
            inode++;
          }

          // Edge 2-0 of the local triangle embedded in tetra face
          // 0-1-3
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_2;
            coord(inode, ETA) = tri_eta_2;
            coord(inode, ZTA) = tri_zta_2 - i * dz;
            inode++;
          }
        }

        tri_ksi_0 += dx;
        // tri_eta_0 remains unchanged
        tri_zta_0 += dz;

        tri_ksi_1 -= 2. * dx;
        // tri_eta_1 remains unchanged
        tri_zta_1 += dz;

        tri_ksi_2 += dx;
        // tri_eta_2 remains unchanged
        tri_zta_2 -= 2. * dz;

        face_poly_order -= 3;
      }

      // ------------------------------
      // Node coordinates of face 0-3-2
      // ------------------------------
      face_poly_order = current_poly_order - 3;

      tri_ksi_0 = ksi_0;
      tri_eta_0 = eta_0 + dy;
      tri_zta_0 = zta_0 + dz;

      tri_ksi_1 = ksi_3;
      tri_eta_1 = eta_3 + dy;
      tri_zta_1 = zta_3 - 2. * dz;

      tri_ksi_2 = ksi_2;
      tri_eta_2 = eta_2 - 2. * dy;
      tri_zta_2 = zta_2 + dz;

      while (face_poly_order >= 0)
      {

        if (face_poly_order == P0)
        {
          coord(inode, KSI) = tri_ksi_0;
          coord(inode, ETA) = tri_eta_0;
          coord(inode, ZTA) = tri_zta_0;
          inode++;
        }

        if (face_poly_order > P0)
        {
          coord(inode, KSI) = tri_ksi_0;
          coord(inode, ETA) = tri_eta_0;
          coord(inode, ZTA) = tri_zta_0;
          inode++;

          coord(inode, KSI) = tri_ksi_1;
          coord(inode, ETA) = tri_eta_1;
          coord(inode, ZTA) = tri_zta_1;
          inode++;

          coord(inode, KSI) = tri_ksi_2;
          coord(inode, ETA) = tri_eta_2;
          coord(inode, ZTA) = tri_zta_2;
          inode++;
        }

        if (face_poly_order > P1)
        {
          // Edge 0-1 of the local triangle embedded in tetra face
          // 0-3-2
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_0;
            coord(inode, ETA) = tri_eta_0;
            coord(inode, ZTA) = tri_zta_0 + i * dz;
            inode++;
          }
          // Edge 1-2 of the local triangle embedded in tetra face
          // 0-3-2
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_1;
            coord(inode, ETA) = tri_eta_1 + i * dy;
            coord(inode, ZTA) = tri_zta_1 - i * dz;
            inode++;
          }

          // Edge 2-0 of the local triangle embedded in tetra face
          // 0-3-2
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_2;
            coord(inode, ETA) = tri_eta_2 - i * dy;
            coord(inode, ZTA) = tri_zta_2;
            inode++;
          }
        }

        // tri_ksi_0 remains unchanged
        tri_eta_0 += dy;
        tri_zta_0 += dz;

        // tri_ksi_1 remains unchanged
        tri_eta_1 += dy;
        tri_zta_1 -= 2. * dz;

        // tri_ksi_2 remains unchanged
        tri_eta_2 -= 2. * dy;
        tri_zta_2 += dz;

        face_poly_order -= 3;
      }

      // ------------------------------
      // Node coordinates of face 3-1-2
      // ------------------------------
      face_poly_order = current_poly_order - 3;

      tri_ksi_0 = ksi_3 + dx;
      tri_eta_0 = eta_3 + dy;
      tri_zta_0 = zta_3 - 2. * dz;

      tri_ksi_1 = ksi_1 - 2. * dx;
      tri_eta_1 = eta_1 + dy;
      tri_zta_1 = zta_1 + dz;

      tri_ksi_2 = ksi_2 + dx;
      tri_eta_2 = eta_2 - 2. * dy;
      tri_zta_2 = zta_2 + dz;

      while (face_poly_order >= 0)
      {

        if (face_poly_order == P0)
        {
          coord(inode, KSI) = tri_ksi_0;
          coord(inode, ETA) = tri_eta_0;
          coord(inode, ZTA) = tri_zta_0;
          inode++;
        }

        if (face_poly_order > P0)
        {
          coord(inode, KSI) = tri_ksi_0;
          coord(inode, ETA) = tri_eta_0;
          coord(inode, ZTA) = tri_zta_0;
          inode++;

          coord(inode, KSI) = tri_ksi_1;
          coord(inode, ETA) = tri_eta_1;
          coord(inode, ZTA) = tri_zta_1;
          inode++;

          coord(inode, KSI) = tri_ksi_2;
          coord(inode, ETA) = tri_eta_2;
          coord(inode, ZTA) = tri_zta_2;
          inode++;
        }

        if (face_poly_order > P1)
        {
          // Edge 0-1 of the local triangle embedded in tetra face
          // 3-1-2
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_0 + i * dx;
            coord(inode, ETA) = tri_eta_0;
            coord(inode, ZTA) = tri_zta_0 - i * dz;
            inode++;
          }
          // Edge 1-2 of the local triangle embedded in tetra face
          // 3-1-2
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_1 - i * dx;
            coord(inode, ETA) = tri_eta_1 + i * dy;
            coord(inode, ZTA) = tri_zta_1;
            inode++;
          }

          // Edge 2-0 of the local triangle embedded in tetra face
          // 0-3-2
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_2;
            coord(inode, ETA) = tri_eta_2 - i * dy;
            coord(inode, ZTA) = tri_zta_2 + i * dz;
            inode++;
          }
        }

        tri_ksi_0 += dx;
        tri_eta_0 += dy;
        tri_zta_0 -= 2. * dz;

        tri_ksi_1 -= 2. * dx;
        tri_eta_1 += dy;
        tri_zta_1 += dz;

        tri_ksi_2 += dx;
        tri_eta_2 -= 2. * dy;
        tri_zta_2 += dz;

        face_poly_order -= 3;
      }

    } // If TETRA polynomial order > 2

    ksi_0 += dx;
    eta_0 += dy;
    zta_0 += dz;

    ksi_1 -= 3. * dx;
    eta_1 += dy;
    zta_1 += dz;

    ksi_2 += dx;
    eta_2 -= 3. * dy;
    zta_2 += dz;

    ksi_3 += dx;
    eta_3 += dy;
    zta_3 -= 3. * dz;

    current_poly_order -= 4;

  } // While TETRA polyorder >= 0
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTetra::fill_facet_normals(math::DenseDMat<Real> &normals)
{
  normals.resize(4, 3);

  // Face 0-2-1
  normals(0, XI0) = 0.0;
  normals(0, XI1) = 0.0;
  normals(0, XI2) = -1.0;

  // Face 0-1-3
  normals(1, XI0) = 0.0;
  normals(1, XI1) = -1.0;
  normals(1, XI2) = 0.0;

  // Face 0-3-2
  normals(2, XI0) = -1.0;
  normals(2, XI1) = 0.0;
  normals(2, XI2) = 0.0;

  // Face 3-1-2
  normals(3, XI0) = 1.0 / std::sqrt(3.);
  normals(3, XI1) = 1.0 / std::sqrt(3.);
  normals(3, XI2) = 1.0 / std::sqrt(3.);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionTetra::fill_permutation(const Uint poly_order,
                                              const EntityRealignCode &permutation_code,
                                              std::vector<Uint> &permutation_vec)
{
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
