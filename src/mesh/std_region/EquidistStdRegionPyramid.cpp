#include <cassert>

#include "mesh/std_region/EquidistStdRegionPyramid.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

/// A Pyramid has 1 volume, 5 faces and 8 edges
const std::array<Uint, 4> EquidistStdRegionPyramid::TopologyStorage = {1, 8, 5, 1};

// ----------------------------------------------------------------------------

/// A PYRAMID OF POLYNOMIAL ORDER P HAS (P+1)*(P+2)*(2*P+3)/6 NODES
// ----------------------------------------------------------------------------

Uint EquidistStdRegionPyramid::nb_dof(const Uint poly_order)
{
  return ((poly_order + 1) * (poly_order + 2) * (2 * poly_order + 3) / 6);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionPyramid::fill_topology_relations(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences)
{
  // edge to node incidences - position (1,0)
  common::BlockArray<SUint, SUint> &edge_to_node = incidences[(_3D + 1) * _1D + _0D];
  // Edge 0-1           Eight entities - the edges of the pyramid
  // Edge 0-3           The values 0,1,2, ... on the rhs are POSITIONS
  // Edge 0-4           of reference entities as they will be stored
  // Edge 1-2           in entity vector for entities of dim 1
  // Edge 1-4
  // Edge 2-3
  // Edge 2-4
  // Edge 3-4
  edge_to_node.build({0, 1, 2, 3, 4, 5, 6, 7}, {8});

  // edge to edge incidences - position (1,1)
  common::BlockArray<SUint, SUint> &edge_to_edge = incidences[(_3D + 1) * _1D + _1D];
  edge_to_edge.build({1, 2, 3, 4, 0, 2, 5, 7, 0, 1, 4, 6, 7, 0, 4, 5, 6, 0,
                      3, 2, 6, 7, 3, 6, 1, 7, 3, 5, 2, 4, 7, 1, 5, 2, 4, 6},
                     {4, 4, 5, 4, 5, 4, 5, 5});

  // edge to face incidences - position (1,2)
  // In a pyramid, each edge is incident either to 4 faces or to all 5 faces
  common::BlockArray<SUint, SUint> &edge_to_face = incidences[(_3D + 1) * _1D + _2D];
  edge_to_face.build({0, 1, 2, 4, 0, 1, 3, 4, 0, 1, 2, 3, 4, 0, 2, 3, 4, 0,
                      1, 2, 3, 4, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4},
                     {4, 4, 5, 4, 5, 4, 5, 5});

  // face-to-node incidence should be stored at position (2,0)
  common::BlockArray<SUint, SUint> &face_to_node = incidences[(_3D + 1) * _2D + _0D];
  // Four triangular faces
  face_to_node.build({0, 1, 2, 3, 4}, {5});

  // face-to-edge incidence - should be stored at position (2,1)
  common::BlockArray<SUint, SUint> &face_to_edge = incidences[(_3D + 1) * _2D + _1D];
  face_to_edge.build({0, 4, 2, 1, 2, 7, 3, 6, 4, 5, 7, 6, 0, 1, 5, 3}, {3, 3, 3, 3, 4});

  // face-to-face incidence - position (2,2)
  // Each face is incident to all remaining faces
  common::BlockArray<SUint, SUint> &face_to_face = incidences[(_3D + 1) * _2D + _2D];
  face_to_face.build({1, 2, 3, 4, 0, 2, 3, 4, 0, 1, 3, 4, 0, 1, 2, 4, 0, 1, 2, 3}, {4, 4, 4, 4, 4});

  // volume-to-node incidence should be stored at position (3,0)
  common::BlockArray<SUint, SUint> &volume_to_node = incidences[(_3D + 1) * _3D + _0D];
  // One entity - volume of the reference element
  volume_to_node.build({0}, {1});

  // volume-to-edge incidence: (3,1)
  common::BlockArray<SUint, SUint> &volume_to_edge = incidences[(_3D + 1) * _3D + _1D];
  // Only one row (we have 1 volume) with 8 columns (for 8 edges)
  volume_to_edge.build({0, 1, 2, 3, 4, 5, 6, 7}, {8});

  // volume-to-face incidence is stored at position (3,2)
  common::BlockArray<SUint, SUint> &volume_to_face = incidences[(_3D + 1) * _3D + _2D];
  volume_to_face.build({0, 1, 2, 3, 4}, {5});

  // volume-to-volume incidence should be stored at position (3,3)
  common::BlockArray<SUint, SUint> &volume_to_volume = incidences[(_3D + 1) * _3D + _3D];
  // Volume 0 (the only volume present) is incident only to itself
  // This incidence has to be filled, otherwise 'local_transform'
  // in MeshEntity will fail when the entity is transformed to
  // itself!
  volume_to_volume.build({0}, {1});
}

// ----------------------------------------------------------------------------

void EquidistStdRegionPyramid::fill_edge_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // A pyramid has 8 edges ...
  // entity_vec.resize(8);

  std::vector<Uint> vert_ids;
  std::vector<Uint> p1_vert_flags;

  // First edge: [ 0 - 4 - 5 - 6 - ....  - 1 ] - 0 and 1 are the end vertices,
  // the number 4 is assigned to the third vertex of the edge,
  // so the numbers 4, 5, 6, ... are the internal nodes of the first edge

  StdRegionEntity &edge0 = *entity_vec[0];
  // edge0.resize(poly_order+1);

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

  vert_ids.resize(Edge03::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge03::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge03::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge03::node_id(poly_order, i);
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

  vert_ids.resize(Edge04::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge04::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge04::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge04::node_id(poly_order, i);
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

  edge3.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 3);

  // ---

  StdRegionEntity &edge4 = *entity_vec[4];

  vert_ids.resize(Edge14::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge14::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge14::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge14::node_id(poly_order, i);
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

  vert_ids.resize(Edge23::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge23::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge23::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge23::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge5.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 5);

  // ---

  StdRegionEntity &edge6 = *entity_vec[6];

  vert_ids.resize(Edge24::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge24::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge24::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge24::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge6.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 6);

  /// ---

  StdRegionEntity &edge7 = *entity_vec[7];

  vert_ids.resize(Edge34::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge34::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge34::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge34::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge7.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 7);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionPyramid::fill_face_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // A pyramid has 5 faces ...
  // entity_vec.resize(5);
  // Set correct ids, dimension and element type of all face entities

  assert(entity_vec.size() == 5);
  // The first four faces are triangles

  std::vector<Uint> vert_ids;
  std::vector<Uint> p1_vert_flags;

  /// **********************************
  /// FACE 0-1-4
  /// **********************************
  StdRegionEntity &face0 = *entity_vec[0];
  // face0.resize( (poly_order+1)*(poly_order+2)/2 );

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.assign(p1_vert_flags.size(), 0);

  vert_ids[0] = 0;
  vert_ids[1] = 1;
  vert_ids[2] = 4;

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;

  Uint entry_pos = 3;
  Uint face_node_id;

  if (poly_order > _1D)
  {

    // Fill interior nodes of edge 0-1
    for (Uint i = 2; i < Edge01::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge01::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 1-4
    for (Uint i = 2; i < Edge14::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge14::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 4-0
    for (Uint i = 2; i < Edge04::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge04::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge04::node_id(poly_order, reverse_idx_on_edge);
    }

    if (poly_order > _2D)
    {
      face_node_id                   = 4 + 8 * (poly_order - 1) + 1;
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
  /// FACE 3-0-4
  /// **********************************
  StdRegionEntity &face1 = *entity_vec[1];

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.assign(p1_vert_flags.size(), 0);

  vert_ids[0] = 3;
  vert_ids[1] = 0;
  vert_ids[2] = 4;

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;

  entry_pos = 3;

  if (poly_order > _1D)
  {
    // Fill interior nodes of edge 3-0
    for (Uint i = 2; i < Edge03::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge03::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge03::node_id(poly_order, reverse_idx_on_edge);
    }

    // Fill interior nodes of edge 0-4
    for (Uint i = 2; i < Edge04::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge04::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 4-3
    for (Uint i = 2; i < Edge34::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge34::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge34::node_id(poly_order, reverse_idx_on_edge);
    }

    if (poly_order > _2D)
    {
      //           vertices           edges         1 previous face
      face_node_id = 4 + 8 * (poly_order - 1) + (poly_order - 2) * (poly_order - 1) / 2 + 1;
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
  /// FACE 1-2-4
  /// **********************************
  StdRegionEntity &face2 = *entity_vec[2];

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.assign(p1_vert_flags.size(), 0);

  vert_ids[0] = 1;
  vert_ids[1] = 2;
  vert_ids[2] = 4;

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;

  entry_pos = 3;

  if (poly_order > _1D)
  {
    // Fill interior nodes of edge 1-2
    for (Uint i = 2; i < Edge12::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge12::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 2-4
    for (Uint i = 2; i < Edge24::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge24::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 4-1
    for (Uint i = 2; i < Edge14::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge14::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge14::node_id(poly_order, reverse_idx_on_edge);
    }

    if (poly_order > _2D)
    {
      //           vertices           edges         2 previous faces
      face_node_id = 4 + 8 * (poly_order - 1) + 2 * (poly_order - 2) * (poly_order - 1) / 2 + 1;
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
  /// FACE 2-3-4
  /// **********************************
  StdRegionEntity &face3 = *entity_vec[3];
  // face3.resize( (poly_order+1)*(poly_order+2)/2 );

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.resize((poly_order + 1) * (poly_order + 2) / 2);
  p1_vert_flags.assign(p1_vert_flags.size(), 0);

  vert_ids[0] = 2;
  vert_ids[1] = 3;
  vert_ids[2] = 4;

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;

  entry_pos = 3;

  if (poly_order > _1D)
  {
    // Fill interior nodes of edge 2-3
    for (Uint i = 2; i < Edge23::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge23::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 3-4
    for (Uint i = 2; i < Edge34::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge34::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 4-2
    for (Uint i = 2; i < Edge24::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge24::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge24::node_id(poly_order, reverse_idx_on_edge);
    }

    if (poly_order > _2D)
    {
      //           vertices           edges         3 previous faces
      face_node_id = 4 + 8 * (poly_order - 1) + 3 * (poly_order - 2) * (poly_order - 1) / 2 + 1;
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

  /// **********************************
  /// FACE 0-3-2-1
  /// **********************************
  StdRegionEntity &face4 = *entity_vec[4];
  // face2.resize( (poly_order+1)*(poly_order+1) );

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids.resize((poly_order + 1) * (poly_order + 1));
  p1_vert_flags.resize((poly_order + 1) * (poly_order + 1));
  p1_vert_flags.assign(p1_vert_flags.size(), 0);

  vert_ids[0] = 0;
  vert_ids[1] = 3;
  vert_ids[2] = 2;
  vert_ids[3] = 1;

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;
  p1_vert_flags[3] = 1;

  entry_pos = 4;

  if (poly_order > _1D)
  {
    // Fill interior nodes of edge 0-3
    for (Uint i = 2; i < Edge03::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge03::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 3-2
    for (Uint i = 2; i < Edge23::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge23::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge23::node_id(poly_order, reverse_idx_on_edge);
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
      //           vertices           edges         4 previous faces
      face_node_id = 4 + 8 * (poly_order - 1) + 4 * (poly_order - 2) * (poly_order - 1) / 2 + 1;
      const Uint nb_edge_to_node_pts = poly_order + 1;
      for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 2)); ++i)
      {
        vert_ids[entry_pos++] = face_node_id++;
      }
    }
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face4.construct(PointSetTag(ElemShape::Quad, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 4);

} // end of method

// ----------------------------------------------------------------------------

void EquidistStdRegionPyramid::fill_volume_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // entity_vec.resize(1);
  StdRegionEntity &volume = *entity_vec[0];

  const Uint nb_nodes = ((poly_order + 1) * (poly_order + 2) * (2 * poly_order + 3) / 6);

  std::vector<Uint> vert_ids(nb_nodes);
  std::iota(vert_ids.begin(), vert_ids.end(), 0);

  std::vector<Uint> p1_vert_flags(nb_nodes);
  p1_vert_flags.assign(nb_nodes, 0);
  for (Uint i = 0; i < 5; ++i)
  {
    p1_vert_flags[i] = 1;
  }

  const auto vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  const auto p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  volume.construct(PointSetTag(ElemShape::Pyramid, poly_order, PointSetID::Equidist), vert_ids_view,
                   p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionPyramid::fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord)
{
  const Uint nb_nodes = ((poly_order + 1) * (poly_order + 2) * (2 * poly_order + 3) / 6);
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

  Real ksi_2 = 1.;
  Real eta_2 = 1.;
  Real zta_2 = -1.;

  Real ksi_3 = -1.;
  Real eta_3 = 1.;
  Real zta_3 = -1.;

  Real ksi_4 = -1.;
  Real eta_4 = -1.;
  Real zta_4 = 1.;

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

      // Fifth corner
      coord(inode, KSI) = ksi_4;
      coord(inode, ETA) = eta_4;
      coord(inode, ZTA) = zta_4;
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

      // Edge 0-3
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_0;
        coord(inode, ETA) = eta_0 + (i + 1) * dy;
        coord(inode, ZTA) = zta_0;
        inode++;
      }

      // Edge 0-4
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_0;
        coord(inode, ETA) = eta_0;
        coord(inode, ZTA) = zta_0 + (i + 1) * dz;
        inode++;
      }

      // Edge 1-2
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_1;
        coord(inode, ETA) = eta_1 + (i + 1) * dy;
        coord(inode, ZTA) = zta_1;
        inode++;
      }

      // Edge 1-4
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_1 - (i + 1) * dx;
        coord(inode, ETA) = eta_1;
        coord(inode, ZTA) = zta_1 + (i + 1) * dz;
        inode++;
      }

      // Edge 2-3
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_2 - (i + 1) * dx;
        coord(inode, ETA) = eta_2;
        coord(inode, ZTA) = zta_2;
        inode++;
      }

      // Edge 2-4
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_2 - (i + 1) * dx;
        coord(inode, ETA) = eta_2 - (i + 1) * dy;
        coord(inode, ZTA) = zta_2 + (i + 1) * dz;
        inode++;
      }

      // Edge 3-4
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_3;
        coord(inode, ETA) = eta_3 - (i + 1) * dy;
        coord(inode, ZTA) = zta_3 + (i + 1) * dz;
        inode++;
      }

    } // If Pyramid polynomial order > 1

    if (current_poly_order > P2)
    {
      Real tri_ksi_0, tri_eta_0, tri_zta_0, tri_ksi_1, tri_eta_1, tri_zta_1, tri_ksi_2, tri_eta_2,
          tri_zta_2;

      Int face_poly_order = current_poly_order - 3;

      // ------------------------------
      // Node coordinates of face 0-1-4
      // ------------------------------
      tri_ksi_0 = ksi_0 + dx;
      tri_eta_0 = eta_0;
      tri_zta_0 = zta_0 + dz;

      tri_ksi_1 = ksi_1 - 2. * dx;
      tri_eta_1 = eta_1;
      tri_zta_1 = zta_1 + dz;

      tri_ksi_2 = ksi_4 + dx;
      tri_eta_2 = eta_4;
      tri_zta_2 = zta_4 - 2 * dz;

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
          // Edge 0-1 of the local triangle embedded in pyramid face
          // 0-1-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_0 + i * dx;
            coord(inode, ETA) = tri_eta_0;
            coord(inode, ZTA) = tri_zta_0;
            inode++;
          }
          // Edge 1-2 of the local triangle embedded in pyramid face
          // 0-1-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_1 - i * dx;
            coord(inode, ETA) = tri_eta_1;
            coord(inode, ZTA) = tri_zta_1 + i * dz;
            inode++;
          }

          // Edge 2-0 of the local triangle embedded in pyramid face
          // 0-1-4
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
      // Node coordinates of face 3-0-4
      // ------------------------------
      face_poly_order = current_poly_order - 3;

      tri_ksi_0 = ksi_3;
      tri_eta_0 = eta_3 - 2. * dy;
      tri_zta_0 = zta_3 + dz;

      tri_ksi_1 = ksi_0;
      tri_eta_1 = eta_0 + dy;
      tri_zta_1 = zta_0 + dz;

      tri_ksi_2 = ksi_4;
      tri_eta_2 = eta_4 + dy;
      tri_zta_2 = zta_4 - 2. * dz;

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
          // Edge 0-1 of the local triangle embedded in pyramid face
          // 3-0-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_0;
            coord(inode, ETA) = tri_eta_0 - i * dy;
            coord(inode, ZTA) = tri_zta_0;
            inode++;
          }
          // Edge 1-2 of the local triangle embedded in pyramid face
          // 3-0-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_1;
            coord(inode, ETA) = tri_eta_1;
            coord(inode, ZTA) = tri_zta_1 + i * dz;
            inode++;
          }

          // Edge 2-0 of the local triangle embedded in pyramid face
          // 3-0-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_2;
            coord(inode, ETA) = tri_eta_2 + i * dy;
            coord(inode, ZTA) = tri_zta_2 - i * dz;
            inode++;
          }
        }

        // tri_ksi_0 remains unchanged
        tri_eta_0 -= 2. * dy;
        tri_zta_0 += dz;

        // tri_ksi_1 remains unchanged
        tri_eta_1 += dy;
        tri_zta_1 += dz;

        // tri_ksi_2 remains unchanged
        tri_eta_2 += dy;
        tri_zta_2 -= 2. * dz;

        face_poly_order -= 3;
      }

      // ------------------------------
      // Node coordinates of face 1-2-4
      // ------------------------------
      face_poly_order = current_poly_order - 3;

      tri_ksi_0 = ksi_1 - dx;
      tri_eta_0 = eta_1 + dy;
      tri_zta_0 = zta_1 + dz;

      tri_ksi_1 = ksi_2 - dx;
      tri_eta_1 = eta_2 - 2. * dy;
      tri_zta_1 = zta_2 + dz;

      tri_ksi_2 = ksi_4 + 2. * dx;
      tri_eta_2 = eta_4 + dy;
      tri_zta_2 = zta_4 - 2. * dz;

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
          // Edge 0-1 of the local triangle embedded in pyramid face
          // 1-2-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_0;
            coord(inode, ETA) = tri_eta_0 + i * dy;
            coord(inode, ZTA) = tri_zta_0;
            inode++;
          }
          // Edge 1-2 of the local triangle embedded in pyramid face
          // 1-2-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_1 - i * dx;
            coord(inode, ETA) = tri_eta_1 - i * dy;
            coord(inode, ZTA) = tri_zta_1 + i * dz;
            inode++;
          }

          // Edge 2-0 of the local triangle embedded in pyramid face
          // 1-2-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_2 + i * dx;
            coord(inode, ETA) = tri_eta_2;
            coord(inode, ZTA) = tri_zta_2 - i * dz;
            inode++;
          }
        }

        tri_ksi_0 -= dx;
        tri_eta_0 += dy;
        tri_zta_0 += dz;

        tri_ksi_1 -= dx;
        tri_eta_1 -= 2. * dy;
        tri_zta_1 += dz;

        tri_ksi_2 += dx;
        tri_eta_2 += dy;
        tri_zta_2 -= 2. * dz;

        face_poly_order -= 3;
      }

      // ------------------------------
      // Node coordinates of face 2-3-4
      // ------------------------------
      face_poly_order = current_poly_order - 3;

      tri_ksi_0 = ksi_2 - 2. * dx;
      tri_eta_0 = eta_2 - dy;
      tri_zta_0 = zta_2 + dz;

      tri_ksi_1 = ksi_3 + dx;
      tri_eta_1 = eta_3 - dy;
      tri_zta_1 = zta_3 + dz;

      tri_ksi_2 = ksi_4 + dx;
      tri_eta_2 = eta_4 + 2. * dy;
      tri_zta_2 = zta_4 - 2. * dz;

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
          // Edge 0-1 of the local triangle embedded in pyramid face
          // 2-3-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_0 - i * dx;
            coord(inode, ETA) = tri_eta_0;
            coord(inode, ZTA) = tri_zta_0;
            inode++;
          }
          // Edge 1-2 of the local triangle embedded in pyramid face
          // 2-3-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_1;
            coord(inode, ETA) = tri_eta_1 - i * dy;
            coord(inode, ZTA) = tri_zta_1 + i * dz;
            inode++;
          }

          // Edge 2-0 of the local triangle embedded in pyramid face
          // 2-3-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = tri_ksi_2 + i * dx;
            coord(inode, ETA) = tri_eta_2 + i * dy;
            coord(inode, ZTA) = tri_zta_2 - i * dz;
            inode++;
          }
        }

        tri_ksi_0 -= 2. * dx;
        tri_eta_0 -= dy;
        tri_zta_0 += dz;

        tri_ksi_1 += dx;
        tri_eta_1 -= dy;
        tri_zta_1 += dz;

        tri_ksi_2 += dx;
        tri_eta_2 += dy;
        tri_zta_2 -= 2. * dz;

        face_poly_order -= 3;
      }

    } // If PYRAMID polynomial order > 2

    // ***********************************************************
    // THE CONDITION FOR BASE QUADRILATERAL FACE HAS TO BE
    // TREATED SEPARATELY: P2 TRIANGLE FACES HAVE NO INTERNAL
    // NODES, BUT P2 QUAD (THE BASE FACE) HAS ONE INTERIOR NODE
    // THIS IS THE REASON WHILE FOR THE OTHER FACES, THE
    // CONDITION WAS 'if (current_poly_order > P2)', BUT HERE
    // IT IS RATHER 'if (current_poly_order >= P2)'
    // ***********************************************************

    // --------------------------------
    // Node coordinates of face 0-3-2-1
    // --------------------------------

    if (current_poly_order >= P2)
    {
      Real face_ksi_0, face_eta_0, face_zta_0, face_ksi_1, face_eta_1, face_zta_1, face_ksi_2,
          face_eta_2, face_zta_2, face_ksi_3, face_eta_3, face_zta_3;

      Int face_poly_order = current_poly_order - 2;
      face_ksi_0          = ksi_0 + dx;
      face_eta_0          = eta_0 + dy;
      face_zta_0          = zta_0;

      face_ksi_1 = ksi_3 + dx;
      face_eta_1 = eta_3 - dy;
      face_zta_1 = zta_3;

      face_ksi_2 = ksi_2 - dx;
      face_eta_2 = eta_2 - dy;
      face_zta_2 = zta_2;

      face_ksi_3 = ksi_1 - dx;
      face_eta_3 = eta_1 + dy;
      face_zta_3 = zta_1;

      while (face_poly_order >= 0)
      {

        if (face_poly_order == P0)
        {
          coord(inode, KSI) = face_ksi_0;
          coord(inode, ETA) = face_eta_0;
          coord(inode, ZTA) = face_zta_0;
          inode++;
        }

        if (face_poly_order > P0)
        {
          coord(inode, KSI) = face_ksi_0;
          coord(inode, ETA) = face_eta_0;
          coord(inode, ZTA) = face_zta_0;
          inode++;

          coord(inode, KSI) = face_ksi_1;
          coord(inode, ETA) = face_eta_1;
          coord(inode, ZTA) = face_zta_1;
          inode++;

          coord(inode, KSI) = face_ksi_2;
          coord(inode, ETA) = face_eta_2;
          coord(inode, ZTA) = face_zta_2;
          inode++;

          coord(inode, KSI) = face_ksi_3;
          coord(inode, ETA) = face_eta_3;
          coord(inode, ZTA) = face_zta_3;
          inode++;
        }

        if (face_poly_order > P1)
        {
          // Edge 0-1 of the local quadrilateral embedded in quad face
          // 0-3-2-1
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_0;
            coord(inode, ETA) = face_eta_0 + i * dy;
            coord(inode, ZTA) = face_zta_0;
            inode++;
          }

          // Edge 1-2 of the local quadrilateral embedded in quad face
          // 0-3-2-1
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_1 + i * dx;
            coord(inode, ETA) = face_eta_1;
            coord(inode, ZTA) = face_zta_1;
            inode++;
          }

          // Edge 2-3 of the local quadrilateral embedded in quad face
          // 0-3-2-1
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_2;
            coord(inode, ETA) = face_eta_2 - i * dy;
            coord(inode, ZTA) = face_zta_2;
            inode++;
          }

          // Edge 3-0 of the local quadrilateral embedded in quad face
          // 0-3-2-1
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_3 - i * dx;
            coord(inode, ETA) = face_eta_3;
            coord(inode, ZTA) = face_zta_3;
            inode++;
          }
        }

        face_ksi_0 += dx;
        face_eta_0 += dy;
        // face_zta_0 remains unchanged

        face_ksi_1 += dx;
        face_eta_1 -= dy;
        // face_zta_1 remains unchanged

        face_ksi_2 -= dx;
        face_eta_2 -= dy;
        // face_zta_2 remains unchanged

        face_ksi_3 -= dx;
        face_eta_3 += dy;
        // face_zta_3 remains unchanged

        face_poly_order -= 2;
      }

    } // If current_poly_order >= 2 (quad base of pyramid)

    ksi_0 += dx;
    eta_0 += dy;
    zta_0 += dz;

    ksi_1 -= 2. * dx;
    eta_1 += dy;
    zta_1 += dz;

    ksi_2 -= 2. * dx;
    eta_2 -= 2. * dy;
    zta_2 += dz;

    ksi_3 += dx;
    eta_3 -= 2. * dy;
    zta_3 += dz;

    ksi_4 += dx;
    eta_4 += dy;
    zta_4 -= 2. * dz;

    current_poly_order -= 3;

  } // While PYRAMID polyorder >= 0
}

// ----------------------------------------------------------------------------

void EquidistStdRegionPyramid::fill_facet_normals(math::DenseDMat<Real> &normals)
{
  normals.resize(5, 3);

  // 0) Face 0-1-4
  normals(0, XI0) = 0.0;
  normals(0, XI1) = -1.0;
  normals(0, XI2) = 0.0;

  // 1) Face 3-0-4
  normals(1, XI0) = -1.0;
  normals(1, XI1) = 0.0;
  normals(1, XI2) = 0.0;

  // 2) Face 1-2-4
  normals(2, XI0) = 1.0;
  normals(2, XI1) = 0.0;
  normals(2, XI2) = 1.0;

  // 3) Face 2-3-4
  normals(3, XI0) = 0.0;
  normals(3, XI1) = 1.0;
  normals(3, XI2) = 1.0;

  // 4) Face 0-3-2-1
  normals(4, XI0) = 0.0;
  normals(4, XI1) = 0.0;
  normals(4, XI2) = -1.0;
}

// ----------------------------------------------------------------------------

void EquidistStdRegionPyramid::fill_permutation(const Uint poly_order,
                                                const EntityRealignCode &permutation_code,
                                                std::vector<Uint> &permutation_vec)
{
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
