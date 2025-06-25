#include <cassert>

#include "mesh/std_region/EquidistStdRegionHexa.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

/// A hexa has 12 edges, 6 quadrilateral faces and 1 volume
const std::array<Uint, 4> EquidistStdRegionHexa::TopologyStorage = {1, 12, 6, 1};

// ----------------------------------------------------------------------------

/// A HEXAHEDRON OF POLYNOMIAL ORDER P HAS (P+1)^3 NODES
Uint EquidistStdRegionHexa::nb_dof(const Uint poly_order)
{
  return ((poly_order + 1) * (poly_order + 1) * (poly_order + 1));
}

// ----------------------------------------------------------------------------

void EquidistStdRegionHexa::fill_topology_relations(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences)
{
  // edge to node incidences - position (1,0)
  common::BlockArray<SUint, SUint> &edge_to_node = incidences[(_3D + 1) * _1D + _0D];
  // Edge 0-1     Twelve entities - the edges of the hexahedron
  // Edge 0-3     The values 0,1,2, ... on the rhs are POSITIONS
  // Edge 0-4     of reference entities as they will be stored
  // Edge 1-2     in entity vector for entities of dim 1
  // Edge 1-5
  // Edge 2-3
  // Edge 2-6     REMARK: shouldn't this rather be
  // Edge 3-7     edge_to_node(0,0) = 0 // Edge 0 is incident to nodes
  // Edge 4-5     edge_to_node(0,1) = 1 // 0 and 1
  // Edge 4-7         ...
  // Edge 5-6     edge_to_node(6,0) = 2 // Edge 6 is incident to nodes
  // Edge 6-7     edge_to_node(6,1) = 6 // 2 and 6
  edge_to_node.build({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, {12});

  // edge to edge incidences - position (1,1)
  common::BlockArray<SUint, SUint> &edge_to_edge = incidences[(_3D + 1) * _1D + _1D];
  edge_to_edge.build({1, 2, 3, 4,  0, 2, 5, 7,  0, 1, 8,  9,  0, 4,  5, 6,
                      0, 3, 8, 10, 3, 6, 1, 7,  3, 5, 10, 11, 1, 5,  9, 11,
                      2, 9, 4, 10, 2, 8, 7, 11, 4, 8, 6,  11, 6, 10, 7, 9},
                     {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4});

  // edge to face incidences - position (1,2)
  // In a hexahedron, each edge is incident to four faces of the element
  common::BlockArray<SUint, SUint> &edge_to_face = incidences[(_3D + 1) * _1D + _2D];
  edge_to_face.build({0, 1, 2, 3, 0, 1, 2, 4, 0, 1, 2, 5, 0, 1, 3, 4, 0, 1, 3, 5, 0, 2, 3, 4,
                      0, 3, 4, 5, 0, 2, 4, 5, 1, 2, 3, 5, 1, 2, 4, 5, 1, 3, 4, 5, 2, 3, 4, 5},
                     {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4});

  // face-to-node incidence should be stored at position (2,0)
  common::BlockArray<SUint, SUint> &face_to_node = incidences[(_3D + 1) * _2D + _0D];
  // Six quadrilateral faces
  face_to_node.build({0, 1, 2, 3, 4, 5}, {6});

  // face-to-edge incidence - should be stored at position (2,1)
  common::BlockArray<SUint, SUint> &face_to_edge = incidences[(_3D + 1) * _2D + _1D];
  face_to_edge.build({1, 5, 3, 0, 0, 4, 8, 2, 2, 9, 7, 1, 3, 6, 10, 4, 5, 7, 11, 6, 8, 10, 11, 9},
                     {4, 4, 4, 4, 4, 4});

  // face-to-face incidence - position (2,2)
  // Each face is incident to all remaining faces except the opposite one
  common::BlockArray<SUint, SUint> &face_to_face = incidences[(_3D + 1) * _2D + _2D];
  face_to_face.build({1, 2, 3, 4, 0, 2, 3, 5, 0, 1, 4, 5, 0, 1, 4, 5, 0, 2, 3, 5, 1, 2, 3, 4},
                     {4, 4, 4, 4, 4, 4});

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
  volume_to_face.build({0, 1, 2, 3, 4, 5}, {6});

  // volume-to-volume incidence should be stored at position (3,3)
  common::BlockArray<SUint, SUint> &volume_to_volume = incidences[(_3D + 1) * _3D + _3D];
  // Volume 0 (the only volume present) is incident only to itself
  // This incidence has to be filled, otherwise 'local_transform'
  // in MeshEntity will fail when the entity is transformed to
  // itself!
  volume_to_volume.build({0}, {1});
}

// ----------------------------------------------------------------------------

void EquidistStdRegionHexa::fill_edge_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // A hexahedron has 12 edges ...

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

  vert_ids.resize(Edge15::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge15::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge15::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge15::node_id(poly_order, i);
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

  vert_ids.resize(Edge26::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge26::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge26::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge26::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge6.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 6);

  // ---

  StdRegionEntity &edge7 = *entity_vec[7];

  vert_ids.resize(Edge37::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge37::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge37::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge37::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge7.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 7);

  // ---

  StdRegionEntity &edge8 = *entity_vec[8];

  vert_ids.resize(Edge45::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge45::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge45::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge45::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge8.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 8);

  // ---

  StdRegionEntity &edge9 = *entity_vec[9];

  vert_ids.resize(Edge47::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge47::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge47::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge47::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge9.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 9);

  // ---

  StdRegionEntity &edge10 = *entity_vec[10];

  vert_ids.resize(Edge56::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge56::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge56::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge56::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge10.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                   p1_vert_flags_view, 10);

  // ---

  StdRegionEntity &edge11 = *entity_vec[11];

  vert_ids.resize(Edge67::nb_nodes(poly_order));
  p1_vert_flags.resize(Edge67::nb_nodes(poly_order));

  for (Uint i = 0; i < Edge67::nb_nodes(poly_order); ++i)
  {
    vert_ids[i]      = Edge67::node_id(poly_order, i);
    p1_vert_flags[i] = 0;
  }

  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  edge11.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist), vert_ids_view,
                   p1_vert_flags_view, 11);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionHexa::fill_face_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // A hexahedron has 6 faces ...
  // Set correct ids, dimension and element type of all face entities

  assert(entity_vec.size() == 6);

  std::vector<Uint> vert_ids;
  std::vector<Uint> p1_vert_flags;

  /// **********************************
  /// FACE 0-3-2-1
  /// **********************************
  StdRegionEntity &face0 = *entity_vec[0];
  // face0.resize( (poly_order+1)*(poly_order+1) );

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

  Uint entry_pos = 4;
  Uint face_node_id;

  if (poly_order > P1)
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

    face_node_id                   = 7 + 12 * (poly_order - 1) + 1;
    const Uint nb_edge_to_node_pts = poly_order + 1;
    for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 2)); ++i)
    {
      vert_ids[entry_pos++] = face_node_id++;
    }
  }

  auto vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  auto p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face0.construct(PointSetTag(ElemShape::Quad, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 0);

  /// **********************************
  /// FACE 0-1-5-4
  /// **********************************
  StdRegionEntity &face1 = *entity_vec[1];
  // face1.resize( (poly_order+1)*(poly_order+1) );

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids[0] = 0;
  vert_ids[1] = 1;
  vert_ids[2] = 5;
  vert_ids[3] = 4;

  p1_vert_flags.assign(p1_vert_flags.size(), 0);
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;
  p1_vert_flags[3] = 1;

  entry_pos = 4;

  if (poly_order > P1)
  {
    // Fill interior nodes of edge 0-1
    for (Uint i = 2; i < Edge01::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge01::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 1-5
    for (Uint i = 2; i < Edge15::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge15::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 5-4
    for (Uint i = 2; i < Edge45::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge45::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge45::node_id(poly_order, reverse_idx_on_edge);
    }

    // Fill interior nodes of edge 4-0
    for (Uint i = 2; i < Edge04::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge04::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge04::node_id(poly_order, reverse_idx_on_edge);
    }

    //           vertices      edges             1 previous face
    face_node_id = 7 + 12 * (poly_order - 1) + ((poly_order - 1) * (poly_order - 1)) + 1;
    const Uint nb_edge_to_node_pts = poly_order + 1;
    for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 2)); ++i)
    {
      vert_ids[entry_pos++] = face_node_id++;
    }
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face1.construct(PointSetTag(ElemShape::Quad, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 1);

  /// **********************************
  /// FACE 0-4-7-3
  /// **********************************
  StdRegionEntity &face2 = *entity_vec[2];
  // face3.resize( (poly_order+1)*(poly_order+1) );

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids[0] = 0;
  vert_ids[1] = 4;
  vert_ids[2] = 7;
  vert_ids[3] = 3;

  p1_vert_flags.assign(p1_vert_flags.size(), 0);
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;
  p1_vert_flags[3] = 1;

  entry_pos = 4;

  if (poly_order > P1)
  {
    // Fill interior nodes of edge 0-4
    for (Uint i = 2; i < Edge04::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge04::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 4-7
    for (Uint i = 2; i < Edge47::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge47::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 7-3
    for (Uint i = 2; i < Edge37::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge37::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge37::node_id(poly_order, reverse_idx_on_edge);
    }

    // Fill interior nodes of edge 3-0
    for (Uint i = 2; i < Edge03::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge03::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge03::node_id(poly_order, reverse_idx_on_edge);
    }

    //           vertices           edges         2 previous faces
    face_node_id = 7 + 12 * (poly_order - 1) + (2 * (poly_order - 1) * (poly_order - 1)) + 1;
    const Uint nb_edge_to_node_pts = poly_order + 1;
    for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 2)); ++i)
    {
      vert_ids[entry_pos++] = face_node_id++;
    }
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face2.construct(PointSetTag(ElemShape::Quad, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 2);

  /// **********************************
  /// FACE 1-2-6-5
  /// **********************************
  StdRegionEntity &face3 = *entity_vec[3];
  // face2.resize( (poly_order+1)*(poly_order+1) );

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes

  vert_ids[0] = 1;
  vert_ids[1] = 2;
  vert_ids[2] = 6;
  vert_ids[3] = 5;

  p1_vert_flags.assign(p1_vert_flags.size(), 0);
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;
  p1_vert_flags[3] = 1;

  entry_pos = 4;

  if (poly_order > P1)
  {
    // Fill interior nodes of edge 1-2
    for (Uint i = 2; i < Edge12::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge12::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 2-6
    for (Uint i = 2; i < Edge26::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge26::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 6-5
    for (Uint i = 2; i < Edge56::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge56::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge56::node_id(poly_order, reverse_idx_on_edge);
    }

    // Fill interior nodes of edge 5-1
    for (Uint i = 2; i < Edge15::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge15::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge15::node_id(poly_order, reverse_idx_on_edge);
    }

    //           vertices           edges         3 previous faces
    face_node_id = 7 + 12 * (poly_order - 1) + 3 * ((poly_order - 1) * (poly_order - 1)) + 1;
    const Uint nb_edge_to_node_pts = poly_order + 1;
    for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 2)); ++i)
    {
      vert_ids[entry_pos++] = face_node_id++;
    }
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face3.construct(PointSetTag(ElemShape::Quad, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 3);

  /// **********************************
  /// FACE 2-3-7-6
  /// **********************************
  StdRegionEntity &face4 = *entity_vec[4];
  // face4.resize( (poly_order+1)*(poly_order+1) );

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes
  vert_ids[0] = 2;
  vert_ids[1] = 3;
  vert_ids[2] = 7;
  vert_ids[3] = 6;

  p1_vert_flags.assign(p1_vert_flags.size(), 0);
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;
  p1_vert_flags[3] = 1;

  entry_pos = 4;

  if (poly_order > P1)
  {
    // Fill interior nodes of edge 2-3
    for (Uint i = 2; i < Edge23::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge23::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 3-7
    for (Uint i = 2; i < Edge37::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge37::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 7-6
    for (Uint i = 2; i < Edge67::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge67::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge67::node_id(poly_order, reverse_idx_on_edge);
    }

    // Fill interior nodes of edge 6-2
    for (Uint i = 2; i < Edge26::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge26::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge26::node_id(poly_order, reverse_idx_on_edge);
    }

    //           vertices           edges         4 previous faces
    face_node_id = 7 + 12 * (poly_order - 1) + 4 * ((poly_order - 1) * (poly_order - 1)) + 1;
    const Uint nb_edge_to_node_pts = poly_order + 1;
    for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 2)); ++i)
    {
      vert_ids[entry_pos++] = face_node_id++;
    }
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face4.construct(PointSetTag(ElemShape::Quad, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 4);

  /// **********************************
  /// FACE 4-5-6-7
  /// **********************************
  StdRegionEntity &face5 = *entity_vec[5];
  // face5.resize( (poly_order+1)*(poly_order+1) );

  // First have to be filled p1 nodes (corners),
  // then edge nodes, then the remaining internal nodes
  vert_ids[0] = 4;
  vert_ids[1] = 5;
  vert_ids[2] = 6;
  vert_ids[3] = 7;

  p1_vert_flags.assign(p1_vert_flags.size(), 0);
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;
  p1_vert_flags[2] = 1;
  p1_vert_flags[3] = 1;

  entry_pos = 4;

  if (poly_order > P1)
  {
    // Fill interior nodes of edge 4-5
    for (Uint i = 2; i < Edge45::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge45::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 5-6
    for (Uint i = 2; i < Edge56::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge56::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 6-7
    for (Uint i = 2; i < Edge67::nb_nodes(poly_order); ++i)
    {
      vert_ids[entry_pos++] = Edge67::node_id(poly_order, i);
    }

    // Fill interior nodes of edge 7-4
    for (Uint i = 2; i < Edge47::nb_nodes(poly_order); ++i)
    {
      const Uint reverse_idx_on_edge = Edge47::nb_nodes(poly_order) - i + 1;
      vert_ids[entry_pos++]          = Edge47::node_id(poly_order, reverse_idx_on_edge);
    }

    //           vertices           edges         5 previous faces
    face_node_id = 7 + 12 * (poly_order - 1) + 5 * ((poly_order - 1) * (poly_order - 1)) + 1;
    const Uint nb_edge_to_node_pts = poly_order + 1;
    for (Uint i = 0; i < ((nb_edge_to_node_pts - 2) * (nb_edge_to_node_pts - 2)); ++i)
    {
      vert_ids[entry_pos++] = face_node_id++;
    }
  }

  vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  face5.construct(PointSetTag(ElemShape::Quad, poly_order, PointSetID::Equidist), vert_ids_view,
                  p1_vert_flags_view, 5);

} // end of method

// ----------------------------------------------------------------------------

void EquidistStdRegionHexa::fill_volume_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // entity_vec.resize(1);
  StdRegionEntity &volume = *entity_vec[0];

  const Uint nb_nodes = ((poly_order + 1) * (poly_order + 1) * (poly_order + 1));
  std::vector<Uint> vert_ids(nb_nodes);
  std::iota(vert_ids.begin(), vert_ids.end(), 0);

  std::vector<Uint> p1_vert_flags(nb_nodes);
  p1_vert_flags.assign(nb_nodes, 0);
  for (Uint n = 0; n < 8; ++n)
  {
    p1_vert_flags[n] = 1;
  }

  const auto vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  const auto p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  volume.construct(PointSetTag(ElemShape::Hexa, poly_order, PointSetID::Equidist), vert_ids_view,
                   p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionHexa::fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord)
{
  const Uint nb_nodes = ((poly_order + 1) * (poly_order + 1) * (poly_order + 1));
  coord.resize(nb_nodes, _3D);

  const Real dx = 2. / poly_order;
  const Real dy = 2. / poly_order;
  const Real dz = 2. / poly_order;

  // Coordinates of local vertex 0
  Real ksi_0 = -1.;
  Real eta_0 = -1.;
  Real zta_0 = -1.;

  // Coordinates of local vertex 1
  Real ksi_1 = 1.;
  Real eta_1 = -1.;
  Real zta_1 = -1.;

  // Coordinates of local vertex 2
  Real ksi_2 = 1.;
  Real eta_2 = 1.;
  Real zta_2 = -1.;

  // Coordinates of local vertex 3
  Real ksi_3 = -1.;
  Real eta_3 = 1.;
  Real zta_3 = -1.;

  // Coordinates of local vertex 4
  Real ksi_4 = -1.;
  Real eta_4 = -1.;
  Real zta_4 = 1.;

  // Coordinates of local vertex 5
  Real ksi_5 = 1.;
  Real eta_5 = -1.;
  Real zta_5 = 1.;

  // Coordinates of local vertex 6
  Real ksi_6 = 1.;
  Real eta_6 = 1.;
  Real zta_6 = 1.;

  // Coordinates of local vertex 7
  Real ksi_7 = -1.;
  Real eta_7 = 1.;
  Real zta_7 = 1.;

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

      // Sixth corner
      coord(inode, KSI) = ksi_5;
      coord(inode, ETA) = eta_5;
      coord(inode, ZTA) = zta_5;
      inode++;

      // Seventh corner
      coord(inode, KSI) = ksi_6;
      coord(inode, ETA) = eta_6;
      coord(inode, ZTA) = zta_6;
      inode++;

      // Eighth corner
      coord(inode, KSI) = ksi_7;
      coord(inode, ETA) = eta_7;
      coord(inode, ZTA) = zta_7;
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

      // Edge 1-5
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_1;
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

      // Edge 2-6
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_2;
        coord(inode, ETA) = eta_2;
        coord(inode, ZTA) = zta_2 + (i + 1) * dz;
        inode++;
      }

      // Edge 3-7
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_3;
        coord(inode, ETA) = eta_3;
        coord(inode, ZTA) = zta_3 + (i + 1) * dz;
        inode++;
      }

      // Edge 4-5
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_4 + (i + 1) * dx;
        coord(inode, ETA) = eta_4;
        coord(inode, ZTA) = zta_4;
        inode++;
      }

      // Edge 4-7
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_4;
        coord(inode, ETA) = eta_4 + (i + 1) * dy;
        coord(inode, ZTA) = zta_4;
        inode++;
      }

      // Edge 5-6
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_5;
        coord(inode, ETA) = eta_5 + (i + 1) * dy;
        coord(inode, ZTA) = zta_5;
        inode++;
      }

      // Edge 6-7
      for (Uint i = 0; i < nb_interior_pts_on_edge; ++i)
      {
        coord(inode, KSI) = ksi_6 - (i + 1) * dx;
        coord(inode, ETA) = eta_6;
        coord(inode, ZTA) = zta_6;
        inode++;
      }

    } // If HEXA polynomial order > 1

    if (current_poly_order >= P2)
    {
      Real face_ksi_0, face_eta_0, face_zta_0, face_ksi_1, face_eta_1, face_zta_1, face_ksi_2,
          face_eta_2, face_zta_2, face_ksi_3, face_eta_3, face_zta_3;

      Int face_poly_order = current_poly_order - 2;

      // --------------------------------
      // Node coordinates of face 0-3-2-1
      // --------------------------------
      face_ksi_0 = ksi_0 + dx;
      face_eta_0 = eta_0 + dy;
      face_zta_0 = zta_0;

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

      // --------------------------------
      // Node coordinates of face 0-1-5-4
      // --------------------------------
      face_poly_order = current_poly_order - 2;

      face_ksi_0 = ksi_0 + dx;
      face_eta_0 = eta_0;
      face_zta_0 = zta_0 + dz;

      face_ksi_1 = ksi_1 - dx;
      face_eta_1 = eta_1;
      face_zta_1 = zta_1 + dz;

      face_ksi_2 = ksi_5 - dx;
      face_eta_2 = eta_5;
      face_zta_2 = zta_5 - dz;

      face_ksi_3 = ksi_4 + dx;
      face_eta_3 = eta_4;
      face_zta_3 = zta_4 - dz;

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
          // 0-1-5-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_0 + i * dx;
            coord(inode, ETA) = face_eta_0;
            coord(inode, ZTA) = face_zta_0;
            inode++;
          }

          // Edge 1-2 of the local quadrilateral embedded in quad face
          // 0-1-5-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_1;
            coord(inode, ETA) = face_eta_1;
            coord(inode, ZTA) = face_zta_1 + i * dz;
            inode++;
          }

          // Edge 2-3 of the local quadrilateral embedded in quad face
          // 0-1-5-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_2 - i * dx;
            coord(inode, ETA) = face_eta_2;
            coord(inode, ZTA) = face_zta_2;
            inode++;
          }

          // Edge 3-0 of the local quadrilateral embedded in quad face
          // 0-1-5-4
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_3;
            coord(inode, ETA) = face_eta_3;
            coord(inode, ZTA) = face_zta_3 - i * dz;
            inode++;
          }
        }

        face_ksi_0 += dx;
        // face_eta_0 remains unchanged
        face_zta_0 += dz;

        face_ksi_1 -= dx;
        // face_eta_1 remains unchanged
        face_zta_1 += dz;

        face_ksi_2 -= dx;
        // face_eta_2 remains unchanged
        face_zta_2 -= dz;

        face_ksi_3 += dx;
        // face_eta_3 remains unchanged
        face_zta_3 -= dz;

        face_poly_order -= 2;
      }

      // --------------------------------
      // Node coordinates of face 0-4-7-3
      // --------------------------------
      face_poly_order = current_poly_order - 2;

      face_ksi_0 = ksi_0;
      face_eta_0 = eta_0 + dy;
      face_zta_0 = zta_0 + dz;

      face_ksi_1 = ksi_4;
      face_eta_1 = eta_4 + dy;
      face_zta_1 = zta_4 - dz;

      face_ksi_2 = ksi_7;
      face_eta_2 = eta_7 - dy;
      face_zta_2 = zta_7 - dz;

      face_ksi_3 = ksi_3;
      face_eta_3 = eta_3 - dy;
      face_zta_3 = zta_3 + dz;

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
          // 0-4-7-3
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_0;
            coord(inode, ETA) = face_eta_0;
            coord(inode, ZTA) = face_zta_0 + i * dz;
            inode++;
          }

          // Edge 1-2 of the local quadrilateral embedded in quad face
          // 0-4-7-3
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_1;
            coord(inode, ETA) = face_eta_1 + i * dy;
            coord(inode, ZTA) = face_zta_1;
            inode++;
          }

          // Edge 2-3 of the local quadrilateral embedded in quad face
          // 0-4-7-3
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_2;
            coord(inode, ETA) = face_eta_2;
            coord(inode, ZTA) = face_zta_2 - i * dz;
            inode++;
          }

          // Edge 3-0 of the local quadrilateral embedded in quad face
          // 0-4-7-3
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_3;
            coord(inode, ETA) = face_eta_3 - i * dy;
            coord(inode, ZTA) = face_zta_3;
            inode++;
          }
        }

        // face_ksi_0 remains unchanged
        face_eta_0 += dy;
        face_zta_0 += dz;

        // face_ksi_1 remains unchanged
        face_eta_1 += dy;
        face_zta_1 -= dz;

        // face_ksi_2 remains unchanged
        face_eta_2 -= dy;
        face_zta_2 -= dz;

        // face_ksi_3 remains unchanged
        face_eta_3 -= dy;
        face_zta_3 += dz;

        face_poly_order -= 2;
      }

      // --------------------------------
      // Node coordinates of face 1-2-6-5
      // --------------------------------
      face_poly_order = current_poly_order - 2;

      face_ksi_0 = ksi_1;
      face_eta_0 = eta_1 + dy;
      face_zta_0 = zta_1 + dz;

      face_ksi_1 = ksi_2;
      face_eta_1 = eta_2 - dy;
      face_zta_1 = zta_2 + dz;

      face_ksi_2 = ksi_6;
      face_eta_2 = eta_6 - dy;
      face_zta_2 = zta_6 - dz;

      face_ksi_3 = ksi_5;
      face_eta_3 = eta_5 + dy;
      face_zta_3 = zta_5 - dz;

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
          // 1-2-6-5
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_0;
            coord(inode, ETA) = face_eta_0 + i * dy;
            coord(inode, ZTA) = face_zta_0;
            inode++;
          }

          // Edge 1-2 of the local quadrilateral embedded in quad face
          // 1-2-6-5
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_1;
            coord(inode, ETA) = face_eta_1;
            coord(inode, ZTA) = face_zta_1 + i * dz;
            inode++;
          }

          // Edge 2-3 of the local quadrilateral embedded in quad face
          // 1-2-6-5
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_2;
            coord(inode, ETA) = face_eta_2 - i * dy;
            coord(inode, ZTA) = face_zta_2;
            inode++;
          }

          // Edge 3-0 of the local quadrilateral embedded in quad face
          // 1-2-6-5
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_3;
            coord(inode, ETA) = face_eta_3;
            coord(inode, ZTA) = face_zta_3 - i * dz;
            inode++;
          }
        }

        // face_ksi_0 remains unchanged
        face_eta_0 += dy;
        face_zta_0 += dz;

        // face_ksi_1 remains unchanged
        face_eta_1 -= dy;
        face_zta_1 += dz;

        // face_ksi_2 remains unchanged
        face_eta_2 -= dy;
        face_zta_2 -= dz;

        // face_ksi_3 remains unchanged
        face_eta_3 += dy;
        face_zta_3 -= dz;

        face_poly_order -= 2;
      }

      // --------------------------------
      // Node coordinates of face 2-3-7-6
      // --------------------------------
      face_poly_order = current_poly_order - 2;

      face_ksi_0 = ksi_2 - dx;
      face_eta_0 = eta_2;
      face_zta_0 = zta_2 + dz;

      face_ksi_1 = ksi_3 + dx;
      face_eta_1 = eta_3;
      face_zta_1 = zta_3 + dz;

      face_ksi_2 = ksi_7 + dx;
      face_eta_2 = eta_7;
      face_zta_2 = zta_7 - dz;

      face_ksi_3 = ksi_6 - dx;
      face_eta_3 = eta_6;
      face_zta_3 = zta_6 - dz;

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
          // 2-3-7-6
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_0 - i * dx;
            coord(inode, ETA) = face_eta_0;
            coord(inode, ZTA) = face_zta_0;
            inode++;
          }

          // Edge 1-2 of the local quadrilateral embedded in quad face
          // 2-3-7-6
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_1;
            coord(inode, ETA) = face_eta_1;
            coord(inode, ZTA) = face_zta_1 + i * dz;
            inode++;
          }

          // Edge 2-3 of the local quadrilateral embedded in quad face
          // 2-3-7-6
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_2 + i * dx;
            coord(inode, ETA) = face_eta_2;
            coord(inode, ZTA) = face_zta_2;
            inode++;
          }

          // Edge 3-0 of the local quadrilateral embedded in quad face
          // 2-3-7-6
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_3;
            coord(inode, ETA) = face_eta_3;
            coord(inode, ZTA) = face_zta_3 - i * dz;
            inode++;
          }
        }

        face_ksi_0 -= dx;
        // face_eta_0 remains unchanged
        face_zta_0 += dz;

        face_ksi_1 += dx;
        // face_eta_1 remains unchanged
        face_zta_1 += dz;

        face_ksi_2 += dx;
        // face_eta_2 remains unchanged
        face_zta_2 -= dz;

        face_ksi_3 -= dx;
        // face_eta_3 remains unchanged
        face_zta_3 -= dz;

        face_poly_order -= 2;
      }

      // --------------------------------
      // Node coordinates of face 4-5-6-7
      // --------------------------------
      face_poly_order = current_poly_order - 2;

      face_ksi_0 = ksi_4 + dx;
      face_eta_0 = eta_4 + dy;
      face_zta_0 = zta_4;

      face_ksi_1 = ksi_5 - dx;
      face_eta_1 = eta_5 + dy;
      face_zta_1 = zta_5;

      face_ksi_2 = ksi_6 - dx;
      face_eta_2 = eta_6 - dy;
      face_zta_2 = zta_6;

      face_ksi_3 = ksi_7 + dx;
      face_eta_3 = eta_7 - dy;
      face_zta_3 = zta_7;

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
          // 4-5-6-7
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_0 + i * dx;
            coord(inode, ETA) = face_eta_0;
            coord(inode, ZTA) = face_zta_0;
            inode++;
          }

          // Edge 1-2 of the local quadrilateral embedded in quad face
          // 4-5-6-7
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_1;
            coord(inode, ETA) = face_eta_1 + i * dy;
            coord(inode, ZTA) = face_zta_1;
            inode++;
          }

          // Edge 2-3 of the local quadrilateral embedded in quad face
          // 4-5-6-7
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_2 - i * dx;
            coord(inode, ETA) = face_eta_2;
            coord(inode, ZTA) = face_zta_2;
            inode++;
          }

          // Edge 3-0 of the local quadrilateral embedded in quad face
          // 4-5-6-7
          for (Int i = 1; i < face_poly_order; ++i)
          {
            coord(inode, KSI) = face_ksi_3;
            coord(inode, ETA) = face_eta_3 - i * dy;
            coord(inode, ZTA) = face_zta_3;
            inode++;
          }
        }

        face_ksi_0 += dx;
        face_eta_0 += dy;
        // face_zta_0 remains unchanged

        face_ksi_1 -= dx;
        face_eta_1 += dy;
        // face_zta_1 remains unchanged

        face_ksi_2 -= dx;
        face_eta_2 -= dy;
        // face_zta_2 remains unchanged

        face_ksi_3 += dx;
        face_zta_3 -= dy;
        // face_zta_3 remains unchanged

        face_poly_order -= 2;
      }

    } // If HEXA polynomial order > 2

    ksi_0 += dx;
    eta_0 += dy;
    zta_0 += dz;

    ksi_1 -= dx;
    eta_1 += dy;
    zta_1 += dz;

    ksi_2 -= dx;
    eta_2 -= dy;
    zta_2 += dz;

    ksi_3 += dx;
    eta_3 -= dy;
    zta_3 += dz;

    ksi_4 += dx;
    eta_4 += dy;
    zta_4 -= dz;

    ksi_5 -= dx;
    eta_5 += dy;
    zta_5 -= dz;

    ksi_6 -= dx;
    eta_6 -= dy;
    zta_6 -= dz;

    ksi_7 += dx;
    eta_7 -= dy;
    zta_7 -= dz;

    current_poly_order -= 2;

  } // While HEXA polyorder >= 0
}

// ----------------------------------------------------------------------------

void EquidistStdRegionHexa::fill_facet_normals(math::DenseDMat<Real> &normals)
{
  normals.resize(6, 3);

  // 0) 0-3-2-1
  normals(0, XI0) = 0.0;
  normals(0, XI1) = 0.0;
  normals(0, XI2) = -1.0;

  // 1) 0-1-5-4
  normals(1, XI0) = 0.0;
  normals(1, XI1) = -1.0;
  normals(1, XI2) = 0.0;

  // 2) 0-4-7-3
  normals(2, XI0) = -1.0;
  normals(2, XI1) = 0.0;
  normals(2, XI2) = 0.0;

  // 3) 1-2-6-5
  normals(3, XI0) = 1.0;
  normals(3, XI1) = 0.0;
  normals(3, XI2) = 0.0;

  // 4) 2-3-7-6
  normals(4, XI0) = 0.0;
  normals(4, XI1) = 1.0;
  normals(4, XI2) = 0.0;

  // 5) 4-5-6-7
  normals(5, XI0) = 0.0;
  normals(5, XI1) = 0.0;
  normals(5, XI2) = 1.0;
}

// ----------------------------------------------------------------------------

void EquidistStdRegionHexa::fill_permutation(const Uint poly_order,
                                             const EntityRealignCode &permutation_code,
                                             std::vector<Uint> &permutation_vec)
{
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
