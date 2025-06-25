#include "mesh/std_region/EquidistStdRegionLine.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

const std::array<Uint, 4> EquidistStdRegionLine::TopologyStorage = {1, 1, 1, 1};

// ----------------------------------------------------------------------------

Uint EquidistStdRegionLine::nb_dof(const Uint poly_order)
{
  return poly_order + 1;
}

// ----------------------------------------------------------------------------

void EquidistStdRegionLine::fill_topology_relations(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences)
{
  // edge to node incidences - position (1,0)
  common::BlockArray<SUint, SUint> &edge_to_node = incidences[(_3D + 1) * _1D + _0D];
  // Only one entity - the actual edge
  edge_to_node.build({0}, {1});

  // face-to-face incidence should be stored at position (1,1)
  common::BlockArray<SUint, SUint> &edge_to_edge = incidences[(_3D + 1) * _1D + _1D];
  // Edge 0 (the only edge present) is incident only to itself
  // This incidence has to be filled, otherwise 'local_transform'
  // in MeshEntity will fail when the entity is transformed to
  // itself!
  edge_to_edge.build({0}, {1});

  // face-to-node incidence should be stored at position (2,0)
  common::BlockArray<SUint, SUint> &face_to_node = incidences[(_3D + 1) * _2D + _0D];
  // One entity - default (invalid) face
  face_to_node.build({0}, {1});

  // volume-to-node incidence should be stored at position (3,0)
  common::BlockArray<SUint, SUint> &volume_to_node = incidences[(_3D + 1) * _3D + _0D];
  // One entity - default (invalid) volume
  volume_to_node.build({0}, {1});
}

// ----------------------------------------------------------------------------

void EquidistStdRegionLine::fill_edge_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{

  // Example: numbering of a P4 equidistant point set

  //  0--2--3--4--1
  //  Vertices which in the reference entity are stored on positions 0 and 1
  //  are the end vertices. The vertices that follow are the interior vertices
  //  Example: a 'real' line will node coordinates distributed as   < 10 - 43
  //  -
  // 52 - 11 - 19 >,
  //  but the node numbers (connectivity) are internally stored as  [ 10,  19,
  // 43,  52,  11 ]

  // entity_vec.resize(1);

  StdRegionEntity &current_entity = *entity_vec[0];
  // current_entity.resize(nb_nodes);

  const Uint nb_nodes = poly_order + 1;
  std::vector<Uint> vert_ids(nb_nodes);
  std::iota(vert_ids.begin(), vert_ids.end(), 0);

  std::vector<Uint> p1_vert_flags(nb_nodes);
  p1_vert_flags.assign(nb_nodes, 0);
  p1_vert_flags[0] = 1;
  p1_vert_flags[1] = 1;

  const auto vert_ids_view      = common::make_view<Uint, Uint>(vert_ids);
  const auto p1_vert_flags_view = common::make_view<Uint, Uint>(p1_vert_flags);

  current_entity.construct(PointSetTag(ElemShape::Line, poly_order, PointSetID::Equidist),
                           vert_ids_view, p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionLine::fill_face_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // entity_vec.resize(1);
  StdRegionEntity &current_entity = *entity_vec[0];
  // This entity will not be resized - it has length 0 by default

  const common::ArrayView<const Uint, _1D, Uint> vert_ids_view;
  const common::ArrayView<const Uint, _1D, Uint> p1_vert_flags_view;

  current_entity.construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined),
                           vert_ids_view, p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionLine::fill_volume_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  // entity_vec.resize(1);
  StdRegionEntity &current_entity = *entity_vec[0];
  // This entity will not be resized - it has length 0 by default

  const common::ArrayView<const Uint, _1D, Uint> vert_ids_view;
  const common::ArrayView<const Uint, _1D, Uint> p1_vert_flags_view;

  current_entity.construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined),
                           vert_ids_view, p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void EquidistStdRegionLine::fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord)
{
  const Uint nb_nodes = poly_order + 1;
  coord.resize(nb_nodes, TopoDim);

  // First vertex  (index 0)
  coord(0, KSI) = -1.0;
  // Second vertex (index 1)
  coord(1, KSI) = 1.0;

  if (poly_order > P1)
  {
    const Real dx = 2. / poly_order;

    // Node coordinates of the first edge
    Uint inode = 2;
    for (Uint i = 1; i < poly_order; ++i)
    {
      coord(inode, KSI) = -1. + i * dx; // x-coordinate
      inode++;
    }

  } // If polynomial order > 1
}

// ----------------------------------------------------------------------------

void EquidistStdRegionLine::fill_facet_normals(math::DenseDMat<Real> &normals)
{
  normals.resize(2, 1);
  normals(0, XI0) = -1.0;
  normals(1, XI0) = 1.0;
}

// ----------------------------------------------------------------------------

void EquidistStdRegionLine::fill_permutation(const Uint poly_order,
                                             const EntityRealignCode &permutation_code,
                                             std::vector<Uint> &permutation_vec)
{
  // First set permutation_vec as identity permutation: p[i] = i;
  const Uint nb_nodes = poly_order + 1;
  permutation_vec.resize(nb_nodes);
  for (Uint i = 0; i < permutation_vec.size(); ++i)
  {
    permutation_vec[i] = i;
  }

  // If this permutation is identity, we are done
  if (permutation_code.is_identity(ElemShape::Line))
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

void EquidistStdRegionLine::fill_rotation_permutation(const Uint poly_order,
                                                      std::vector<Uint> &permutation_vec)
{
  const Uint nb_nodes = poly_order + 1;
  permutation_vec.resize(nb_nodes);
  permutation_vec[0] = 1;
  permutation_vec[1] = 0;

  for (Uint index = 2; index < nb_nodes; ++index)
  {
    permutation_vec[index] = nb_nodes + 1 - index;
  }
}

// ----------------------------------------------------------------------------

void EquidistStdRegionLine::fill_flip_permutation(const Uint poly_order,
                                                  std::vector<Uint> &permutation_vec)
{
  const Uint nb_nodes = poly_order + 1;
  permutation_vec.resize(nb_nodes);
  permutation_vec[0] = 1;
  permutation_vec[1] = 0;

  for (Uint index = 2; index < nb_nodes; ++index)
  {
    permutation_vec[index] = nb_nodes + 1 - index;
  }
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
