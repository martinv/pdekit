#include "mesh/TopologyPredicates.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

TopologyPredicates::TopologyPredicates()
{
}

// ----------------------------------------------------------------------------

TopologyPredicates::~TopologyPredicates()
{
}

// ----------------------------------------------------------------------------

bool TopologyPredicates::entity_has_vertex(const MeshEntity &entity, const Uint v)
{
  for (Uint i = 0; i < entity.nb_vert(); ++i)
  {
    if (entity.vertex(i) == v)
    {
      return true;
    }
  }

  return false;
}

// ----------------------------------------------------------------------------

Uint TopologyPredicates::max_entity_vertex(const MeshEntity &entity)
{
  Uint max_vert = entity.vertex(0);
  for (Uint i = 1; i < entity.nb_vert(); ++i)
  {
    max_vert = std::max(max_vert, entity.vertex(i));
  }
  return max_vert;
}

// ----------------------------------------------------------------------------

Uint TopologyPredicates::min_entity_vertex(const MeshEntity &entity)
{
  Uint min_vert = entity.vertex(0);
  for (Uint i = 1; i < entity.nb_vert(); ++i)
  {
    min_vert = std::min(min_vert, entity.vertex(i));
  }
  return min_vert;
}

// ----------------------------------------------------------------------------

Uint TopologyPredicates::sum_min_max_entity_dof(const MeshEntity &entity)
{
  Uint min_vert = entity.vertex(0);
  Uint max_vert = entity.vertex(0);
  for (Uint i = 1; i < entity.nb_vert(); ++i)
  {
    min_vert = std::min(min_vert, entity.vertex(i));
    max_vert = std::max(max_vert, entity.vertex(i));
  }
  return min_vert + max_vert;
}

// ----------------------------------------------------------------------------

bool TopologyPredicates::entities_match(const MeshEntity &entity_left,
                                        const MeshEntity &entity_right,
                                        EntityDofRealign &permutation)
{
  const Uint size_L = entity_left.nb_vert();
  const Uint size_R = entity_right.nb_vert();

  // If the two entities that we are comparing have different number of
  // vertices, they cannot possibly match
  if (size_L != size_R)
  {
    return false;
  }

  Uint value_min_L = entity_left.vertex(0);
  Uint value_max_L = entity_right.vertex(0);

  Uint value_min_R = entity_right.vertex(0);
  Uint value_max_R = entity_right.vertex(0);

  bool entities_match = true;

  // Check if the entities match as they are (without any rotation or flip)
  for (Uint i = 0; i < size_L; ++i)
  {
    if (entity_left.vertex(i) != entity_right.vertex(i))
    {
      entities_match = false;
    }
    value_min_L = std::min(value_min_L, entity_left.vertex(i));
    value_max_L = std::max(value_max_L, entity_left.vertex(i));
    value_min_R = std::min(value_min_R, entity_right.vertex(i));
    value_max_R = std::max(value_max_R, entity_right.vertex(i));
  }

  if (entities_match)
  {
    permutation.change_type(entity_right.pt_set_id(),
                            EntityRealignCode::identity(entity_right.pt_set_id().elem_shape()));
    return true;
  }

  // For 1-dimensional entities (edges), we do not try any rotation,
  // because rotation and flip are identical and if edges match after
  // flip, they should be identified by 'entities match reverse',
  // but not here.

  if (entity_left.topo_dim() == _1D)
  {
    return false;
  }

  // If the minimum number of vertex in the left and right vertex is not the
  // same, the entities are not matching
  if ((value_min_L != value_min_R) || (value_max_L != value_max_R))
  {
    return false;
  }

  EntityRealignCode permutation_code =
      EntityRealignCode::single_rotation(entity_right.pt_set_id().elem_shape());
  permutation.change_type(std::make_pair(entity_right.pt_set_id(), permutation_code));

  entities_match = false;
  while ((entity_right.vertex(0) != entity_right.vertex(permutation.get().vertex(0))) &&
         !entities_match)
  {
    entities_match = true;
    for (Uint i = 0; i < size_L; ++i)
    {
      if (entity_left.vertex(i) != entity_right.vertex(permutation.get().vertex(i)))
        entities_match = false;
    }

    if (entities_match)
      return true;

    permutation_code.add_rotation();
    permutation.change_type(entity_right.pt_set_id(), permutation_code);
  }

  return false;
}

// ----------------------------------------------------------------------------

bool TopologyPredicates::entities_match_reverse(const MeshEntity &entity_left,
                                                const MeshEntity &entity_right,
                                                EntityDofRealign &permutation)
{
  const Uint size_L = entity_left.nb_vert();
  const Uint size_R = entity_right.nb_vert();

  // If the two entities that we are comparing have different number of
  // vertices, they cannot possibly match
  if (size_L != size_R)
  {
    return false;
  }

  Uint value_min_L = entity_left.vertex(0);
  Uint value_max_L = entity_right.vertex(0);

  Uint value_min_R = entity_right.vertex(0);
  Uint value_max_R = entity_right.vertex(0);

  bool entities_match = true;

  // Try to see if it is sufficient to flip the right face to have a match
  // between the left and right faces
  EntityRealignCode permutation_code =
      EntityRealignCode::single_flip(entity_right.pt_set_id().elem_shape());
  permutation.change_type(std::make_pair(entity_right.pt_set_id(), permutation_code));

  // Check if the entities match as they are (without any rotation, using just
  // a flip)
  for (Uint i = 0; i < size_L; ++i)
  {
    if (entity_left.vertex(i) != entity_right.vertex(permutation.get().vertex(i)))
      entities_match = false;
    value_min_L = std::min(value_min_L, entity_left.vertex(i));
    value_max_L = std::max(value_max_L, entity_left.vertex(i));
    value_min_R = std::min(value_min_R, entity_right.vertex(i));
    value_max_R = std::max(value_max_R, entity_right.vertex(i));
  }

  if (entities_match)
    return true;

  // For 1-dimensional entities (edges), we do not try any rotation,
  // because rotation and flip are identical and if edges do not match after
  // flip, they should NOT be identified by performing rotations as they
  // are applied below.
  // The reason is that this function could then detect two edges
  // as matching after applying flip AND rotation (hence turning
  // one edge twice and putting it into original position), which
  // is undesired behaviour

  if (entity_left.topo_dim() == _1D)
  {
    return false;
  }

  // If the minimum number of vertex in the left and right vertex is not the
  // same, the entities are not matching
  if ((value_min_L != value_min_R) || (value_max_L != value_max_R))
  {
    return false;
  }

  // Now we will rotate the right entity trying to match it with the left one
  // The reference vertex is the first permuted vertex after a plain flip(no
  // rotations)
  // We will repeat rotations until the first permuted vertex is equal to this
  // reference vertex again
  const Uint ref_vertex_right = entity_right.vertex(permutation.get().vertex(0));
  // Add one rotation
  permutation_code.add_rotation();
  permutation.change_type(entity_right.pt_set_id(), permutation_code);

  entities_match = false;
  while ((ref_vertex_right != entity_right.vertex(permutation.get().vertex(0))) && !entities_match)
  {
    entities_match = true;
    for (Uint i = 0; i < size_L; ++i)
    {
      if (entity_left.vertex(i) != entity_right.vertex(permutation.get().vertex(i)))
        entities_match = false;
    }

    if (entities_match)
      return true;

    permutation_code.add_rotation();
    permutation.change_type(entity_right.pt_set_id(), permutation_code);
  }

  return false;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
