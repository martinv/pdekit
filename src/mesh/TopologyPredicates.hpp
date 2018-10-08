#ifndef PDEKIT_Mesh_TopologyPredicates_hpp
#define PDEKIT_Mesh_TopologyPredicates_hpp

#include "mesh/EntityDofRealign.hpp"
#include "mesh/MeshEntity.hpp"

namespace pdekit
{

namespace mesh
{

class TopologyPredicates
{

  public:
  /// Default constructor
  TopologyPredicates();

  /// Destructor
  ~TopologyPredicates();

  /// Return true if given entity contains a vertex with number v
  /// @param entity - entity whose vertices are we checking
  /// @param v      - vertex number which we are looking for in 'entity'
  static bool entity_has_vertex(const MeshEntity &entity, const Uint v);

  /// Return the biggest number of a vertex contained in this entity
  /// @param entity - entity whose vertices are we checking
  /// @return         the biggest vertex number in this entity
  static Uint max_entity_vertex(const MeshEntity &entity);

  /// Return the smallest number of a vertex contained in this entity
  /// @param entity - entity whose vertices are we checking
  /// @return         the smallest vertex number in this entity
  static Uint min_entity_vertex(const MeshEntity &entity);

  /// Return the smallest number of a vertex contained in this entity
  /// @param entity - entity whose vertices are we checking
  /// @return         the smallest vertex number in this entity
  static Uint sum_min_max_entity_dof(const MeshEntity &entity);

  /// Return true if the two entities match (they contain the same nodes).
  /// This method supposes that the entities have the same nodes. Example:
  /// entity_left is a quad   4 75 3 99
  /// entity_right is a quad  4 75 3 99
  /// Then the function shoud consider these two entities as matching and
  /// return 'true' In addition, we allow a situation where the indexing is
  /// 'shifted' - the first index of the right entity is not necessarily the
  /// last index of the left entity. Example: entity_left is a quad   4  6 39
  /// 28 entity_right is a quad 39 28  4  6
  static bool entities_match(const MeshEntity &entity_left, const MeshEntity &entity_right,
                             EntityDofRealign &permutation);

  /// Return true if the two entities match (they contain the same nodes).
  /// This method supposes that the entities have the same nodes, but listed
  /// in reverse order. Example: entity_left is a quad   4 75 3 99
  /// entity_right is a quad  99 3 75 4
  /// Then the function shoud consider these two entities as matching and
  /// return 'true' In addition, we allow a situation where the indexing is
  /// 'shifted' - the first index of the right entity is not necessarily the
  /// last index of the left entity, but the indexes in second entity are
  /// still listed in reverse order. Example: entity_left is a quad  4 6 39 28
  /// entity_right is a quad 39 6 4 28
  static bool entities_match_reverse(const MeshEntity &entity_left, const MeshEntity &entity_right,
                                     EntityDofRealign &permutation);

  private:
};

} // namespace mesh

} // namespace pdekit

#endif
