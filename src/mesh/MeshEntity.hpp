#ifndef PDEKIT_Mesh_Mesh_Entity_hpp
#define PDEKIT_Mesh_Mesh_Entity_hpp

#include <iostream>
#include <vector>

#include "mesh/MeshEntity.hpp"
#include "mesh/MeshEntityIterator.hpp"
#include "mesh/std_region/StdRegion.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

class MeshEntity
{
  public:
  /// TYPEDEFS
  // using const_data_ptr = const Uint *;
  /// using data_ptr       = Uint *;
  using const_iterator = mesh::MeshEntityIterator;

  using incidence_list_type = StdRegionInstance::incidence_list_type;
  using coordinates_type    = StdRegionEntity::coordinates_type;

  /// Constructor
  MeshEntity();

  /// Constructor
  MeshEntity(common::ArrayView<const Uint, _1D, Uint> dofs, const Uint idx,
             const PointSetTag type_id);

  /// Constructor
  MeshEntity(common::ArrayView<const Uint, _1D, Uint> dofs, const Uint idx, const StdRegion &type);

  /// Copy constructor
  MeshEntity(const MeshEntity &other);

  /// Assignment operator
  MeshEntity &operator=(const MeshEntity &other);

  /// Destructor
  ~MeshEntity();

  /// Return the index (number) of this entity
  Uint idx() const;

  /// Indexing operator, const version
  Uint vertex(const Uint i) const;

  /// Return true if this vertex is a p1 (corner) vertex
  /// @param  the local index of the vertex in this element
  bool vert_is_p1(const Uint i) const;

  //    /// Indexing operator  which returns the LOCAL number of the vertex
  //    inline const Uint& local_vertex(const Uint i) const
  //    {
  //      return (*m_ref_sub_entity).vertex(i);
  //    }

  /// Size of this entity - number of nodes
  Uint nb_vert() const;

  /// Return the reference coordinates of this entity
  const coordinates_type &ref_coordinates() const;

  /// Return the dimension of this entity
  Uint topo_dim() const;

  /// Return the element type of this entity
  PointSetTag pt_set_id() const;

  /// Return one sub-entity of this entity. For example, if this
  /// entity is a hexahedron, sub_entity(2,0) would return the first
  /// sub-entity of dimension 2, i.e. the first quadrilateral face
  /// @param dim - the topological dimension of the entity for which we are
  ///              asking
  /// @param idx - the index of this sub-entity (counted from 0)
  MeshEntity sub_entity(const Uint dim, const Uint idx) const;

  /// Transform self into another sub-entity of this entity. For example, if
  /// this entity is a hexahedron, local_transform(2,0) would change this
  /// entity into the first sub-entity of dimension 2, i.e. the first
  /// quadrilateral face
  /// @param dim - the topological dimension of the entity for which we are
  ///              asking
  /// @param idx - the index of this sub-entity (counted from 0)
  MeshEntity &local_transform(const Uint dim, const Uint idx);

  /// Return the number of sub-entities of given dimension
  Uint nb_sub_elements(const Uint dim) const;

  /// Return a list of incident entities of dimension dim
  incidence_list_type local_incident_entities(const Uint dim) const;

  /// Set the index of this entity
  void set_idx(const Uint idx);

  /// Reset the data of this entity
  void reinit(common::ArrayView<const Uint, _1D, Uint> dofs, const Uint idx,
              const PointSetTag type_id);

  /// Reset the data of this entity
  void reinit(common::ArrayView<const Uint, _1D, Uint> dofs, const Uint idx, const StdRegion &type);

  const_iterator cbegin() const;

  const_iterator cend() const;

  /// Print the reference topology (i.e. connectivity information
  /// of reference element) for this entity type
  void print_reference_topology() const;

  private:
  /// FRIENDS
  template <Uint GeoDim>
  friend class CellGeometry;

  template <Uint GeoDim>
  friend class DofCoordinates;

  /// DATA

  /// Array view which can access the entity vertex numbers
  common::ArrayView<const Uint, _1D, Uint> m_dofs;

  /// Index (number) of this entity. If this entity represents a cell, for
  /// example, then m_idx is the number of that cell
  Uint m_idx;

  /// Type of this entity (element shape and type of point distribution in
  /// this shape)
  StdRegion m_std_region_type;

  /// Pointer to the active entity of the reference element. For example, a
  /// sub-entity of hexa element can be THE HEXAHEDRAL ELEMENT ITSELF
  /// (default) or one of its quad faces or one of its edges (i.e. a line)
  /// THIS DETERMINES TO WHAT ELEMENT TYPE IS THIS MeshEntity CURRENTLY SET
  std::shared_ptr<StdRegionEntity const> m_active_entity;
};

// ----------------------------------------------------------------------------

inline Uint MeshEntity::idx() const
{
  return m_idx;
}

// ----------------------------------------------------------------------------

inline Uint MeshEntity::vertex(const Uint i) const
{
  return m_dofs[(*m_active_entity).vertex(i)];
}

// ----------------------------------------------------------------------------

inline bool MeshEntity::vert_is_p1(const Uint i) const
{
  return (*m_active_entity).is_p1_vert(i);
}

// ----------------------------------------------------------------------------

inline Uint MeshEntity::nb_vert() const
{
  return (*m_active_entity).nb_vert();
}

// ----------------------------------------------------------------------------

inline MeshEntity::coordinates_type const &MeshEntity::ref_coordinates() const
{
  return (*m_active_entity).coordinates();
}

// ----------------------------------------------------------------------------

inline Uint MeshEntity::topo_dim() const
{
  return (*m_active_entity).topo_dim();
}

// ----------------------------------------------------------------------------

inline PointSetTag MeshEntity::pt_set_id() const
{
  return (*m_active_entity).pt_set_id();
}

// ----------------------------------------------------------------------------

inline MeshEntity &MeshEntity::local_transform(const Uint dim, const Uint idx)
{
  /*
  m_active_entity = m_ref_elem_type.get().elem_entity(dim, idx);
  return *this;
  */

  // Get a list of all reference entities in the master element (i.e. in
  // m_ref_elem_type) which have dimension 'dim' and are incident to the
  // active local entity
  incidence_list_type const ilist(m_std_region_type.get().incident_entities(
      (*m_active_entity).topo_dim(), (*m_active_entity).idx(), dim));

  // Pick the entity with index 'idx' from the list and make it active
  m_active_entity = m_std_region_type.get().elem_entity(dim, ilist[idx]);
  return *this;
}

// ----------------------------------------------------------------------------

inline MeshEntity::const_iterator MeshEntity::cbegin() const
{
  return const_iterator(*this, 0);
}

// ----------------------------------------------------------------------------

inline MeshEntity::const_iterator MeshEntity::cend() const
{
  return const_iterator(*this, nb_vert());
}

// ----------------------------------------------------------------------------

struct MeshEntityVertHasher
{
  inline static Uint entity_xorrer(Uint const &val_l, Uint const &val_r)
  {
    // return val ^ 0x9e3779b9;
    const Uint ret_val = val_l + 0x9e3779b9 + (val_r << 6) + (val_r >> 2);
    return ret_val;
  }

  // ----------------------------------------------------------------------------

  inline static Uint hash(const MeshEntity &entity)
  {
    // return std::accumulate(entity.cbegin(), entity.cend(), 0u, entity_xorrer);

    Uint min_vert = entity.vertex(0);
    Uint max_vert = entity.vertex(0);
    Uint vert_sum = min_vert;

    for (Uint v = 1; v < entity.nb_vert(); ++v)
    {
      vert_sum += entity.vertex(v);
      min_vert = std::min(min_vert, entity.vertex(v));
      max_vert = std::max(max_vert, entity.vertex(v));
    }

    // return min_vert * vert_sum + max_vert;
    return (min_vert * vert_sum) ^ 0x9e377b9 + ((max_vert - min_vert) << 2);
  }
};

// ----------------------------------------------------------------------------

// Output the contents of this mesh entity
std::ostream &operator<<(std::ostream &os, const MeshEntity &entity);

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
