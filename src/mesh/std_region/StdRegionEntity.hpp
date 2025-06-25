#ifndef PDEKIT_Mesh_Std_Region_Entity_hpp
#define PDEKIT_Mesh_Std_Region_Entity_hpp

#include <iostream>
#include <vector>

#include "common/ArrayView.hpp"
#include "common/PtrHandle.hpp"
#include "math/DenseDMat.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

class StdRegionEntity
{

  public:
  /// TYPEDEFS

  using wrapped_ptr       = common::PtrHandle<StdRegionEntity>;
  using const_wrapped_ptr = common::PtrHandle<const StdRegionEntity>;
  using coordinates_type  = math::DenseDMat<Real>;

  /// METHODS

  /// Constructor
  StdRegionEntity();

  /// Constructor
  // ReferenceEntity( const Uint size, const Uint idx, const Uint topo_dim,
  // const Uint type_id );

  /// Copy constructor
  StdRegionEntity(const StdRegionEntity &other);

  /// Assignment operator
  StdRegionEntity &operator=(const StdRegionEntity &other);

  /// Destructor
  ~StdRegionEntity();

  /// Set the element type of this entity
  /// This method:
  /// 1) sets the correct interpolation point set type for the entity
  /// 2) sets its topological dimension
  /// 3) resizes its storage for the number of vertices
  /// The method set_vertex() can only be called AFTER the the reference
  /// entity was set to some valid type, because change_std_region() destroys
  /// the local connectivity!
  void construct(const PointSetTag pt_set_id,
                 const common::ArrayView<const Uint, _1D, Uint> &vertices,
                 const common::ArrayView<const Uint, _1D, Uint> &p1_flags, const SUint entity_id);

  void change_point_set_type(const PointSetID pt_set_id);

  /// Resize this entity
  // void resize(const Uint nb_vert);

  /// Indexing operator, const version
  Uint vertex(const Uint i) const;

  /// Return true if this vertex is a p1 (corner) vertex
  /// @param  the local index of the vertex in this element
  bool is_p1_vert(const Uint i) const;

  /// Return the index of this entity in its reference cell
  SUint idx() const;

  /// Size of this entity - number of incident entities
  Uint nb_vert() const;

  /// Return the dimension of this entity
  SUint topo_dim() const;

  /// Return the element type of this entity
  const PointSetTag pt_set_id() const;

  /// Return pointer to the coordinates storage
  const coordinates_type &coordinates() const;

  /*
  /// Set one vertex of this reference entity
  void set_vertex(const Uint vertex_pos, const Uint vertex_value);

  /// Set which of the vertices are p1 vertices
  void set_p1_vert_flag(const Uint p1_vertex_pos, const bool flag);

  /// Set the index of this entity
  void set_id(const SUint my_id);
  */

  /// Set the dimension of this entity
  // void set_dim(const SUint topo_dim);

  private:
  /// Index of this entity (among other reference entities of the same
  /// dimension in the same element. For example, if this entity is a face of
  /// a tetrahedron, then m_idx can be 0,1,2 or 3
  SUint m_idx;

  /// Topological dimension of this entity
  SUint m_dim;

  /// Type of this entity
  PointSetTag m_std_region_id;

  /// Pointer to the coordinates of this reference entity
  coordinates_type m_coordinates;

  /// List of vertices. These indices are the vertex numbers of the reference
  /// element
  std::vector<Uint> m_vert;

  /// A vector of booleans to determine which vertices are p1 vertices
  /// m_p1_vert_flag[i] is true when m_vert[i] is a p1 vertex of this entity
  std::vector<bool> m_p1_vert_flag;
};

// ----------------------------------------------------------------------------

inline Uint StdRegionEntity::vertex(const Uint i) const
{
  return m_vert[i];
}

// ----------------------------------------------------------------------------

inline bool StdRegionEntity::is_p1_vert(const Uint i) const
{
  return m_p1_vert_flag[i];
}

// ----------------------------------------------------------------------------

inline SUint StdRegionEntity::idx() const
{
  return m_idx;
}

// ----------------------------------------------------------------------------

inline Uint StdRegionEntity::nb_vert() const
{
  return m_vert.size();
}

// ----------------------------------------------------------------------------

inline SUint StdRegionEntity::topo_dim() const
{
  return m_dim;
}

// ----------------------------------------------------------------------------

inline const PointSetTag StdRegionEntity::pt_set_id() const
{
  return m_std_region_id;
}

// ----------------------------------------------------------------------------

inline const StdRegionEntity::coordinates_type &StdRegionEntity::coordinates() const
{
  return m_coordinates;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &, const StdRegionEntity &entity);

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
