#include "mesh/MeshEntity.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

MeshEntity::MeshEntity() : m_dofs(), m_idx(0)
{
  m_std_region_type.change_type(StdRegionInstance::undefined);
  // m_active_entity.reset(nullptr);
}

// ----------------------------------------------------------------------------

MeshEntity::MeshEntity(common::ArrayView<const Uint, _1D, Uint> dofs, const Uint idx,
                       const PointSetTag type_id)
    : m_dofs(dofs), m_idx(idx), m_std_region_type(type_id)
{
  m_std_region_type.change_type(type_id);
  m_active_entity = m_std_region_type.get().elem_entity(m_std_region_type.get().topo_dim(), 0);
}

// ----------------------------------------------------------------------------

MeshEntity::MeshEntity(common::ArrayView<const Uint, _1D, Uint> dofs, const Uint idx,
                       const StdRegion &type)
    : m_dofs(dofs), m_idx(idx), m_std_region_type(type)
{
  m_active_entity = m_std_region_type.get().elem_entity(m_std_region_type.get().topo_dim(), 0);
}

// ----------------------------------------------------------------------------

MeshEntity::MeshEntity(const MeshEntity &other)
{
  m_dofs            = other.m_dofs;
  m_idx             = other.m_idx;
  m_std_region_type = other.m_std_region_type;
  m_active_entity   = other.m_active_entity;
}

// ----------------------------------------------------------------------------

MeshEntity &MeshEntity::operator=(const MeshEntity &other)
{
  m_dofs            = other.m_dofs;
  m_idx             = other.m_idx;
  m_std_region_type = other.m_std_region_type;
  m_active_entity   = other.m_active_entity;
  return *this;
}

// ----------------------------------------------------------------------------

MeshEntity::~MeshEntity()
{
}

// ----------------------------------------------------------------------------

MeshEntity MeshEntity::sub_entity(const Uint dim, const Uint idx) const
{
  MeshEntity sub(*this);
  sub.local_transform(dim, idx);
  // sub.m_active_entity = m_ref_elem_type.get().elem_entity(dim, idx);
  return sub;
}

// ----------------------------------------------------------------------------

Uint MeshEntity::nb_sub_elements(const Uint dim) const
{
  // This doesn not work if the current mesh entity is just a subcell
  // of the original entity (i.e. for example m_ref_elem_type is tetra and
  // m_active_entity is triangle, so one of its subfaces). Then the number
  // of sub-elements of dim 1 would return 6 (edges of tetra), instead of
  // 3 (edges of the triangular face)
  // return m_ref_elem_type.get().nb_entities(dim);

  // This is hopefully more correct
  incidence_list_type const ilist(m_std_region_type.get().incident_entities(
      (*m_active_entity).topo_dim(), (*m_active_entity).idx(), dim));
  return ilist.size();
}

// ----------------------------------------------------------------------------

MeshEntity::incidence_list_type MeshEntity::local_incident_entities(const Uint dim) const
{
  incidence_list_type ilist(m_std_region_type.get().incident_entities(
      (*m_active_entity).topo_dim(), (*m_active_entity).idx(), dim));
  return ilist;
}

// ----------------------------------------------------------------------------

void MeshEntity::set_idx(const Uint idx)
{
  m_idx = idx;
}

// ----------------------------------------------------------------------------

void MeshEntity::reinit(common::ArrayView<const Uint, _1D, Uint> dofs, const Uint idx,
                        const PointSetTag type_id)
{
  m_dofs = dofs;
  m_idx  = idx;
  m_std_region_type.change_type(type_id);
  m_active_entity = m_std_region_type.get().elem_entity(m_std_region_type.get().topo_dim(), 0);
}

// ----------------------------------------------------------------------------

void MeshEntity::reinit(common::ArrayView<const Uint, _1D, Uint> dofs, const Uint idx,
                        const StdRegion &type)
{
  m_dofs            = dofs;
  m_idx             = idx;
  m_std_region_type = type;
  m_active_entity   = m_std_region_type.get().elem_entity(m_std_region_type.get().topo_dim(), 0);
}

// ----------------------------------------------------------------------------

void MeshEntity::print_reference_topology() const
{
  m_std_region_type.get().print_complete_topology();
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const MeshEntity &entity)
{
  //   os.clear();
  os << entity.vertex(0);

  for (Uint i = 1; i < entity.nb_vert(); ++i)
  {
    os << " " << entity.vertex(i);
  }
  return os;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
