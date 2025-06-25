#ifndef PDEKIT_Mesh_Cell_Geometry_hpp
#define PDEKIT_Mesh_Cell_Geometry_hpp

#include "common/PDEKit.hpp"
#include "math/DenseVecView.hpp"
#include "mesh/MeshEntity.hpp"

namespace pdekit
{

namespace mesh
{

/// This class should be a lightweight wrapper providing access to all
/// coordinates of
/// one cell

template <Uint GeoDim>
class CellGeometry
{
  public:
  using node_view_t = math::DenseVecView<const Real>;

  enum
  {
    GDIM = GeoDim
  };

  /// Default constructor
  CellGeometry();

  /// Copy constructor
  CellGeometry(const CellGeometry &other_cell_geom);

  /// Assignement operator
  CellGeometry &operator=(const CellGeometry &other_cell_geom);

  /// Constructor from cell and global mesh coordinate array
  CellGeometry(const math::DenseVecView<const Real> &cell_node_coords, const MeshEntity &entity);

  /// Construct cell coordinates by providing
  /// 1) (all) cell dof ids
  /// 2)  access to ids of active dofs
  /// 3) access to coordinate storage
  CellGeometry(const math::DenseVecView<const Real> &cell_node_coords,
               std::shared_ptr<StdRegionEntity const> const &active_entity);

  /// Destructor
  ~CellGeometry();

  /// Return coordinates of i-th node of this cell
  const node_view_t const_node_view(const Uint i) const;

  /// Return coordinates of i-th node of this cell
  // node_view_t node_view(const Uint i);

  constexpr Uint dim() const;

  /// Return the number of coordinates in this cell
  inline Uint size() const;

  private:
  /// Point to the array of indexes of coordinate nodes
  math::DenseVecView<const Real> m_cell_node_coords;

  /// Pointer to active entity so that the correct
  /// local nodes can be selected
  std::shared_ptr<StdRegionEntity const> m_active_entity;
};

// ----------------------------------------------------------------------------

template <Uint GeoDim>
CellGeometry<GeoDim>::CellGeometry() : m_cell_node_coords(), m_active_entity(nullptr)
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
CellGeometry<GeoDim>::CellGeometry(const CellGeometry &other_cell_geom)
    : m_cell_node_coords(other_cell_geom.m_cell_node_coords),
      m_active_entity(other_cell_geom.m_active_entity)
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
CellGeometry<GeoDim> &CellGeometry<GeoDim>::operator=(const CellGeometry &other_cell_geom)
{
  m_cell_node_coords = other_cell_geom.m_cell_node_coords;
  m_active_entity    = other_cell_geom.m_active_entity;
  return *this;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
CellGeometry<GeoDim>::CellGeometry(const math::DenseVecView<const Real> &cell_node_coords,
                                   const MeshEntity &entity)
    : m_cell_node_coords(cell_node_coords), m_active_entity(entity.m_active_entity)
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
CellGeometry<GeoDim>::CellGeometry(const math::DenseVecView<const Real> &cell_node_coords,
                                   std::shared_ptr<StdRegionEntity const> const &active_entity)
    : m_cell_node_coords(cell_node_coords), m_active_entity(active_entity)
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
CellGeometry<GeoDim>::~CellGeometry()
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
inline const typename CellGeometry<GeoDim>::node_view_t CellGeometry<GeoDim>::const_node_view(
    const Uint i) const
{
  const node_view_t nv = m_cell_node_coords.slice(GDIM * (*m_active_entity).vertex(i), GDIM);
  return nv;
}

// ----------------------------------------------------------------------------

/*
template <Uint GeoDim>
inline typename CellGeometry<GeoDim>::node_view_t
CellGeometry<GeoDim>::node_view(const Uint i)
{
  const node_view_t nv(&m_cell_node_coords[(*m_active_entity).vertex(i)], GDIM,
1); return nv;
}
*/

// ----------------------------------------------------------------------------

template <Uint GeoDim>
constexpr Uint CellGeometry<GeoDim>::dim() const
{
  return GeoDim;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
inline Uint CellGeometry<GeoDim>::size() const
{
  return (*m_active_entity).nb_vert();
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
std::ostream &operator<<(std::ostream &os, const CellGeometry<GeoDim> &coords)
{
  for (Uint n = 0; n < coords.size(); ++n)
  {
    os << coords.const_node_view(n) << std::endl;
  }
  return os;
}

// ----------------------------------------------------------------------------=

} // Namespace mesh

} // Namespace pdekit

#endif
