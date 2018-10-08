#ifndef PDEKIT_Mesh_Dof_Coordinates_hpp
#define PDEKIT_Mesh_Dof_Coordinates_hpp

#include "common/PDEKit.hpp"
#include "math/DenseConstVecView.hpp"
#include "mesh/MeshEntity.hpp"

namespace pdekit
{

namespace mesh
{

/// This class should be a lightweight wrapper providing access to all
/// coordinates of
/// one cell

template <Uint GeoDim>
class DofCoordinates
{
  public:
  using node_coords = math::DenseConstVecView<Real>;

  enum
  {
    GDIM = GeoDim
  };

  /// Default constructor
  DofCoordinates();

  /// Copy constructor
  DofCoordinates(const DofCoordinates &other_cell_coords);

  /// Assignement operator
  DofCoordinates &operator=(const DofCoordinates &other_cell_coords);

  /// Constructor from cell and global mesh coordinate array
  DofCoordinates(const MeshEntity &entity, const Real *global_coordinates);

  /// Construct cell coordinates by providing
  /// 1) (all) cell dof ids
  /// 2)  access to ids of active dofs
  /// 3) access to coordinate storage
  DofCoordinates(const common::ArrayView<const Uint, _1D, Uint> &cell_dof_ids,
                 std::shared_ptr<StdRegionEntity const> const &active_entity,
                 const Real *global_coordinates);

  /// Destructor
  ~DofCoordinates();

  /// Return coordinates of i-th node of this cell
  const node_coords c(const Uint i) const;

  /// Return the number of coordinates in this cell
  inline Uint size() const;

  private:
  /// Point to the array of indexes of coordinate nodes
  common::ArrayView<const Uint, _1D, Uint> m_cell_dof_ids;

  /// Pointer to active entity so that the correct
  /// local nodes can be selected
  std::shared_ptr<StdRegionEntity const> m_active_entity;

  /// Pointer to the array containing the actual coordinates
  const Real *m_coords;
};

// ----------------------------------------------------------------------------

template <Uint GeoDim>
DofCoordinates<GeoDim>::DofCoordinates()
    : m_cell_dof_ids(), m_active_entity(nullptr), m_coords(nullptr)
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
DofCoordinates<GeoDim>::DofCoordinates(const DofCoordinates &other_cell_coords)
    : m_cell_dof_ids(other_cell_coords.m_cell_dof_ids),
      m_active_entity(other_cell_coords.m_active_entity), m_coords(other_cell_coords.m_coords)
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
DofCoordinates<GeoDim> &DofCoordinates<GeoDim>::operator=(const DofCoordinates &other_cell_coords)
{
  m_cell_dof_ids  = other_cell_coords.m_cell_dof_ids;
  m_active_entity = other_cell_coords.m_active_entity;
  m_coords        = other_cell_coords.m_coords;
  return *this;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
DofCoordinates<GeoDim>::DofCoordinates(const MeshEntity &entity, const Real *global_coordinates)
    : m_cell_dof_ids(entity.m_dofs), m_active_entity(entity.m_active_entity),
      m_coords(global_coordinates)
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
DofCoordinates<GeoDim>::DofCoordinates(const common::ArrayView<const Uint, _1D, Uint> &cell_dof_ids,
                                       std::shared_ptr<StdRegionEntity const> const &active_entity,
                                       const Real *global_coordinates)
    : m_cell_dof_ids(cell_dof_ids), m_active_entity(active_entity), m_coords(global_coordinates)
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
DofCoordinates<GeoDim>::~DofCoordinates()
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
inline const typename DofCoordinates<GeoDim>::node_coords DofCoordinates<GeoDim>::c(
    const Uint i) const
{
  // The index of the node whose coordinates we want to retrieve is
  // idx = m_first[(*m_active_entity).vertex(i)]
  // return m_coords->node(m_first_index[(*m_active_entity).vertex(i)]);
  const node_coords nc(m_coords + GDIM * m_cell_dof_ids[(*m_active_entity).vertex(i)], GDIM);
  return nc;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
inline Uint DofCoordinates<GeoDim>::size() const
{
  return (*m_active_entity).nb_vert();
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
std::ostream &operator<<(std::ostream &os, const DofCoordinates<GeoDim> &coords)
{
  for (Uint n = 0; n < coords.size(); ++n)
  {
    os << coords.c(n) << std::endl;
  }
  return os;
}

// ----------------------------------------------------------------------------=

} // Namespace mesh

} // Namespace pdekit

#endif
