#ifndef PDEKIT_Mesh_View_Cell_Dof_View_hpp
#define PDEKIT_Mesh_View_Cell_Dof_View_hpp

#include "mesh/CellGeometry.hpp"
#include "mesh/MeshEntity.hpp"
#include "mesh/MeshIndex.hpp"
#include "mesh/view/ViewTraits.hpp"

namespace pdekit
{

namespace mesh
{

template <typename MeshConfig>
class CellTopologyView;

// This is a proxy class for accessing information related
// to degrees of freedom stored in DofMap and related to
// one element

template <typename DofContainer, typename ConstTrait>
class CellDofView
{
  public:
  static_assert(std::is_same<ConstTrait, ViewIsConst>::value ||
                    std::is_same<ConstTrait, ViewIsNotConst>::value,
                "Use ViewIsConst or ViewIsNotConst as constness trait.");

  static const bool is_immutable_proxy = std::is_same<ConstTrait, ViewIsConst>::value;
  using container_type =
      typename std::conditional<std::is_same<ConstTrait, ViewIsConst>::value,
                                const typename std::remove_const<DofContainer>::type,
                                typename std::remove_const<DofContainer>::type>::type;

  using container_ptr_type = container_type *;
  using container_ref_type = container_type &;

  enum
  {
    GDIM = DofContainer::GDIM
  };

  enum
  {
    TDIM = DofContainer::TDIM
  };

  CellDofView() = default;

  CellDofView(container_ref_type dof_container, const ActiveIdx active_idx);

  CellDofView(const CellDofView &other_view);

  ~CellDofView();

  CellDofView &operator=(const CellDofView &other_view);

  CellDofView &increment();

  CellDofView &decrement();

  static std::ptrdiff_t distance(const CellDofView &first, const CellDofView &last);

  bool idx_less(const ActiveIdx idx) const;

  bool idx_greater(const ActiveIdx idx) const;

  bool idx_less_or_equal(const ActiveIdx idx) const;

  bool idx_greater_or_equal(const ActiveIdx idx) const;

  bool equal(const CellDofView &other_view) const;

  const mesh::CellTopologyView<typename DofContainer::config_t> tcell() const;

  void synchronize_with(const CellDofView<DofContainer, ConstTrait> &other_view);

  const MeshEntity mesh_entity() const;

  const CellGeometry<GDIM> cell_geometry() const;

  PointSetTag geo_pt_set_id() const;

  PointSetTag pt_set_id() const;

  Uint nb_active_cells() const;

  /// Return the tag of underlying cell
  Uint cell_tag() const;

  void print() const;

  private:
  /// Pointer to the group of cells over which we iterate
  container_ptr_type m_dof_container;

  /// Position of this view
  ActiveIdx m_active_idx;
};

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
CellDofView<DofContainer, ConstTrait>::CellDofView(container_ref_type dof_container,
                                                   const ActiveIdx idx)
    : m_dof_container(&dof_container), m_active_idx(idx)
{
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
CellDofView<DofContainer, ConstTrait>::CellDofView(const CellDofView &other_view)
    : m_dof_container(other_view.m_dof_container), m_active_idx(other_view.m_active_idx)
{
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
CellDofView<DofContainer, ConstTrait>::~CellDofView()
{
  m_dof_container = nullptr;
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
CellDofView<DofContainer, ConstTrait> &CellDofView<DofContainer, ConstTrait>::operator=(
    const CellDofView &other_view)
{
  m_dof_container = other_view.m_dof_container;
  m_active_idx    = other_view.m_active_idx;
  return *this;
}
// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
CellDofView<DofContainer, ConstTrait> &CellDofView<DofContainer, ConstTrait>::increment()
{
  m_active_idx++;
  return *this;
}
// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
CellDofView<DofContainer, ConstTrait> &CellDofView<DofContainer, ConstTrait>::decrement()
{
  m_active_idx--;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline std::ptrdiff_t CellDofView<DofContainer, ConstTrait>::distance(const CellDofView &first,
                                                                      const CellDofView &last)
{
  return last.m_active_idx.id() - first.m_active_idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline bool CellDofView<DofContainer, ConstTrait>::idx_less(const ActiveIdx idx) const
{
  return m_active_idx.id() < idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline bool CellDofView<DofContainer, ConstTrait>::idx_greater(const ActiveIdx idx) const
{
  return m_active_idx.id() > idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline bool CellDofView<DofContainer, ConstTrait>::idx_less_or_equal(const ActiveIdx idx) const
{
  return m_active_idx.id() <= idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline bool CellDofView<DofContainer, ConstTrait>::idx_greater_or_equal(const ActiveIdx idx) const
{
  return m_active_idx.id() >= idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
bool CellDofView<DofContainer, ConstTrait>::equal(const CellDofView &other_view) const
{
  return (m_dof_container == other_view.m_dof_container) &&
         (m_active_idx == other_view.m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
const mesh::CellTopologyView<typename DofContainer::config_t> CellDofView<DofContainer,
                                                                          ConstTrait>::tcell() const
{
  return m_dof_container->tcell(m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
void CellDofView<DofContainer, ConstTrait>::synchronize_with(
    const CellDofView<DofContainer, ConstTrait> &other_view)
{
  /// @todo Check that both views point to the same container instance
  m_active_idx = other_view.m_active_idx;
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
const MeshEntity CellDofView<DofContainer, ConstTrait>::mesh_entity() const
{
  return m_dof_container->active_cell(m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
const CellGeometry<CellDofView<DofContainer, ConstTrait>::GDIM> CellDofView<
    DofContainer, ConstTrait>::cell_geometry() const
{
  const mesh::CellTopologyView<typename DofContainer::config_t> tcell =
      m_dof_container->tcell(m_active_idx);
  return tcell.coordinates();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
PointSetTag CellDofView<DofContainer, ConstTrait>::geo_pt_set_id() const
{
  const mesh::CellTopologyView<typename DofContainer::config_t> tcell =
      m_dof_container->tcell(m_active_idx);
  return tcell.pt_set_id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
PointSetTag CellDofView<DofContainer, ConstTrait>::pt_set_id() const
{
  return m_dof_container->active_cell_std_region_id(m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
Uint CellDofView<DofContainer, ConstTrait>::cell_tag() const
{
  return m_dof_container->active_cell_tag(m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline Uint CellDofView<DofContainer, ConstTrait>::nb_active_cells() const
{
  return m_dof_container->nb_active_cells();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
void CellDofView<DofContainer, ConstTrait>::print() const
{
  std::cout << "DofView: curr. pos/nb. entries in container: [" << m_active_idx.id() << "/"
            << m_dof_container->nb_active_cells() << "]" << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
