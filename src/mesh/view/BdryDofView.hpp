#ifndef PDEKIT_Mesh_Containers_Bdry_Dof_View_hpp
#define PDEKIT_Mesh_Containers_Bdry_Dof_View_hpp

#include "mesh/CellGeometry.hpp"
#include "mesh/MeshEntity.hpp"
#include "mesh/MeshIndex.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"
#include "mesh/view/ViewTraits.hpp"

namespace pdekit
{

namespace mesh
{

template <typename MeshConfig>
class CellTopologyView;

// ----------------------------------------------------------------------------
// DofView that does not distinguish between element types
// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
class BdryDofView
{
  public:
  static_assert(std::is_same<ConstTrait, ViewIsConst>::value ||
                    std::is_same<ConstTrait, ViewIsNotConst>::value,
                "Use ViewIsConst or ViewIsNotConst as constness trait.");

  static const bool is_immutable_proxy = std::is_same<ConstTrait, ViewIsConst>::value;
  using dof_container_t =
      typename std::conditional<std::is_same<ConstTrait, ViewIsConst>::value,
                                const typename std::remove_const<Container>::type,
                                typename std::remove_const<Container>::type>::type;

  using trace_container_t =
      typename std::conditional<std::is_same<ConstTrait, ViewIsConst>::value,
                                const typename std::remove_const<BdryContainer>::type,
                                typename std::remove_const<BdryContainer>::type>::type;

  using dof_container_ptr_t = dof_container_t *;
  using dof_container_ref_t = dof_container_t &;

  using trace_container_ptr_t = trace_container_t *;
  using trace_container_ref_t = trace_container_t &;

  enum
  {
    GDIM = BdryContainer::GDIM
  };

  enum
  {
    TDIM = BdryContainer::TDIM
  };

  enum
  {
    BcDIM = BdryContainer::BcDIM
  };

  BdryDofView() = default;

  BdryDofView(dof_container_ref_t dof_container, trace_container_ref_t bdry_container,
              const ActiveIdx idx);

  BdryDofView(const BdryDofView &other_view);

  ~BdryDofView();

  BdryDofView &operator=(const BdryDofView &other_view);

  BdryDofView &increment();

  BdryDofView &decrement();

  static std::ptrdiff_t distance(const BdryDofView &first, const BdryDofView &last);

  bool idx_less(const ActiveIdx idx) const;

  bool idx_greater(const ActiveIdx idx) const;

  bool idx_less_or_equal(const ActiveIdx idx) const;

  bool idx_greater_or_equal(const ActiveIdx idx) const;

  bool equal(const BdryDofView &other_view) const;

  const mesh::CellTopologyView<typename Container::config_t> tcell() const;

  void synchronize_with(const BdryDofView &other_view);

  const MeshEntity mesh_entity() const;

  Uint local_id() const;

  const CellGeometry<GDIM> cell_geometry() const;

  PointSetTag geo_pt_set_id() const;

  PointSetTag pt_set_id() const;

  Uint nb_active_cells() const;

  /// Return the tag of underlying cell
  Uint cell_tag() const;

  void print() const;

  private:
  /// Pointer to the dof map which we view
  dof_container_ptr_t m_dof_container;

  /// Pointer to the boundary information (cells and their faces on boundary)
  trace_container_ptr_t m_bdry_container;

  /// Position of this view
  ActiveIdx m_active_idx;
};

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
BdryDofView<Container, BdryContainer, ConstTrait>::BdryDofView(dof_container_ref_t dof_container,
                                                               trace_container_ref_t bdry_container,
                                                               const ActiveIdx idx)
    : m_dof_container(&dof_container), m_bdry_container(&bdry_container), m_active_idx(idx)
{
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
BdryDofView<Container, BdryContainer, ConstTrait>::BdryDofView(const BdryDofView &other_view)
    : m_dof_container(other_view.m_dof_container), m_bdry_container(other_view.m_bdry_container),
      m_active_idx(other_view.m_active_idx)
{
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
BdryDofView<Container, BdryContainer, ConstTrait>::~BdryDofView()
{
  m_dof_container  = nullptr;
  m_bdry_container = nullptr;
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
BdryDofView<Container, BdryContainer, ConstTrait>
    &BdryDofView<Container, BdryContainer, ConstTrait>::operator=(const BdryDofView &other_view)
{
  m_dof_container  = other_view.m_dof_container;
  m_bdry_container = other_view.m_bdry_container;
  m_active_idx     = other_view.m_active_idx;
  return *this;
}
// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
BdryDofView<Container, BdryContainer, ConstTrait>
    &BdryDofView<Container, BdryContainer, ConstTrait>::increment()
{
  m_active_idx++;
  return *this;
}
// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
BdryDofView<Container, BdryContainer, ConstTrait>
    &BdryDofView<Container, BdryContainer, ConstTrait>::decrement()
{
  m_active_idx--;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
inline std::ptrdiff_t BdryDofView<Container, BdryContainer, ConstTrait>::distance(
    const BdryDofView &first, const BdryDofView &last)
{
  return last.m_active_idx.id() - first.m_active_idx.id();
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
inline bool BdryDofView<Container, BdryContainer, ConstTrait>::idx_less(const ActiveIdx idx) const
{
  return m_active_idx.id() < idx.id();
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
inline bool BdryDofView<Container, BdryContainer, ConstTrait>::idx_greater(
    const ActiveIdx idx) const
{
  return m_active_idx.id() > idx.id();
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
inline bool BdryDofView<Container, BdryContainer, ConstTrait>::idx_less_or_equal(
    const ActiveIdx idx) const
{
  return m_active_idx.id() <= idx.id();
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
inline bool BdryDofView<Container, BdryContainer, ConstTrait>::idx_greater_or_equal(
    const ActiveIdx idx) const
{
  return m_active_idx.id() >= idx.id();
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
bool BdryDofView<Container, BdryContainer, ConstTrait>::equal(const BdryDofView &other_view) const
{
  return (m_dof_container == other_view.m_dof_container) &&
         (m_bdry_container == other_view.m_bdry_container) &&
         (m_active_idx == other_view.m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
const mesh::CellTopologyView<typename Container::config_t> BdryDofView<Container, BdryContainer,
                                                                       ConstTrait>::tcell() const
{
  const FlatIdx adjacent_cell_id = m_bdry_container->parent_cell_id(m_active_idx);
  const CellTopologyView<typename Container::config_t> tcell =
      m_dof_container->tcell(adjacent_cell_id);

  return tcell;
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
void BdryDofView<Container, BdryContainer, ConstTrait>::synchronize_with(
    const BdryDofView<Container, BdryContainer, ConstTrait> &other_view)
{
  /// @todo Check that both views point to the same container instance
  m_active_idx = other_view.m_active_idx;
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
const MeshEntity BdryDofView<Container, BdryContainer, ConstTrait>::mesh_entity() const
{
  return m_bdry_container->active_cell(*m_dof_container, m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
Uint BdryDofView<Container, BdryContainer, ConstTrait>::local_id() const
{
  return m_bdry_container->local_id(m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
const CellGeometry<BdryDofView<Container, BdryContainer, ConstTrait>::GDIM> BdryDofView<
    Container, BdryContainer, ConstTrait>::cell_geometry() const
{
  const FlatIdx adjacent_cell_id = m_bdry_container->parent_cell_id(m_active_idx);
  const Uint local_id            = m_bdry_container->local_id(m_active_idx);
  const CellTopologyView<typename Container::config_t> tcell =
      m_dof_container->tcell(adjacent_cell_id);

  return tcell.coordinates(BcDIM, local_id);
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
PointSetTag BdryDofView<Container, BdryContainer, ConstTrait>::geo_pt_set_id() const
{
  const FlatIdx adjacent_cell_id = m_bdry_container->parent_cell_id(m_active_idx);
  const Uint local_id            = m_bdry_container->local_id(m_active_idx);
  const CellTopologyView<typename Container::config_t> tcell =
      m_dof_container->tcell(adjacent_cell_id);
  return tcell.sub_entity(BcDIM, local_id)->pt_set_id();
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
PointSetTag BdryDofView<Container, BdryContainer, ConstTrait>::pt_set_id() const
{
  const mesh::MeshEntity bdry_entity =
      m_bdry_container->active_cell(*m_dof_container, m_active_idx);
  return bdry_entity.pt_set_id();
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
inline Uint BdryDofView<Container, BdryContainer, ConstTrait>::nb_active_cells() const
{
  return m_bdry_container->nb_active_cells();
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
Uint BdryDofView<Container, BdryContainer, ConstTrait>::cell_tag() const
{
  return m_bdry_container->material_id();
}

// ----------------------------------------------------------------------------

template <typename Container, typename BdryContainer, typename ConstTrait>
void BdryDofView<Container, BdryContainer, ConstTrait>::print() const
{
  std::cout << "BdryDofView: curr. pos/nb. entries in container: [" << m_active_idx.id() << "/"
            << m_bdry_container->nb_active_cells() << "]" << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
