#ifndef PDEKIT_Mesh_View_Trace_Dof_View_hpp
#define PDEKIT_Mesh_View_Trace_Dof_View_hpp

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

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
class TraceDofView
{
  public:
  static_assert(std::is_same<ConstTrait, ViewIsConst>::value ||
                    std::is_same<ConstTrait, ViewIsNotConst>::value,
                "Use ViewIsConst or ViewIsNotConst as constness trait.");

  static const bool is_immutable_proxy = std::is_same<ConstTrait, ViewIsConst>::value;
  using dof_container_t =
      typename std::conditional<std::is_same<ConstTrait, ViewIsConst>::value,
                                const typename std::remove_const<DofContainer>::type,
                                typename std::remove_const<DofContainer>::type>::type;

  using trace_container_t =
      typename std::conditional<std::is_same<ConstTrait, ViewIsConst>::value,
                                const typename std::remove_const<TraceTopoContainer>::type,
                                typename std::remove_const<TraceTopoContainer>::type>::type;

  using dof_container_ptr_t = dof_container_t *;
  using dof_container_ref_t = dof_container_t &;

  using trace_container_ptr_t = trace_container_t *;
  using trace_container_ref_t = trace_container_t &;

  enum
  {
    GDIM = TraceTopoContainer::GDIM
  };

  enum
  {
    TDIM = TraceTopoContainer::TDIM
  };

  enum
  {
    FacetDIM = TraceTopoContainer::FacetDIM
  };

  TraceDofView() = default;

  TraceDofView(dof_container_ref_t dof_container, trace_container_ref_t bdry_container,
               const ActiveIdx idx);

  TraceDofView(const TraceDofView &other_view);

  ~TraceDofView();

  TraceDofView &operator=(const TraceDofView &other_view) = default;

  TraceDofView &increment();

  TraceDofView &decrement();

  static std::ptrdiff_t distance(const TraceDofView &first, const TraceDofView &last);

  bool idx_less(const ActiveIdx idx) const;

  bool idx_greater(const ActiveIdx idx) const;

  bool idx_less_or_equal(const ActiveIdx idx) const;

  bool idx_greater_or_equal(const ActiveIdx idx) const;

  bool equal(const TraceDofView &other_view) const;

  void synchronize_with(const TraceDofView &other_view);

  Uint size() const;

  const MeshEntity mesh_entity(const Uint i) const;

  const CellGeometry<GDIM> cell_geometry(const Uint i) const;

  PointSetTag geo_cell_type(const Uint i) const;

  PointSetTag cell_type(const Uint i) const;

  Uint nb_active_cells() const;

  void print() const;

  private:
  /// Reset trace views to default values (empty)
  void reset_trace_views();

  /// Pointer to the dof map which we view
  dof_container_ptr_t m_dof_container;

  /// Pointer to the boundary information (cells and their faces on boundary)
  trace_container_ptr_t m_trace_container;

  /// Position of this view
  ActiveIdx m_active_idx;

  /// Array views into incidence information and entity permutations
  std::tuple<common::ArrayView<const IncidenceEntry, _1D, Uint>,
             common::ArrayView<const EntityDofRealign, _1D, Uint>>
      m_trace_views;
};

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::TraceDofView(
    dof_container_ref_t dof_container, trace_container_ref_t bdry_container, const ActiveIdx idx)
    : m_dof_container(&dof_container), m_trace_container(&bdry_container),
      m_active_idx(idx), m_trace_views{m_trace_container->active_facet_view(idx)}
{
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::TraceDofView(
    const TraceDofView &other_view)
    : m_dof_container(other_view.m_dof_container), m_trace_container(other_view.m_trace_container),
      m_active_idx(other_view.m_active_idx), m_trace_views(other_view.m_trace_views)
{
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::~TraceDofView()
{
  m_dof_container   = nullptr;
  m_trace_container = nullptr;
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>
    &TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::increment()
{
  m_active_idx++;
  if (m_active_idx.id() < m_trace_container->nb_active_facets())
  {
    m_trace_views = m_trace_container->active_facet_view(m_active_idx);
  }
  else
  {
    reset_trace_views();
  }
  return *this;
}
// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>
    &TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::decrement()
{
  m_active_idx--;
  if (m_active_idx.id() >= 0)
  {
    m_trace_views = m_trace_container->active_facet_view(m_active_idx);
  }
  else
  {
    reset_trace_views();
  }
  return *this;
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
inline std::ptrdiff_t TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::distance(
    const TraceDofView &first, const TraceDofView &last)
{
  return last.m_active_idx.id() - first.m_active_idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
inline bool TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::idx_less(
    const ActiveIdx idx) const
{
  return m_active_idx.id() < idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
inline bool TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::idx_greater(
    const ActiveIdx idx) const
{
  return m_active_idx.id() > idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
inline bool TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::idx_less_or_equal(
    const ActiveIdx idx) const
{
  return m_active_idx.id() <= idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
inline bool TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::idx_greater_or_equal(
    const ActiveIdx idx) const
{
  return m_active_idx.id() >= idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
bool TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::equal(
    const TraceDofView &other_view) const
{
  return (m_dof_container == other_view.m_dof_container) &&
         (m_trace_container == other_view.m_trace_container) &&
         (m_active_idx == other_view.m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
void TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::synchronize_with(
    const TraceDofView<DofContainer, TraceTopoContainer, ConstTrait> &other_view)
{
  /// @todo Check that both views point to the same container instance
  m_active_idx  = other_view.m_active_idx;
  m_trace_views = m_trace_container->active_facet_view(m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
Uint TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::size() const
{
  return std::get<0>(m_trace_views).size();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
const MeshEntity TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::mesh_entity(
    const Uint i) const
{
  const IncidenceEntry incidence = std::get<0>(m_trace_views)[i];
  const ActiveIdx cell_active_idx{static_cast<int>(incidence.cell_idx)};
  MeshEntity trace_entity{m_dof_container->active_cell(cell_active_idx)};
  trace_entity.local_transform(FacetDIM, incidence.local_id);
  return trace_entity;
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
const CellGeometry<TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::GDIM> TraceDofView<
    DofContainer, TraceTopoContainer, ConstTrait>::cell_geometry(const Uint i) const
{
  const IncidenceEntry incidence = std::get<0>(m_trace_views)[i];
  const ActiveIdx cell_active_idx{static_cast<int>(incidence.cell_idx)};
  CellTopologyView<typename DofContainer::config_t> topology_cell{
      m_dof_container->tcell(cell_active_idx)};

  return topology_cell.coordinates(FacetDIM, incidence.local_id);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
PointSetTag TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::geo_cell_type(
    const Uint i) const
{
  const IncidenceEntry incidence = std::get<0>(m_trace_views)[i];
  const ActiveIdx cell_active_idx{static_cast<int>(incidence.cell_idx)};
  CellTopologyView<typename DofContainer::config_t> topology_cell{
      m_dof_container->tcell(cell_active_idx)};

  return topology_cell.sub_entity(FacetDIM, incidence.local_id)->pt_set_id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
PointSetTag TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::cell_type(
    const Uint i) const
{
  const mesh::MeshEntity bdry_entity = this->mesh_entity(i);
  return bdry_entity.pt_set_id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
inline Uint TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::nb_active_cells() const
{
  // Needs to be fixed
  return m_trace_container->nb_active_cells();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
void TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::print() const
{
  std::cout << "TraceDofView: curr. pos/nb. entries in container: [" << m_active_idx.id() << "/"
            << m_trace_container->nb_active_facets() << "]" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename TraceTopoContainer, typename ConstTrait>
void TraceDofView<DofContainer, TraceTopoContainer, ConstTrait>::reset_trace_views()
{
  std::get<0>(m_trace_views) = common::ArrayView<const IncidenceEntry, _1D, Uint>();
  std::get<1>(m_trace_views) = common::ArrayView<const EntityDofRealign, _1D, Uint>();
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
