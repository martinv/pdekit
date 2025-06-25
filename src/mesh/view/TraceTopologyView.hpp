#ifndef PDEKIT_Mesh_View_Trace_Topology_View_hpp
#define PDEKIT_Mesh_View_Trace_Topology_View_hpp

#include "mesh/MeshIndex.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"
#include "mesh/view/ViewTraits.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// TraceTopologyView that does not distinguish between element types
// ----------------------------------------------------------------------------

template <typename Container, typename ConstTrait>
class TraceTopologyView
{
  public:
  static_assert(std::is_same<ConstTrait, ViewIsConst>::value ||
                    std::is_same<ConstTrait, ViewIsNotConst>::value,
                "Use ViewIsConst or ViewIsNotConst as constness trait.");

  static const bool is_immutable_proxy = std::is_same<ConstTrait, ViewIsConst>::value;
  using container_type =
      typename std::conditional<std::is_same<ConstTrait, ViewIsConst>::value,
                                const typename std::remove_const<Container>::type,
                                typename std::remove_const<Container>::type>::type;

  using container_ptr_type = container_type *;
  using container_ref_type = container_type &;

  enum
  {
    GDIM = Container::GDIM
  };

  enum
  {
    TDIM = Container::TDIM
  };

  enum
  {
    FacetDIM = Container::FacetDIM
  };

  TraceTopologyView() = default;

  TraceTopologyView(container_ref_type trace_container, const int position);

  TraceTopologyView(const TraceTopologyView &other_view);

  ~TraceTopologyView();

  TraceTopologyView &operator=(const TraceTopologyView &other_view);

  TraceTopologyView &increment();

  TraceTopologyView &decrement();

  bool equal(const TraceTopologyView &other_view) const;

  const TraceIncidences incidences() const;

  private:
  /// Pointer to the group of cells over which we iterate
  container_ptr_type m_trace_container;

  /// Position of this view
  int m_position;
};

// ----------------------------------------------------------------------------

template <typename Container, typename ConstTrait>
TraceTopologyView<Container, ConstTrait>::TraceTopologyView(container_ref_type trace_container,
                                                            const int position)
    : m_trace_container(&trace_container), m_position(position)
{
}

// ----------------------------------------------------------------------------

template <typename Container, typename ConstTrait>
TraceTopologyView<Container, ConstTrait>::TraceTopologyView(const TraceTopologyView &other_view)
    : m_trace_container(other_view.m_trace_container), m_position(other_view.m_position)
{
}

// ----------------------------------------------------------------------------

template <typename Container, typename ConstTrait>
TraceTopologyView<Container, ConstTrait>::~TraceTopologyView()
{
  m_trace_container = nullptr;
}

// ----------------------------------------------------------------------------

template <typename Container, typename ConstTrait>
TraceTopologyView<Container, ConstTrait> &TraceTopologyView<Container, ConstTrait>::operator=(
    const TraceTopologyView &other_view)
{
  m_trace_container = other_view.m_trace_container;
  m_position        = other_view.m_position;
  return *this;
}
// ----------------------------------------------------------------------------

template <typename Container, typename ConstTrait>
TraceTopologyView<Container, ConstTrait> &TraceTopologyView<Container, ConstTrait>::increment()
{
  m_position++;
  return *this;
}
// ----------------------------------------------------------------------------

template <typename Container, typename ConstTrait>
TraceTopologyView<Container, ConstTrait> &TraceTopologyView<Container, ConstTrait>::decrement()
{
  m_position--;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename Container, typename ConstTrait>
bool TraceTopologyView<Container, ConstTrait>::equal(const TraceTopologyView &other_view) const
{
  return (m_trace_container == other_view.m_trace_container) &&
         (m_position == other_view.m_position);
}

// ----------------------------------------------------------------------------

template <typename Container, typename ConstTrait>
const TraceIncidences TraceTopologyView<Container, ConstTrait>::incidences() const
{
  return m_trace_container->facet_data(FlatIdx(m_position));
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
