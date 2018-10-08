#ifndef PDEKIT_Mesh_Iterators_Trace_Topology_Iterator_hpp
#define PDEKIT_Mesh_Iterators_Trace_Topology_Iterator_hpp

#include <iostream>

#include "mesh/iterators/IteratorInterface.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// TraceTopologyIterator: can loop over the traces of mesh elements
// ----------------------------------------------------------------------------

template <typename ViewType>
class TraceTopologyIterator : public IteratorInterface<TraceTopologyIterator<ViewType>, ViewType>
{
  public:
  /// TYPEDEFS
  using value_type =
      typename std::conditional<ViewType::is_immutable_proxy, const ViewType, ViewType>::type;
  using reference_type = value_type &;
  using ptr_type       = value_type *;

  private:
  using base_t = IteratorInterface<TraceTopologyIterator<ViewType>, ViewType>;

  public:
  /// Default constructor
  TraceTopologyIterator();

  /// Constructor
  TraceTopologyIterator(typename ViewType::container_ref_type skeleton, const int position = 0);

  /// Copy constructor
  TraceTopologyIterator(const TraceTopologyIterator &other_iterator);

  /// Assignement operator
  TraceTopologyIterator &operator=(const TraceTopologyIterator &other_iterator);

  static constexpr Uint geo_dim();

  static constexpr Uint topo_dim();

  static constexpr Uint traversal_dim();

  private:
  /// Step to next entity (prefix increment)
  inline TraceTopologyIterator &prefix_increment();

  /// Step to the previous entity (prefix decrease)
  inline TraceTopologyIterator &prefix_decrement();

  /// Dereference operator
  reference_type dereference();

  /// Member access operator
  ptr_type arrow_operator();

  /// Comparison operator
  bool equal(const TraceTopologyIterator &other_it) const;

  /// Comparison operator
  bool not_equal(const TraceTopologyIterator &other_it) const;

  private:
  /// FRIENDS

  friend class IteratorAccess;

  /// MEMBER VARIABLES

  /// Pointer to the group of entities over which we iterate
  ViewType m_view;
};

// ----------------------------------------------------------------------------

template <typename ViewType>
TraceTopologyIterator<ViewType>::TraceTopologyIterator() : base_t(), m_view()
{
}

// ----------------------------------------------------------------------------

template <typename ViewType>
TraceTopologyIterator<ViewType>::TraceTopologyIterator(
    typename ViewType::container_ref_type skeleton, const int position)
    : base_t(), m_view(skeleton, position)
{
}

// ----------------------------------------------------------------------------

template <typename ViewType>
TraceTopologyIterator<ViewType>::TraceTopologyIterator(const TraceTopologyIterator &other_iterator)
    : m_view(other_iterator.m_view)
{
}

// ----------------------------------------------------------------------------

template <typename ViewType>
TraceTopologyIterator<ViewType> &TraceTopologyIterator<ViewType>::operator=(
    const TraceTopologyIterator &other_iterator)
{
  m_view = other_iterator.m_view;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename ViewType>
constexpr Uint TraceTopologyIterator<ViewType>::geo_dim()
{
  return ViewType::GDIM;
}

// ----------------------------------------------------------------------------

template <typename ViewType>
constexpr Uint TraceTopologyIterator<ViewType>::topo_dim()
{
  return ViewType::TDIM;
}

// ----------------------------------------------------------------------------

template <typename ViewType>
constexpr Uint TraceTopologyIterator<ViewType>::traversal_dim()
{
  return ViewType::FacetDIM;
}

// ----------------------------------------------------------------------------

template <typename ViewType>
TraceTopologyIterator<ViewType> &TraceTopologyIterator<ViewType>::prefix_increment()
{
  m_view.increment();
  return *this;
}

// ----------------------------------------------------------------------------

template <typename ViewType>
TraceTopologyIterator<ViewType> &TraceTopologyIterator<ViewType>::prefix_decrement()
{
  m_view.decrement();
  return *this;
}

// ----------------------------------------------------------------------------

template <typename ViewType>
typename TraceTopologyIterator<ViewType>::reference_type TraceTopologyIterator<
    ViewType>::dereference()
{
  return *arrow_operator();
}

// ----------------------------------------------------------------------------

template <typename ViewType>
typename TraceTopologyIterator<ViewType>::ptr_type TraceTopologyIterator<ViewType>::arrow_operator()
{
  return &m_view;
  ;
}

// ----------------------------------------------------------------------------

template <typename ViewType>
bool TraceTopologyIterator<ViewType>::equal(const TraceTopologyIterator &other_it) const
{
  // Use const_cast to use operator* inside comparison, which automatically
  // updates the entity index corresponding to pos *before* comparison (since
  // update of entity delays until request for entity)

  //       return ((const_cast<CellIterator *>(this))->operator*()
  //             == (const_cast<CellIterator *>(&other_it))->operator*()
  //             && m_connectivity == other_it.m_connectivity && m_pos ==
  // other_it.m_pos );

  return m_view.equal(other_it.m_view);
}

// ----------------------------------------------------------------------------

template <typename ViewType>
bool TraceTopologyIterator<ViewType>::not_equal(const TraceTopologyIterator &other_it) const
{
  return !equal(other_it);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
