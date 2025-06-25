#ifndef PDEKIT_Mesh_Iterators_Cell_Topology_Iterator_hpp
#define PDEKIT_Mesh_Iterators_Cell_Topology_Iterator_hpp

#include <iostream>

#include "mesh/MeshIndex.hpp"
#include "mesh/iterators/CellTopologyIteratorFilter.hpp"
#include "mesh/iterators/IteratorInterface.hpp"
#include "mesh/view/CellTopologyView.hpp"

namespace pdekit
{

namespace mesh
{

// This is a parameter of a friend function of DofIterator (see below), hence
// it must be declared here

// ----------------------------------------------------------------------------

namespace detail
{
template <typename Filter>
struct CellTopologyIteratorInternals
{
  template <typename TopoCellType>
  inline static void cell_topo_iterator_increment_impl(TopoCellType &cell, const Filter &filter)
  {
    cell.increment();

    for (; cell.position() < static_cast<int>(cell.nb_all_cells()); cell.increment())
    {
      if (filter.filter_pass(cell))
        break;
    }
  }

  template <typename TopoCellType>
  inline static void cell_topo_iterator_decrement_impl(TopoCellType &cell, const Filter &filter)
  {
    cell.decrement();

    for (; cell.position() >= 0; cell.decrement())
    {
      if (filter.filter_pass(cell))
        break;
    }
  }

  // Distance returns distance in terms of active indices
  template <typename TopoCellType>
  static std::ptrdiff_t distance(const TopoCellType &tcell_begin, const TopoCellType &tcell_end)
  {
    // tcell_end is 'past end' by 1 increment - we have to decrement it once
    TopoCellType tcell_tmp = tcell_end;
    tcell_tmp.decrement();
    const std::ptrdiff_t dist = tcell_tmp.active_idx().id() - tcell_begin.active_idx().id();
    return dist;
  }
};

// Specialize for defailt filter (which is 'all pass' filter)
// In this case, increment and decrement does not need any loop
// and repeated checking using the filter

template <>
struct CellTopologyIteratorInternals<CellTopologyIterFilterDefault>
{
  template <typename TopoCellType>
  inline static void cell_topo_iterator_increment_impl(TopoCellType &cell,
                                                       const CellTopologyIterFilterDefault &filter)
  {
    cell.increment();
  }

  template <typename TopoCellType>
  inline static void cell_topo_iterator_decrement_impl(TopoCellType &cell,
                                                       const CellTopologyIterFilterDefault &filter)
  {
    cell.decrement();
  }

  // Distance returns distance in terms of active indices
  template <typename TopoCellType>
  static std::ptrdiff_t distance(const TopoCellType &tcell_begin, const TopoCellType &tcell_end)
  {
    // tcell_end is 'past end' by 1 increment - we have to decrement it once
    TopoCellType tcell_tmp = tcell_end;
    tcell_tmp.decrement();
    const std::ptrdiff_t dist = tcell_tmp.active_idx().id() - tcell_begin.active_idx().id();
    return dist;
  }
};

} // namespace detail

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter = CellTopologyIterFilterDefault>
class CellTopologyIterator
    : public IteratorInterface<CellTopologyIterator<TopoCellType, Filter>, TopoCellType>
{
  public:
  /// TYPEDEFS
  using value_type = typename std::conditional<TopoCellType::is_immutable_proxy, const TopoCellType,
                                               TopoCellType>::type;
  using reference_type = value_type &;
  using ptr_type       = value_type *;

  using container_type     = typename TopoCellType::container_type;
  using container_ptr_type = typename TopoCellType::container_ptr_type;
  using container_ref_type = typename TopoCellType::container_ref_type;

  private:
  using base_t = IteratorInterface<CellTopologyIterator<TopoCellType, Filter>, TopoCellType>;

  public:
  /// Default constructor
  CellTopologyIterator();

  /// Constructor
  template <typename... FilterArgs>
  CellTopologyIterator(typename TopoCellType::container_ref_type connectivity, const int position,
                       FilterArgs &&... args);

  /// Copy constructor
  CellTopologyIterator(const CellTopologyIterator &other_iterator);

  /// Assignement operator
  CellTopologyIterator &operator=(const CellTopologyIterator &other_iterator);

  /// Return distance between two iterators
  static std::ptrdiff_t distance(const CellTopologyIterator &begin,
                                 const CellTopologyIterator &end);

  private:
  /// Step to next entity (prefix increment)
  inline CellTopologyIterator &prefix_increment();

  /// Step to the previous entity (prefix decrease)
  inline CellTopologyIterator &prefix_decrement();

  /// Dereference operator
  reference_type dereference();

  /// Member access operator
  ptr_type arrow_operator();

  /// Comparison operator
  bool equal(const CellTopologyIterator &other_it) const;

  /// Comparison operator
  bool not_equal(const CellTopologyIterator &other_it) const;

  static constexpr Uint geo_dim();

  static constexpr Uint topo_dim();

  static constexpr Uint traversal_dim();

  private:
  /// FRIENDS

  friend class IteratorAccess;

  /// MEMBER VARIABLES

  /// Current entity
  TopoCellType m_topo_cell;

  /// Filter - checks whether TopologyCell is valid (i.e. satisfies filter)
  Filter m_filter;
};

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
CellTopologyIterator<TopoCellType, Filter>::CellTopologyIterator() : base_t(), m_topo_cell()
{
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
template <typename... FilterArgs>
CellTopologyIterator<TopoCellType, Filter>::CellTopologyIterator(
    typename TopoCellType::container_ref_type connectivity, const int position,
    FilterArgs &&... args)
    : base_t(), m_topo_cell(&connectivity, FlatIdx(position)),
      m_filter(std::forward<FilterArgs>(args)...)
{
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
CellTopologyIterator<TopoCellType, Filter>::CellTopologyIterator(
    const CellTopologyIterator &other_iterator)
    : m_topo_cell(other_iterator.m_topo_cell), m_filter(other_iterator.m_filter)
{
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
CellTopologyIterator<TopoCellType, Filter> &CellTopologyIterator<TopoCellType, Filter>::operator=(
    const CellTopologyIterator &other_iterator)
{
  m_topo_cell = other_iterator.m_topo_cell;
  m_filter    = other_iterator.m_filter;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
std::ptrdiff_t CellTopologyIterator<TopoCellType, Filter>::distance(
    const CellTopologyIterator &begin, const CellTopologyIterator &end)
{
  return detail::CellTopologyIteratorInternals<Filter>::distance(begin.m_topo_cell,
                                                                 end.m_topo_cell);
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
CellTopologyIterator<TopoCellType, Filter>
    &CellTopologyIterator<TopoCellType, Filter>::prefix_increment()
{
  detail::CellTopologyIteratorInternals<Filter>::cell_topo_iterator_increment_impl(m_topo_cell,
                                                                                   m_filter);
  return *this;
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
CellTopologyIterator<TopoCellType, Filter>
    &CellTopologyIterator<TopoCellType, Filter>::prefix_decrement()
{
  detail::CellTopologyIteratorInternals<Filter>::cell_topo_iterator_increment_impl(m_topo_cell,
                                                                                   m_filter);
  return *this;
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
typename CellTopologyIterator<TopoCellType, Filter>::reference_type CellTopologyIterator<
    TopoCellType, Filter>::dereference()
{
  return *arrow_operator();
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
typename CellTopologyIterator<TopoCellType, Filter>::ptr_type CellTopologyIterator<
    TopoCellType, Filter>::arrow_operator()
{
  return &m_topo_cell;
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
bool CellTopologyIterator<TopoCellType, Filter>::equal(const CellTopologyIterator &other_it) const
{
  // Use const_cast to use operator* inside comparison, which automatically
  // updates the entity index corresponding to pos *before* comparison (since
  // update of entity delays until request for entity)

  //       return ((const_cast<CellIterator *>(this))->operator*()
  //             == (const_cast<CellIterator *>(&other_it))->operator*()
  //             && m_connectivity == other_it.m_connectivity && m_pos ==
  // other_it.m_pos );

  return m_topo_cell.equal(other_it.m_topo_cell);
  /// @todo CHECK THAT FILTERS ARE EQUAL TOO
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
bool CellTopologyIterator<TopoCellType, Filter>::not_equal(
    const CellTopologyIterator &other_it) const
{
  return !equal(other_it);
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
constexpr Uint CellTopologyIterator<TopoCellType, Filter>::geo_dim()
{
  return TopoCellType::GDIM;
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
constexpr Uint CellTopologyIterator<TopoCellType, Filter>::topo_dim()
{
  return TopoCellType::TDIM;
}

// ----------------------------------------------------------------------------

template <typename TopoCellType, typename Filter>
constexpr Uint CellTopologyIterator<TopoCellType, Filter>::traversal_dim()
{
  return TopoCellType::TDIM;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
