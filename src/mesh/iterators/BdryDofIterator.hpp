#ifndef PDEKIT_Mesh_Iterators_Bdry_Dof_Iterator_hpp
#define PDEKIT_Mesh_Iterators_Bdry_Dof_Iterator_hpp

#include <iostream>

#include "mesh/MeshIndex.hpp"
#include "mesh/iterators/DofIteratorFilter.hpp"
#include "mesh/iterators/IteratorInterface.hpp"

namespace pdekit
{

namespace mesh
{

// This is a parameter of a friend function of BdryDofIterator (see below),
// hence it must be declared here

// ----------------------------------------------------------------------------

namespace detail
{
template <typename Filter>
struct BdryDofIteratorInternals
{
  template <typename ViewType>
  inline static void dof_iterator_increment_impl(ViewType &view, const Filter &filter)
  {
    view.increment();

    for (; view.idx_less(ActiveIdx(view.nb_active_cells())); view.increment())
    {
      if (filter.filter_pass(view))
        break;
    }
  }

  template <typename ViewType>
  inline static void dof_iterator_decrement_impl(ViewType &view, const Filter &filter)
  {
    view.decrement();

    for (; view.idx_greater_or_equal(ActiveIdx(0)); view.decrement())
    {
      if (filter.filter_pass(view))
        break;
    }
  }
};

// Specialize for defailt filter (which is 'all pass' filter)
// In this case, increment and decrement does not need any loop
// and repeated checking using the filter

template <>
struct BdryDofIteratorInternals<DofIterFilterDefault>
{
  template <typename ViewType>
  inline static void dof_iterator_increment_impl(ViewType &view, const DofIterFilterDefault &filter)
  {
    view.increment();
  }

  template <typename ViewType>
  inline static void dof_iterator_decrement_impl(ViewType &view, const DofIterFilterDefault &filter)
  {
    view.decrement();
  }
};

} // namespace detail

// ----------------------------------------------------------------------------
// BdryDofsIterator loops over boundary of DofMap
// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter = DofIterFilterDefault>
class BdryDofIterator : public IteratorInterface<BdryDofIterator<ViewType, Filter>, ViewType>
{
  public:
  /// TYPEDEFS
  using value_type =
      typename std::conditional<ViewType::is_immutable_proxy, const ViewType, ViewType>::type;
  using reference_type = value_type &;
  using ptr_type       = value_type *;

  private:
  using base_t = IteratorInterface<BdryDofIterator<ViewType, Filter>, ViewType>;

  public:
  /// Default constructor
  BdryDofIterator();

  /// Constructor
  BdryDofIterator(typename ViewType::dof_container_ref_t container,
                  typename ViewType::trace_container_ref_t bdry_container, const ActiveIdx idx);

  /// Constructor
  template <typename... FilterArgs>
  BdryDofIterator(typename ViewType::dof_container_ref_t container,
                  typename ViewType::trace_container_ref_t bdry_container, const ActiveIdx idx,
                  FilterArgs &&... args);

  /// Copy constructor
  BdryDofIterator(const BdryDofIterator &other_iterator);

  /// Assignement operator
  BdryDofIterator &operator=(const BdryDofIterator &other_iterator);

  /// Static function to compute distance of two iterators
  static std::ptrdiff_t distance(const BdryDofIterator &begin, const BdryDofIterator &end);

  static constexpr Uint geo_dim();

  static constexpr Uint topo_dim();

  static constexpr Uint traversal_dim();

  /// Print the position of the iterator
  void print() const;

  private:
  template <typename ViewType1, typename Filter1, typename ViewType2, typename Filter2>
  friend void synchronize_dof_iterators(BdryDofIterator<ViewType1, Filter1> const &it1,
                                        BdryDofIterator<ViewType2, Filter2> &it2);

  /// Step to next cell (prefix increment)
  inline BdryDofIterator &prefix_increment();

  /// Step to the previous cell (prefix decrease)
  inline BdryDofIterator &prefix_decrement();

  /// Dereference operator
  inline reference_type dereference();

  /// Member access operator
  inline ptr_type arrow_operator();

  /// Comparison operator
  inline bool equal(const BdryDofIterator &other_it) const;

  /// Comparison operator
  inline bool not_equal(const BdryDofIterator &other_it) const;

  private:
  /// FRIENDS
  friend class IteratorAccess;
  /// MEMBER VARIABLES

  /// Pointer to the group of cells over which we iterate
  ViewType m_view;

  /// Filter - checks whether m_view is valid (i.e. satisfies filter)
  Filter m_filter;
};

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
BdryDofIterator<ViewType, Filter>::BdryDofIterator() : base_t(), m_view(), m_filter()
{
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
BdryDofIterator<ViewType, Filter>::BdryDofIterator(
    typename ViewType::dof_container_ref_t container,
    typename ViewType::trace_container_ref_t bdry_container, const ActiveIdx idx)
    : base_t(), m_view(container, bdry_container, idx), m_filter()
{
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
template <typename... FilterArgs>
BdryDofIterator<ViewType, Filter>::BdryDofIterator(
    typename ViewType::dof_container_ref_t container,
    typename ViewType::trace_container_ref_t bdry_container, const ActiveIdx idx,
    FilterArgs &&... args)
    : base_t(), m_view(container, bdry_container, idx), m_filter(std::forward<FilterArgs>(args)...)
{
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
BdryDofIterator<ViewType, Filter>::BdryDofIterator(const BdryDofIterator &other_iterator)
    : m_view(other_iterator.m_view), m_filter(other_iterator.m_filter)
{
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
BdryDofIterator<ViewType, Filter> &BdryDofIterator<ViewType, Filter>::operator=(
    const BdryDofIterator &other_iterator)
{
  m_view   = other_iterator.m_view;
  m_filter = other_iterator.m_filter;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
std::ptrdiff_t BdryDofIterator<ViewType, Filter>::distance(const BdryDofIterator &begin,
                                                           const BdryDofIterator &end)
{
  return ViewType::distance(begin.m_view, end.m_view);
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
constexpr Uint BdryDofIterator<ViewType, Filter>::geo_dim()
{
  return ViewType::GDIM;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
constexpr Uint BdryDofIterator<ViewType, Filter>::topo_dim()
{
  return ViewType::TDIM;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
constexpr Uint BdryDofIterator<ViewType, Filter>::traversal_dim()
{
  return ViewType::BcDIM;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
void BdryDofIterator<ViewType, Filter>::print() const
{
  m_view.print();
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
BdryDofIterator<ViewType, Filter> &BdryDofIterator<ViewType, Filter>::prefix_increment()
{
  detail::BdryDofIteratorInternals<Filter>::dof_iterator_increment_impl(m_view, m_filter);
  return *this;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
BdryDofIterator<ViewType, Filter> &BdryDofIterator<ViewType, Filter>::prefix_decrement()
{
  detail::BdryDofIteratorInternals<Filter>::dof_iterator_decrement_impl(m_view, m_filter);
  return *this;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
typename BdryDofIterator<ViewType, Filter>::reference_type BdryDofIterator<ViewType,
                                                                           Filter>::dereference()
{
  return *arrow_operator();
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
typename BdryDofIterator<ViewType, Filter>::ptr_type BdryDofIterator<ViewType,
                                                                     Filter>::arrow_operator()
{
  return &m_view;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
bool BdryDofIterator<ViewType, Filter>::equal(const BdryDofIterator &other_it) const
{
  // Use const_cast to use operator* inside comparison, which automatically
  // updates the entity index corresponding to pos *before* comparison (since
  // update of entity delays until request for entity)

  //       return ((const_cast<AllCellIterator *>(this))->operator*()
  //             == (const_cast<AllCellIterator *>(&other_it))->operator*()
  //             && m_connectivity == other_it.m_connectivity && m_pos ==
  // other_it.m_pos );

  return m_view.equal(other_it.m_view);
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
bool BdryDofIterator<ViewType, Filter>::not_equal(const BdryDofIterator &other_it) const
{
  return !equal(other_it);
}

// ----------------------------------------------------------------------------

template <typename ViewType1, typename Filter1, typename ViewType2, typename Filter2>
inline void synchronize_dof_iterators(BdryDofIterator<ViewType1, Filter1> const &it1,
                                      BdryDofIterator<ViewType2, Filter2> &it2)
{
  it2.m_view.synchronize_with(it1.m_view);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
