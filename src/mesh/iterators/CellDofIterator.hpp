#ifndef PDEKIT_Mesh_Iterators_Cell_Dof_Iterator_hpp
#define PDEKIT_Mesh_Iterators_Cell_Dof_Iterator_hpp

#include <cstddef>
#include <iostream>

#include "mesh/MeshIndex.hpp"
#include "mesh/iterators/DofIteratorFilter.hpp"
#include "mesh/iterators/IteratorInterface.hpp"

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
struct CellDofIteratorInternals
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
struct CellDofIteratorInternals<DofIterFilterDefault>
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

template <typename ViewType, typename Filter = DofIterFilterDefault>
class CellDofIterator : public IteratorInterface<CellDofIterator<ViewType, Filter>, ViewType>
{
  public:
  /// TYPEDEFS
  using value_type =
      typename std::conditional<ViewType::is_immutable_proxy, const ViewType, ViewType>::type;
  using reference_type = value_type &;
  using ptr_type       = value_type *;

  private:
  using base_t = IteratorInterface<CellDofIterator<ViewType, Filter>, ViewType>;

  public:
  /// Default constructor
  CellDofIterator();

  /// Constructor
  template <typename... FilterArgs>
  CellDofIterator(typename ViewType::container_ref_type dof_handler, const ActiveIdx idx,
                  FilterArgs &&... args);

  /// Copy constructor
  CellDofIterator(const CellDofIterator &other_iterator);

  /// Assignement operator
  CellDofIterator &operator=(const CellDofIterator &other_iterator);

  /// Static function to compute distance of two iterators
  static std::ptrdiff_t distance(const CellDofIterator &begin, const CellDofIterator &end);

  static constexpr Uint geo_dim();

  static constexpr Uint topo_dim();

  static constexpr Uint traversal_dim();

  /// Print the position of the iterator
  void print() const;

  private:
  /*
  template <typename ViewType1, typename ViewType2>
  friend void synchronize_dof_iterators(DofIteratorTyped<ViewType1> const
  &it1, DofIterator<ViewType2> &it2);
  */

  template <typename ViewType1, typename Filter1, typename ViewType2, typename Filter2>
  friend void synchronize_dof_iterators(CellDofIterator<ViewType1, Filter1> const &it1,
                                        CellDofIterator<ViewType2, Filter2> &it2);

  /// Step to next cell (prefix increment)
  CellDofIterator &prefix_increment();

  /// Step to the previous cell (prefix decrease)
  CellDofIterator &prefix_decrement();

  /// Dereference operator
  reference_type dereference();

  /// Member access operator
  ptr_type arrow_operator();

  /// Comparison operator
  bool equal(const CellDofIterator &other_it) const;

  /// Comparison operator
  bool not_equal(const CellDofIterator &other_it) const;

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
CellDofIterator<ViewType, Filter>::CellDofIterator() : base_t(), m_view(), m_filter()
{
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
template <typename... FilterArgs>
CellDofIterator<ViewType, Filter>::CellDofIterator(
    typename ViewType::container_ref_type dof_handler, const ActiveIdx idx, FilterArgs &&... args)
    : base_t(), m_view(dof_handler, idx), m_filter(std::forward<FilterArgs>(args)...)
{
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
CellDofIterator<ViewType, Filter>::CellDofIterator(const CellDofIterator &other_iterator)
    : m_view(other_iterator.m_view), m_filter(other_iterator.m_filter)
{
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
CellDofIterator<ViewType, Filter> &CellDofIterator<ViewType, Filter>::operator=(
    const CellDofIterator &other_iterator)
{
  m_view   = other_iterator.m_view;
  m_filter = other_iterator.m_filter;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
std::ptrdiff_t CellDofIterator<ViewType, Filter>::distance(const CellDofIterator &begin,
                                                           const CellDofIterator &end)
{
  return ViewType::distance(begin.m_view, end.m_view);
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
constexpr Uint CellDofIterator<ViewType, Filter>::geo_dim()
{
  return ViewType::GDIM;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
constexpr Uint CellDofIterator<ViewType, Filter>::topo_dim()
{
  return ViewType::TDIM;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
constexpr Uint CellDofIterator<ViewType, Filter>::traversal_dim()
{
  return ViewType::TDIM;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
void CellDofIterator<ViewType, Filter>::print() const
{
  m_view.print();
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
inline CellDofIterator<ViewType, Filter> &CellDofIterator<ViewType, Filter>::prefix_increment()
{
  // m_view.increment();
  /*
  m_filter.view_increment(m_view);
  return *this;
  */

  detail::CellDofIteratorInternals<Filter>::dof_iterator_increment_impl(m_view, m_filter);
  return *this;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
inline CellDofIterator<ViewType, Filter> &CellDofIterator<ViewType, Filter>::prefix_decrement()
{
  // m_view.decrement();
  /*
  m_filter.view_decrement(m_view);
  return *this;
  */

  detail::CellDofIteratorInternals<Filter>::dof_iterator_decrement_impl(m_view, m_filter);
  return *this;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
inline typename CellDofIterator<ViewType, Filter>::reference_type CellDofIterator<
    ViewType, Filter>::dereference()
{
  return *arrow_operator();
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
inline typename CellDofIterator<ViewType, Filter>::ptr_type CellDofIterator<
    ViewType, Filter>::arrow_operator()
{
  return &m_view;
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
inline bool CellDofIterator<ViewType, Filter>::equal(const CellDofIterator &other_it) const
{
  // Use const_cast to use operator* inside comparison, which automatically
  // updates the entity index corresponding to pos *before* comparison (since
  // update of entity delays until request for entity)

  //       return ((const_cast<AllCellIterator *>(this))->operator*()
  //             == (const_cast<AllCellIterator *>(&other_it))->operator*()
  //             && m_connectivity == other_it.m_connectivity && m_pos ==
  // other_it.m_pos );

  return m_view.equal(other_it.m_view);
  /// @todo CHECK THAT FILTERS ARE EQUAL TOO
}

// ----------------------------------------------------------------------------

template <typename ViewType, typename Filter>
inline bool CellDofIterator<ViewType, Filter>::not_equal(const CellDofIterator &other_it) const
{
  return !equal(other_it);
}

// ----------------------------------------------------------------------------

template <typename ViewType1, typename Filter1, typename ViewType2, typename Filter2>
void synchronize_dof_iterators(CellDofIterator<ViewType1, Filter1> const &it1,
                               CellDofIterator<ViewType2, Filter2> &it2)
{
  it2.m_view.synchronize_with(it1.m_view);
}

// ----------------------------------------------------------------------------

/*
template <typename C>
auto begin(C &c) -> decltype((c.begin()));

template <typename C>
auto end(C &c) -> decltype(c.end());
*/

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
