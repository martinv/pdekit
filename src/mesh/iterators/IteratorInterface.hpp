#ifndef PDEKIT_Mesh_Cell_Iterator_Interface_hpp
#define PDEKIT_Mesh_Cell_Iterator_Interface_hpp

#include <type_traits>

namespace pdekit
{

namespace mesh
{

/// ---------------------------------------------------------------------------
/// Helper class to perform operations on the iterator
/// This class is a friend of CellIterator and AllCellIterator and therefore
/// has access to their private methods (see CellIterator.hpp)
/// ---------------------------------------------------------------------------

class IteratorAccess
{
  public:
  /// Perform prefix increment (++it)
  template <typename IterType>
  inline static void prefix_increment(IterType &it)
  {
    it.prefix_increment();
  }

  /// Perform prefix decrement (--it)
  template <typename IterType>
  inline static void prefix_decrement(IterType &it)
  {
    it.prefix_decrement();
  }

  template <typename IterType>
  inline static typename IterType::reference_type dereference(IterType &it)
  {
    return it.dereference();
  }

  template <typename IterType>
  inline static typename IterType::ptr_type arrow_operator(IterType &it)
  {
    return it.arrow_operator();
  }

  template <typename IterType>
  inline static bool equal(const IterType &it_a, const IterType &it_b)
  {
    return it_a.equal(it_b);
  }

  template <typename IterType>
  inline static bool not_equal(const IterType &it_a, const IterType &it_b)
  {
    return it_a.not_equal(it_b);
  }
};

/// ---------------------------------------------------------------------------
/// CellIteratorInterface: use CRTP to define the interface of the iterator
/// The actual iterators that inherit from this class are in CellIterator.hpp
/// ---------------------------------------------------------------------------

template <typename DerivedIt, typename ViewType>
class IteratorInterface
{
  private:
  using value_type =
      typename std::conditional<ViewType::is_immutable_proxy, const ViewType, ViewType>::type;

  public:
  /// Increment the iterator
  inline DerivedIt &operator++();

  /// Decrease the operator
  inline DerivedIt &operator--();

  /// Get the value the iterator points to
  inline value_type &operator*();

  /// Get the address (pointer) the iterator points to
  inline value_type *operator->();

  private:
  inline DerivedIt &derived()
  {
    return static_cast<DerivedIt &>(*this);
  }

  inline DerivedIt const &derived() const
  {
    return static_cast<DerivedIt const &>(*this);
  }
};

// ----------------------------------------------------------------------------

template <typename DerivedIt, typename ViewType>
DerivedIt &IteratorInterface<DerivedIt, ViewType>::operator++()
{
  IteratorAccess::prefix_increment(derived());
  return derived();
}

// ----------------------------------------------------------------------------

template <typename DerivedIt, typename ViewType>
DerivedIt &IteratorInterface<DerivedIt, ViewType>::operator--()
{
  IteratorAccess::prefix_decrement(derived());
  return derived();
}

// ----------------------------------------------------------------------------

template <typename DerivedIt, typename ViewType>
typename IteratorInterface<DerivedIt, ViewType>::value_type
    &IteratorInterface<DerivedIt, ViewType>::operator*()
{
  return IteratorAccess::dereference(derived());
}

// ----------------------------------------------------------------------------

template <typename DerivedIt, typename ViewType>
typename IteratorInterface<DerivedIt, ViewType>::value_type
    *IteratorInterface<DerivedIt, ViewType>::operator->()
{
  return IteratorAccess::arrow_operator(derived());
}

// ----------------------------------------------------------------------------
// Free functions to determine if two iterators are equal or different
// ----------------------------------------------------------------------------

template <typename DerivedIt, typename ViewType>
inline bool operator==(const IteratorInterface<DerivedIt, ViewType> &it_lhs,
                       const IteratorInterface<DerivedIt, ViewType> &it_rhs)
{
  return IteratorAccess::equal(static_cast<DerivedIt const &>(it_lhs),
                               static_cast<DerivedIt const &>(it_rhs));
}

// ----------------------------------------------------------------------------

template <typename DerivedIt, typename ViewType>
inline bool operator!=(const IteratorInterface<DerivedIt, ViewType> &it_lhs,
                       const IteratorInterface<DerivedIt, ViewType> &it_rhs)
{
  return IteratorAccess::not_equal(static_cast<DerivedIt const &>(it_lhs),
                                   static_cast<DerivedIt const &>(it_rhs));
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
