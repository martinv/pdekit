
#ifndef PDEKIT_Common_TupleOf_hpp
#define PDEKIT_Common_TupleOf_hpp

#include <tuple>

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

template <size_t I, typename T>
struct tuple_n
{
  // This concatenates T and the pack 'Args' into 'type'
  template <typename... Args>
  using type = typename tuple_n<I - 1, T>::template type<T, Args...>;
};

template <typename T>
struct tuple_n<0, T>
{
  template <typename... Args>
  using type = std::tuple<Args...>;
};

template <size_t I, typename T>
using tuple_of = typename tuple_n<I, T>::template type<>;

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
