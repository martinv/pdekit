#ifndef PDEKIT_Common_Tuple_Meta_hpp
#define PDEKIT_Common_Tuple_Meta_hpp

#include <vector>

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace common
{

namespace meta
{

// ----------------------------------------------------------------------------

template <Uint N, typename... T>
struct tuple_nth_elem_type;

template <typename T0, typename... T>
struct tuple_nth_elem_type<0, T0, T...>
{
  using type = T0;
};

template <Uint N, typename T0, typename... T>
struct tuple_nth_elem_type<N, T0, T...>
{
  using type = typename tuple_nth_elem_type<N - 1, T...>::type;
};

// ----------------------------------------------------------------------------

template <typename... T>
struct TransformToVector
{
  using type = std::tuple<std::vector<T>...>;
};

// ----------------------------------------------------------------------------

} // Namespace meta

} // Namespace common

} // Namespace pdekit

#endif
