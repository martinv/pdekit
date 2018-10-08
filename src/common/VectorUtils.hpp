#ifndef PDEKIT_Common_VECTOR_UTILS_hpp
#define PDEKIT_Common_VECTOR_UTILS_hpp

#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

#include "common/PDEKit.hpp"

namespace pdekit
{

// -----------------------------------------------------------------------------

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec)
{
  for (typename std::vector<T>::size_type i = 0; (i + 1) < vec.size(); ++i)
  {
    os << vec[i] << " ";
  }
  os << vec.back();
  return os;
}

namespace common
{

// -----------------------------------------------------------------------------
// VectorPack is a class that holds a tuple of reference wrappers to std::vector
// -----------------------------------------------------------------------------

// Every argument in the parameter pack 'Types' should be of type std::vector< ... >

template <typename... Types>
class VectorPack
{
  public:
  using pack_t = std::tuple<std::reference_wrapper<Types>...>;

  VectorPack(Types &... pack) : m_pack(std::tuple<std::reference_wrapper<Types>...>(pack...))
  {
  }

  template <Uint Pos>
  typename std::tuple_element<Pos, pack_t>::type get()
  {
    return typename std::tuple_element<Pos, pack_t>::type(std::get<Pos>(m_pack));
  }

  private:
  pack_t m_pack;
};

// -----------------------------------------------------------------------------

template <Uint Pos, typename... Types>
struct VectorPackEntryType
{
  using pack_t = std::tuple<std::reference_wrapper<Types>...>;
  // the first ::type returns std::reference_wrapper< ... >
  // the second ::type returns the underlying type, which should be std::vector< ... >
  using type = typename std::tuple_element<Pos, pack_t>::type::type;
};

// -----------------------------------------------------------------------------

template <typename... Types>
VectorPack<Types...> pack_vectors(Types &... args)
{
  return VectorPack<Types...>(args...);
}

// -----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
