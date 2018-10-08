#ifndef PDEKIT_Mesh_Point_hpp
#define PDEKIT_Mesh_Point_hpp

#include "common/PDEKit.hpp"
#include "math/DenseSVec.hpp"

namespace pdekit
{

namespace mesh
{

/// Class representing one mesh point
/// Internally, it uses a fixed-length vector for storage

template <Uint Dim, typename T = Real>
class Point
{
  public:
  /// Default constructor
  Point();

  /// Copy constructor
  Point(const Point &other_point);

  /// Assignement operator
  Point &operator=(const Point &other_point);

  /// Return the dimension of the point
  constexpr Uint dim() const;

  /// Constant indexing operator
  const T operator[](const Uint idx) const;

  /// Indexing operator, mutable
  T &operator[](const Uint idx);

  private:
  math::DenseSVec<T, Dim> m_data;
};

// ----------------------------------------------------------------------------

template <Uint Dim, typename T>
Point<Dim, T>::Point() : m_data()
{
}

// ----------------------------------------------------------------------------

template <Uint Dim, typename T>
Point<Dim, T>::Point(const Point &other_point) : m_data(other_point.m_data)
{
}

// ----------------------------------------------------------------------------

template <Uint Dim, typename T>
Point<Dim, T> &Point<Dim, T>::operator=(const Point &other_point)
{
  m_data = other_point.m_data;
  return *this;
}

// ----------------------------------------------------------------------------

template <Uint Dim, typename T>
constexpr Uint Point<Dim, T>::dim() const
{
  return Dim;
}

// ----------------------------------------------------------------------------

template <Uint Dim, typename T>
inline const T Point<Dim, T>::operator[](const Uint idx) const
{
  return m_data[idx];
}

// ----------------------------------------------------------------------------

template <Uint Dim, typename T>
inline T &Point<Dim, T>::operator[](const Uint idx)
{
  return m_data[idx];
}

// ----------------------------------------------------------------------------

template <Uint Dim, typename T>
std::ostream &operator<<(std::ostream &os, const Point<Dim, T> &point)
{
  for (Uint d = 0; (d + 1) < point.dim(); ++d)
  {
    os << point[d] << " ";
  }
  os << point[point.dim() - 1];
  return os;
}

// ----------------------------------------------------------------------------=

} // Namespace mesh

} // Namespace pdekit

#endif
