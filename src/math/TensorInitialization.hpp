#ifndef PDEKIT_Math_Tensor_Initialization_hpp
#define PDEKIT_Math_Tensor_Initialization_hpp

#include <vector>

#include "common/PDEKit.hpp"
#include "math/MathMeta.hpp"

namespace pdekit
{

namespace math
{

template <typename T>
class TensorInitializer
{
  public:
  // Constructor
  TensorInitializer(const T val);

  // Bracket operator
  TensorInitializer &operator()(const T val);

  // Access operator
  const T &operator[](const Uint i) const;

  // Return the size of m_data
  Uint size() const;

  private:
  std::vector<T> m_data;
};

/// HELPER TEMPLATE FUNCTION

template <typename T>
TensorInitializer<T> values_list(const T value)
{
  return TensorInitializer<T>(value);
}

// ----------------------------------------------------------------------------

template <typename T>
TensorInitializer<T>::TensorInitializer(const T val)
{
  m_data.push_back(val);
}

// ----------------------------------------------------------------------------

template <typename T>
TensorInitializer<T> &TensorInitializer<T>::operator()(const T val)
{
  m_data.push_back(val);
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T>
const T &TensorInitializer<T>::operator[](const Uint i) const
{
  return m_data[i];
}

// ----------------------------------------------------------------------------

template <typename T>
Uint TensorInitializer<T>::size() const
{
  return m_data.size();
}

// ----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif
