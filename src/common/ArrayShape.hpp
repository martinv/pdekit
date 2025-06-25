#ifndef PDEKIT_Common_Array_Shape_hpp
#define PDEKIT_Common_Array_Shape_hpp

#include "common/PDEKit.hpp"
#include <ostream>

namespace pdekit
{

namespace common
{

template <Uint DIM, typename SizeType>
class ArrayShape;

template <Uint DIM, typename SizeType>
bool operator==(const ArrayShape<DIM, SizeType> &lhs, const ArrayShape<DIM, SizeType> &rhs);

template <Uint DIM, typename SizeType>
bool operator!=(const ArrayShape<DIM, SizeType> &lhs, const ArrayShape<DIM, SizeType> &rhs);

template <Uint DIM, typename SizeType>
class ArrayShape
{
  public:
  /// TYPEDEFS
  using size_type = SizeType;

  /// Default constructor
  ArrayShape();

  /// Construct from sizes in each dimension
  template <typename... Ts>
  ArrayShape(Ts... args);

  /// Copy constructor
  ArrayShape(const ArrayShape &other);

  /// Destructor
  ~ArrayShape();

  /// Assigment operator
  ArrayShape &operator=(const ArrayShape &rhs);

  /// Get the size in one dimension
  SizeType size(const Uint dim) const;

  /// Return the number of dimensions in shape
  constexpr Uint nb_dims() const;

  /// Check if two shapes are the same
  friend bool operator==
      <DIM, SizeType>(const ArrayShape<DIM, SizeType> &lhs, const ArrayShape<DIM, SizeType> &rhs);

  /// Check if two shapes are the same
  friend bool operator!=
      <DIM, SizeType>(const ArrayShape<DIM, SizeType> &lhs, const ArrayShape<DIM, SizeType> &rhs);

  private:
  /// TYPES
  struct PosIdx
  {
    Uint value;
    SizeType *dim_storage;
  };

  /// FUNCTIONS
  template <typename LastSize>
  void process_params(PosIdx &pos_idx, LastSize size)
  {
    pos_idx.dim_storage[pos_idx.value] = size;
  }

  template <typename Size, typename... OtherSizes>
  void process_params(PosIdx &pos_idx, Size size, OtherSizes... sizes)
  {
    pos_idx.dim_storage[pos_idx.value] = size;
    pos_idx.value++;
    process_params(pos_idx, sizes...);
  }

  /// DATA
  SizeType m_sizes[DIM];
};

// ----------------------------------------------------------------------------

template <Uint DIM, typename SizeType>
ArrayShape<DIM, SizeType>::ArrayShape()
{
  for (Uint dim = 0; dim < DIM; ++dim)
  {
    m_sizes[dim] = SizeType{};
  }
}

// ----------------------------------------------------------------------------

template <Uint DIM, typename SizeType>
template <typename... Ts>
ArrayShape<DIM, SizeType>::ArrayShape(Ts... args)
{
  static_assert(sizeof...(args) <= DIM, "ArrayShape constructor: too many arguments");
  for (Uint dim = 0; dim < DIM; ++dim)
  {
    m_sizes[dim] = 0;
  }

  PosIdx pos_idx;
  pos_idx.value       = 0;
  pos_idx.dim_storage = &m_sizes[0];
  process_params(pos_idx, args...);

  /*
  std::cout << "Constructed values of array shape:" << std::endl;
  for (Uint dim = 0; dim < DIM; ++dim)
  {
    std::cout << "[" << dim << "] = " << m_sizes[dim] << std::endl;
  }
  */
}

// ----------------------------------------------------------------------------

template <Uint DIM, typename SizeType>
ArrayShape<DIM, SizeType>::ArrayShape(const ArrayShape &other)
{
  for (Uint dim = 0; dim < DIM; ++dim)
  {
    m_sizes[dim] = other.m_sizes[dim];
  }
}

// ----------------------------------------------------------------------------

template <Uint DIM, typename SizeType>
ArrayShape<DIM, SizeType>::~ArrayShape()
{
}

// ----------------------------------------------------------------------------

template <Uint DIM, typename SizeType>
ArrayShape<DIM, SizeType> &ArrayShape<DIM, SizeType>::operator=(const ArrayShape &rhs)
{
  for (Uint dim = 0; dim < DIM; ++dim)
  {
    m_sizes[dim] = rhs.m_sizes[dim];
  }
  return *this;
}

// ----------------------------------------------------------------------------

template <Uint DIM, typename SizeType>
inline SizeType ArrayShape<DIM, SizeType>::size(const Uint dim) const
{
  return m_sizes[dim];
}

// ----------------------------------------------------------------------------

template <Uint DIM, typename SizeType>
inline constexpr Uint ArrayShape<DIM, SizeType>::nb_dims() const
{
  return DIM;
}

// ----------------------------------------------------------------------------

template <Uint DIM, typename SizeType>
std::ostream &operator<<(std::ostream &os, const ArrayShape<DIM, SizeType> &arr_shape)
{
  os << "[";
  for (Uint d = 0; d < DIM; ++d)
  {
    os << " " << arr_shape.size(d);
  }

  os << " ]";
  return os;
}

// ----------------------------------------------------------------------------

template <Uint DIM, typename SizeType>
bool operator==(const ArrayShape<DIM, SizeType> &lhs, const ArrayShape<DIM, SizeType> &rhs)
{
  for (Uint d = 0; d < DIM; ++d)
  {
    if (lhs.m_sizes[d] != rhs.m_sizes[d])
    {
      return false;
    }
  }
  return true;
}

// ----------------------------------------------------------------------------

template <Uint DIM, typename SizeType>
bool operator!=(const ArrayShape<DIM, SizeType> &lhs, const ArrayShape<DIM, SizeType> &rhs)
{
  return !operator==(lhs, rhs);
}

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit

#endif // PDEKIT_Common_Array_Shape_hpp
