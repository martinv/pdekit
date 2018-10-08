#ifndef PDEKIT_Math_Dense_S_Vec_hpp
#define PDEKIT_Math_Dense_S_Vec_hpp

#include <cmath>
#include <iostream>

#include "math/MathForward.hpp"
#include "math/OperationEvalTime.hpp"
#include "math/TensorInitialization.hpp"
#include "math/TensorRank.hpp"
#include "math/binary_ops/DVecEvalExpr.hpp"
#include "math/binary_ops/VectorOps.hpp"
#include "math/traits/AddTrait.hpp"
#include "math/traits/DivTrait.hpp"
#include "math/traits/MultTrait.hpp"
#include "math/traits/SubTrait.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------
//    Forward declarations to enable friend functions are in MathForward.hpp
// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
std::ostream &operator<<(std::ostream &os, const DenseSVec<T, N, TF> &vec);

template <typename T, Uint N, bool TF>
std::istream &operator>>(std::istream &, DenseSVec<T, N, TF> &vec);

// ----------------------------------------------------------------------------
//             Implementation of class DenseSVec
// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF = DefaultVectorTransposeFlag>
class DenseSVec : public DenseVector<DenseSVec<T, N, TF>, TF>
{

  public:
  using tensor_type    = DenseSVec<T, N, TF>;
  using composite_type = DenseSVec<T, N, TF> const &;
  using value_type     = T;

  enum
  {
    is_expression = 0
  };
  enum
  {
    evaluates_fast = 1
  };
  enum
  {
    owns_data = 1
  };

  /// Default constructor
  explicit DenseSVec();

  // /// Copy constructor
  // DenseSVec ( const DenseSVec<T,TF>& rhs );

  /// Construct the DenseSVec from an initializing list
  /// @param init - list of values to build the DenseSVec from
  DenseSVec(const TensorInitializer<T> &init);

  /// Construct from another dense vector
  template <typename VT>
  explicit DenseSVec(const DenseVector<VT, TF> &init);

  /// Destructor
  ~DenseSVec();

  /// Indexing operator
  /// @param i - index of the element which should be returned
  /// @return the i-th element
  inline T &operator[](const Uint i);

  /// Indexing operator, const version
  /// @param i - index of the element which should be returned
  /// @return the i-th element
  inline const T &operator[](const Uint i) const;

  template <typename VT2>
  inline void assign(const math::DenseVector<VT2, TF> &vec_rhs);

  /// Operator returning the length of the vector
  inline Uint size() const
  {
    return N;
  }

  /// Fill the vector with the same scalar value
  void fill(const T val);

  /// Return a constant subvector (block)
  DenseConstVecView<T, TF> const_block(const Uint first_idx, const Uint size) const;

  /// Return a mutable subvector (block)
  DenseVecView<T, TF> block(const Uint first_idx, const Uint size) const;

  /// Overloaded output operator <<
  friend std::ostream &operator<<<T, N, TF>(std::ostream &os, const DenseSVec<T, N, TF> &vec);

  /// Overloaded input operator >>
  friend std::istream &operator>><T, N, TF>(std::istream &, DenseSVec<T, N, TF> &vec);

  // ----------------------------------------------------------------------------
  //             Definition of Vector operators
  // ----------------------------------------------------------------------------

  /// Vector assignement operator
  DenseSVec<T, N, TF> &operator=(const DenseSVec<T, N, TF> &vec_r);

  /// Assignment operator from a vector expression
  /// @param rhs - a vector expression (for example, a sum of two vectors or
  ///              a product matrix * vector)
  template <typename VT>
  inline DenseSVec<T, N, TF> &operator=(const DenseVector<VT, TF> &rhs);

  /// Accumulation from matrix expression
  template <typename VT>
  inline DenseSVec<T, N, TF> &operator+=(const DenseVector<VT, TF> &rhs);

  /// Subtraction of a matrix expression
  template <typename VT>
  inline DenseSVec<T, N, TF> &operator-=(const DenseVector<VT, TF> &rhs);

  private:
  /// DATA:
  T m_data[N];
};

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
DenseSVec<T, N, TF>::DenseSVec() : DenseVector<DenseSVec<T, N, TF>, TF>(*this)
{
  for (Uint i = 0; i < N; ++i)
  {
    m_data[i] = T();
  }
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
DenseSVec<T, N, TF>::DenseSVec(const TensorInitializer<T> &init)
{
  for (Uint i = 0; i < std::min(N, init.size()); ++i)
  {
    m_data[i] = init[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
template <typename VT>
DenseSVec<T, N, TF>::DenseSVec(const DenseVector<VT, TF> &rhs)
{
  free_assign(*this, rhs.wrapped_type());
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
DenseSVec<T, N, TF>::~DenseSVec()
{
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
T &DenseSVec<T, N, TF>::operator[](const Uint i)
{
  return m_data[i];
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
const T &DenseSVec<T, N, TF>::operator[](const Uint i) const
{
  return m_data[i];
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
template <typename VT2>
void DenseSVec<T, N, TF>::assign(const math::DenseVector<VT2, TF> &vec_rhs)
{
  VT2 const &wrapped_vec_rhs = vec_rhs.wrapped_type();

  for (Uint i = 0; i < N; ++i)
  {
    m_data[i] = wrapped_vec_rhs[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
void DenseSVec<T, N, TF>::fill(const T val)
{
  for (Uint i = 0; i < N; ++i)
  {
    m_data[i] = val;
  }
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
DenseConstVecView<T, TF> DenseSVec<T, N, TF>::const_block(const Uint first_idx,
                                                          const Uint size) const
{
  DenseConstVecView<T, TF> res(m_data + first_idx, size, 1);
  return res;
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
DenseVecView<T, TF> DenseSVec<T, N, TF>::block(const Uint first_idx, const Uint size) const
{
  DenseVecView<T, TF> res(m_data + first_idx, size, 1);
  return res;
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
std::ostream &operator<<(std::ostream &os, const DenseSVec<T, N, TF> &vec)
{
  for (Uint i = 0; i < (N - 1); ++i)
    os << vec[i] << " ";
  os << vec[N - 1];

  return os;
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
std::istream &operator>>(std::istream &is, DenseSVec<T, N, TF> &vec)
{
  for (Uint i = 0; i < N; ++i)
    is >> vec[i];
  return is;
}

// ----------------------------------------------------------------------------
//             Implementation of Vector operators
// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
DenseSVec<T, N, TF> &DenseSVec<T, N, TF>::operator=(const DenseSVec<T, N, TF> &vec_r)
{
  for (Uint i = 0; i < N; ++i)
  {
    m_data[i] = vec_r.m_data[i];
  }
  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
template <typename VT>
DenseSVec<T, N, TF> &DenseSVec<T, N, TF>::operator=(const DenseVector<VT, TF> &rhs)
{
  free_assign(*this, rhs.wrapped_type());

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
template <typename VT>
DenseSVec<T, N, TF> &DenseSVec<T, N, TF>::operator+=(const DenseVector<VT, TF> &rhs)
{
  const VT &rhs_expr = rhs.wrapped_type();

  for (Uint i = 0; i < N; ++i)
  {
    m_data[i] += rhs_expr[i];
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
template <typename VT>
DenseSVec<T, N, TF> &DenseSVec<T, N, TF>::operator-=(const DenseVector<VT, TF> &rhs)
{
  const VT &rhs_expr = rhs.wrapped_type();

  for (Uint i = 0; i < N; ++i)
  {
    m_data[i] -= rhs_expr[i];
  }

  return (*this);
}

// ----------------------------------------------------------------------------
// AddTrait specializations
// ----------------------------------------------------------------------------
template <typename T1, typename T2, Uint N, bool TF>
struct AddTrait<DenseSVec<T1, N, TF>, DenseSVec<T2, N, TF>>
{
  using type = DenseSVec<typename AddTrait<T1, T2>::type, N, TF>;
};

// ----------------------------------------------------------------------------
// SubTrait specializations
// ----------------------------------------------------------------------------

template <typename T1, typename T2, Uint N, bool TF>
struct SubTrait<DenseSVec<T1, N, TF>, DenseSVec<T2, N, TF>>
{
  using type = DenseSVec<typename SubTrait<T1, T2>::type, N, TF>;
};

// ----------------------------------------------------------------------------
// MultTrait specializations
// ----------------------------------------------------------------------------

// scalar * DenseSVec
template <typename T, Uint N, bool TF>
struct MultTrait<Real, DenseSVec<T, N, TF>>
{
  using type = DenseSVec<typename MultTrait<Real, T>::type, N, TF>;
};

// ----------------------------------------------------------------------------
// DivTrait specializations
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TENSOR RANK OF VECTOR = 1
// ----------------------------------------------------------------------------

template <typename T, Uint N, bool TF>
struct TensorRank<DenseSVec<T, N, TF>>
{
  enum
  {
    value = tensor_rank_1
  };
};

} // Namespace math

} // Namespace pdekit

#endif
