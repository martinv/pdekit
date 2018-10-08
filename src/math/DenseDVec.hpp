#ifndef PDEKIT_Math_Dense_D_Vec_hpp
#define PDEKIT_Math_Dense_D_Vec_hpp

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
//             Forward declarations to enable friend functions
// ----------------------------------------------------------------------------

// template<typename T, bool TF> class DenseDVec;

template <typename T, bool TF>
std::ostream &operator<<(std::ostream &os, const DenseDVec<T, TF> &vec);

template <typename T, bool TF>
std::istream &operator>>(std::istream &, DenseDVec<T, TF> &vec);

// ----------------------------------------------------------------------------
//             Implementation of class DenseDVec
// ----------------------------------------------------------------------------

template <typename T, bool TF = DefaultVectorTransposeFlag>
class DenseDVec : public DenseVector<DenseDVec<T, TF>, TF>
{

  public:
  using tensor_type    = DenseDVec<T, TF>;
  using composite_type = DenseDVec<T, TF> const &;
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
  explicit DenseDVec();

  /// Constructor, set the number of rows
  explicit DenseDVec(Uint m);

  /// The copy constructor is explicitly defined due to the required
  /// dynamic memory management and in order to enable/facilitate NRV
  /// optimization.
  DenseDVec(const DenseDVec &rhs);

  /// Construct the DenseDVec from an initializing list
  /// @param init - list of values to build the DenseDVec from
  DenseDVec(const TensorInitializer<T> &init);

  /// Construct from another dense vector
  template <typename VT>
  explicit DenseDVec(const DenseVector<VT, TF> &init);

  /// Destructor
  ~DenseDVec();

  /// Resize the vector
  /// @param new_size - the new number of rows of the vector
  void resize(const Uint new_size);

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
    return m_size;
  }

  /// Fill the vector with the same scalar value
  void fill(const T val);

  /// Return a constant subvector (block)
  DenseConstVecView<T, TF> const_block(const Uint first_idx, const Uint size) const;

  /// Return a mutable subvector (block)
  DenseVecView<T, TF> block(const Uint first_idx, const Uint size) const;

  /// Overloaded output operator <<
  friend std::ostream &operator<<<T, TF>(std::ostream &os, const DenseDVec<T, TF> &vec);

  /// Overloaded input operator >>
  friend std::istream &operator>><T, TF>(std::istream &, DenseDVec<T, TF> &vec);

  // --------------------------------------------------------------------------
  //             Definition of Vector operators
  // --------------------------------------------------------------------------

  /// Vector assignement operator
  DenseDVec<T, TF> &operator=(const DenseDVec<T, TF> &vec_r);

  /// Assignment operator from a vector expression
  /// @param rhs - a vector expression (for example, a sum of two vectors or
  ///              a product matrix * vector)
  template <typename VT>
  inline DenseDVec<T, TF> &operator=(const DenseVector<VT, TF> &rhs);

  /// Accumulation from matrix expression
  template <typename VT>
  inline DenseDVec<T, TF> &operator+=(const DenseVector<VT, TF> &rhs);

  /// Subtraction of a matrix expression
  template <typename VT>
  inline DenseDVec<T, TF> &operator-=(const DenseVector<VT, TF> &rhs);

  private:
  /// DATA:
  T *m_data;

  Uint m_size;
};

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseDVec<T, TF>::DenseDVec() : DenseVector<DenseDVec<T, TF>, TF>(*this), m_data(0), m_size(0)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseDVec<T, TF>::DenseDVec(Uint m) : DenseVector<DenseDVec<T, TF>, TF>(*this)
{
  m_data = new T[m];
  m_size = m;
  for (Uint i = 0; i < m; ++i)
  {
    m_data[i] = T();
  }
  // resize(m);
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseDVec<T, TF>::DenseDVec(const DenseDVec &rhs)
{
  m_data = new T[rhs.size()];
  m_size = rhs.size();

  for (Uint i = 0; i < m_size; ++i)
  {
    m_data[i] = rhs.m_data[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseDVec<T, TF>::DenseDVec(const TensorInitializer<T> &init)
{
  // resize(init.size());

  /*
  // This is wrong: the variable m_size has to be properly initialized
  // first. m_size is not necessarily 0!
  // The if condition might randomly fail ...
  if (m_size == 0u)
  {
    m_data = new T[init.size()];
    m_size = init.size();
  }
  */

  m_data = new T[init.size()];
  m_size = init.size();

  for (Uint i = 0; i < std::min(m_size, init.size()); ++i)
  {
    m_data[i] = init[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
template <typename VT>
DenseDVec<T, TF>::DenseDVec(const DenseVector<VT, TF> &rhs)
{
  const VT &rhs_expr = rhs.wrapped_type();

  // resize(rhs_expr.size());

  m_data = new T[rhs_expr.size()];
  m_size = rhs_expr.size();

  free_assign(*this, rhs.wrapped_type());
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseDVec<T, TF>::~DenseDVec()
{
  delete[] m_data;
  m_data = 0;
  m_size = 0;
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
void DenseDVec<T, TF>::resize(const Uint new_size)
{
  if (new_size == 0)
  {
    delete[] m_data;
    m_data = 0;
    m_size = 0;
    return;
  }

  if (m_size != new_size)
  {
    delete[] m_data;
    m_data = new T[new_size];
    m_size = new_size;
  }

  fill(T());
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
T &DenseDVec<T, TF>::operator[](const Uint i)
{
  return m_data[i];
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
const T &DenseDVec<T, TF>::operator[](const Uint i) const
{
  return m_data[i];
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
template <typename VT2>
void DenseDVec<T, TF>::assign(const math::DenseVector<VT2, TF> &vec_rhs)
{
  VT2 const &wrapped_vec_rhs = vec_rhs.wrapped_type();

  for (Uint i = 0; i < m_size; ++i)
  {
    m_data[i] = wrapped_vec_rhs[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
void DenseDVec<T, TF>::fill(const T val)
{
  for (Uint i = 0; i < m_size; ++i)
  {
    m_data[i] = val;
  }
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseConstVecView<T, TF> DenseDVec<T, TF>::const_block(const Uint first_idx, const Uint size) const
{
  DenseConstVecView<T, TF> res(m_data + first_idx, size, 1);
  return res;
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseVecView<T, TF> DenseDVec<T, TF>::block(const Uint first_idx, const Uint size) const
{
  DenseVecView<T, TF> res(m_data + first_idx, size, 1);
  return res;
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
std::ostream &operator<<(std::ostream &os, const DenseDVec<T, TF> &vec)
{
  for (Uint i = 0; i < (vec.m_size - 1); ++i)
    os << vec[i] << " ";
  os << vec[vec.m_size - 1];

  return os;
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
std::istream &operator>>(std::istream &is, DenseDVec<T, TF> &vec)
{
  for (Uint i = 0; i < vec.m_size; ++i)
    is >> vec[i];
  return is;
}

// ----------------------------------------------------------------------------
//             Implementation of Vector operators
// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseDVec<T, TF> &DenseDVec<T, TF>::operator=(const DenseDVec<T, TF> &vec_r)
{
  for (Uint i = 0; i < m_size; ++i)
  {
    m_data[i] = vec_r.m_data[i];
  }
  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
template <typename VT>
DenseDVec<T, TF> &DenseDVec<T, TF>::operator=(const DenseVector<VT, TF> &rhs)
{
  const VT &rhs_expr = rhs.wrapped_type();

  if (m_size != rhs_expr.size())
  {
    resize(rhs_expr.size());
  }

  free_assign(*this, rhs_expr);
  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
template <typename VT>
DenseDVec<T, TF> &DenseDVec<T, TF>::operator+=(const DenseVector<VT, TF> &rhs)
{
  const VT &rhs_expr = rhs.wrapped_type();

  // Do not use 'resize', because 'resize' sets all values to 0, but we want
  // to ACCUMULATE here - all previously accumulated values would be
  // destroyed! resize(rhs_expr.size());

  for (Uint i = 0; i < m_size; ++i)
  {
    m_data[i] += rhs_expr[i];
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
template <typename VT>
DenseDVec<T, TF> &DenseDVec<T, TF>::operator-=(const DenseVector<VT, TF> &rhs)
{
  const VT &rhs_expr = rhs.wrapped_type();

  // Do not use 'resize', because 'resize' sets all values to 0, but we want
  // to ACCUMULATE here - all previously accumulated values would be
  // destroyed! resize(rhs_expr.size());

  for (Uint i = 0; i < m_size; ++i)
  {
    m_data[i] -= rhs_expr[i];
  }

  return (*this);
}

// ----------------------------------------------------------------------------
// AddTrait specializations
// ----------------------------------------------------------------------------

template <typename T1, typename T2, bool TF>
struct AddTrait<DenseDVec<T1, TF>, DenseDVec<T2, TF>>
{
  using type = DenseDVec<typename AddTrait<T1, T2>::type, TF>;
};

template <typename T1, bool TF, typename T2, Uint N>
struct AddTrait<DenseDVec<T1, TF>, DenseSVec<T2, N, TF>>
{
  using type = DenseDVec<typename AddTrait<T1, T2>::type, TF>;
};

template <typename T1, Uint N, bool TF, typename T2>
struct AddTrait<DenseSVec<T1, N, TF>, DenseDVec<T2, TF>>
{
  using type = DenseDVec<typename AddTrait<T1, T2>::type, TF>;
};

// ----------------------------------------------------------------------------
// SubTrait specializations
// ----------------------------------------------------------------------------

template <typename T1, typename T2, bool TF>
struct SubTrait<DenseDVec<T1, TF>, DenseDVec<T2, TF>>
{
  using type = DenseDVec<typename SubTrait<T1, T2>::type, TF>;
};

template <typename T1, bool TF, typename T2, Uint N>
struct SubTrait<DenseDVec<T1, TF>, DenseSVec<T2, N, TF>>
{
  using type = DenseDVec<typename SubTrait<T1, T2>::type, TF>;
};

template <typename T1, Uint N, bool TF, typename T2>
struct SubTrait<DenseSVec<T1, N, TF>, DenseDVec<T2, TF>>
{
  using type = DenseDVec<typename SubTrait<T1, T2>::type, TF>;
};

// ----------------------------------------------------------------------------
// MultTrait specializations
// ----------------------------------------------------------------------------

// scalar * DenseDVec
template <typename T, bool TF>
struct MultTrait<Real, DenseDVec<T, TF>>
{
  using type = DenseDVec<typename MultTrait<Real, T>::type, TF>;
};

// ----------------------------------------------------------------------------
// DivTrait specializations
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TENSOR RANK OF VECTOR = 1
// ----------------------------------------------------------------------------

template <typename T, bool TF>
struct TensorRank<DenseDVec<T, TF>>
{
  enum
  {
    value = tensor_rank_1
  };
};

} // Namespace math

} // Namespace pdekit

#endif
