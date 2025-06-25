#ifndef Math_Dense_Vec_View_hpp
#define Math_Dense_Vec_View_hpp

#include <cmath>
#include <iostream>

#include "math/MathForward.hpp"
#include "math/OperationEvalTime.hpp"
#include "math/TensorInitialization.hpp"
#include "math/TensorRank.hpp"
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
//             Implementation of class DenseVecView
// ----------------------------------------------------------------------------

template <typename T, bool TF = DefaultVectorTransposeFlag>
class DenseVecView : public DenseVector<DenseVecView<T, TF>, TF>
{

  public:
  using tensor_type    = DenseVecView<T, TF>;
  using composite_type = DenseVecView<T, TF> const &;
  using value_type     = T;

  using ref_type =
      typename std::conditional<std::is_const<T>::value, const value_type &, value_type &>::type;
  using ptr_type =
      typename std::conditional<std::is_const<T>::value, const value_type *, value_type *>::type;

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
    owns_data = 0
  };

  /// Default constructor
  explicit DenseVecView();

  /// Constructor
  explicit DenseVecView(ptr_type first, ptr_type last, const Uint stride = 1);

  /// Constructor
  explicit DenseVecView(ptr_type first, const Uint n, const Uint stride = 1);

  /// Copy constructor
  DenseVecView(const DenseVecView<T, TF> &slice_r);

  /// Destructor
  ~DenseVecView();

  /// Indexing operator
  ref_type operator[](const Uint i);

  /// Indexing operator, const version
  const value_type &operator[](const Uint i) const;

  /// Operator returning the length of the vector
  Uint size() const;

  /// Return a slice of the current view
  const DenseVecView slice(const Uint idx_first, const Uint length) const;

  /// Fill the vector with the same scalar value
  void fill(const T val);

  // ----------------------------------------------------------------------------
  //                         DenseVecView operators
  // ----------------------------------------------------------------------------

  /// DenseVecView assignement operator
  DenseVecView &operator=(const DenseVecView<T, TF> &view_r);

  /// Assignment operator from a vector expression
  /// @param rhs - a vector expression (for example, a sum of two vectors or
  ///              a product matrix * vector)
  template <typename VT>
  inline DenseVecView<T, TF> &operator=(const DenseVector<VT, TF> &rhs);

  /// Accumulation from matrix expression
  template <typename VT>
  inline DenseVecView<T, TF> &operator+=(const DenseVector<VT, TF> &rhs);

  /// Subtraction of a matrix expression
  template <typename VT>
  inline DenseVecView<T, TF> &operator-=(const DenseVector<VT, TF> &rhs);

  private:
  /// Pointer to raw data, first entry
  ptr_type m_first;

  /// Pointer to last entry of raw data
  ptr_type m_last;

  /// Stride - how many entries to skip between entries [i] and [i+1]
  Uint m_stride;
};

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseVecView<T, TF>::DenseVecView()
    : DenseVector<DenseVecView<T, TF>, TF>(*this), m_first(nullptr), m_last(nullptr), m_stride(1)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseVecView<T, TF>::DenseVecView(ptr_type first, ptr_type last, const Uint stride)
    : DenseVector<DenseVecView<T, TF>, TF>(*this), m_first(first), m_last(last), m_stride(stride)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseVecView<T, TF>::DenseVecView(ptr_type first, const Uint n, const Uint stride)
    : DenseVector<DenseVecView<T, TF>, TF>(*this), m_first(first), m_last(first + (n - 1) * stride),
      m_stride(stride)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseVecView<T, TF>::DenseVecView(const DenseVecView<T, TF> &slice_r)
{
  m_first  = slice_r.m_first;
  m_last   = slice_r.m_last;
  m_stride = slice_r.m_stride;
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseVecView<T, TF>::~DenseVecView()
{
  m_first  = 0;
  m_last   = 0;
  m_stride = 1;
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
inline typename DenseVecView<T, TF>::ref_type DenseVecView<T, TF>::operator[](const Uint i)
{
  return m_first[i * m_stride];
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
inline const typename DenseVecView<T, TF>::value_type &DenseVecView<T, TF>::operator[](
    const Uint i) const
{
  return m_first[i * m_stride];
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
inline Uint DenseVecView<T, TF>::size() const
{
  return (m_last - m_first) / m_stride + 1;
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
const DenseVecView<T, TF> DenseVecView<T, TF>::slice(const Uint idx_first, const Uint length) const
{
  const Uint max_len = idx_first + length > size() ? size() - idx_first : length;
  DenseVecView<T, TF> result(m_first + idx_first * m_stride,
                             m_first + (idx_first + max_len - 1) * m_stride, m_stride);
  return result;
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
void DenseVecView<T, TF>::fill(const T val)
{
  for (T *it = m_first; it <= m_last; ++it)
  {
    *it = val;
  }
}

// ----------------------------------------------------------------------------
//             Implementation of DenseVewView operators
// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseVecView<T, TF> &DenseVecView<T, TF>::operator=(const DenseVecView<T, TF> &view_r)
{
  m_first  = view_r.m_first;
  m_last   = view_r.m_last;
  m_stride = view_r.m_stride;
  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
template <typename VT>
DenseVecView<T, TF> &DenseVecView<T, TF>::operator=(const DenseVector<VT, TF> &rhs)
{
  static_assert((std::is_const<T>::value == std::is_const<typename VT::value_type>::value) ||
                    (std::is_const<T>::value && !std::is_const<typename VT::value_type>::value),
                "Can't create non-const view from const vector.");

  const VT &rhs_expr = rhs.wrapped_type();

  const Uint size = (m_last - m_first) / m_stride + 1;
  for (Uint i = 0; i < size; ++i)
  {
    m_first[i * m_stride] = rhs_expr[i];
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
template <typename VT>
DenseVecView<T, TF> &DenseVecView<T, TF>::operator+=(const DenseVector<VT, TF> &rhs)
{
  const VT &rhs_expr = rhs.wrapped_type();

  const Uint size = (m_last - m_first) / m_stride + 1;
  for (Uint i = 0; i < size; ++i)
  {
    m_first[i * m_stride] += rhs_expr[i];
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
template <typename VT>
DenseVecView<T, TF> &DenseVecView<T, TF>::operator-=(const DenseVector<VT, TF> &rhs)
{
  const VT &rhs_expr = rhs.wrapped_type();

  const Uint size = (m_last - m_first) / m_stride + 1;
  for (Uint i = 0; i < size; ++i)
  {
    m_first[i * m_stride] -= rhs_expr[i];
  }

  return (*this);
}

// ----------------------------------------------------------------------------
// AddTrait specializations
// ----------------------------------------------------------------------------

template <typename T1, typename T2, bool TF>
struct AddTrait<DenseVecView<T1, TF>, DenseVecView<T2, TF>>
{
  using type = DenseDVec<typename AddTrait<T1, T2>::type, TF>;
};

// DenseDVec + DenseVecView
template <typename T1, typename T2, bool TF>
struct AddTrait<DenseDVec<T1, TF>, DenseVecView<T2, TF>>
{
  using type = DenseDVec<typename AddTrait<T1, T2>::type, TF>;
};

// DenseVecView + DenseDVec
template <typename T1, typename T2, bool TF>
struct AddTrait<DenseVecView<T1, TF>, DenseDVec<T2, TF>>
{
  using type = DenseDVec<typename AddTrait<T1, T2>::type, TF>;
};

// ----------------------------------------------------------------------------
// SubTrait specializations
// ----------------------------------------------------------------------------

template <typename T1, typename T2, bool TF>
struct SubTrait<DenseVecView<T1, TF>, DenseVecView<T2, TF>>
{
  using type = DenseDVec<typename SubTrait<T1, T2>::type, TF>;
};

template <typename T1, typename T2, bool TF>
struct SubTrait<DenseDVec<T1, TF>, DenseVecView<T2, TF>>
{
  using type = DenseDVec<typename SubTrait<T1, T2>::type, TF>;
};

template <typename T1, typename T2, bool TF>
struct SubTrait<DenseVecView<T1, TF>, DenseDVec<T2, TF>>
{
  using type = DenseDVec<typename SubTrait<T1, T2>::type, TF>;
};

template <typename T1, Uint N, bool TF, typename T2>
struct SubTrait<DenseSVec<T1, N, TF>, DenseVecView<T2, TF>>
{
  using type = DenseSVec<typename SubTrait<T1, T2>::type, N, TF>;
};

template <typename T1, bool TF, typename T2, Uint N>
struct SubTrait<DenseVecView<T1, TF>, DenseSVec<T2, N, TF>>
{
  using type = DenseSVec<typename SubTrait<T1, T2>::type, N, TF>;
};

// ----------------------------------------------------------------------------
// MultTrait specializations
// ----------------------------------------------------------------------------

// scalar * DenseVecView
template <typename SValueT, typename T, bool TF>
struct MultTrait<SValueT, DenseVecView<T, TF>>
{
  using type = DenseDVec<typename MultTrait<SValueT, T>::type, TF>;
};

// ----------------------------------------------------------------------------
// DivTrait specializations
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TENSOR RANK OF VECTOR = 1
// ----------------------------------------------------------------------------

template <typename T, bool TF>
struct TensorRank<DenseVecView<T, TF>>
{
  enum
  {
    value = tensor_rank_1
  };
};

// ----------------------------------------------------------------------------

template <typename T, bool TF>
std::ostream &operator<<(std::ostream &os, const DenseVecView<T, TF> &slice)
{
  //   os.clear();
  for (Uint i = 0; (i + 1) < slice.size(); ++i)
  {
    os << slice[i] << " ";
  }
  if (slice.size() > 0)
  {
    os << slice[slice.size() - 1];
  }

  return os;
}

// ----------------------------------------------------------------------------

} // Namespace math

} // Namespace pdekit

#endif
