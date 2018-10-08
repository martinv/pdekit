#ifndef Math_Dense_Const_Vec_View_hpp
#define Math_Dense_Const_Vec_View_hpp

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
//             Implementation of class DenseConstVecView
// ----------------------------------------------------------------------------

template <typename T, bool TF = DefaultVectorTransposeFlag>
class DenseConstVecView : public DenseVector<DenseConstVecView<T, TF>, TF>
{

  public:
  using tensor_type    = DenseConstVecView<T, TF>;
  using composite_type = DenseConstVecView<T, TF> const &;
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
    owns_data = 0
  };

  /// Default constructor
  explicit DenseConstVecView();

  /// Constructor
  explicit DenseConstVecView(const T *first, const T *last, const Uint stride = 1);

  /// Constructor
  explicit DenseConstVecView(const T *first, const Uint n, const Uint stride = 1);

  /// Copy constructor
  DenseConstVecView(const DenseConstVecView<T, TF> &slice_r);

  /// Destructor
  ~DenseConstVecView();

  /// Indexing operator
  //  T& operator[] ( const Uint i );

  /// Indexing operator, const version
  const T &operator[](const Uint i) const;

  /// Operator returning the length of the vector
  Uint size() const;

  /// Return a slice of the current view
  const DenseConstVecView slice(const Uint idx_first, const Uint length) const;

  // ----------------------------------------------------------------------------
  //                      DenseConstVecView operators
  // ----------------------------------------------------------------------------

  /// DenseConstVecView assignment operator
  DenseConstVecView &operator=(const DenseConstVecView<T, TF> &slice_r);

  private:
  /// Pointer to raw data, first entry
  T const *m_first; // Pointer to a constant object

  /// Pointer to last entry of raw data
  T const *m_last;

  /// Stride - how many entries to skip between entries [i] and [i+1]
  Uint m_stride;
};

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseConstVecView<T, TF>::DenseConstVecView()
    : DenseVector<DenseConstVecView<T, TF>, TF>(*this), m_first(nullptr), m_last(nullptr),
      m_stride(1)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseConstVecView<T, TF>::DenseConstVecView(const T *first, const T *last, const Uint stride)
    : DenseVector<DenseConstVecView<T, TF>, TF>(*this), m_first(first), m_last(last),
      m_stride(stride)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseConstVecView<T, TF>::DenseConstVecView(const T *first, const Uint n, const Uint stride)
    : DenseVector<DenseConstVecView<T, TF>, TF>(*this), m_first(first),
      m_last(first + (n - 1) * stride), m_stride(stride)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseConstVecView<T, TF>::DenseConstVecView(const DenseConstVecView<T, TF> &slice_r)
{
  m_first  = slice_r.m_first;
  m_last   = slice_r.m_last;
  m_stride = slice_r.m_stride;
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseConstVecView<T, TF>::~DenseConstVecView()
{
  m_first  = 0;
  m_last   = 0;
  m_stride = 1;
}

// ----------------------------------------------------------------------------

// template<typename T, bool TF>
// inline T& DenseConstVecView<T,TF>::operator[] ( const Uint i )
//{
//    return m_first[i*m_stride];
//}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
inline const T &DenseConstVecView<T, TF>::operator[](const Uint i) const
{
  return m_first[i * m_stride];
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
inline Uint DenseConstVecView<T, TF>::size() const
{
  return (m_last - m_first) / m_stride + 1;
}

// ----------------------------------------------------------------------------

template <typename T, bool TF>
const DenseConstVecView<T, TF> DenseConstVecView<T, TF>::slice(const Uint idx_first,
                                                               const Uint length) const
{
  const Uint max_len = idx_first + length > size() ? size() - idx_first : length;
  DenseConstVecView<T, TF> result(m_first + idx_first * m_stride,
                                  m_first + (idx_first + max_len - 1) * m_stride, m_stride);
  return result;
}

// ----------------------------------------------------------------------------
//             Implementation of DenseConstVecView operators
// ----------------------------------------------------------------------------

template <typename T, bool TF>
DenseConstVecView<T, TF> &DenseConstVecView<T, TF>::operator=(
    const DenseConstVecView<T, TF> &slice_r)
{
  m_first  = slice_r.m_first;
  m_last   = slice_r.m_last;
  m_stride = slice_r.m_stride;
  return (*this);
}

// ----------------------------------------------------------------------------
// AddTrait specializations
// ----------------------------------------------------------------------------

template <typename T1, typename T2, bool TF>
struct AddTrait<DenseConstVecView<T1, TF>, DenseConstVecView<T2, TF>>
{
  using type = DenseDVec<typename AddTrait<T1, T2>::type, TF>;
};

// DenseDVec + DenseConstVecView
template <typename T1, typename T2, bool TF>
struct AddTrait<DenseDVec<T1, TF>, DenseConstVecView<T2, TF>>
{
  using type = DenseDVec<typename AddTrait<T1, T2>::type, TF>;
};

// DenseConstVecView + DenseDVec
template <typename T1, typename T2, bool TF>
struct AddTrait<DenseConstVecView<T1, TF>, DenseDVec<T2, TF>>
{
  using type = DenseDVec<typename AddTrait<T1, T2>::type, TF>;
};

// ----------------------------------------------------------------------------
// SubTrait specializations
// ----------------------------------------------------------------------------

template <typename T1, typename T2, bool TF>
struct SubTrait<DenseConstVecView<T1, TF>, DenseConstVecView<T2, TF>>
{
  using type = DenseDVec<typename SubTrait<T1, T2>::type, TF>;
};

template <typename T1, typename T2, bool TF>
struct SubTrait<DenseDVec<T1, TF>, DenseConstVecView<T2, TF>>
{
  using type = DenseDVec<typename SubTrait<T1, T2>::type, TF>;
};

template <typename T1, typename T2, bool TF>
struct SubTrait<DenseConstVecView<T1, TF>, DenseDVec<T2, TF>>
{
  using type = DenseDVec<typename SubTrait<T1, T2>::type, TF>;
};

template <typename T1, Uint N, bool TF, typename T2>
struct SubTrait<DenseSVec<T1, N, TF>, DenseConstVecView<T2, TF>>
{
  using type = DenseSVec<typename SubTrait<T1, T2>::type, N, TF>;
};

template <typename T1, bool TF, typename T2, Uint N>
struct SubTrait<DenseConstVecView<T1, TF>, DenseSVec<T2, N, TF>>
{
  using type = DenseSVec<typename SubTrait<T1, T2>::type, N, TF>;
};

// ----------------------------------------------------------------------------
// MultTrait specializations
// ----------------------------------------------------------------------------

// scalar * DenseConstVecView
template <typename SValueT, typename T, bool TF>
struct MultTrait<SValueT, DenseConstVecView<T, TF>>
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
struct TensorRank<DenseConstVecView<T, TF>>
{
  enum
  {
    value = tensor_rank_1
  };
};

// ----------------------------------------------------------------------------

template <typename T, bool TF>
std::ostream &operator<<(std::ostream &os, const DenseConstVecView<T, TF> &view)
{
  //   os.clear();
  for (Uint i = 0; (i + 1) < view.size(); ++i)
  {
    os << view[i] << " ";
  }
  if (view.size() > 0)
  {
    os << view[view.size() - 1];
  }

  return os;
}

// ----------------------------------------------------------------------------

} // Namespace math

} // Namespace pdekit

#endif
