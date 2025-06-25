#ifndef PDEKIT_Math_Dense_Const_Mat_View_hpp
#define PDEKIT_Math_Dense_Const_Mat_View_hpp

#include <cmath>
#include <iostream>

#include "math/MathForward.hpp"
#include "math/MatrixStorageOrder.hpp"
#include "math/OperationEvalTime.hpp"
#include "math/TensorInitialization.hpp"
#include "math/TensorRank.hpp"
#include "math/binary_ops/MatrixOps.hpp"
#include "math/traits/AddTrait.hpp"
#include "math/traits/DivTrait.hpp"
#include "math/traits/MultTrait.hpp"
#include "math/traits/SubTrait.hpp"
#include "math/unary_ops/MatrixUnaryOps.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------
//    Forward declarations to enable friend functions are in MathForward.hpp
// ----------------------------------------------------------------------------

template <typename T, bool SO>
std::ostream &operator<<(std::ostream &, const DenseConstMatView<T, SO> &view);

// ----------------------------------------------------------------------------
//             Implementation of class MatrixBlock
// ----------------------------------------------------------------------------

template <typename T, bool SO = DefaultMatrixStorageOrder>
class DenseConstMatView : public DenseMatrix<DenseConstMatView<T, SO>, SO>
{

  public:
  using tensor_type    = DenseConstMatView<T, SO>;
  using composite_type = DenseConstMatView<T, SO> const &;
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
  explicit DenseConstMatView();

  /// Constructor
  explicit DenseConstMatView(const T *first, const Uint source_ld, const Uint m, const Uint n);

  /// Copy constructor
  DenseConstMatView(const DenseConstMatView<T, SO> &other_block);

  /// Destructor
  ~DenseConstMatView();

  /// Indexing operator
  // inline T& operator() (const Uint i, const Uint j);

  /// Indexing operator, const version
  const T &operator()(const Uint i, const Uint j) const;

  /// Return just one row of the matrix
  DenseConstVecView<T, RowVector> row(const Uint i) const;

  /// Return one column of the matrix
  DenseConstVecView<T, ColumnVector> col(const Uint i) const;

  /// Return the transpose of one row of the matrix, i.e. a column vector
  DenseConstVecView<T, ColumnVector> row_transpose(const Uint i) const;

  /// Operator returning the number of rows
  Uint rows() const;

  /// Operator returning the number of columns
  Uint cols() const;

  /// Overloaded output operator <<
  friend std::ostream &operator<<<T>(std::ostream &, const DenseConstMatView<T, SO> &view);

  // --------------------------------------------------------------------------
  //                      DenseConstMatView operators
  // --------------------------------------------------------------------------

  /// ConstMatrixBlock assignement operator
  DenseConstMatView &operator=(const DenseConstMatView<T, SO> &other_block);

  private:
  /// Pointer to raw data (first entry)
  T const *m_first; // Pointer to a constant object

  /// Length of row/column (depending on storage type) of the
  /// source matrix from which we're taking this block
  Uint m_source_ld;

  /// Number of rows
  Uint m_rows;

  /// Number of columns
  Uint m_cols;
};

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstMatView<T, SO>::DenseConstMatView()
    : DenseMatrix<DenseConstMatView<T, SO>, SO>(*this), m_first(nullptr), m_source_ld(0), m_rows(0),
      m_cols(0)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstMatView<T, SO>::DenseConstMatView(const T *first, const Uint source_ld, const Uint m,
                                            const Uint n)
    : DenseMatrix<DenseConstMatView<T, SO>, SO>(*this), m_first(first), m_source_ld(source_ld),
      m_rows(m), m_cols(n)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstMatView<T, SO>::DenseConstMatView(const DenseConstMatView<T, SO> &other_block)
{
  m_first     = other_block.m_first;
  m_source_ld = other_block.m_source_ld;
  m_rows      = other_block.m_rows;
  m_cols      = other_block.m_cols;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstMatView<T, SO>::~DenseConstMatView()
{
  m_first     = nullptr;
  m_source_ld = 0;
  m_rows      = 0;
  m_cols      = 0;
}

// ----------------------------------------------------------------------------

// template<typename T,bool SO>
// T& MatrixBlock<T,SO>::operator() ( const Uint i, const Uint j )
//{
//  return m_first[i*m_source_ld+j];
//}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
inline const T &DenseConstMatView<T, SO>::operator()(const Uint i, const Uint j) const
{
  return m_first[i * m_source_ld + j];
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstVecView<T, RowVector> DenseConstMatView<T, SO>::row(const Uint i) const
{
  //  return VectorView<T>(m_first+i*m_cols,m_first+(i+1)*m_cols-1);
  return DenseConstVecView<T, RowVector>(m_first + i * m_source_ld, m_cols);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstVecView<T, ColumnVector> DenseConstMatView<T, SO>::col(const Uint i) const
{
  DenseConstVecView<T, ColumnVector> column(m_first + i, m_rows, m_source_ld);
  return column;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstVecView<T, ColumnVector> DenseConstMatView<T, SO>::row_transpose(const Uint i) const
{
  //  return VectorView<T>(m_first+i*m_cols,m_first+(i+1)*m_cols-1);
  return DenseConstVecView<T, ColumnVector>(m_first + i * m_source_ld, m_cols);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
inline Uint DenseConstMatView<T, SO>::rows() const
{
  return m_rows;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
inline Uint DenseConstMatView<T, SO>::cols() const
{
  return m_cols;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
std::ostream &operator<<(std::ostream &os, const DenseConstMatView<T, SO> &view)
{
  for (Uint i = 0; i < view.m_rows; ++i)
  {
    for (Uint j = 0; j < view.m_cols - 1; ++j)
    {
      os << view.m_first[i * view.m_source_ld + j] << " ";
    }
    os << view.m_first[i * view.m_source_ld + view.m_cols - 1] << std::endl;
  }
  return os;
}

// ----------------------------------------------------------------------------
//         Implementation of ConstMatrixBlock operators
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstMatView<T, SO> &DenseConstMatView<T, SO>::operator=(
    const DenseConstMatView<T, SO> &other_block)
{
  m_first     = other_block.m_first;
  m_source_ld = other_block.m_source_ld;
  m_rows      = other_block.m_rows;
  m_cols      = other_block.m_cols;
  return (*this);
}

// ----------------------------------------------------------------------------
// AddTrait specializations
// ----------------------------------------------------------------------------

// MatrixBlock + MatrixBlock

// T1 - value type in first MatrixBlock
// T2 - value type in second MatrixBlock
// SO - storage order (common to both matrix blocks)

template <typename T1, typename T2, bool SO>
struct AddTrait<DenseConstMatView<T1, SO>, DenseConstMatView<T2, SO>>
{
  using type = DenseDMat<typename AddTrait<T1, T2>::type, SO>;
};

// ----------------------------------------------------------------------------
// SubTrait specializations
// ----------------------------------------------------------------------------

// MatrixBlock - MatrixBlock

// T1 - value type in first MatrixBlock
// T2 - value type in second MatrixBlock
// SO - storage order (common to both matrix blocks)

template <typename T1, typename T2, bool SO>
struct SubTrait<DenseConstMatView<T1, SO>, DenseConstMatView<T2, SO>>
{
  using type = DenseDMat<typename SubTrait<T1, T2>::type, SO>;
};

// ----------------------------------------------------------------------------
// MultTrait specializations
// ----------------------------------------------------------------------------

// scalar * MatrixBlock
template <typename SValueT, typename T, bool SO>
struct MultTrait<SValueT, DenseConstMatView<T, SO>>
{
  using type = DenseDMat<typename MultTrait<SValueT, T>::type, SO>;
};

// MatrixBlock * DenseDVec
template <typename T1, typename T2, bool SO, bool TF>
struct MultTrait<DenseConstMatView<T1, SO>, DenseDVec<T2, TF>>
{
  using type = DenseDVec<typename MultTrait<T1, T2>::type, TF>;
};

// MatrixBlock * DenseSVec
template <typename T1, bool SO, typename T2, Uint N, bool TF>
struct MultTrait<DenseConstMatView<T1, SO>, DenseSVec<T2, N, TF>>
{
  using type = DenseDVec<typename MultTrait<T1, T2>::type, TF>;
};

// MatrixBlock * DenseConstVecView
template <typename T1, typename T2, bool SO, bool TF>
struct MultTrait<DenseConstMatView<T1, SO>, DenseConstVecView<T2, TF>>
{
  using type = DenseDVec<typename MultTrait<T1, T2>::type, ColumnVector>;
};

// DenseConstMatView * DenseConstMatView

// T1 - value type in first MatrixBlock
// T2 - value type in second MatrixBlock
// S0 - storage order in both matrix blocks

template <typename T1, typename T2, bool SO>
struct MultTrait<DenseConstMatView<T1, SO>, DenseConstMatView<T2, SO>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO>;
};

// DenseConstMatView * DenseMatView

// T1 - value type in first MatrixBlock
// T2 - value type in second MatrixBlock
// S0 - storage order in both matrix blocks

template <typename T1, typename T2, bool SO>
struct MultTrait<DenseConstMatView<T1, SO>, DenseMatView<T2, SO>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO>;
};

// DenseMatView * DenseConstMatView

// T1 - value type in first MatrixBlock
// T2 - value type in second MatrixBlock
// S0 - storage order in both matrix blocks

template <typename T1, typename T2, bool SO>
struct MultTrait<DenseMatView<T1, SO>, DenseConstMatView<T2, SO>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO>;
};

// ----------------------------------------------------------------------------
// DivTrait specializations
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TENSOR RANK OF MATRIX = 2
// ----------------------------------------------------------------------------

template <typename T, bool SO>
struct TensorRank<DenseConstMatView<T, SO>>
{
  enum
  {
    value = tensor_rank_2
  };
};

} // Namespace math

} // Namespace pdekit

#endif
