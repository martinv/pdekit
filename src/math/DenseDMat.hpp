#ifndef PDEKIT_Math_Dense_Dynamic_Matrix_hpp
#define PDEKIT_Math_Dense_Dynamic_Matrix_hpp

#include <cmath>
#include <iostream>

#include "math/MathForward.hpp"
#include "math/MatrixStorageOrder.hpp"
#include "math/OperationEvalTime.hpp"
#include "math/TensorInitialization.hpp"
#include "math/TensorRank.hpp"
#include "math/binary_ops/DMatEvalExpr.hpp"
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
//             Forward declarations to enable friend functions
// ----------------------------------------------------------------------------

// template<typename T, bool SO> class DenseDMat;

template <typename T, bool SO>
std::ostream &operator<<(std::ostream &os, const DenseDMat<T, SO> &mat);

template <typename T, bool SO>
std::istream &operator>>(std::istream &, DenseDMat<T, SO> &mat);

// ----------------------------------------------------------------------------
//             Implementation of class Matrix
// ----------------------------------------------------------------------------

template <typename T, bool SO = DefaultMatrixStorageOrder>
class DenseDMat : public DenseMatrix<DenseDMat<T, SO>, SO>
{

  public:
  using tensor_type    = DenseDMat<T, SO>;
  using composite_type = DenseDMat<T, SO> const &;
  using value_type     = T;

  using row_slice_type    = DenseConstVecView<T, RowVector>;
  using column_slice_type = DenseConstVecView<T, ColumnVector>;

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
  explicit DenseDMat();

  /// Constructor, takes the number of rows and columns
  explicit DenseDMat(Uint m, Uint n);

  /// The copy constructor is explicitly defined due to the required
  /// dynamic memory management and in order to enable/facilitate NRV
  /// optimization.
  DenseDMat(const DenseDMat &rhs);

  template <typename MT>
  explicit DenseDMat(const DenseMatrix<MT, SO> &init);

  /// Intentionally disabled - DenseDMat cannot be constructed
  /// from initializer list, because the shape of the matrix would be
  /// undefined!
  /// explicit DenseDMat( const TensorInitializer<T>& init );

  /// Destructor
  ~DenseDMat();

  /// Resize the matrix
  void resize(const Uint m, const Uint n);

  /// Indexing operator
  /// @param the row and column index of the element to be returned
  inline T &operator()(const Uint i, const Uint j);

  /// Indexing operator, const version
  inline const T &operator()(const Uint i, const Uint j) const;

  template <typename MT2>
  inline void assign(const math::DenseMatrix<MT2, SO> &mat_rhs);

  /*
  /// Give access to the raw data:
  inline const T* const data() const
  {
    return m_data;
  }
  */

  /// Operator returning the number of the rows of the matrix
  inline Uint rows() const
  {
    return m_rows;
  }

  /// Operator returning the number of the columns of the matrix
  inline Uint cols() const
  {
    return m_cols;
  }

  /// Matrix assignement operator
  DenseDMat<T, SO> &operator=(const DenseDMat<T, SO> &mat_r);

  /// Fill the matrix with the same value
  void fill(const T val);

  /// Return one row (a row vector)
  /// @param  i index of the row
  /// @return the matrix row
  DenseConstVecView<T, RowVector> const_row(const Uint i) const;

  /// Return one row and transpose it so that it is a column vector
  /// @param  i index of the row
  /// @return the matrix row
  DenseConstVecView<T, ColumnVector> const_row_transp(const Uint i) const;

  /// Return one row (a row vector)
  /// @param  i index of the row
  /// @return the matrix row
  DenseVecView<T, RowVector> row(const Uint i);

  /// Return one row and transpose it so that it is a column vector
  /// @param  i index of the row
  /// @return the matrix row
  DenseVecView<T, ColumnVector> row_transp(const Uint i);

  /// Return one constant column
  /// @param  i index of the column
  /// @return the i-th column of the matrix
  DenseConstVecView<T, ColumnVector> const_col(const Uint i) const;

  /// Return one column
  /// @param  i index of the column
  /// @return the i-th column of the matrix
  DenseVecView<T, ColumnVector> col(const Uint i);

  /// Return a constant submatrix (block)
  DenseConstMatView<T, SO> const_block(const Uint row_start, const Uint col_start,
                                       const Uint nb_rows, const Uint nb_cols) const;

  /// Return a mutable submatrix (block)
  DenseMatView<T, SO> block(const Uint row_start, const Uint col_start, const Uint nb_rows,
                            const Uint nb_cols) const;

  /// Insert one row into the matrix
  /// @param row - the number of row which should be insterted
  /// @param row_vector - data to be inserted
  template <typename VectorType>
  void insert_row(const Uint const_row, const VectorType &row_vector);

  /// Insert one column into the matrix
  /// @param col - the number of column which should be inserted
  /// @param col_vector - column values to be inserted
  template <typename VectorType>
  void insert_col(const Uint const_col, const VectorType &col_vector);

  /// Overloaded output operator <<
  friend std::ostream &operator<<<T, SO>(std::ostream &os, const DenseDMat<T, SO> &mat);

  /// Overloaded input operator >>
  friend std::istream &operator>><T, SO>(std::istream &, DenseDMat<T, SO> &mat);

  /// Computation of determinant
  const T det() const;

  /// Compute matrix inverse
  template <typename OtherMatrixType>
  void inv(DenseMatrix<OtherMatrixType, SO> &mat) const;

  /// Transpose the matrix
  void transpose_in_place();

  // ----------------------------------------------------------------------------
  //             Definition of Matrix operators
  // ----------------------------------------------------------------------------

  /// Assignment operator from a matrix expression
  /// @param rhs - a matrix expression (for example, product of two matrices)
  template <typename MT>
  inline DenseDMat<T, SO> &operator=(const DenseMatrix<MT, SO> &rhs);

  /// Assign values from an initializing list
  DenseDMat<T, SO> &operator=(const TensorInitializer<T> &init);

  /// Accumulation from matrix expression
  template <typename MT>
  inline DenseDMat<T, SO> &operator+=(const DenseMatrix<MT, SO> &rhs);

  /// Subtraction of matrix expression
  template <typename MT>
  inline DenseDMat<T, SO> &operator-=(const DenseMatrix<MT, SO> &rhs);

  private:
  /// DATA:
  T *m_data;

  Uint m_rows;
  Uint m_cols;
};

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMat<T, SO>::DenseDMat()
    : DenseMatrix<DenseDMat<T, SO>, SO>(*this), m_data(nullptr), m_rows(0), m_cols(0)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMat<T, SO>::DenseDMat(Uint m, Uint n) : DenseMatrix<DenseDMat<T, SO>, SO>(*this)
{
  // resize(m,n);
  if (m * n != 0u)
  {
    m_data = new T[m * n];
    m_rows = m;
    m_cols = n;
  }
  else
  {
    m_data = nullptr;
    m_rows = 0u;
    m_cols = 0u;
  }
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMat<T, SO>::DenseDMat(const DenseDMat &rhs)
{
  m_data = new T[rhs.rows() * rhs.cols()];
  m_rows = rhs.rows();
  m_cols = rhs.cols();

  for (Uint i = 0; i < m_rows * m_cols; ++i)
  {
    m_data[i] = rhs.m_data[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
template <typename MT>
DenseDMat<T, SO>::DenseDMat(const DenseMatrix<MT, SO> &init)
{
  const MT &rhs_expr = init.wrapped_type();

  // resize(rhs_expr.rows(),rhs_expr.cols());

  m_data = new T[rhs_expr.rows() * rhs_expr.cols()];
  m_rows = rhs_expr.rows();
  m_cols = rhs_expr.cols();

  free_assign(*this, init.wrapped_type());
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMat<T, SO>::~DenseDMat()
{
  delete[] m_data;
  m_data = nullptr;
  m_rows = 0;
  m_cols = 0;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void DenseDMat<T, SO>::resize(const Uint m, const Uint n)
{
  if (m * n == 0)
  {
    if (m_data)
    {
      delete[] m_data;
      m_data = nullptr;
    }
    m_rows = 0;
    m_cols = 0;
    return;
  }

  if (m_rows * m_cols != m * n)
  {
    if (m_data != nullptr)
    {
      delete[] m_data;
    }
    m_data = new T[m * n];
  }
  m_rows = m;
  m_cols = n;

  fill(T());
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
inline T &DenseDMat<T, SO>::operator()(const Uint i, const Uint j)
{
  return m_data[i * m_cols + j];
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
inline const T &DenseDMat<T, SO>::operator()(const Uint i, const Uint j) const
{
  return m_data[i * m_cols + j];
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
template <typename MT2>
void DenseDMat<T, SO>::assign(const math::DenseMatrix<MT2, SO> &mat_rhs)
{
  MT2 const &wrapped_mat_rhs = mat_rhs.wrapped_type();

  for (Uint i = 0; i < m_rows; ++i)
  {
    for (Uint j = 0; j < m_cols; ++j)
    {
      m_data[i * m_cols + j] = wrapped_mat_rhs(i, j);
    }
  }
}

// ----------------------------------------------------------------------------
// NOTE THAT THE SHAPE OF THE MATRIX TO WHICH WE ASSIGN IS NOT CHANGED!

template <typename T, bool SO>
DenseDMat<T, SO> &DenseDMat<T, SO>::operator=(const DenseDMat<T, SO> &mat_r)
{
  for (Uint i = 0; i < m_rows * m_cols; i++)
  {
    m_data[i] = mat_r.m_data[i];
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void DenseDMat<T, SO>::fill(const T val)
{
  for (Uint i = 0; i < m_rows * m_cols; ++i)
  {
    m_data[i] = val;
  }
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstVecView<T, RowVector> DenseDMat<T, SO>::const_row(const Uint i) const
{
  return DenseConstVecView<T, RowVector>(m_data + i * m_cols, m_cols);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstVecView<T, ColumnVector> DenseDMat<T, SO>::const_row_transp(const Uint i) const
{
  return DenseConstVecView<T, ColumnVector>(m_data + i * m_cols, m_cols);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseVecView<T, RowVector> DenseDMat<T, SO>::row(const Uint i)
{
  return DenseVecView<T, RowVector>(m_data + i * m_cols, m_cols);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseVecView<T, ColumnVector> DenseDMat<T, SO>::row_transp(const Uint i)
{
  return DenseVecView<T, ColumnVector>(m_data + i * m_cols, m_cols);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstVecView<T, ColumnVector> DenseDMat<T, SO>::const_col(const Uint i) const
{
  return DenseConstVecView<T, ColumnVector>(m_data + i, m_rows, m_cols);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseVecView<T, ColumnVector> DenseDMat<T, SO>::col(const Uint i)
{
  return DenseVecView<T, ColumnVector>(m_data + i, m_rows, m_cols);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstMatView<T, SO> DenseDMat<T, SO>::const_block(const Uint row_start, const Uint col_start,
                                                       const Uint nb_rows, const Uint nb_cols) const
{
  DenseConstMatView<T, SO> block(m_data + row_start * m_cols + col_start, m_cols, nb_rows, nb_cols);
  return block;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseMatView<T, SO> DenseDMat<T, SO>::block(const Uint row_start, const Uint col_start,
                                            const Uint nb_rows, const Uint nb_cols) const
{
  DenseMatView<T, SO> block(m_data + row_start * m_cols + col_start, m_cols, nb_rows, nb_cols);
  return block;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
template <typename VectorType>
void DenseDMat<T, SO>::insert_row(const Uint row, const VectorType &row_vector)
{
  const Uint offset = m_cols * row;

  for (Uint i = 0; i < m_cols; ++i)
  {
    m_data[offset + i] = row_vector[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
template <typename VectorType>
void DenseDMat<T, SO>::insert_col(const Uint col, const VectorType &col_vector)
{
  for (Uint i = 0; i < m_rows; ++i)
  {
    m_data[col + i * m_cols] = col_vector[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
std::ostream &operator<<(std::ostream &os, const DenseDMat<T, SO> &mat)
{
  //   os.clear();

  for (Uint i = 0; i < mat.m_rows; ++i)
  {
    for (Uint j = 0; j < mat.m_cols - 1; ++j)
    {
      os << mat(i, j) << " ";
    }
    os << mat(i, mat.m_cols - 1) << std::endl;
  }

  return os;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
std::istream &operator>>(std::istream &is, DenseDMat<T, SO> &mat)
{
  for (Uint i = 0; i < mat.size(); ++i)
    is >> mat[i];
  return is;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
const T DenseDMat<T, SO>::det() const
{
  return DynamicMatrixDeterminant<T>::compute(m_data, m_rows);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
template <typename OtherMatrixType>
void DenseDMat<T, SO>::inv(DenseMatrix<OtherMatrixType, SO> &mat) const
{
  DynamicMatrixInverter<T>::invert(m_data, mat.wrapped_type().m_data, m_rows);
  return;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void DenseDMat<T, SO>::transpose_in_place()
{
  if (m_rows == m_cols) // The matrix is square
  {
    for (Uint i = 0; i < (m_rows - 1); ++i)
      for (Uint j = i + 1; j < m_cols; ++j)
      {
        // Swap the values on position (i,j) and (j,i)
        std::swap(m_data[i * m_cols + j], m_data[j * m_cols + i]);
      }
  }
  else
  {
    std::vector<T> tmp(m_rows * m_cols, T());
    Uint pos = 0;
    for (Uint j = 0; j < m_cols; ++j)
    {
      for (Uint i = 0; i < m_rows; ++i)
      {
        tmp[pos++] = m_data[i * m_cols + j];
      }
    }
    for (Uint i = 0; i < m_rows * m_cols; ++i)
    {
      m_data[i] = tmp[i];
    }
    std::swap(m_rows, m_cols);
  }
}

// ----------------------------------------------------------------------------
//             Implementation of Matrix operators
// ----------------------------------------------------------------------------

template <typename T, bool SO>
template <typename MT>
DenseDMat<T, SO> &DenseDMat<T, SO>::operator=(const DenseMatrix<MT, SO> &rhs)
{
  const MT &rhs_expr = rhs.wrapped_type();

  resize(rhs_expr.rows(), rhs_expr.cols());

  free_assign(*this, rhs.wrapped_type());

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMat<T, SO> &DenseDMat<T, SO>::operator=(const TensorInitializer<T> &init)
{
  for (Uint i = 0; i < std::min(m_rows * m_cols, init.size()); ++i)
  {
    m_data[i] = init[i];
  }
  for (Uint i = std::min(m_rows * m_cols, init.size()); i < m_rows * m_cols; ++i)
  {
    m_data[i] = T();
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
template <typename MT>
DenseDMat<T, SO> &DenseDMat<T, SO>::operator+=(const DenseMatrix<MT, SO> &rhs)
{
  const MT &rhs_expr = rhs.wrapped_type();

  // Do not use 'resize', because 'resize' sets all values to 0, but we want
  // to ACCUMULATE here - all previously accumulated values would be
  // destroyed! resize(rhs_expr.rows(), rhs_expr.cols());

  for (Uint i = 0; i < m_rows; ++i)
  {
    for (Uint j = 0; j < m_cols; ++j)
    {
      m_data[i * m_cols + j] += rhs_expr(i, j);
    }
  }
  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
template <typename MT>
DenseDMat<T, SO> &DenseDMat<T, SO>::operator-=(const DenseMatrix<MT, SO> &rhs)
{
  const MT &rhs_expr = rhs.wrapped_type();

  // Do not use 'resize', because 'resize' sets all values to 0, but we want
  // to ACCUMULATE here - all previously accumulated values would be
  // destroyed! resize(rhs_expr.rows(), rhs_expr.cols());

  for (Uint i = 0; i < m_rows; ++i)
  {
    for (Uint j = 0; j < m_cols; ++j)
    {
      m_data[i * m_cols + j] -= rhs_expr(i, j);
    }
  }
  return (*this);
}

// ----------------------------------------------------------------------------
// AddTrait specializations
// ----------------------------------------------------------------------------

// DenseDMat + DenseDMat

// T1 - value type in first DenseDMat
// T2 - value type in second DenseDMat
// M0S0 - size in dimension '0' in DenseDMat 0 ( = nb rows in M0 )
// M0S1 - number of columns in M0
// M1S0 - number of rows in M1
// M1S1 - number of columns in M1

template <typename T1, typename T2, bool SO>
struct AddTrait<DenseDMat<T1, SO>, DenseDMat<T2, SO>>
{
  using type = DenseDMat<typename AddTrait<T1, T2>::type, SO>;
};

// ----------------------------------------------------------------------------
// SubTrait specializations
// ----------------------------------------------------------------------------

// DenseDMat - DenseDMat

// T1 - value type in first DenseDMat
// T2 - value type in second DenseDMat
// M0S0 - size in dimension '0' in DenseDMat 0 ( = nb rows in M0 )
// M0S1 - number of columns in M0
// M1S0 - number of rows in M1
// M1S1 - number of columns in M1

template <typename T1, typename T2, bool SO>
struct SubTrait<DenseDMat<T1, SO>, DenseDMat<T2, SO>>
{
  using type = DenseDMat<typename SubTrait<T1, T2>::type, SO>;
};

// ----------------------------------------------------------------------------
// MultTrait specializations
// ----------------------------------------------------------------------------

// scalar * DenseDMat
template <typename SValueT, typename T, bool SO>
struct MultTrait<SValueT, DenseDMat<T, SO>>
{
  using type = DenseDMat<typename MultTrait<SValueT, T>::type, SO>;
};

// DenseDMat * DenseDVec
template <typename T1, typename T2, bool SO, bool TF>
struct MultTrait<DenseDMat<T1, SO>, DenseDVec<T2, TF>>
{
  using type = DenseDVec<typename MultTrait<T1, T2>::type, TF>;
};

// DenseDMat * DenseSVec
template <typename T1, bool SO, typename T2, Uint N, bool TF>
struct MultTrait<DenseDMat<T1, SO>, DenseSVec<T2, N, TF>>
{
  using type = DenseDVec<typename MultTrait<T1, T2>::type, TF>;
};

// DenseDMat * DenseConstVecView
template <typename T1, typename T2, bool SO, bool TF>
struct MultTrait<DenseDMat<T1, SO>, DenseConstVecView<T2, TF>>
{
  using type = DenseDVec<typename MultTrait<T1, T2>::type, ColumnVector>;
};

// DenseDMat * DenseDMat

// T1 - value type in first DenseDMat
// T2 - value type in second DenseDMat
// M0S0 - size in dimension '0' in DenseDMat 0 ( = nb rows in M0 )
// M0S1 - number of columns in M0
// M1S0 - number of rows in M1
// M1S1 - number of columns in M1

template <typename T1, typename T2, bool SO>
struct MultTrait<DenseDMat<T1, SO>, DenseDMat<T2, SO>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO>;
};

// DenseDMat * DenseSMat

// T1  - value type in DenseDMat
// SO1 - storage order in DenseDMat
// T2  - value type in DenseSMat
// M   - number of rows in DenseSMat
// N   - number of columns in DenseSMat
// SO2 - storage order in DenseSMat

template <typename T1, bool SO1, typename T2, Uint M, Uint N, bool SO2>
struct MultTrait<DenseDMat<T1, SO1>, DenseSMat<T2, M, N, SO2>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO1>;
};

// DenseSMat * DenseDMat

// T1  - value type in DenseSMat
// M   - number of rows in DenseSMat
// N   - number of columns in DenseSMat
// SO1 - storage order in DenseSMat
// T2  - value type in DenseDMat
// SO2 - storage order in DenseDMat

template <typename T1, Uint M, Uint N, bool SO1, typename T2, bool SO2>
struct MultTrait<DenseSMat<T1, M, N, SO1>, DenseDMat<T2, SO2>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO1>;
};

// DenseDMat * DenseConstMatView

// T1  - value type in DenseDMat
// SO1 - storage order in DenseDMat
// T2  - value type in DenseConstMatView
// SO2 - storage order in DenseConstMatView

template <typename T1, bool SO1, typename T2, bool SO2>
struct MultTrait<DenseDMat<T1, SO1>, DenseConstMatView<T2, SO2>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO1>;
};

// DenseConstMatView * DenseDMat

// T1  - value type in DenseConstMatView
// SO1 - storage order in DenseConstMatView
// T2  - value type in DenseDMat
// SO2 - storage order in DenseDMat

template <typename T1, bool SO1, typename T2, bool SO2>
struct MultTrait<DenseConstMatView<T1, SO1>, DenseDMat<T2, SO2>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO1>;
};

// DenseDMat * DenseMatView

// T1  - value type in DenseDMat
// SO1 - storage order in DenseDMat
// T2  - value type in DenseMatView
// SO2 - storage order in DenseMatView

template <typename T1, bool SO1, typename T2, bool SO2>
struct MultTrait<DenseDMat<T1, SO1>, DenseMatView<T2, SO2>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO1>;
};

// DenseMatView * DenseDMat

// T1  - value type in DenseMatView
// SO1 - storage order in DenseMatView
// T2  - value type in DenseDMat
// SO2 - storage order in DenseDMat

template <typename T1, bool SO1, typename T2, bool SO2>
struct MultTrait<DenseMatView<T1, SO1>, DenseDMat<T2, SO2>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO1>;
};

// ----------------------------------------------------------------------------
// DivTrait specializations
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TENSOR RANK OF MATRIX = 2
// ----------------------------------------------------------------------------

template <typename T, bool SO>
struct TensorRank<DenseDMat<T, SO>>
{
  enum
  {
    value = tensor_rank_2
  };
};

} // Namespace math

} // Namespace pdekit

#endif // PDEKIT_Math_Matrix_hpp
