#ifndef PDEKIT_Math_Dense_Static_Matrix_hpp
#define PDEKIT_Math_Dense_Static_Matrix_hpp

#include <cmath>
#include <iostream>

#include "common/Meta.hpp"

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
//    Forward declarations to enable friend functions are in MathForward.hpp
// ----------------------------------------------------------------------------

// template<typename T, Uint M, Uint N, bool SO> class DenseSMat;

template <typename T, Uint M, Uint N, bool SO>
std::ostream &operator<<(std::ostream &os, const DenseSMat<T, M, N, SO> &mat);

template <typename T, Uint M, Uint N, bool SO>
std::istream &operator>>(std::istream &, DenseSMat<T, M, N, SO> &mat);

// ----------------------------------------------------------------------------
//             Implementation of class Matrix
// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO = DefaultMatrixStorageOrder>
class DenseSMat : public DenseMatrix<DenseSMat<T, M, N, SO>, SO>
{

  public:
  using tensor_type    = DenseSMat<T, M, N, SO>;
  using composite_type = DenseSMat<T, M, N, SO> const &;
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
  explicit DenseSMat();

  // /// Copy constructor
  // DenseSMat ( const Matrix<T,SO>& rhs );

  template <typename MT>
  explicit DenseSMat(const DenseMatrix<MT, SO> &init);

  /// Construct from tensor initializer
  /// This is allowed ONLY FOR STATIC MATRICES, because for
  /// dynamic matrices, we don't know what the shape of the matrix should be!
  DenseSMat(const TensorInitializer<T> &init);

  /// Destructor
  ~DenseSMat();

  /// Indexing operator
  /// @param the row and column index of the element to be returned
  ///        This is different from the square bracket operator!
  inline T &operator()(const Uint i, const Uint j);

  /// Indexing operator, const version
  inline const T &operator()(const Uint i, const Uint j) const;

  template <typename MT2>
  inline void assign(const math::DenseMatrix<MT2, SO> &mat_rhs);

  /// Operator returning the number of the rows of the matrix
  inline Uint rows() const
  {
    return M;
  }

  /// Operator returning the number of the columns of the matrix
  inline Uint cols() const
  {
    return N;
  }

  /// Matrix assignement operator
  DenseSMat<T, M, N, SO> &operator=(const DenseSMat<T, M, N, SO> &mat_r);

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
  DenseVecView<T, ColumnVector> row_transp(const Uint i) const;

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
  friend std::ostream &operator<<<T, M, N, SO>(std::ostream &os, const DenseSMat<T, M, N, SO> &mat);

  /// Overloaded input operator >>
  friend std::istream &operator>><T, M, N, SO>(std::istream &, DenseSMat<T, M, N, SO> &mat);

  /// Computation of determinant
  const T det() const;

  /// Compute matrix inverse
  template <typename OtherMatrixType>
  void inv(DenseMatrix<OtherMatrixType, SO> &mat) const;

  /// Transpose the matrix in place
  void transpose_in_place();

  // ----------------------------------------------------------------------------
  //             Definition of Matrix operators
  // ----------------------------------------------------------------------------

  /// Assignment operator from a matrix expression
  /// @param rhs - a matrix expression (for example, product of two matrices)
  template <typename MT>
  inline DenseSMat<T, M, N, SO> &operator=(const DenseMatrix<MT, SO> &rhs);

  /// Assign the matrix values from an initializing list
  DenseSMat<T, M, N, SO> &operator=(const TensorInitializer<T> &init);

  /// Accumulation from matrix expression
  template <typename MT>
  inline DenseSMat<T, M, N, SO> &operator+=(const DenseMatrix<MT, SO> &rhs);

  /// Subtraction of matrix expression
  template <typename MT>
  inline DenseSMat<T, M, N, SO> &operator-=(const DenseMatrix<MT, SO> &rhs);

  private:
  /// DATA:
  T m_data[M * N];
};

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseSMat<T, M, N, SO>::DenseSMat() : DenseMatrix<DenseSMat<T, M, N, SO>, SO>(*this)
{
  for (Uint i = 0; i < M * N; ++i)
  {
    m_data[i] = T();
  }
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseSMat<T, M, N, SO>::DenseSMat(const TensorInitializer<T> &init)
{
  for (Uint i = 0; i < std::min(M * N, init.size()); ++i)
  {
    m_data[i] = init[i];
  }
  for (Uint i = std::min(M * N, init.size()); i < M * N; ++i)
  {
    m_data[i] = T();
  }
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
template <typename MT>
DenseSMat<T, M, N, SO>::DenseSMat(const DenseMatrix<MT, SO> &init)
{
  /*
  const MT &rhs_expr = init.wrapped_type();

  for(Uint i = 0; i < M*N; ++i)
  {
    m_data[i] = rhs_expr[i];
  }
  */

  free_assign(*this, init.wrapped_type());
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseSMat<T, M, N, SO>::~DenseSMat()
{
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
inline T &DenseSMat<T, M, N, SO>::operator()(const Uint i, const Uint j)
{
  return m_data[i * N + j];
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
inline const T &DenseSMat<T, M, N, SO>::operator()(const Uint i, const Uint j) const
{
  return m_data[i * N + j];
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
template <typename MT2>
void DenseSMat<T, M, N, SO>::assign(const math::DenseMatrix<MT2, SO> &mat_rhs)
{
  MT2 const &wrapped_mat_rhs = mat_rhs.wrapped_type();

  for (Uint i = 0; i < M; ++i)
  {
    for (Uint j = 0; j < M; ++j)
    {
      m_data[i * N + j] = wrapped_mat_rhs(i, j);
    }
  }
}

// ----------------------------------------------------------------------------
// NOTE THAT THE SHAPE OF THE MATRIX TO WHICH WE ASSIGN IS NOT CHANGED!

template <typename T, Uint M, Uint N, bool SO>
DenseSMat<T, M, N, SO> &DenseSMat<T, M, N, SO>::operator=(const DenseSMat<T, M, N, SO> &mat_r)
{
  for (Uint i = 0; i < M * N; ++i)
  {
    m_data[i] = mat_r.m_data[i];
  }
  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
void DenseSMat<T, M, N, SO>::fill(const T val)
{
  for (Uint i = 0; i < M * N; ++i)
  {
    m_data[i] = val;
  }
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseConstVecView<T, RowVector> DenseSMat<T, M, N, SO>::const_row(const Uint i) const
{
  return DenseConstVecView<T, RowVector>(m_data + i * N, N);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseConstVecView<T, ColumnVector> DenseSMat<T, M, N, SO>::const_row_transp(const Uint i) const
{
  return DenseConstVecView<T, ColumnVector>(m_data + i * N, N);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseVecView<T, RowVector> DenseSMat<T, M, N, SO>::row(const Uint i)
{
  return DenseVecView<T, RowVector>(m_data + i * N, N);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseVecView<T, ColumnVector> DenseSMat<T, M, N, SO>::row_transp(const Uint i) const
{
  return DenseVecView<T, ColumnVector>(m_data + i * N, N);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseConstVecView<T, ColumnVector> DenseSMat<T, M, N, SO>::const_col(const Uint i) const
{
  return DenseConstVecView<T, ColumnVector>(m_data + i, M, N);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseVecView<T, ColumnVector> DenseSMat<T, M, N, SO>::col(const Uint i)
{
  return DenseVecView<T, ColumnVector>(m_data + i, M, N);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseConstMatView<T, SO> DenseSMat<T, M, N, SO>::const_block(const Uint row_start,
                                                             const Uint col_start,
                                                             const Uint nb_rows,
                                                             const Uint nb_cols) const
{
  DenseConstMatView<T, SO> block(m_data + row_start * N + col_start, N, nb_rows, nb_cols);
  return block;
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseMatView<T, SO> DenseSMat<T, M, N, SO>::block(const Uint row_start, const Uint col_start,
                                                  const Uint nb_rows, const Uint nb_cols) const
{
  DenseMatView<T, SO> block(m_data + row_start * N + col_start, N, nb_rows, nb_cols);
  return block;
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
template <typename VectorType>
void DenseSMat<T, M, N, SO>::insert_row(const Uint row, const VectorType &row_vector)
{
  const Uint offset = N * row;

  for (Uint i = 0; i < N; ++i)
  {
    m_data[offset + i] = row_vector[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
template <typename VectorType>
void DenseSMat<T, M, N, SO>::insert_col(const Uint col, const VectorType &col_vector)
{
  for (Uint i = 0; i < M; ++i)
  {
    m_data[col + i * N] = col_vector[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
std::ostream &operator<<(std::ostream &os, const DenseSMat<T, M, N, SO> &mat)
{
  //   os.clear();

  for (Uint i = 0; i < M; ++i)
  {
    for (Uint j = 0; j < N - 1; ++j)
    {
      os << mat(i, j) << " ";
    }
    os << mat(i, N - 1) << std::endl;
  }
  return os;
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
std::istream &operator>>(std::istream &is, DenseSMat<T, M, N, SO> &mat)
{
  for (Uint i = 0; i < M * N; ++i)
    is >> mat[i];
  return is;
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
const T DenseSMat<T, M, N, SO>::det() const
{
  using DetType = typename common::SelectType<(M < 5), MatrixDeterminant<T, M>,
                                              DynamicMatrixDeterminant<T>>::type;
  return DetType::compute(m_data, M);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
template <typename OtherMatrixType>
void DenseSMat<T, M, N, SO>::inv(DenseMatrix<OtherMatrixType, SO> &mat) const
{
  using inverter_type =
      typename common::SelectType<(M < 5), DirectMatrixInverse<T, M>, LapackMatrixInverse<T>>::type;
  inverter_type::compute(m_data, mat.wrapped_type().m_data, M);
  return;
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
void DenseSMat<T, M, N, SO>::transpose_in_place()
{
  if (M == N) // The matrix is square
  {
    for (Uint i = 0; i < (M - 1); ++i)
      for (Uint j = i + 1; j < N; ++j)
      {
        // Swap the values on position (i,j) and (j,i)
        std::swap(m_data[i * N + j], m_data[j * N + i]);
      }
  }
  else
  {
    std::cerr << "Error:Cannot transpose DenseSMat[" << M << "," << N
              << "] in place because it is not square." << std::endl;
  }
}

// ----------------------------------------------------------------------------
//             Implementation of Matrix operators
// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
template <typename MT>
DenseSMat<T, M, N, SO> &DenseSMat<T, M, N, SO>::operator=(const DenseMatrix<MT, SO> &rhs)
{
  free_assign(*this, rhs.wrapped_type());
  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
DenseSMat<T, M, N, SO> &DenseSMat<T, M, N, SO>::operator=(const TensorInitializer<T> &init)
{
  for (Uint i = 0; i < std::min(M * N, init.size()); ++i)
  {
    m_data[i] = init[i];
  }
  for (Uint i = std::min(M * N, init.size()); i < (M * N); ++i)
  {
    m_data[i] = T();
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
template <typename MT>
DenseSMat<T, M, N, SO> &DenseSMat<T, M, N, SO>::operator+=(const DenseMatrix<MT, SO> &rhs)
{
  const MT &rhs_expr = rhs.wrapped_type();

  // resize(rhs_expr.rows(),rhs_expr.cols());

  for (Uint i = 0; i < M; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      m_data[i * N + j] += rhs_expr(i, j);
    }
  }
  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
template <typename MT>
DenseSMat<T, M, N, SO> &DenseSMat<T, M, N, SO>::operator-=(const DenseMatrix<MT, SO> &rhs)
{
  const MT &rhs_expr = rhs.wrapped_type();

  // resize(rhs_expr.rows(),rhs_expr.cols());

  for (Uint i = 0; i < M; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      m_data[i * N + j] -= rhs_expr(i, j);
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

// DenseSMat + DenseSMat

/// Binary operation OP(DenseSMat,DenseSMat) traits:
// T1 - value type in first DenseDMat
// T2 - value type in second DenseDMat
// M0S0 - size in dimension '0' in DenseDMat 0 ( = nb rows in M0 )
// M0S1 - number of columns in M0
// M1S0 - number of rows in M1
// M1S1 - number of columns in M1

template <typename T1, typename T2, Uint N, bool SO>
struct AddTrait<DenseSMat<T1, N, N, SO>, DenseSMat<T2, N, N, SO>>
{
  using type = DenseSMat<typename AddTrait<T1, T2>::type, N, N, SO>;
};

// ----------------------------------------------------------------------------
// SubTrait specializations
// ----------------------------------------------------------------------------

// DenseSMat - DenseSMat

/// Binary operation OP(DenseSMat,DenseSMat) traits:
// T1 - value type in first DenseDMat
// T2 - value type in second DenseDMat
// M0S0 - size in dimension '0' in DenseDMat 0 ( = nb rows in M0 )
// M0S1 - number of columns in M0
// M1S0 - number of rows in M1
// M1S1 - number of columns in M1

template <typename T1, typename T2, Uint N, bool SO>
struct SubTrait<DenseSMat<T1, N, N, SO>, DenseSMat<T2, N, N, SO>>
{
  using type = DenseSMat<typename SubTrait<T1, T2>::type, N, N, SO>;
};

// ----------------------------------------------------------------------------
// MultTrait specializations
// ----------------------------------------------------------------------------

// scalar * DenseSMat
template <typename SValueT, typename T, Uint M, Uint N, bool SO>
struct MultTrait<SValueT, DenseSMat<T, M, N, SO>>
{
  using type = DenseSMat<typename MultTrait<SValueT, T>::type, M, N, SO>;
};

// DenseSMat * DenseSVec
template <typename T1, typename T2, Uint M, Uint N, bool SO, bool TF>
struct MultTrait<DenseSMat<T1, M, N, SO>, DenseSVec<T2, N, TF>>
{
  using type = DenseSVec<typename MultTrait<T1, T2>::type, M, TF>;
};

// DenseSMat * DenseDVec
template <typename T1, typename T2, Uint M, Uint N, bool SO, bool TF>
struct MultTrait<DenseSMat<T1, M, N, SO>, DenseDVec<T2, TF>>
{
  using type = DenseSVec<typename MultTrait<T1, T2>::type, M, TF>;
};

// DenseSMat * DenseConstVecView
template <typename T1, typename T2, bool SO, Uint M, Uint N, bool TF>
struct MultTrait<DenseSMat<T1, M, N, SO>, DenseConstVecView<T2, TF>>
{
  using type = DenseSVec<typename MultTrait<T1, T2>::type, M, ColumnVector>;
};

// DenseSMat * DenseSMat

// T1 - value type in first DenseSMat
// T2 - value type in second DenseSMat
// M0S0 - size in dimension '0' in DenseSMat 0 ( = nb rows in M0 )
// M0S1 - number of columns in M0
// M1S0 - number of rows in M1
// M1S1 - number of columns in M1

template <typename T1, typename T2, Uint N, bool SO>
struct MultTrait<DenseSMat<T1, N, N, SO>, DenseSMat<T2, N, N, SO>>
{
  using type = DenseSMat<typename MultTrait<T1, T2>::type, N, N, SO>;
};

// MatrixBlock * DenseSMat

// T1  - value type in MatrixBlock
// SO1 - storage order in MatrixBlock
// T2  - value type in DenseSMat
// M   - number of rows in DenseSMat
// N   - number of columns in DenseSMat
// SO2 - storage order in DenseSMat

template <typename T1, bool SO1, typename T2, Uint M, Uint N, bool SO2>
struct MultTrait<DenseConstMatView<T1, SO1>, DenseSMat<T2, M, N, SO2>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO1>;
};

// DenseSMat * MatrixBlock

// T1  - value type in DenseSMat
// M   - number of rows in DenseSMat
// N   - number of columns in DenseSMat
// SO1 - storage order in DenseSMat
// T2  - value type in MatrixBlock
// SO2 - storage order in MatrixBlock

template <typename T1, Uint M, Uint N, bool SO1, typename T2, bool SO2>
struct MultTrait<DenseSMat<T1, M, N, SO1>, DenseConstMatView<T2, SO2>>
{
  using type = DenseDMat<typename MultTrait<T1, T2>::type, SO1>;
};

// ----------------------------------------------------------------------------
// TENSOR RANK OF MATRIX = 2
// ----------------------------------------------------------------------------

template <typename T, Uint M, Uint N, bool SO>
struct TensorRank<DenseSMat<T, M, N, SO>>
{
  enum
  {
    value = tensor_rank_2
  };
};

} // Namespace math

} // Namespace pdekit

#endif // PDEKIT_Math_Matrix_hpp
