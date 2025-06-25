#ifndef PDEKIT_Math_Matrix_Cond_Number_hpp
#define PDEKIT_Math_Matrix_Cond_Number_hpp

#include "math/LapackInterface.hpp"
#include "math/Matrix.hpp"
#include "math/unary_ops/matrix_unary_ops/MatrixExpressionCondNumber.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

namespace detail
{

// ----------------------------------------------------------------------------

template <typename MT, bool SO>
struct Norm1CondNumberEvaluator
{
  /// Resulting value type
  using value_type = typename MT::value_type;

  /// Normally we would choose the MatrixArg type to be either
  /// a) const MT&
  /// or
  /// b) MT::tensor_type (which is StaticMatrix or DynamicMatrix)
  /// not to have to allocate an extra matrix. However, here we will perform
  /// LU decomposition of the matrix (Lapack dgetrf) to estimate its
  /// condition number. Therefore we create a copy of the matrix:
  /// The MatrixArg is __ALWAYS__ MT::tensor_type
  // using MatrixArg = typename common::SelectType<MT::evaluates_fast, const
  // MT
  // &,
  //                                    const typename MT::tensor_type>::type;

  using MatrixArg = typename MT::tensor_type;

  inline static value_type const evaluate(const MT &mat)
  {
    MatrixArg matrix(mat);

    const int m          = static_cast<int>(matrix.rows());
    const lapack_int lda = static_cast<lapack_int>(matrix.rows());

    // double LAPACKE_dlange( int matrix_order, char norm, lapack_int m,
    //                            lapack_int n, const double* a, lapack_int
    //                            lda
    //                            );
    const value_type norm1 = Lapack::dlange(LAPACK_ROW_MAJOR, '1', m, m, &matrix(0, 0), lda);

    /*
    value_type norm1 = value_type();

    for (Uint j = 0; j < matrix.cols(); ++j)
    {
      value_type colsum = value_type();

      for (Uint i = 0; i < matrix.rows(); ++i)
      {
        colsum += std::abs(matrix(i, j));
      }
      norm1 = std::max(norm1, colsum);
    }
    */

    lapack_int *ipiv = new lapack_int[matrix.rows()];

    // lapack_int LAPACKE_dgetrf( int matrix_order, lapack_int m, lapack_int
    // n,
    //                            double* a, lapack_int lda, lapack_int*
    //                            ipiv );
    Lapack::dgetrf(LAPACK_ROW_MAJOR, m, m, &matrix(0, 0), lda, ipiv);

    delete[] ipiv;

    value_type cond_number = value_type();

    // lapack_int LAPACKE_dgecon( int matrix_order, char norm, lapack_int n,
    //                            const double* a, lapack_int lda, double
    //                            anorm, double* rcond );
    Lapack::dgecon(LAPACK_ROW_MAJOR, '1', static_cast<lapack_int>(matrix.rows()), &matrix(0, 0),
                   lda, norm1, &cond_number);

    return 1. / cond_number;
  }
};

// ----------------------------------------------------------------------------

template <typename MT, bool SO>
struct InfNormCondNumberEvaluator
{
  /// Resulting value type
  using value_type = typename MT::value_type;

  /// Normally we would choose the MatrixArg type to be either
  /// a) const MT&
  /// or
  /// b) MT::tensor_type (which is StaticMatrix or DynamicMatrix)
  /// not to have to allocate an extra matrix. However, here we will perform
  /// LU decomposition of the matrix (Lapack dgetrf) to estimate its
  /// condition number. Therefore we create a copy of the matrix:
  /// The MatrixArg is __ALWAYS__ MT::tensor_type
  // using MatrixArg = typename common::SelectType<MT::evaluates_fast, const
  // MT
  // &,
  //                                    const typename MT::tensor_type>::type;

  using MatrixArg = typename MT::tensor_type;

  inline static value_type const evaluate(const MT &mat)
  {
    MatrixArg matrix(mat);

    const int m          = static_cast<int>(matrix.rows());
    const lapack_int lda = static_cast<lapack_int>(matrix.rows());

    // double LAPACKE_dlange( int matrix_order, char norm, lapack_int m,
    //                            lapack_int n, const double* a, lapack_int
    //                            lda
    //                            );
    // const value_type norm_inf = LAPACKE_dlange(LAPACK_ROW_MAJOR, 'I', m,
    // m, &matrix(0, 0), lda);

    value_type norm_inf = value_type();

    for (Uint i = 0; i < matrix.rows(); ++i)
    {
      value_type rowsum = value_type();

      for (Uint j = 0; j < matrix.cols(); ++j)
      {
        rowsum += std::abs(matrix(i, j));
      }
      norm_inf = std::max(norm_inf, rowsum);
    }

    lapack_int *ipiv = new lapack_int[matrix.rows()];

    // lapack_int LAPACKE_dgetrf( int matrix_order, lapack_int m, lapack_int
    // n,
    //                            double* a, lapack_int lda, lapack_int*
    //                            ipiv );
    Lapack::dgetrf(LAPACK_ROW_MAJOR, m, m, &matrix(0, 0), lda, ipiv);

    delete[] ipiv;

    value_type cond_number = value_type();

    // lapack_int LAPACKE_dgecon( int matrix_order, char norm, lapack_int n,
    //                            const double* a, lapack_int lda, double
    //                            anorm, double* rcond );
    Lapack::dgecon(LAPACK_ROW_MAJOR, 'I', static_cast<lapack_int>(matrix.rows()), &matrix(0, 0),
                   lda, norm_inf, &cond_number);

    return 1. / cond_number;
  }
};

// ----------------------------------------------------------------------------

} // namespace detail

// ----------------------------------------------------------------------------

template <typename MT, bool SO>
typename MatrixExpressionCondNumber<MT, SO, detail::Norm1CondNumberEvaluator<MT, SO>>::value_type
cond_number_1_norm(math::DenseMatrix<MT, SO> const &mat)
{
  MatrixExpressionCondNumber<MT, SO, detail::Norm1CondNumberEvaluator<MT, SO>> cond_number(
      mat.wrapped_type());
  return cond_number.value();
}

// ----------------------------------------------------------------------------

template <typename MT, bool SO>
typename MatrixExpressionCondNumber<MT, SO, detail::InfNormCondNumberEvaluator<MT, SO>>::value_type
cond_number_inf_norm(math::DenseMatrix<MT, SO> const &mat)
{
  MatrixExpressionCondNumber<MT, SO, detail::InfNormCondNumberEvaluator<MT, SO>> cond_number(
      mat.wrapped_type());
  return cond_number.value();
}

// ----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif // Matrix_Cond_Number_hpp
