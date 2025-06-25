#ifndef PDEKIT_Math_Matrix_Inverse_hpp
#define PDEKIT_Math_Matrix_Inverse_hpp

#include <typeinfo>

#include "common/PDEKit.hpp"
#include "math/LapackInterface.hpp"
#include "math/unary_ops/matrix_unary_ops/MatrixDeterminant.hpp"

namespace pdekit
{

namespace math
{

// GENERAL STRATEGY FOR INVERTING DENSE MATRICES:
// For static and dynamic matrices of size <= 4, use DirectMatrixInverse,
// which is hardcoded matrix inverse for small matrices
// For static and dynamic matrix of size > 4, use LapackMatrixInverse, which
// should work for any size, but it is slower

// ----------------------------------------------------------------------------
// MATRIX INVERSE COMPUTATION WITH LAPACK
// ----------------------------------------------------------------------------

template <typename T>
class LapackMatrixInverse
{
  public:
  /// TYPEDEFS
  inline static void compute(const T *const m1, T *m2, const Uint size);

  private:
};

template <typename T>
void LapackMatrixInverse<T>::compute(const T *const m1, T *m2, const Uint size)
{
  std::cerr << "LapackMatrixInverse<T>: sorry, the inverse is not supported ";
  std::cerr << "for type \"" << typeid(T).name() << "\"" << std::endl;
}

// ----------------------------------------------------------------------------
// Lapack matrix inverse specialization for single precision real numbers
// ----------------------------------------------------------------------------

template <>
class LapackMatrixInverse<Float>
{
  public:
  /// TYPEDEFS
  inline static void compute(const Float *const m1, Float *m2, const Uint size);

  private:
};

// This function is inline, we define it directly in the header
void LapackMatrixInverse<Float>::compute(const Float *const m1, Float *m2, const Uint size)
{
  for (Uint i = 0; i < size * size; ++i)
  {
    m2[i] = m1[i];
  }

  Int *IPIV = new Int[size + 1];

  Lapack::sgetrf(LAPACK_ROW_MAJOR, size, size, m2, size, IPIV);
  Lapack::sgetri(LAPACK_ROW_MAJOR, size, m2, size, IPIV);

  delete[] IPIV;
}

// ----------------------------------------------------------------------------
// Lapack matrix inverse specialization for double precision real numbers
// ----------------------------------------------------------------------------

template <>
class LapackMatrixInverse<Real>
{
  public:
  inline static void compute(const Real *const m1, Real *m2, const Uint size);

  private:
};

// This function is inline, we define it directly in the header
void LapackMatrixInverse<Real>::compute(const Real *const m1, Real *m2, const Uint size)
{
  for (Uint i = 0; i < size * size; ++i)
  {
    m2[i] = m1[i];
  }

  Int *IPIV = new Int[size + 1];

  Lapack::dgetrf(LAPACK_ROW_MAJOR, size, size, m2, size, IPIV);
  Lapack::dgetri(LAPACK_ROW_MAJOR, size, m2, size, IPIV);

  delete[] IPIV;
}

// ----------------------------------------------------------------------------
// DIRECT MATRIX INVERSE COMPUTATION
// ----------------------------------------------------------------------------

template <typename T, Uint M>
class DirectMatrixInverse;

// -----------------------------
// SPECIALIZATION FOR 1x1 MATRIX
// -----------------------------

template <typename T>
class DirectMatrixInverse<T, 1>
{
  public:
  inline static void compute(const T *const m1, T *m2, const Uint size = 2);

  private:
};

template <typename T>
void DirectMatrixInverse<T, 1>::compute(const T *const m1, T *m2, const Uint size)
{
  m2[0] = 1.0 / m1[0];
  return;
}

// -----------------------------
// SPECIALIZATION FOR 2x2 MATRIX
// -----------------------------

template <typename T>
class DirectMatrixInverse<T, 2>
{
  public:
  inline static void compute(const T *const m1, T *m2, const Uint size = 2);

  private:
};

template <typename T>
void DirectMatrixInverse<T, 2>::compute(const T *const m1, T *m2, const Uint size)
{
  const Real det    = m1[0] * m1[3] - m1[1] * m1[2];
  const Real invDet = 1. / det;

  m2[0] = m1[3] * invDet;
  m2[1] = -m1[1] * invDet;
  m2[2] = -m1[2] * invDet;
  m2[3] = m1[0] * invDet;
}

// -----------------------------
// SPECIALIZATION FOR 3x3 MATRIX
// -----------------------------

template <typename T>
class DirectMatrixInverse<T, 3>
{
  public:
  inline static void compute(const T *const m1, T *m2, const Uint size = 3);

  private:
};

template <typename T>
void DirectMatrixInverse<T, 3>::compute(const T *const m1, T *m2, const Uint size)
{
  const Real det = MatrixDeterminant<T, 3>::compute(m1);
  //  const Real det = m1.det();
  //  cf_assert(MathChecks::isNotZero(det));
  const Real invDet = 1. / det;

  m2[0] = (m1[4] * m1[8] - m1[5] * m1[7]) * invDet;
  m2[1] = -(m1[1] * m1[8] - m1[2] * m1[7]) * invDet;
  m2[2] = (m1[1] * m1[5] - m1[4] * m1[2]) * invDet;
  m2[3] = -(m1[3] * m1[8] - m1[5] * m1[6]) * invDet;
  m2[4] = (m1[0] * m1[8] - m1[2] * m1[6]) * invDet;
  m2[5] = -(m1[0] * m1[5] - m1[2] * m1[3]) * invDet;
  m2[6] = (m1[3] * m1[7] - m1[4] * m1[6]) * invDet;
  m2[7] = -(m1[0] * m1[7] - m1[1] * m1[6]) * invDet;
  m2[8] = (m1[0] * m1[4] - m1[1] * m1[3]) * invDet;
}

// -----------------------------
// SPECIALIZATION FOR 4x4 MATRIX
// -----------------------------

template <typename T>
class DirectMatrixInverse<T, 4>
{
  public:
  inline static void compute(const T *const m1, T *m2, const Uint size = 4);

  private:
};

template <typename T>
void DirectMatrixInverse<T, 4>::compute(const T *const m1, T *m2, const Uint size)
{
  const Real t14 = m1[0] * m1[5];
  const Real t15 = m1[10] * m1[15];
  const Real t17 = m1[11] * m1[14];
  const Real t19 = m1[0] * m1[9];
  const Real t20 = m1[6] * m1[15];
  const Real t22 = m1[7] * m1[14];
  const Real t24 = m1[0] * m1[13];
  const Real t25 = m1[6] * m1[11];
  const Real t27 = m1[7] * m1[10];
  const Real t29 = m1[4] * m1[1];
  const Real t32 = m1[4] * m1[9];
  const Real t33 = m1[2] * m1[15];
  const Real t35 = m1[3] * m1[14];
  const Real t37 = m1[4] * m1[13];
  const Real t38 = m1[2] * m1[11];
  const Real t40 = m1[3] * m1[10];
  const Real t42 = t14 * t15 - t14 * t17 - t19 * t20 + t19 * t22 + t24 * t25 - t24 * t27 -
                   t29 * t15 + t29 * t17 + t32 * t33 - t32 * t35 - t37 * t38 + t37 * t40;
  const Real t43 = m1[8] * m1[1];
  const Real t46 = m1[8] * m1[5];
  const Real t49 = m1[8] * m1[13];
  const Real t50 = m1[2] * m1[7];
  const Real t52 = m1[3] * m1[6];
  const Real t54 = m1[12] * m1[1];
  const Real t57 = m1[12] * m1[5];
  const Real t60 = m1[12] * m1[9];
  const Real t63 = t43 * t20 - t43 * t22 - t46 * t33 + t46 * t35 + t49 * t50 - t49 * t52 -
                   t54 * t25 + t54 * t27 + t57 * t38 - t57 * t40 - t60 * t50 + t60 * t52;

  const Real ddet = t42 + t63;
  //  cf_assert(MathChecks::isNotZero(ddet));
  const Real deter = 1. / ddet;

  const Real t71  = m1[9] * m1[2];
  const Real t73  = m1[9] * m1[3];
  const Real t75  = m1[13] * m1[2];
  const Real t77  = m1[13] * m1[3];
  const Real t81  = m1[1] * m1[6];
  const Real t83  = m1[1] * m1[7];
  const Real t85  = m1[5] * m1[2];
  const Real t87  = m1[5] * m1[3];
  const Real t119 = m1[8] * m1[2];
  const Real t121 = m1[8] * m1[3];
  const Real t123 = m1[12] * m1[2];
  const Real t125 = m1[12] * m1[3];
  const Real t129 = m1[0] * m1[6];
  const Real t131 = m1[0] * m1[7];
  const Real t133 = m1[4] * m1[2];
  const Real t135 = m1[4] * m1[3];

  m2[0] = (m1[5] * m1[10] * m1[15] - m1[5] * m1[11] * m1[14] - m1[9] * m1[6] * m1[15] +
           m1[9] * m1[7] * m1[14] + m1[13] * m1[6] * m1[11] - m1[13] * m1[7] * m1[10]) *
          deter;
  m2[1] = -(m1[1] * m1[10] * m1[15] - m1[1] * m1[11] * m1[14] - t71 * m1[15] + t73 * m1[14] +
            t75 * m1[11] - t77 * m1[10]) *
          deter;
  m2[2] = (t81 * m1[15] - t83 * m1[14] - t85 * m1[15] + t87 * m1[14] + t75 * m1[7] - t77 * m1[6]) *
          deter;
  m2[3] = -(t81 * m1[11] - t83 * m1[10] - t85 * m1[11] + t87 * m1[10] + t71 * m1[7] - t73 * m1[6]) *
          deter;

  m2[4] = (-m1[4] * m1[10] * m1[15] + m1[4] * m1[11] * m1[14] + m1[8] * m1[6] * m1[15] -
           m1[8] * m1[7] * m1[14] - m1[12] * m1[6] * m1[11] + m1[12] * m1[7] * m1[10]) *
          deter;
  m2[5] = (m1[0] * m1[10] * m1[15] - m1[0] * m1[11] * m1[14] - t119 * m1[15] + t121 * m1[14] +
           t123 * m1[11] - t125 * m1[10]) *
          deter;
  m2[6] = -(t129 * m1[15] - t131 * m1[14] - t133 * m1[15] + t135 * m1[14] + t123 * m1[7] -
            t125 * m1[6]) *
          deter;
  m2[7] = (t129 * m1[11] - t131 * m1[10] - t133 * m1[11] + t135 * m1[10] + t119 * m1[7] -
           t121 * m1[6]) *
          deter;

  m2[8] =
      -(-t32 * m1[15] + t37 * m1[11] + t46 * m1[15] - t49 * m1[7] - t57 * m1[11] + t60 * m1[7]) *
      deter;
  m2[9] =
      -(t19 * m1[15] - t24 * m1[11] - t43 * m1[15] + t121 * m1[13] + t54 * m1[11] - t125 * m1[9]) *
      deter;
  m2[10] =
      (t14 * m1[15] - t24 * m1[7] - t29 * m1[15] + t135 * m1[13] + t54 * m1[7] - t125 * m1[5]) *
      deter;
  m2[11] =
      -(t14 * m1[11] - t19 * m1[7] - t29 * m1[11] + t135 * m1[9] + t43 * m1[7] - t121 * m1[5]) *
      deter;

  m2[12] =
      -(t32 * m1[14] - t37 * m1[10] - t46 * m1[14] + t49 * m1[6] + t57 * m1[10] - t60 * m1[6]) *
      deter;
  m2[13] =
      (t19 * m1[14] - t24 * m1[10] - t43 * m1[14] + t119 * m1[13] + t54 * m1[10] - t123 * m1[9]) *
      deter;
  m2[14] =
      -(t14 * m1[14] - t24 * m1[6] - t29 * m1[14] + t133 * m1[13] + t54 * m1[6] - t123 * m1[5]) *
      deter;
  m2[15] = (t14 * m1[10] - t19 * m1[6] - t29 * m1[10] + t133 * m1[9] + t43 * m1[6] - t119 * m1[5]) *
           deter;
}

// ----------------------------------------------------
// Convenience wrapper for DynamicMatrices. It picks
// DirectMatrixInverse of LapackMatrixInverse depending
// on the size of the matrix
// ----------------------------------------------------

template <typename T>
class DynamicMatrixInverter
{
  public:
  /// TYPEDEFS:

  // Compute the matrix inverse
  inline static void invert(const T *const m1, T *m2, const Uint size)
  {
    switch (size)
    {
      case 1:
      {
        m2[0] = 1.0 / m1[0];
        return;
      }
      case 2:
      {
        DirectMatrixInverse<T, 2>::compute(m1, m2, 2);
        return;
      }
      case 3:
      {
        DirectMatrixInverse<T, 3>::compute(m1, m2, 3);
        return;
      }
      case 4:
      {
        DirectMatrixInverse<T, 4>::compute(m1, m2, 4);
        return;
      }
      default:
        LapackMatrixInverse<T>::compute(m1, m2, size);
    };
  }

  private:
};

} // Namespace math

} // Namespace pdekit

#endif // Matrix_Inverse_hpp
