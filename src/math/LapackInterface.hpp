#ifndef PDEKIT_Math_Lapack_Interface_hpp
#define PDEKIT_Math_Lapack_Interface_hpp

#include "PDEKit_Config.hpp"

#if PDEKIT_BLAS_LAPACK_PROVIDER_IS_NETLIB

// This is the C++ header for complex numbers
#include <complex>

// First define lapack_complex_float and lapack_complex_double
// This will prevent the preprocessor from pulling in complex.h
// after lapacke.h is included
// The C-header complex.h has compilation issues in c++11 !
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include "lapacke.h"

#elif PDEKIT_BLAS_LAPACK_PROVIDER_IS_OPENBLAS

// This is the C++ header for complex numbers
#include <complex>

// First define lapack_complex_float and lapack_complex_double
// This will prevent the preprocessor from pulling in complex.h
// after lapacke.h is included
// The C-header complex.h has compilation issues in c++11 !
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include "lapacke.h"

#endif

namespace pdekit
{
namespace math
{
class Lapack
{
  public:
  static float slange(const int matrix_order, const char norm, const int m, const int n,
                      const float *a, const int lda)
  {
    return LAPACKE_slange(matrix_order, norm, m, n, a, lda);
  }

  static double dlange(const int matrix_order, const char norm, const int m, const int n,
                       const double *a, const int lda)
  {
    return LAPACKE_dlange(matrix_order, norm, m, n, a, lda);
  }

  static inline void sgetri(const int matrix_order, const int n, float *a, const int lda,
                            const int *ipiv)
  {
    LAPACKE_sgetri(matrix_order, n, a, lda, ipiv);
  }

  static inline void dgetri(const int matrix_order, const int n, double *a, const int lda,
                            const int *ipiv)
  {
    LAPACKE_dgetri(matrix_order, n, a, lda, ipiv);
  }

  static inline void sgetrf(const int matrix_order, const int m, const int n, float *a,
                            const int lda, int *ipiv)
  {
    LAPACKE_sgetrf(matrix_order, m, n, a, lda, ipiv);
  }

  static inline void dgetrf(const int matrix_order, const int m, const int n, double *a,
                            const int lda, int *ipiv)
  {
    LAPACKE_dgetrf(matrix_order, m, n, a, lda, ipiv);
  }

  static inline int sgecon(const int matrix_order, const char norm, const int n, const float *a,
                           const int lda, float anorm, float *rcond)
  {
    return LAPACKE_sgecon(matrix_order, norm, n, a, lda, anorm, rcond);
  }

  static inline int dgecon(const int matrix_order, const char norm, const int n, const double *a,
                           const int lda, double anorm, double *rcond)
  {
    return LAPACKE_dgecon(matrix_order, norm, n, a, lda, anorm, rcond);
  }

  static inline int sgeev(const int matrix_order, const char jobvl, const char jobvr, const int n,
                          float *a, const int lda, float *wr, float *wi, float *vl, const int ldvl,
                          float *vr, const int ldvr)
  {
    return LAPACKE_sgeev(matrix_order, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr);
  }

  static inline int dgeev(const int matrix_order, const char jobvl, const char jobvr, const int n,
                          double *a, const int lda, double *wr, double *wi, double *vl,
                          const int ldvl, double *vr, const int ldvr)
  {
    return LAPACKE_dgeev(matrix_order, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr);
  }
};
} // namespace math
} // namespace pdekit

#endif
