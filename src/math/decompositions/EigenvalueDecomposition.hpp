#ifndef PDEKIT_Math_Decompositions_Eigenvalue_Decomposition_hpp
#define PDEKIT_Math_Decompositions_Eigenvalue_Decomposition_hpp

#include "common/Meta.hpp"
#include "math/LapackInterface.hpp"
#include "math/Matrix.hpp"
#include "math/TensorRank.hpp"

namespace pdekit
{

namespace math
{

// ============================================================================

template <typename T>
class EigenvalueDecomposition
{
  public:
  /// Default constructor
  EigenvalueDecomposition();

  /// Default destructor
  ~EigenvalueDecomposition();

  /// Set the rank of matrix to decompose
  void set_rank(const Uint N);

  /// Decompose a matrix, compute the matrix of eigenvalues and right
  /// eigenvectors
  template <typename MT1, typename MT2, typename MT3, bool SO>
  void factorize(DenseMatrix<MT1, SO> &M, DenseMatrix<MT2, SO> &Eig, DenseMatrix<MT3, SO> &R);

  /// Decompose a matrix, compute the matrix of eigenvalues, matrix of right
  /// eigenvectors and also its inverse
  template <typename MT1, typename MT2, typename MT3, typename MT4, bool SO>
  void factorize(DenseMatrix<MT1, SO> &M, DenseMatrix<MT2, SO> &Eig, DenseMatrix<MT3, SO> &R,
                 DenseMatrix<MT4, SO> &Rinv);

  private:
  /// Rank of matrix to decompose
  Uint m_rank;

  /// Real parts of the computed eigenvalues
  T *m_wr;

  /// Imaginary parts of the computed eigenvalues
  T *m_wi;

  /// Array containing the left eigenvectors of A
  T *m_vl;

  /// Leading dimension of the array m_vl
  Int m_ldvl;

  /// Array containing the right eigenvectors of A
  T *m_vr;

  /// Leading dimension of the array m_vr
  Int m_ldvr;

  /// Temporary work array
  T *m_work;

  /// The dimension of the array WORK.  LWORK >= max(1,3*N), and
  /// if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
  /// performance, LWORK must generally be larger.
  Int m_lwork;
  // Int* m_pivot;
};

// ============================================================================

template <typename T>
EigenvalueDecomposition<T>::EigenvalueDecomposition()
    : m_rank(0), m_wr(nullptr), m_wi(nullptr), m_vl(nullptr), m_ldvl(0), m_vr(nullptr), m_ldvr(0),
      m_work(nullptr), m_lwork(0)

{
}

// ============================================================================

template <typename T>
EigenvalueDecomposition<T>::~EigenvalueDecomposition()
{
  delete[] m_wr;
  delete[] m_wi;
  delete[] m_vl;
  delete[] m_vr;
  delete[] m_work;
}

// ============================================================================

template <typename T>
void EigenvalueDecomposition<T>::set_rank(const Uint N)
{

  if (m_rank != N)
  {
    m_rank = N;

    delete[] m_wr;
    m_wr = new T[m_rank];

    delete[] m_wi;
    m_wi = new T[m_rank];

    delete[] m_vl;
    m_vl   = new T[m_rank * m_rank];
    m_ldvl = m_rank;

    delete[] m_vr;
    m_vr   = new T[m_rank * m_rank];
    m_ldvr = m_rank;

    delete[] m_work;
    m_lwork = 10 * m_rank;
    m_work  = new T[m_lwork];
  }
}

// ============================================================================

template <typename T>
template <typename MT1, typename MT2, typename MT3, bool SO>
void EigenvalueDecomposition<T>::factorize(DenseMatrix<MT1, SO> &M, DenseMatrix<MT2, SO> &Eig,
                                           DenseMatrix<MT3, SO> &R)
{
  static_assert(common::TypesAreIdentical<T, typename MT1::value_type>::value,
                "Matrix to decompose ('M') has wrong value type");

  static_assert(common::TypesAreIdentical<T, typename MT2::value_type>::value,
                "Matrix of eigenvalues ('Eig') has wrong value type");

  static_assert(common::TypesAreIdentical<T, typename MT3::value_type>::value,
                "Matrix of right eigenvectors ('R') has wrong value type");

  MT1 &Mw   = M.wrapped_type();
  MT2 &Eigw = Eig.wrapped_type();
  MT3 &Rw   = R.wrapped_type();

  if (Mw.rows() != m_rank)
  {
    set_rank(Mw.rows());
  }

  /*
  int matrix_order, char jobvl, char jobvr,
  lapack_int n, double* a, lapack_int lda, double* wr,
  double* wi, double* vl, lapack_int ldvl, double* vr,
  lapack_int ldvr
   */

  // jobvl = 'N' => do not compute left eigenvectors
  // jobvr = 'V' => compute right eigenvectors
  // #if PDEKIT_HAVE_LAPACK
  Lapack::dgeev(LAPACK_ROW_MAJOR, 'N', 'V', static_cast<lapack_int>(m_rank), &Mw(0, 0),
                static_cast<lapack_int>(m_rank), m_wr, m_wi, m_vl, m_ldvl, m_vr, m_ldvr);
  // #endif

  for (Uint i = 0; i < m_rank; ++i)
  {
    Eigw(i, i) = m_wr[i];
  }

  for (Uint i = 0; i < m_rank; ++i)
  {
    for (Uint j = 0; j < m_rank; ++j)
    {
      Rw(i, j) = m_vr[i * m_rank + j];
    }
  }
}

// ============================================================================

template <typename T>
template <typename MT1, typename MT2, typename MT3, typename MT4, bool SO>
void EigenvalueDecomposition<T>::factorize(DenseMatrix<MT1, SO> &M, DenseMatrix<MT2, SO> &Eig,
                                           DenseMatrix<MT3, SO> &R, DenseMatrix<MT4, SO> &Rinv)
{
  factorize(M, Eig, R);

  MT3 &Rw    = R.wrapped_type();
  MT4 &Rwinv = Rinv.wrapped_type();

  Rw.inv(Rwinv);
}

// ============================================================================

} // Namespace math

} // Namespace pdekit

#endif
