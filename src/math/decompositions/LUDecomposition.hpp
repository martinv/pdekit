#ifndef PDEKIT_Math_Decompositions_LU_Decomposition_hpp
#define PDEKIT_Math_Decompositions_LU_Decomposition_hpp

#include "common/Meta.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/LapackInterface.hpp"
#include "math/Matrix.hpp"
#include "math/TensorRank.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename T>
class LUDecomposition
{
  public:
  /// Default constructor
  LUDecomposition();

  /// Default destructor
  ~LUDecomposition();

  /// Set the rank of matrix to decompose
  void set_rank(const Uint N);

  /// Decompose a matrix, compute the matrix of eigenvalues and right
  /// eigenvectors
  template <typename MT1, typename MT2, typename MT3, bool SO>
  void factorize(DenseMatrix<MT1, SO> &M, DenseMatrix<MT2, SO> &L, DenseMatrix<MT3, SO> &U);

  /// Return the reordering
  const math::DenseConstVecView<Int> reordering() const;

  private:
  /// Rank of matrix to decompose
  Uint m_rank;

  /// The pivot indices; for 1 <= i <= min(M,N), row i of the
  /// *          matrix was interchanged with row IPIV(i).
  Int *m_ipiv;
  // Int* m_pivot;
};

// ----------------------------------------------------------------------------

template <typename T>
LUDecomposition<T>::LUDecomposition() : m_rank(0), m_ipiv(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename T>
LUDecomposition<T>::~LUDecomposition()
{
  delete[] m_ipiv;
}

// ----------------------------------------------------------------------------

template <typename T>
void LUDecomposition<T>::set_rank(const Uint N)
{

  if (m_rank != N)
  {
    m_rank = N;

    delete[] m_ipiv;
    m_ipiv = new Int[m_rank];
  }
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MT1, typename MT2, typename MT3, bool SO>
void LUDecomposition<T>::factorize(DenseMatrix<MT1, SO> &M, DenseMatrix<MT2, SO> &L,
                                   DenseMatrix<MT3, SO> &U)
{
  static_assert(common::TypesAreIdentical<T, typename MT1::value_type>::value,
                "Matrix to decompose ('M') has wrong value type");

  static_assert(common::TypesAreIdentical<T, typename MT2::value_type>::value,
                "Lower diagonal matrix ('L') has wrong value type");

  static_assert(common::TypesAreIdentical<T, typename MT3::value_type>::value,
                "Uppwer diagonal matrix ('U') has wrong value type");

  MT1 &Mw = M.wrapped_type();
  MT2 &Lw = L.wrapped_type();
  MT3 &Uw = U.wrapped_type();

  if (Mw.rows() != m_rank)
  {
    set_rank(Mw.rows());
  }

  // The decomposition is done in-place, so we copy M to L, and dempose memory
  // held by L The result is stored completely in L, the digonal values are
  // equal to 1 and are assumed implicitly
  for (Uint i = 0; i < m_rank; ++i)
  {
    for (Uint j = 0; j < m_rank; ++j)
    {
      Lw(i, j) = Mw(i, j);
    }
  }

  // lapack_int LAPACKE_dgetrf( int matrix_order, lapack_int m, lapack_int n,
  //                            double* a, lapack_int lda, lapack_int* ipiv );

  Lapack::dgetrf(LAPACK_ROW_MAJOR, static_cast<lapack_int>(m_rank), static_cast<lapack_int>(m_rank),
                 &Lw(0, 0), static_cast<lapack_int>(Mw.cols()), m_ipiv);

  Uw.fill(0.0);

  for (Uint i = 0; i < m_rank; ++i)
  {
    std::cout << " " << m_ipiv[i];
  }
  std::cout << std::endl;

  // Copy upper part from Lw to Uw and set corresponding entries in Lw to zero
  for (Uint i = 0; i < m_rank; ++i)
  {
    for (Uint j = i; j < m_rank; ++j)
    {
      Uw(i, j) = Lw(i, j);
      Lw(i, j) = 0.0;
    }
  }

  for (Uint i = 0; i < m_rank; ++i)
  {
    Lw(i, i) = 1.0;
  }

  /*
  // This is the ROW reordering (aka permutation matrix P) we would have to
  apply
  // TO THE ORIGINAL MATRIX M in order to have L * U = P * M
  for(Uint i = 0; i < m_rank; ++i)
  {
    if ( static_cast<Uint>(m_ipiv[i]-1) != i )
    {
      for(Uint j = 0; j < m_rank; ++j)
      {
        std::swap(Lw(i,j), Lw(m_ipiv[i]-1,j));
      }
    }
  }
  */

  // Note that the reordering coming from dgetrf is in Fortran indexing, i.e.
  // starting from 1 Therefore we subtract one from each value
  for (Uint i = 0; i < m_rank; ++i)
  {
    m_ipiv[i]--;
  }
}

// ----------------------------------------------------------------------------

template <typename T>
const math::DenseConstVecView<Int> LUDecomposition<T>::reordering() const
{
  return math::DenseConstVecView<Int>(m_ipiv, m_rank, 1u);
}

// ----------------------------------------------------------------------------

} // Namespace math

} // Namespace pdekit

#endif
