#ifndef PDEKIT_Math_PolyLib_hpp
#define PDEKIT_Math_PolyLib_hpp

#include "common/PDEKit.hpp"
#include "math/DenseDVec.hpp"

/*
LIBRARY ROUTINES FOR ORTHOGONAL POLYNOMIAL CALCULUS AND INTERPOLATION

  Spencer Sherwin
  Aeronautics, Imperial College London

  Based on codes by Einar Ronquist and Ron Henderson

  Abbreviations
  - z    -   Set of collocation/quadrature points
  - w    -   Set of quadrature weights
  - D    -   Derivative matrix
  - h    -   Lagrange Interpolant
  - I    -   Interpolation matrix
  - g    -   Gauss
  - gr   -   Gauss-Radau
  - gl   -   Gauss-Lobatto
  - j    -   Jacobi
  - l    -   Legendre  (Jacobi with alpha = beta =  0.0)
  - c    -   Chebychev (Jacobi with alpha = beta = -0.5)
  - m    -   point at minus 1 in Radau rules
  - p    -   point at plus  1 in Radau rules

  -----------------------------------------------------------------------
                         M A I N     R O U T I N E S
  -----------------------------------------------------------------------

  Points and Weights:

  zwgj        Compute Gauss-Jacobi         points and weights
  zwgrjm      Compute Gauss-Radau-Jacobi   points and weights (z=-1)
  zwgrjp      Compute Gauss-Radau-Jacobi   points and weights (z= 1)
  zwglj       Compute Gauss-Lobatto-Jacobi points and weights

  Derivative Matrices:

  Dgj         Compute Gauss-Jacobi         derivative matrix
  Dgrjm       Compute Gauss-Radau-Jacobi   derivative matrix (z=-1)
  Dgrjp       Compute Gauss-Radau-Jacobi   derivative matrix (z= 1)
  Dglj        Compute Gauss-Lobatto-Jacobi derivative matrix

  Lagrange Interpolants:

  hgj         Compute Gauss-Jacobi         Lagrange interpolants
  hgrjm       Compute Gauss-Radau-Jacobi   Lagrange interpolants (z=-1)
  hgrjp       Compute Gauss-Radau-Jacobi   Lagrange interpolants (z= 1)
  hglj        Compute Gauss-Lobatto-Jacobi Lagrange interpolants

  Interpolation Operators:

  Imgj        Compute interpolation operator gj->m
  Imgrjm      Compute interpolation operator grj->m (z=-1)
  Imgrjp      Compute interpolation operator grj->m (z= 1)
  Imglj       Compute interpolation operator glj->m

  Polynomial Evaluation:

  jacobfd     Returns value and derivative of Jacobi poly. at point z
  jacobd      Returns derivative of Jacobi poly. at point z (valid at z=-1,1)

  -----------------------------------------------------------------------
                     L O C A L      R O U T I N E S
  -----------------------------------------------------------------------

  jacobz      Returns Jacobi polynomial zeros
  gammaf      Gamma function for integer values and halves

  -----------------------------------------------------------------------
                         M A C R O S
  -----------------------------------------------------------------------

  Legendre  polynomial alpha = beta = 0
  Chebychev polynomial alpha = beta = -0.5

  Points and Weights:

  zwgl        Compute Gauss-Legendre          points and weights
  zwgrlm      Compute Gauss-Radau-Legendre    points and weights (z=-1)
  zwgrlp      Compute Gauss-Radau-Legendre    points and weights (z=+1)
  zwgll       Compute Gauss-Lobatto-Legendre  points and weights

  zwgc        Compute Gauss-Chebychev         points and weights
  zwgrcm      Compute Gauss-Radau-Chebychev   points and weights (z=-1)
  zwgrcp      Compute Gauss-Radau-Chebychev   points and weights (z=+1)
  zwglc       Compute Gauss-Lobatto-Chebychev points and weights

  Derivative Operators:

  Dgl         Compute Gauss-Legendre          derivative matrix
  Dgrlm       Compute Gauss-Radau-Legendre    derivative matrix (z=-1)
  Dgrlp       Compute Gauss-Radau-Legendre    derivative matrix (z=+1)
  Dgll        Compute Gauss-Lobatto-Legendre  derivative matrix

  Dgc         Compute Gauss-Chebychev         derivative matrix
  Dgrcm       Compute Gauss-Radau-Chebychev   derivative matrix (z=-1)
  Dgrcp       Compute Gauss-Radau-Chebychev   derivative matrix (z=+1)
  Dglc        Compute Gauss-Lobatto-Chebychev derivative matrix

  Lagrangian Interpolants:

  hgl         Compute Gauss-Legendre          Lagrange interpolants
  hgrlm       Compute Gauss-Radau-Legendre    Lagrange interpolants (z=-1)
  hgrlp       Compute Gauss-Radau-Legendre    Lagrange interpolants (z=+1)
  hgll        Compute Gauss-Lobatto-Legendre  Lagrange interpolants

  hgc         Compute Gauss-Chebychev         Lagrange interpolants
  hgrcm       Compute Gauss-Radau-Chebychev   Lagrange interpolants (z=-1)
  hgrcp       Compute Gauss-Radau-Chebychev   Lagrange interpolants (z=+1)
  hglc        Compute Gauss-Lobatto-Chebychev Lagrange interpolants

  Interpolation Operators:

  Imgl        Compute interpolation operator gl->m
  Imgrlm      Compute interpolation operator grl->m (z=-1)
  Imgrlp      Compute interpolation operator grl->m (z=+1)
  Imgll       Compute interpolation operator gll->m

  Imgc        Compute interpolation operator gc->m
  Imgrcm      Compute interpolation operator grc->m (z=-1)
  Imgrcp      Compute interpolation operator grc->m (z=+1)
  Imglc       Compute interpolation operator glc->m

  ------------------------------------------------------------------------

  Useful references:

  - [1] Gabor Szego: Orthogonal Polynomials, American Mathematical Society,
      Providence, Rhode Island, 1939.
  - [2] Abramowitz \& Stegun: Handbook of Mathematical Functions,
      Dover, New York, 1972.
  - [3] Canuto, Hussaini, Quarteroni \& Zang: Spectral Methods in Fluid
      Dynamics, Springer-Verlag, 1988.
  - [4] Ghizzetti \& Ossicini: Quadrature Formulae, Academic Press, 1970.
  - [5] Karniadakis \& Sherwin: Spectral/hp element methods for CFD, 1999


  NOTES
  -----
  (1) All routines are double precision.
  (2) All array subscripts start from zero, i.e. vector[0..N-1]
*/

namespace pdekit
{

namespace math
{

// ============================================================================

/**
   \brief  Gauss-Jacobi zeros and weights.

   \li Generate \a np Gauss Jacobi zeros, \a z, and weights,\a w,
   associated with the Jacobi polynomial \f$ P^{\alpha,\beta}_{np}(z)
   \f$,

   \li Exact for polynomials of order \a 2np-1 or less
*/

void zwgj(DenseDVec<Real> &z, DenseDVec<Real> &w, const Uint np, const Real alpha, const Real beta);

// ============================================================================

/**
    \brief Routine to calculate Jacobi polynomials, \f$
    P^{\alpha,\beta}_n(z) \f$, and their first derivative, \f$
    \frac{d}{dz} P^{\alpha,\beta}_n(z) \f$.

    \li This function returns the vectors \a poly_in and \a poly_d
    containing the value of the \f$ n^th \f$ order Jacobi polynomial
    \f$ P^{\alpha,\beta}_n(z) \alpha > -1, \beta > -1 \f$ and its
    derivative at the \a np points in \a z[i]

    - If \a poly_in = NULL then only calculate derivatice

    - If \a polyd   = NULL then only calculate polynomial

    - To calculate the polynomial this routine uses the recursion
    relationship (see appendix A ref [4]) :
    \f$ \begin{array}{rcl}
    P^{\alpha,\beta}_0(z) &=& 1 \\
    P^{\alpha,\beta}_1(z) &=& \frac{1}{2} [ \alpha-\beta+(\alpha+\beta+2)z] \\
    a^1_n P^{\alpha,\beta}_{n+1}(z) &=& (a^2_n + a^3_n z)
    P^{\alpha,\beta}_n(z) - a^4_n P^{\alpha,\beta}_{n-1}(z) \\
    a^1_n &=& 2(n+1)(n+\alpha + \beta + 1)(2n + \alpha + \beta) \\
    a^2_n &=& (2n + \alpha + \beta + 1)(\alpha^2 - \beta^2)  \\
    a^3_n &=& (2n + \alpha + \beta)(2n + \alpha + \beta + 1)
    (2n + \alpha + \beta + 2)  \\
    a^4_n &=& 2(n+\alpha)(n+\beta)(2n + \alpha + \beta + 2)
    \end{array} \f$

    - To calculate the derivative of the polynomial this routine uses
    the relationship (see appendix A ref [4]) :
    \f$ \begin{array}{rcl}
    b^1_n(z)\frac{d}{dz} P^{\alpha,\beta}_n(z)&=&b^2_n(z)P^{\alpha,\beta}_n(z)
    + b^3_n(z) P^{\alpha,\beta}_{n-1}(z) \hspace{2.2cm} \\
    b^1_n(z) &=& (2n+\alpha + \beta)(1-z^2) \\
    b^2_n(z) &=& n[\alpha - \beta - (2n+\alpha + \beta)z]\\
    b^3_n(z) &=& 2(n+\alpha)(n+\beta)
    \end{array} \f$

    - Note the derivative from this routine is only valid for -1 < \a z < 1.
**/

void jacobfd(const Uint np, DenseDVec<Real> &z, DenseDVec<Real> &poly_in, DenseDVec<Real> &polyd,
             const Uint n, const Real alpha, const Real beta);

// ============================================================================

/**
   \brief Calculate the  derivative of Jacobi polynomials

   \li Generates a vector \a poly of values of the derivative of the
   \a n th order Jacobi polynomial \f$ P^(\alpha,\beta)_n(z)\f$ at the
   \a np points \a z.

   \li To do this we have used the relation
   \n
   \f$ \frac{d}{dz} P^{\alpha,\beta}_n(z)
   = \frac{1}{2} (\alpha + \beta + n + 1)  P^{\alpha,\beta}_n(z) \f$

   \li This formulation is valid for \f$ -1 \leq z \leq 1 \f$

**/

void jacobd(const Uint np, DenseDVec<Real> &z, DenseDVec<Real> &polyd, const Uint n,
            const Real alpha, const Real beta);

// ============================================================================

/**
   \brief Calculate the Gamma function , \f$ \Gamma(n)\f$, for integer
   values and halves.

   Determine the value of \f$\Gamma(n)\f$ using:

   \f$ \Gamma(n) = (n-1)!  \mbox{ or  }  \Gamma(n+1/2) = (n-1/2)\Gamma(n-1/2)\f$

   where \f$ \Gamma(1/2) = \sqrt(\pi)\f$
**/

Real gammaF(const Real x);

// ============================================================================

/**
   \brief Zero determination through the eigenvalues of a tridiagonal
   matrix from teh three term recursion relationship.

   Set up a symmetric tridiagonal matrix

   \f$ \left [  \begin{array}{ccccc}
   a[0] & b[0]   &        &        & \\
   b[0] & a[1]   & b[1]   &        & \\
    0   & \ddots & \ddots & \ddots &  \\
        &        & \ddots & \ddots & b[n-2] \\
        &        &        & b[n-2] & a[n-1] \end{array} \right ] \f$

   Where the coefficients a[n], b[n] come from the  recurrence relation

   \f$  b_j p_j(z) = (z - a_j ) p_{j-1}(z) - b_{j-1}   p_{j-2}(z) \f$

   where \f$ j=n+1\f$ and \f$p_j(z)\f$ are the Jacobi (normalized)
   orthogonal polynomials \f$ \alpha,\beta > -1\f$( integer values and
   halves). Since the polynomials are orthonormalized, the tridiagonal
   matrix is guaranteed to be symmetric. The eigenvalues of this
   matrix are the zeros of the Jacobi polynomial.
**/

void jac_zeros(const Uint n, DenseDVec<Real> &a, const Real alpha, const Real beta);

// ============================================================================

/** \brief QL algorithm for symmetric tridiagonal matrix

   This subroutine is a translation of an algol procedure,
   num. math. \b 12, 377-383(1968) by Martin and Wilkinson, as modified
   in num. math. \b 15, 450(1970) bu Dubrulle. Handbook for
   auto. comp., vol.ii-linear algebra, 241-248(1971). This is a
   modified version from numerical recipes.

   This subroutine finds the eigenvalues and first components of the
   eigenvectors of a symmetric tridiagonal matrix by the implicit QL
   method.

   on input:
   - n is the order of the matrix;
   - d contains the diagonal elements of the input matrix;
   - e contains the subdiagonal elements of the input matrix
       in its first n-1 positions. e(n) is arbitrary

   on output:
   - d contains the eigenvalues in ascending order
   - e has been destroyed
 **/

void tri_QL(const Uint n, DenseDVec<Real> &d, DenseDVec<Real> &e);

// ============================================================================

} // namespace math

} // namespace pdekit

#endif
