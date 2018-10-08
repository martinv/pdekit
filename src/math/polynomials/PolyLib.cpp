#include "math/polynomials/PolyLib.hpp"

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

void zwgj(DenseDVec<Real> &z, DenseDVec<Real> &w, const Uint np, const Real alpha, const Real beta)
{
  Real fac;
  const Real one = 1.0;
  const Real two = 2.0;
  const Real apb = alpha + beta;

  jac_zeros(np, z, alpha, beta);
  jacobd(np, z, w, np, alpha, beta);

  fac = std::pow(two, apb + one) * gammaF(alpha + np + one) * gammaF(beta + np + one);
  fac /= gammaF(np + one) * gammaF(apb + np + one);

  for (Uint i = 0; i < np; ++i)
  {
    w[i] = fac / (w[i] * w[i] * (one - z[i] * z[i]));
  }

  return;
}

// ============================================================================

void jacobfd(const Uint np, DenseDVec<Real> &z, DenseDVec<Real> &poly_in, DenseDVec<Real> &polyd,
             const Uint n, const Real alpha, const Real beta)
{
  // register int i;
  const Real zero = 0.0;
  const Real one  = 1.0;
  const Real two  = 2.0;

  if (n == 0)
  {
    if (poly_in.size() > 0)
    {
      for (Uint i = 0; i < np; ++i)
      {
        poly_in[i] = one;
      }
    }
    if (polyd.size() > 0)
    {
      for (Uint i = 0; i < np; ++i)
      {
        polyd[i] = zero;
      }
    }
  }
  else if (n == 1)
  {
    if (poly_in.size() > 0)
    {
      for (Uint i = 0; i < np; ++i)
      {
        poly_in[i] = 0.5 * (alpha - beta + (alpha + beta + two) * z[i]);
      }
    }

    if (polyd.size() > 0)
    {
      for (Uint i = 0; i < np; ++i)
      {
        polyd[i] = 0.5 * (alpha + beta + two);
      }
    }
  }
  else
  {
    Real a1, a2, a3, a4;
    const Real apb = alpha + beta;
    DenseDVec<Real> poly(np);
    DenseDVec<Real> polyn1(np);
    DenseDVec<Real> polyn2(np);

    if (poly_in.size() > 0)
    { // switch for case of no polynomial function return
      for (Uint i = 0; i < np; ++i)
      {
        poly[i] = poly_in[i];
      }
    }

    for (Uint i = 0; i < np; ++i)
    {
      polyn2[i] = one;
      polyn1[i] = 0.5 * (alpha - beta + (alpha + beta + two) * z[i]);
    }

    for (Uint k = 2; k <= n; ++k)
    {
      a1 = two * k * (k + apb) * (two * k + apb - two);
      a2 = (two * k + apb - one) * (alpha * alpha - beta * beta);
      a3 = (two * k + apb - two) * (two * k + apb - one) * (two * k + apb);
      a4 = two * (k + alpha - one) * (k + beta - one) * (two * k + apb);

      a2 /= a1;
      a3 /= a1;
      a4 /= a1;

      for (Uint i = 0; i < np; ++i)
      {
        poly[i]   = (a2 + a3 * z[i]) * polyn1[i] - a4 * polyn2[i];
        polyn2[i] = polyn1[i];
        polyn1[i] = poly[i];
      }
    }

    if (polyd.size() > 0)
    {
      a1 = n * (alpha - beta);
      a2 = n * (two * n + alpha + beta);
      a3 = two * (n + alpha) * (n + beta);
      a4 = (two * n + alpha + beta);
      a1 /= a4;
      a2 /= a4;
      a3 /= a4;

      // note polyn2 points to polyn1 at end of poly iterations
      for (Uint i = 0; i < np; ++i)
      {
        polyd[i] = (a1 - a2 * z[i]) * poly[i] + a3 * polyn2[i];
        polyd[i] /= (one - z[i] * z[i]);
      }
    }

    // Copy back into poly_in
    if (poly_in.size() > 0)
    { // switch for case of no polynomial function return
      for (Uint i = 0; i < np; ++i)
      {
        poly_in[i] = poly[i];
      }
    }
  }

  return;
}

// ============================================================================

void jacobd(const Uint np, DenseDVec<Real> &z, DenseDVec<Real> &polyd, const Uint n,
            const Real alpha, const Real beta)
{
  const Real one = 1.0;
  if (n == 0)
  {
    for (Uint i = 0; i < np; ++i)
    {
      polyd[i] = 0.0;
    }
  }
  else
  {
    DenseDVec<Real> dummy(0U);

    // jacobf(np,z,polyd,n-1,alpha+one,beta+one);
    jacobfd(np, z, polyd, dummy, n - 1, alpha + one, beta + one);
    for (Uint i = 0; i < np; ++i)
    {
      polyd[i] *= 0.5 * (alpha + beta + (Real)n + one);
    }
  }
  return;
}

// ============================================================================

Real gammaF(const Real x)
{
  Real gamma = 1.0;

  if (x == -0.5)
  {
    gamma = -2.0 * std::sqrt(PI);
  }
  else if (!x)
  {
    return gamma;
  }
  else if ((x - (Int)x) == 0.5)
  {
    Int n    = (Int)x;
    Real tmp = x;

    gamma = std::sqrt(PI);
    while (n--)
    {
      tmp -= 1.0;
      gamma *= tmp;
    }
  }
  else if ((x - (Int)x) == 0.0)
  {
    Int n    = (Int)x;
    Real tmp = x;

    while (--n)
    {
      tmp -= 1.0;
      gamma *= tmp;
    }
  }
  else
  {
    std::cerr << x << " is not of integer or half order" << std::endl;
  }
  return gamma;
}

// ============================================================================

void jac_zeros(const Uint n, DenseDVec<Real> &a, const Real alpha, const Real beta)
{
  Real apbi;
  DenseDVec<Real> b(n);

  // generate normalised terms
  const Real apb  = alpha + beta;
  const Real a2b2 = beta * beta - alpha * alpha;

  apbi = 2.0 + apb;

  b[n - 1] = std::pow(2.0, apb + 1.0) * gammaF(alpha + 1.0) * gammaF(beta + 1.0) / gammaF(apbi);
  a[0]     = (beta - alpha) / apbi;
  b[0]     = std::sqrt(4.0 * (1.0 + alpha) * (1.0 + beta) / ((apbi + 1.0) * apbi * apbi));

  for (Uint i = 1; i < n - 1; ++i)
  {
    apbi = 2.0 * (i + 1) + apb;
    a[i] = a2b2 / ((apbi - 2.0) * apbi);
    b[i] = std::sqrt(4.0 * (i + 1) * (i + 1 + alpha) * (i + 1 + beta) * (i + 1 + apb) /
                     ((apbi * apbi - 1) * apbi * apbi));
  }

  apbi     = 2.0 * n + apb;
  a[n - 1] = a2b2 / ((apbi - 2.0) * apbi);

  // find eigenvalues
  tri_QL(n, a, b);

  return;
}

// ============================================================================

void tri_QL(const Uint n, DenseDVec<Real> &d, DenseDVec<Real> &e)
{
  const Uint STOP = 30;

  Int m, l, iter, i, k; // These have to be (signed) integers !!!
                        // Otherwise overflow will occur !!!
  const Int N = static_cast<Int>(n);

  Real s, r, p, g, f, dd, c, b;

  for (l = 0; l < N; l++)
  {
    iter = 0;
    do
    {
      for (m = l; m < N - 1; m++)
      {
        dd = std::abs(d[m]) + std::abs(d[m + 1]);
        if (std::abs(e[m]) + dd == dd)
          break;
      }

      if (m != l)
      {
        if (iter++ == STOP)
        {
          std::cerr << "triQL: Too many iterations in TQLI" << std::endl;
          return;
        }

        g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        r = std::sqrt((g * g) + 1.0);

        // sign(a,b) ((b)<0 ? -fabs(a) : fabs(a))
        const Real real_sign = g < 0.0 ? -std::abs(r) : std::abs(r);

        g = d[m] - d[l] + e[l] / (g + real_sign);
        s = c = 1.0;
        p     = 0.0;
        for (i = m - 1; i >= l; i--)
        {
          f = s * e[i];
          b = c * e[i];
          if (std::abs(f) >= std::abs(g))
          {
            c        = g / f;
            r        = std::sqrt((c * c) + 1.0);
            e[i + 1] = f * r;
            c *= (s = 1.0 / r);
          }
          else
          {
            s        = f / g;
            r        = std::sqrt((s * s) + 1.0);
            e[i + 1] = g * r;
            s *= (c = 1.0 / r);
          }

          g        = d[i + 1] - p;
          r        = (d[i] - g) * s + 2.0 * c * b;
          p        = s * r;
          d[i + 1] = g + p;
          g        = c * r - b;
        }
        d[l] = d[l] - p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }

  // order eigenvalues
  for (i = 0; i < N - 1; ++i)
  {
    k = i;
    p = d[i];
    for (l = i + 1; l < N; ++l)
    {
      if (d[l] < p)
      {
        k = l;
        p = d[l];
      }
    }
    d[k] = d[i];
    d[i] = p;
  }
}

// ============================================================================

} // namespace math

} // namespace pdekit
