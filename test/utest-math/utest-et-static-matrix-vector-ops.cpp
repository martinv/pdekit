#include <cmath>
#include <ctime>
#include <iostream>

#include <boost/assign.hpp>

#include "math/DenseConstMatView.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"
#include "math/ScalarConstant.hpp"
#include "math/unary_ops/VectorNorm.hpp"

using namespace pdekit;
using namespace pdekit::math;

int main()
{

  clock_t start, end;
  Real elapsed;
  const Uint N       = 10;
  const Uint n_loops = 1e7;

  // --------------------------------------------------------------------------

  std::cout << std::endl
            << "************ RUNNING VECTOR BINARY OPERATIONS TEST ************" << std::endl;

  /// TESTS FOR ARITHMETIC OPERATIONS WITH VECTORS:

  DenseSVec<Real, N> v1, v2, v3, v4, ref1, diff;

  const Real PI = 3.14159265358979323846;

  for (Uint i = 0; i < N; ++i)
  {
    v1[i] = 10.0 * std::sin(0.1 * PI * i + 1.0);
    v2[i] = 5.0 * std::sin(0.2 * PI * i + 0.5);
    v3[i] = 7.0 * std::cos(0.2 * PI * i + 0.5);
  }

  start = clock();
  for (Uint iter = 0; iter < n_loops; iter++)
  {
    for (Uint i = 0; i < N; ++i)
    {
      ref1[i] = v1[i] + v2[i] + v3[i] + v1[i];
    }
  }
  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "CPU time (vector operations [reference]) = " << elapsed << " s" << std::endl;

  start = clock();

  for (Uint i = 0; i < n_loops; i++)
  {
    v4 = v1 + v2 + v3 + v1;
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (vector operations) = " << elapsed << " s" << std::endl;

  for (Uint i = 0; i < N; ++i)
    diff[i] = v4[i] - ref1[i];

  std::cout << "The norm ||v4 - v_ref|| = " << norm_max(diff) << std::endl;
  std::cout << std::endl;

  // --------------------------------------------------------------------------

  std::cout << "************ RUNNING MATRIX BINARY OPERATIONS TEST ************" << std::endl;

  /// TESTS FOR MATRIX-VECTOR AND MATRIX-MATRIX MULTIPLICATION:

  DenseSMat<Real, N, N> A, B, C;

  for (Uint i = 0; i < N; ++i)
    for (Uint j = 0; j < N; ++j)
    {
      A(i, j) = 10.0 * std::sin(0.1 * PI * i + 0.2 * PI * j + 1.0);
      B(i, j) = 5.0 * std::sin(0.2 * PI * i - 0.3 * PI * j + 0.5);
      C(i, j) = 7.0 * std::cos(0.2 * PI * i + 0.2 * PI * j + 0.5);
    }

  // Check for matrix-vector multiplication:
  ref1.fill(0.0);

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      ref1[i] = ref1[i] + A(i, j) * v1[j] + B(i, j) * (v2[j] + v3[j]);
    }
    ref1[i] += v2[i];
  }

  start = clock();
  for (Uint i = 0; i < n_loops; ++i)
  {
    v4 = A * v1 + v2 + B * (v2 + v3);
  }
  end = clock();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (matrix operations) = " << elapsed << " s" << std::endl;

  for (Uint i = 0; i < N; ++i)
    diff[i] = v4[i] - ref1[i];

  std::cout << "The norm ||v4 - v_ref|| = " << norm_max(diff) << std::endl;

  std::cout << std::endl;

  return 0;
}

// ----------------------------------------------------------------------------
