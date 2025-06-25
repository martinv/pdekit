#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>

#include <boost/assign.hpp>

#include "math/DenseConstMatView.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "math/ScalarConstant.hpp"
#include "math/unary_ops/VectorNorm.hpp"

using namespace pdekit;
using namespace pdekit::math;

int main()
{
  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::milliseconds elapsed;

  const Uint N       = 10;
  const Uint n_loops = 1e7;

  // --------------------------------------------------------------------------

  std::cout << std::endl
            << "************ RUNNING VECTOR BINARY OPERATIONS TEST ************" << std::endl;

  /// TESTS FOR ARITHMETIC OPERATIONS WITH VECTORS:

  DenseDVec<Real> v1(N), v2(N), v3(N), v4(N), ref1(N), diff(N);

  const Real PI = 3.14159265358979323846;

  for (Uint i = 0; i < N; ++i)
  {
    v1[i] = 10.0 * std::sin(0.1 * PI * i + 1.0);
    v2[i] = 5.0 * std::sin(0.2 * PI * i + 0.5);
    v3[i] = 7.0 * std::cos(0.2 * PI * i + 0.5);
  }

  start = std::chrono::high_resolution_clock::now();
  for (Uint iter = 0; iter < n_loops; iter++)
  {
    for (Uint i = 0; i < N; ++i)
    {
      ref1[i] = v1[i] - v2[i] + v3[i] + 1.5 * v1[i];
    }
  }
  end     = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  std::cout.setf(std::ios::fixed);
  std::cout.precision(7);
  std::cout << "CPU time (vector operations [reference]) = " << elapsed.count() << " ms"
            << std::endl;

  start = std::chrono::high_resolution_clock::now();

  for (Uint i = 0; i < n_loops; i++)
  {
    v4 = v1 - v2 + v3 + 1.5 * v1;
  }

  end = std::chrono::high_resolution_clock::now();

  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "CPU time (vector operations) = " << elapsed.count() << " ms" << std::endl;

  for (Uint i = 0; i < N; ++i)
    diff[i] = v4[i] - ref1[i];

  std::cout << "The norm ||v4 - v_ref|| = " << norm_max(diff) << std::endl;
  std::cout << std::endl;

  // --------------------------------------------------------------------------

  std::cout << "************ RUNNING MATRIX BINARY OPERATIONS TEST ************" << std::endl;

  /// TESTS FOR MATRIX-VECTOR AND MATRIX-MATRIX MULTIPLICATION:

  DenseDMat<Real> A(N, N), B(N, N), C(N, N);

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
      ref1[i] = ref1[i] + 2.0 * A(i, j) * v1[j] + B(i, j) * (v2[j] + v3[j]);
    }
    ref1[i] += v2[i];
  }

  start = std::chrono::high_resolution_clock::now();

  for (Uint i = 0; i < n_loops; ++i)
  {
    v4 = 2.0 * A * v1 + v2 + B * (v2 + v3);
  }

  end = std::chrono::high_resolution_clock::now();

  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "CPU time (matrix operations) = " << elapsed.count() << " ms" << std::endl;

  for (Uint i = 0; i < N; ++i)
    diff[i] = v4[i] - ref1[i];

  std::cout << "The norm ||v4 - v_ref|| = " << norm_max(diff) << std::endl;

  std::cout << std::endl;

  // --------------------------------------------------------------------------

  std::cout << "************ RUNNING SLICE BINARY OPERATIONS TEST *************" << std::endl;

  DenseConstVecView<Real> v1_slice(&v1[0], &v1[N - 1]);
  DenseConstVecView<Real> v2_slice(&v2[0], &v2[N - 1]);
  DenseConstVecView<Real> v3_slice(&v3[0], &v3[N - 1]);

  // VectorView<Real> v1_slice(v1.data(), v1.data()+N);
  // VectorView<Real> v2_slice(v2.data(), v2.data()+N);
  // VectorView<Real> v3_slice(v3.data(), v3.data()+N);

  for (Uint i = 0; i < N; ++i)
  {
    ref1[i] = v1[i] + v2[i] + v3[i] + v1[i];
  }

  start = std::chrono::high_resolution_clock::now();

  for (Uint i = 0; i < n_loops; i++)
  {
    v4 = v1_slice + v2_slice + v3_slice + v1_slice;
  }

  end = std::chrono::high_resolution_clock::now();

  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "CPU time (vector slice operations) = " << elapsed.count() << " ms" << std::endl;

  for (Uint i = 0; i < N; ++i)
    diff[i] = v4[i] - ref1[i];

  std::cout << "The norm ||v4 - v_ref|| = " << norm_max(diff) << std::endl;
  std::cout << std::endl;

  /// -------------------------------------------------------------------------
  /// CHECK THE TIMING DIFFERENCE BETWEEN REPEATED MULTIPLICATION OF SMALL
  /// MATRICES
  /// AND ONE MULTIPLICATION OF BIG MATRICES
  /// -------------------------------------------------------------------------

  std::cout << "************ RUNNING TIMINGS FOR EFFICIENT ASSEMBLY *************" << std::endl;

  DenseDMat<Real> V, B_small, C_small, B_big, C_big;

  const Uint nb_qd_pt        = 7;
  const Uint nb_eqns         = 5;
  const Uint nb_dof_per_elem = 15;
  const Uint nb_elem         = 1500;

  V.resize(nb_qd_pt, nb_dof_per_elem);
  V.fill(1.0);

  B_small.resize(nb_dof_per_elem, nb_eqns);
  C_small.resize(nb_qd_pt, nb_eqns);

  B_small.fill(0.1);

  B_big.resize(nb_dof_per_elem, nb_eqns * nb_elem);
  C_big.resize(nb_qd_pt, nb_eqns * nb_elem);

  B_big.fill(0.1);

  start = std::chrono::high_resolution_clock::now();

  for (Uint e = 0; e < nb_elem; ++e)
  {
    C_small = V * B_small;
  }

  end     = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "CPU time (assembly element by element) = " << elapsed.count() << " ms "
            << std::endl;

  start = std::chrono::high_resolution_clock::now();

  C_big = V * B_big;

  end     = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "CPU time (assembly with multiple elements) = " << elapsed.count() << " ms"
            << std::endl;

  return 0;
}

// ----------------------------------------------------------------------------
