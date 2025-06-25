/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE expression_templates_corectness_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iomanip>
#include <iostream>

/// PDEKIT headers
#include "math/DenseConstMatView.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseMatView.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"
#include "math/DenseVecView.hpp"
#include "math/TensorInitialization.hpp"
#include "math/unary_ops/MatrixNorm.hpp"
#include "math/unary_ops/VectorNorm.hpp"

using namespace pdekit;
using namespace boost::unit_test;

// ----------------------------------------------------------------------------

struct ETCorectnessUtestFixture
{

  enum
  {
    N = 10
  };

  static const Real PI;

  // Common setup for each test case
  ETCorectnessUtestFixture()
  {
    // Resize vectors
    v1_d.resize(N);
    v2_d.resize(N);
    v3_d.resize(N);
    vec_result_d.resize(N);
    vec_ref.resize(N);
    vec_diff.resize(N);

    // Assign initial values to the vectors
    for (Uint i = 0; i < N; ++i)
    {
      v1_d[i] = 10.0 * std::sin(0.1 * PI * i + 1.0);
      v2_d[i] = 5.0 * std::sin(0.2 * PI * i + 0.5);
      v3_d[i] = 7.0 * std::cos(0.2 * PI * i + 0.5);

      v1_s[i] = 10.0 * std::sin(0.1 * PI * i + 1.0);
      v2_s[i] = 5.0 * std::sin(0.2 * PI * i + 0.5);
      v3_s[i] = 7.0 * std::cos(0.2 * PI * i + 0.5);
    }

    // Resize matrices
    M1_d.resize(N, N);
    M2_d.resize(N, N);
    M3_d.resize(N, N);
    mat_result_d.resize(N, N);
    mat_ref.resize(N, N);
    mat_diff.resize(N, N);

    // Assign initial values to matrices
    for (Uint i = 0; i < N; ++i)
      for (Uint j = 0; j < N; ++j)
      {
        M1_d(i, j) = 10.0 * std::sin(0.1 * PI * i + 0.2 * PI * j + 1.0);
        M2_d(i, j) = 5.0 * std::sin(0.2 * PI * i - 0.3 * PI * j + 0.5);
        M3_d(i, j) = 7.0 * std::cos(0.2 * PI * i + 0.2 * PI * j + 0.5);

        M1_s(i, j) = 10.0 * std::sin(0.1 * PI * i + 0.2 * PI * j + 1.0);
        M2_s(i, j) = 5.0 * std::sin(0.2 * PI * i - 0.3 * PI * j + 0.5);
        M3_s(i, j) = 7.0 * std::cos(0.2 * PI * i + 0.2 * PI * j + 0.5);
      }
  }

  // Common tear-down for each test case
  ~ETCorectnessUtestFixture()
  {
  }

  // Member variables
  math::DenseDVec<Real> v1_d, v2_d, v3_d, vec_result_d, vec_ref, vec_diff;
  math::DenseSVec<Real, N> v1_s, v2_s, v3_s, vec_result_s;

  math::DenseDMat<Real> M1_d, M2_d, M3_d, mat_result_d, mat_ref, mat_diff;
  math::DenseSMat<Real, N, N> M1_s, M2_s, M3_s, mat_result_s;
};

const Real ETCorectnessUtestFixture::PI = 3.14159265358979323846;

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(ETCorectness_TestSuite, ETCorectnessUtestFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(check_tensor_ranks)
{
  // We cannot compare enums directly (compiler warning), but instead
  // we need to explicitely cast them to integer values
  const Uint rank1      = math::TensorRank<math::DenseSVec<Real, 10>>::value;
  const Uint reference1 = math::tensor_rank_1;
  static_assert(rank1 == reference1, "Static vector should have rank 1");

  const Uint rank2      = math::TensorRank<math::DenseDVec<Uint>>::value;
  const Uint reference2 = math::tensor_rank_1;
  static_assert(rank2 == reference2, "Dynamic vector should have rank 1");

  const Uint rank3      = math::TensorRank<math::DenseConstVecView<Real>>::value;
  const Uint reference3 = math::tensor_rank_1;
  static_assert(rank3 == reference3, "Vector block should have rank 1");

  const Uint rank4      = math::TensorRank<math::DenseSMat<Real, 2, 3>>::value;
  const Uint reference4 = math::tensor_rank_2;
  static_assert(rank4 == reference4, "Static matrix should have rank 2");

  const Uint rank5      = math::TensorRank<math::DenseDMat<Real>>::value;
  const Uint reference5 = math::tensor_rank_2;
  static_assert(rank5 == reference5, "Dynamic matrix should have rank 2");

  const Uint rank6      = math::TensorRank<math::DenseConstMatView<Real>>::value;
  const Uint reference6 = math::tensor_rank_2;
  static_assert(rank6 == reference6, "Matrix block should have rank 2");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(initialize_vector)
{
  const math::DenseSVec<Real, 9> vs = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(9.);

  math::DenseDVec<Real> vd(9);
  vd = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(9.);

  for (Uint i = 0; i < 9; ++i)
  {
    BOOST_CHECK_EQUAL(vs[i], i + 1);
    BOOST_CHECK_EQUAL(vd[i], i + 1);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(initialize_matrix)
{
  math::DenseSMat<Real, 3, 3> ms = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(9.);

  math::DenseDMat<Real> md(3, 3);
  md = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(9.);

  Uint idx = 1;
  for (Uint i = 0; i < 3; ++i)
  {
    for (Uint j = 0; j < 3; ++j)
    {
      BOOST_CHECK_EQUAL(ms(i, j), idx);
      BOOST_CHECK_EQUAL(md(i, j), idx);
      idx++;
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(const_vector_block_creation)
{
  const math::DenseSMat<Real, 3, 3> ms = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(9.);

  math::DenseDMat<Real> md(3, 3);
  md = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(9.);

  // Get the first row of matrix 'ms'
  math::DenseConstVecView<Real, math::RowVector> vb1 = ms.const_row(0);

  // Get the second column of matrix 'md'
  math::DenseConstVecView<Real, math::ColumnVector> vb2 = md.const_col(1);

  BOOST_CHECK_EQUAL(vb1.size(), 3U);
  BOOST_CHECK_EQUAL(vb1[0], 1);
  BOOST_CHECK_EQUAL(vb1[1], 2);
  BOOST_CHECK_EQUAL(vb1[2], 3);

  BOOST_CHECK_EQUAL(vb2.size(), 3U);
  BOOST_CHECK_EQUAL(vb2[0], 2);
  BOOST_CHECK_EQUAL(vb2[1], 5);
  BOOST_CHECK_EQUAL(vb2[2], 8);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mutable_vector_block_creation)
{
  math::DenseSMat<Real, 3, 3> ms = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(9.);

  math::DenseDMat<Real> md(3, 3);
  md = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(9.);

  // Get the first row of matrix 'ms'
  math::DenseVecView<Real, math::RowVector> vb1 = ms.row(0);

  vb1[0] = 21.;
  // vb1[1] = 22.;
  ms.row(0)[1] = 22.;
  vb1[2]       = 23.;

  BOOST_CHECK_EQUAL(vb1.size(), 3u);
  BOOST_CHECK_EQUAL(vb1[0], 21.);
  BOOST_CHECK_EQUAL(vb1[1], 22.);
  BOOST_CHECK_EQUAL(vb1[2], 23.);

  // Check the first row of matrix ms
  BOOST_CHECK_EQUAL(ms(0, 0), 21.);
  BOOST_CHECK_EQUAL(ms(0, 1), 22.);
  BOOST_CHECK_EQUAL(ms(0, 2), 23.);

  // Check the second row of matrix ms
  BOOST_CHECK_EQUAL(ms(1, 0), 4.);
  BOOST_CHECK_EQUAL(ms(1, 1), 5.);
  BOOST_CHECK_EQUAL(ms(1, 2), 6.);

  // Check the third row of matrix ms
  BOOST_CHECK_EQUAL(ms(2, 0), 7.);
  BOOST_CHECK_EQUAL(ms(2, 1), 8.);
  BOOST_CHECK_EQUAL(ms(2, 2), 9.);

  // Get the second column of matrix 'md'
  math::DenseVecView<Real, math::ColumnVector> vb2 = md.col(1);

  vb2[0] = 31.;
  // vb2[1] = 32.;
  md.col(1)[1] = 32.;
  vb2[2]       = 33.;

  BOOST_CHECK_EQUAL(vb2.size(), 3u);
  BOOST_CHECK_EQUAL(vb2[0], 31.);
  BOOST_CHECK_EQUAL(vb2[1], 32.);
  BOOST_CHECK_EQUAL(vb2[2], 33.);

  // Check the first column of matrix md
  BOOST_CHECK_EQUAL(md(0, 0), 1.);
  BOOST_CHECK_EQUAL(md(1, 0), 4.);
  BOOST_CHECK_EQUAL(md(2, 0), 7.);

  // Check the second column of matrix md
  BOOST_CHECK_EQUAL(md(0, 1), 31.);
  BOOST_CHECK_EQUAL(md(1, 1), 32.);
  BOOST_CHECK_EQUAL(md(2, 1), 33.);

  // Check the third column of matrix md
  BOOST_CHECK_EQUAL(md(0, 2), 3.);
  BOOST_CHECK_EQUAL(md(1, 2), 6.);
  BOOST_CHECK_EQUAL(md(2, 2), 9.);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(const_matrix_block_creation)
{
  math::DenseDMat<Real> md(5, 5);
  md = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(9.)(10.)(11.)(12.)(13.)(14.)(15.)(16.)(
      17.)(18.)(19.)(20.)(21.)(22.)(23.)(24.)(25.);

  math::DenseConstMatView<Real> block = md.const_block(1, 2, 4, 3);

  // Check the dimensions of the block
  BOOST_CHECK_EQUAL(block.rows(), 4U);
  BOOST_CHECK_EQUAL(block.cols(), 3U);

  // Check the entries in the first row of block
  BOOST_CHECK_EQUAL(block(0, 0), 8.);
  BOOST_CHECK_EQUAL(block(0, 1), 9.);
  BOOST_CHECK_EQUAL(block(0, 2), 10.);

  // Check the entries in the second row of block
  BOOST_CHECK_EQUAL(block(1, 0), 13.);
  BOOST_CHECK_EQUAL(block(1, 1), 14.);
  BOOST_CHECK_EQUAL(block(1, 2), 15.);

  // Check the entries in the third row of block
  BOOST_CHECK_EQUAL(block(2, 0), 18.);
  BOOST_CHECK_EQUAL(block(2, 1), 19.);
  BOOST_CHECK_EQUAL(block(2, 2), 20.);

  // Check the entries in the fourth row of block
  BOOST_CHECK_EQUAL(block(3, 0), 23.);
  BOOST_CHECK_EQUAL(block(3, 1), 24.);
  BOOST_CHECK_EQUAL(block(3, 2), 25.);

  math::DenseConstVecView<Real, math::RowVector> row0 = block.row(0);
  BOOST_CHECK_EQUAL(row0[0], 8.0);
  BOOST_CHECK_EQUAL(row0[1], 9.0);
  BOOST_CHECK_EQUAL(row0[2], 10.0);

  math::DenseConstVecView<Real, math::RowVector> row1 = block.row(1);
  BOOST_CHECK_EQUAL(row1[0], 13.0);
  BOOST_CHECK_EQUAL(row1[1], 14.0);
  BOOST_CHECK_EQUAL(row1[2], 15.0);

  math::DenseConstVecView<Real, math::RowVector> row2 = block.row(2);
  BOOST_CHECK_EQUAL(row2[0], 18.0);
  BOOST_CHECK_EQUAL(row2[1], 19.0);
  BOOST_CHECK_EQUAL(row2[2], 20.0);

  math::DenseConstVecView<Real, math::RowVector> row3 = block.row(3);
  BOOST_CHECK_EQUAL(row3[0], 23.0);
  BOOST_CHECK_EQUAL(row3[1], 24.0);
  BOOST_CHECK_EQUAL(row3[2], 25.0);

  math::DenseConstVecView<Real, math::ColumnVector> col0 = block.col(0);
  BOOST_CHECK_EQUAL(col0[0], 8.0);
  BOOST_CHECK_EQUAL(col0[1], 13.0);
  BOOST_CHECK_EQUAL(col0[2], 18.0);
  BOOST_CHECK_EQUAL(col0[3], 23.0);

  math::DenseConstVecView<Real, math::ColumnVector> col1 = block.col(1);
  BOOST_CHECK_EQUAL(col1[0], 9.0);
  BOOST_CHECK_EQUAL(col1[1], 14.0);
  BOOST_CHECK_EQUAL(col1[2], 19.0);
  BOOST_CHECK_EQUAL(col1[3], 24.0);

  math::DenseConstVecView<Real, math::ColumnVector> col2 = block.col(2);
  BOOST_CHECK_EQUAL(col2[0], 10.0);
  BOOST_CHECK_EQUAL(col2[1], 15.0);
  BOOST_CHECK_EQUAL(col2[2], 20.0);
  BOOST_CHECK_EQUAL(col2[3], 25.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mutable_matrix_block_creation)
{
  math::DenseDMat<Real> md(5, 5);
  md = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(9.)(10.)(11.)(12.)(13.)(14.)(15.)(16.)(
      17.)(18.)(19.)(20.)(21.)(22.)(23.)(24.)(25.);

  math::DenseMatView<Real> block = md.block(1, 2, 4, 3);

  block(0, 0) = 108.;
  block(0, 1) = 109.;
  block(0, 2) = 110.;
  block(1, 0) = 113.;
  block(1, 1) = 114.;
  block(1, 2) = 115.;
  block(2, 0) = 118.;
  block(2, 1) = 119.;
  block(2, 2) = 120.;
  block(3, 0) = 123.;
  block(3, 1) = 124.;
  block(3, 2) = 125.;

  // The matrix should look as follows:
  // |   1   2   3    4    5   |
  // |   6   7  108  109  110  |
  // |  11  12  113  114  115  |
  // |  16  17  118  119  120  |
  // |  21  22  123  124  125  |

  // Check the dimensions of the block
  BOOST_CHECK_EQUAL(block.rows(), 4U);
  BOOST_CHECK_EQUAL(block.cols(), 3U);

  // Check the entries in the first row of the original matrix
  BOOST_CHECK_EQUAL(md(0, 0), 1.);
  BOOST_CHECK_EQUAL(md(0, 1), 2.);
  BOOST_CHECK_EQUAL(md(0, 2), 3.);
  BOOST_CHECK_EQUAL(md(0, 3), 4.);
  BOOST_CHECK_EQUAL(md(0, 4), 5.);

  // Check the entries in the second row of the original matrix
  BOOST_CHECK_EQUAL(md(1, 0), 6.);
  BOOST_CHECK_EQUAL(md(1, 1), 7.);
  BOOST_CHECK_EQUAL(md(1, 2), 108.);
  BOOST_CHECK_EQUAL(md(1, 3), 109.);
  BOOST_CHECK_EQUAL(md(1, 4), 110.);

  // Check the entries in the third row of the original matrix
  BOOST_CHECK_EQUAL(md(2, 0), 11.);
  BOOST_CHECK_EQUAL(md(2, 1), 12.);
  BOOST_CHECK_EQUAL(md(2, 2), 113.);
  BOOST_CHECK_EQUAL(md(2, 3), 114.);
  BOOST_CHECK_EQUAL(md(2, 4), 115.);

  // Check the entries in the fourth row of the original matrix
  BOOST_CHECK_EQUAL(md(3, 0), 16.);
  BOOST_CHECK_EQUAL(md(3, 1), 17.);
  BOOST_CHECK_EQUAL(md(3, 2), 118.);
  BOOST_CHECK_EQUAL(md(3, 3), 119.);
  BOOST_CHECK_EQUAL(md(3, 4), 120.);

  // Check the entries in the fifth row of the original matrix
  BOOST_CHECK_EQUAL(md(4, 0), 21.);
  BOOST_CHECK_EQUAL(md(4, 1), 22.);
  BOOST_CHECK_EQUAL(md(4, 2), 123.);
  BOOST_CHECK_EQUAL(md(4, 3), 124.);
  BOOST_CHECK_EQUAL(md(4, 4), 125.);

  math::DenseConstVecView<Real, math::RowVector> row0 = block.row(0);
  BOOST_CHECK_EQUAL(row0[0], 108.0);
  BOOST_CHECK_EQUAL(row0[1], 109.0);
  BOOST_CHECK_EQUAL(row0[2], 110.0);

  math::DenseConstVecView<Real, math::RowVector> row1 = block.row(1);
  BOOST_CHECK_EQUAL(row1[0], 113.0);
  BOOST_CHECK_EQUAL(row1[1], 114.0);
  BOOST_CHECK_EQUAL(row1[2], 115.0);

  math::DenseConstVecView<Real, math::RowVector> row2 = block.row(2);
  BOOST_CHECK_EQUAL(row2[0], 118.0);
  BOOST_CHECK_EQUAL(row2[1], 119.0);
  BOOST_CHECK_EQUAL(row2[2], 120.0);

  math::DenseConstVecView<Real, math::RowVector> row3 = block.row(3);
  BOOST_CHECK_EQUAL(row3[0], 123.0);
  BOOST_CHECK_EQUAL(row3[1], 124.0);
  BOOST_CHECK_EQUAL(row3[2], 125.0);

  math::DenseConstVecView<Real, math::ColumnVector> col0 = block.col(0);
  BOOST_CHECK_EQUAL(col0[0], 108.0);
  BOOST_CHECK_EQUAL(col0[1], 113.0);
  BOOST_CHECK_EQUAL(col0[2], 118.0);
  BOOST_CHECK_EQUAL(col0[3], 123.0);

  math::DenseConstVecView<Real, math::ColumnVector> col1 = block.col(1);
  BOOST_CHECK_EQUAL(col1[0], 109.0);
  BOOST_CHECK_EQUAL(col1[1], 114.0);
  BOOST_CHECK_EQUAL(col1[2], 119.0);
  BOOST_CHECK_EQUAL(col1[3], 124.0);

  math::DenseConstVecView<Real, math::ColumnVector> col2 = block.col(2);
  BOOST_CHECK_EQUAL(col2[0], 110.0);
  BOOST_CHECK_EQUAL(col2[1], 115.0);
  BOOST_CHECK_EQUAL(col2[2], 120.0);
  BOOST_CHECK_EQUAL(col2[3], 125.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(scalar_vector_mult)
{
  for (Uint i = 0; i < N; ++i)
  {
    vec_ref[i] = 3.1 * v1_d[i];
  }

  vec_result_d = 3.1 * v1_d;
  vec_result_s = 3.1 * v1_s;

  for (Uint i = 0; i < N; ++i)
  {
    BOOST_CHECK_EQUAL(vec_result_d[i], vec_ref[i]);
    BOOST_CHECK_EQUAL(vec_result_s[i], vec_ref[i]);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(vector_lin_combination)
{
  for (Uint i = 0; i < N; ++i)
  {
    vec_ref[i] = 3.1 * v1_d[i] - 4.6 * v2_d[i] + 0.01 * v3_d[i];
  }

  vec_result_d = 3.1 * v1_d - 4.6 * v2_d + 0.01 * v3_d;
  vec_result_s = 3.1 * v1_s - 4.6 * v2_s + 0.01 * v3_s;

  for (Uint i = 0; i < N; ++i)
  {
    BOOST_CHECK_EQUAL(vec_result_d[i], vec_ref[i]);
    BOOST_CHECK_EQUAL(vec_result_s[i], vec_ref[i]);
  }

  math::DenseDMat<Real> md(N, N);
  for (Uint r = 0; r < N; ++r)
  {
    for (Uint c = 0; c < N; ++c)
    {
      md(r, c) = r + c;
    }
  }

  // Create a matrix block
  math::DenseMatView<Real> block = md.block(0, 0, N, N);

  // Take two rows slices (ConstVectorBlocks) from the matrix block
  const math::DenseConstVecView<Real> first_row = block.row_transpose(0);
  const math::DenseConstVecView<Real> last_row  = block.row_transpose(N - 1);

  // Accumulate the rows to a dynamic vector
  vec_result_d.fill(0.0);
  vec_result_d += first_row;
  vec_result_d += last_row;

  // Check
  for (Uint i = 0; i < N; ++i)
  {
    BOOST_CHECK_EQUAL(vec_result_d[i], first_row[i] + last_row[i]);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(scalar_matrix_mult)
{
  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      mat_ref(i, j) = 3.1 * M1_d(i, j);
    }
  }

  mat_result_d = 3.1 * M1_d;
  mat_result_s = 3.1 * M1_s;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_d(i, j), mat_ref(i, j));
      BOOST_CHECK_EQUAL(mat_result_s(i, j), mat_ref(i, j));
    }
  }

  math::DenseConstMatView<Real> mat_block = M1_d.const_block(0, 0, N, N);

  mat_result_d = 3.1 * mat_block;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_d(i, j), mat_ref(i, j));
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(matrix_lin_combination)
{
  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      mat_ref(i, j) = 3.1 * M1_d(i, j) + 20.43 * M2_d(i, j) - 50.0 * M3_d(i, j);
    }
  }

  mat_result_d = 3.1 * M1_d + 20.43 * M2_d - 50.0 * M3_d;
  mat_result_s = 3.1 * M1_s + 20.43 * M2_s - 50.0 * M3_s;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_d(i, j), mat_ref(i, j));
      BOOST_CHECK_EQUAL(mat_result_s(i, j), mat_ref(i, j));
    }
  }

  math::DenseConstMatView<Real> mat_block1 = M1_d.const_block(0, 0, N, N);
  math::DenseConstMatView<Real> mat_block2 = M2_d.const_block(0, 0, N, N);
  math::DenseConstMatView<Real> mat_block3 = M3_d.const_block(0, 0, N, N);

  mat_result_d = 3.1 * mat_block1 + 20.43 * mat_block2 - 50.0 * mat_block3;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_d(i, j), mat_ref(i, j));
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(matrix_vector_mult)
{
  vec_ref.fill(0.0);

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      vec_ref[i] += M1_d(i, j) * v2_d[j];
    }
  }

  vec_result_d = M1_d * v2_d;
  vec_result_s = M1_s * v2_s;

  for (Uint i = 0; i < N; ++i)
  {
    BOOST_CHECK_EQUAL(vec_result_d[i], vec_ref[i]);
    BOOST_CHECK_EQUAL(vec_result_s[i], vec_ref[i]);
  }

  math::DenseConstMatView<Real> mat_block = M1_d.const_block(0, 0, N, N);

  // Matrix block x static vector, assigned to dynamic vector
  vec_result_d = mat_block * v2_s;

  for (Uint i = 0; i < N; ++i)
  {
    BOOST_CHECK_EQUAL(vec_result_d[i], vec_ref[i]);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(matrix_matrix_mult)
{
  mat_ref.fill(0.0);

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      for (Uint k = 0; k < N; ++k)
      {
        mat_ref(i, j) += M1_d(i, k) * M3_d(k, j);
      }
    }
  }

  // Multiply two dynamic matrices
  mat_result_d.fill(0.0);
  mat_result_d = M1_d * M3_d;

  // Multiply two static matrices
  mat_result_s.fill(0.0);
  mat_result_s = M1_s * M3_s;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_d(i, j), mat_ref(i, j));
      BOOST_CHECK_EQUAL(mat_result_s(i, j), mat_ref(i, j));
    }
  }

  // Multiply dynamic and static matrix
  mat_result_s.fill(0.0);
  mat_result_s = M1_s * M3_d;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_s(i, j), mat_ref(i, j));
    }
  }

  // Multiply two matrix blocks
  math::DenseConstMatView<Real> M1_block = M1_s.const_block(0, 0, N, N);
  math::DenseConstMatView<Real> M3_block = M3_s.const_block(0, 0, N, N);

  mat_result_s.fill(0.0);
  mat_result_s = M1_block * M3_block;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_s(i, j), mat_ref(i, j));
    }
  }

  // Multiply two matrix blocks, one is constant, the other is not
  math::DenseConstMatView<Real> M1_block_const = M1_d.const_block(0, 0, N, N);
  math::DenseMatView<Real> M3_block_mutable    = M3_d.block(0, 0, N, N);

  mat_result_s.fill(0.0);
  mat_result_s = M1_block_const * M3_block_mutable;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_s(i, j), mat_ref(i, j));
    }
  }

  // Multiply matrix block and dynamic matrix
  mat_result_d.fill(0.0);
  mat_result_d = M1_block * M3_d;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_d(i, j), mat_ref(i, j));
    }
  }

  // Multiply dynamic matrix and matrix block
  mat_result_d.fill(0.0);
  mat_result_d = M1_d * M3_block;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_d(i, j), mat_ref(i, j));
    }
  }

  // Multiply matrix block and static matrix
  mat_result_s.fill(0.0);
  mat_result_s = M1_block * M3_s;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_s(i, j), mat_ref(i, j));
    }
  }

  // Multiply static matrix and matrix block
  mat_result_s.fill(0.0);
  mat_result_s = M1_s * M3_block;

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      BOOST_CHECK_EQUAL(mat_result_s(i, j), mat_ref(i, j));
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(matrix_vector_lin_combination)
{
  vec_ref.fill(0.0);

  for (Uint i = 0; i < N; ++i)
  {
    for (Uint j = 0; j < N; ++j)
    {
      vec_ref[i] = vec_ref[i] + 2.0 * M1_d(i, j) * v1_d[j] + M2_d(i, j) * (v2_d[j] + v3_d[j]);
    }
    vec_ref[i] += v2_d[i];
  }

  // Linear combination of dynamic vectors/matrices
  vec_result_d = 2.0 * M1_d * v1_d + v2_d + M2_d * (v2_d + v3_d);

  // Linear combination of static vectors/matrices
  vec_result_s = 2.0 * M1_s * v1_s + v2_s + M2_s * (v2_s + v3_s);

  for (Uint i = 0; i < N; ++i)
  {
    BOOST_CHECK_CLOSE(vec_result_d[i], vec_ref[i], 1.e-13);
    BOOST_CHECK_CLOSE(vec_result_s[i], vec_ref[i], 1.e-13);
  }

  // Blocks cut from static matrices
  math::DenseConstMatView<Real> M1_block = M1_s.const_block(0, 0, N, N);
  math::DenseConstMatView<Real> M2_block = M2_s.const_block(0, 0, N, N);

  // Linear combination of static vectors with matrices blocks
  vec_result_s.fill(0.0);
  vec_result_s = 2.0 * M1_block * v1_s + v2_s + M2_block * (v2_s + v3_s);

  for (Uint i = 0; i < N; ++i)
  {
    BOOST_CHECK_CLOSE(vec_result_s[i], vec_ref[i], 1.e-12);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(vector_norms)
{

  const Real nmax = norm_max(M1_s * (v1_d + v2_s));
  BOOST_CHECK_CLOSE(nmax, 571.643216410270270, 1.e-13);

  const Real ne2 = norm_e2(M1_d * (v1_s + v2_d));
  BOOST_CHECK_CLOSE(ne2, 1283.856857346584093, 1.e-13);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(matrix_norms)
{

  const Real nmax = norm_max(M1_s * M2_d);
  BOOST_CHECK_CLOSE(nmax, 247.944616925085484, 1.e-13);

  const Real ne2 = norm_e2(M1_d * M2_s);
  BOOST_CHECK_CLOSE(ne2, 1767.766952966368535, 1.e-13);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
