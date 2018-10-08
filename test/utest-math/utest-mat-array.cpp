/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE matrix_array_utest
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>

/// PDEKIT headers
#include "math/DenseDMatArray.hpp"
#include "math/DenseDMatArrayUniform.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dense_dynamic_matrix_array_utest)
{
  math::DenseDMatArray<Real> mat_array;
  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> mat_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>);
  mat_shapes->resize(2);
  (*mat_shapes)[0] = common::ArrayShape<_2D, SUint>(3, 2);
  (*mat_shapes)[1] = common::ArrayShape<_2D, SUint>(3, 3);

  mat_array.allocate(std::move(mat_shapes));

  // mat_array.debug_print();

  math::DenseMatView<Real> mat0a = mat_array.mat_view(0);
  math::DenseMatView<Real> mat1a = mat_array.mat_view(1);

  mat0a(0, 0) = 101.0;
  mat0a(0, 1) = 102.0;
  mat0a(1, 0) = 103.0;
  mat0a(1, 1) = 104.0;
  mat0a(2, 0) = 105.0;
  mat0a(2, 1) = 106.0;

  mat1a(0, 0) = 1.0;
  mat1a(0, 1) = 2.0;
  mat1a(0, 2) = 3.0;
  mat1a(1, 0) = 4.0;
  mat1a(1, 1) = 5.0;
  mat1a(1, 2) = 6.0;
  mat1a(2, 0) = 7.0;
  mat1a(2, 1) = 8.0;
  mat1a(2, 2) = 9.0;

  math::DenseConstMatView<Real> mat0b = mat_array.const_mat_view(0);
  math::DenseConstMatView<Real> mat1b = mat_array.const_mat_view(1);

  // std::cout << mat0b << std::endl;
  // std::cout << mat1b << std::endl;

  // Test copy constructor
  math::DenseDMatArray<Real> mat_array2 = mat_array;

  // Test assignment operator
  math::DenseDMatArray<Real> mat_array3;
  mat_array3 = mat_array;

  math::DenseConstMatView<Real> mat2_0 = mat_array2.const_mat_view(0);
  math::DenseConstMatView<Real> mat2_1 = mat_array2.const_mat_view(1);

  math::DenseConstMatView<Real> mat3_0 = mat_array3.const_mat_view(0);
  math::DenseConstMatView<Real> mat3_1 = mat_array3.const_mat_view(1);

  BOOST_CHECK_EQUAL(mat2_0.rows(), 3u);
  BOOST_CHECK_EQUAL(mat2_0.cols(), 2u);
  BOOST_CHECK_EQUAL(mat2_1.rows(), 3u);
  BOOST_CHECK_EQUAL(mat2_1.cols(), 3u);

  BOOST_CHECK_EQUAL(mat3_0.rows(), 3u);
  BOOST_CHECK_EQUAL(mat3_0.cols(), 2u);
  BOOST_CHECK_EQUAL(mat3_1.rows(), 3u);
  BOOST_CHECK_EQUAL(mat3_1.cols(), 3u);

  // Check the values in matrices
  // First matrix
  BOOST_CHECK_EQUAL(mat2_0(0, 0), 101.0);
  BOOST_CHECK_EQUAL(mat2_0(0, 1), 102.0);
  BOOST_CHECK_EQUAL(mat2_0(1, 0), 103.0);
  BOOST_CHECK_EQUAL(mat2_0(1, 1), 104.0);
  BOOST_CHECK_EQUAL(mat2_0(2, 0), 105.0);
  BOOST_CHECK_EQUAL(mat2_0(2, 1), 106.0);

  BOOST_CHECK_EQUAL(mat3_0(0, 0), 101.0);
  BOOST_CHECK_EQUAL(mat3_0(0, 1), 102.0);
  BOOST_CHECK_EQUAL(mat3_0(1, 0), 103.0);
  BOOST_CHECK_EQUAL(mat3_0(1, 1), 104.0);
  BOOST_CHECK_EQUAL(mat3_0(2, 0), 105.0);
  BOOST_CHECK_EQUAL(mat3_0(2, 1), 106.0);

  // Second matrix
  BOOST_CHECK_EQUAL(mat2_1(0, 0), 1.0);
  BOOST_CHECK_EQUAL(mat2_1(0, 1), 2.0);
  BOOST_CHECK_EQUAL(mat2_1(0, 2), 3.0);
  BOOST_CHECK_EQUAL(mat2_1(1, 0), 4.0);
  BOOST_CHECK_EQUAL(mat2_1(1, 1), 5.0);
  BOOST_CHECK_EQUAL(mat2_1(1, 2), 6.0);
  BOOST_CHECK_EQUAL(mat2_1(2, 0), 7.0);
  BOOST_CHECK_EQUAL(mat2_1(2, 1), 8.0);
  BOOST_CHECK_EQUAL(mat2_1(2, 2), 9.0);

  BOOST_CHECK_EQUAL(mat3_1(0, 0), 1.0);
  BOOST_CHECK_EQUAL(mat3_1(0, 1), 2.0);
  BOOST_CHECK_EQUAL(mat3_1(0, 2), 3.0);
  BOOST_CHECK_EQUAL(mat3_1(1, 0), 4.0);
  BOOST_CHECK_EQUAL(mat3_1(1, 1), 5.0);
  BOOST_CHECK_EQUAL(mat3_1(1, 2), 6.0);
  BOOST_CHECK_EQUAL(mat3_1(2, 0), 7.0);
  BOOST_CHECK_EQUAL(mat3_1(2, 1), 8.0);
  BOOST_CHECK_EQUAL(mat3_1(2, 2), 9.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(uniform_dynamic_matrix_array_utest)
{
  math::DenseDMatArrayUniform<Real> mat_array;
  common::ArrayShape<_2D, SUint> mat_shape(2, 2);

  mat_array.allocate(mat_shape, 2);

  // mat_array.debug_print();

  math::DenseMatView<Real> mat0a = mat_array.mat_view(0);
  math::DenseMatView<Real> mat1a = mat_array.mat_view(1);

  mat0a(0, 0) = 101.0;
  mat0a(0, 1) = 102.0;
  mat0a(1, 0) = 103.0;
  mat0a(1, 1) = 104.0;

  mat1a(0, 0) = 1.0;
  mat1a(0, 1) = 2.0;
  mat1a(1, 0) = 4.0;
  mat1a(1, 1) = 5.0;

  math::DenseConstMatView<Real> mat0b = mat_array.const_mat_view(0);
  math::DenseConstMatView<Real> mat1b = mat_array.const_mat_view(1);

  // std::cout << mat0b << std::endl;
  // std::cout << mat1b << std::endl;

  // Test copy constructor
  math::DenseDMatArrayUniform<Real> mat_array2 = mat_array;

  // Test assignment operator
  math::DenseDMatArrayUniform<Real> mat_array3;
  mat_array3 = mat_array;

  math::DenseConstMatView<Real> mat2_0 = mat_array2.const_mat_view(0);
  math::DenseConstMatView<Real> mat2_1 = mat_array2.const_mat_view(1);

  math::DenseConstMatView<Real> mat3_0 = mat_array3.const_mat_view(0);
  math::DenseConstMatView<Real> mat3_1 = mat_array3.const_mat_view(1);

  BOOST_CHECK_EQUAL(mat2_0.rows(), 2u);
  BOOST_CHECK_EQUAL(mat2_0.cols(), 2u);
  BOOST_CHECK_EQUAL(mat2_1.rows(), 2u);
  BOOST_CHECK_EQUAL(mat2_1.cols(), 2u);

  BOOST_CHECK_EQUAL(mat3_0.rows(), 2u);
  BOOST_CHECK_EQUAL(mat3_0.cols(), 2u);
  BOOST_CHECK_EQUAL(mat3_1.rows(), 2u);
  BOOST_CHECK_EQUAL(mat3_1.cols(), 2u);

  // Check the values in matrices
  // First matrix
  BOOST_CHECK_EQUAL(mat2_0(0, 0), 101.0);
  BOOST_CHECK_EQUAL(mat2_0(0, 1), 102.0);
  BOOST_CHECK_EQUAL(mat2_0(1, 0), 103.0);
  BOOST_CHECK_EQUAL(mat2_0(1, 1), 104.0);

  BOOST_CHECK_EQUAL(mat3_0(0, 0), 101.0);
  BOOST_CHECK_EQUAL(mat3_0(0, 1), 102.0);
  BOOST_CHECK_EQUAL(mat3_0(1, 0), 103.0);
  BOOST_CHECK_EQUAL(mat3_0(1, 1), 104.0);

  // Second matrix
  BOOST_CHECK_EQUAL(mat2_1(0, 0), 1.0);
  BOOST_CHECK_EQUAL(mat2_1(0, 1), 2.0);
  BOOST_CHECK_EQUAL(mat2_1(1, 0), 4.0);
  BOOST_CHECK_EQUAL(mat2_1(1, 1), 5.0);

  BOOST_CHECK_EQUAL(mat3_1(0, 0), 1.0);
  BOOST_CHECK_EQUAL(mat3_1(0, 1), 2.0);
  BOOST_CHECK_EQUAL(mat3_1(1, 0), 4.0);
  BOOST_CHECK_EQUAL(mat3_1(1, 1), 5.0);
}

// ----------------------------------------------------------------------------
