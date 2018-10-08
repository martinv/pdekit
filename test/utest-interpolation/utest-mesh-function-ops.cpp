/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE mesh_function_algebraic_ops_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(scalar_mesh_function_utest)
{
  interpolation::ScalarMeshFunction<Real> f1("space", "f1"), f2("space", "f2"), f3("space", "f3");

  const Uint M = 10;

  f1.resize(M);
  f2.resize(M);
  f3.resize(M);

  for (Uint i = 0; i < M; ++i)
  {
    f1[i] = 2. * i;
    f2[i] = 3. * i;
  }

  // --------------
  // Unary addition
  // --------------
  f3 = -f1;

  for (Uint i = 0; i < M; ++i)
  {
    BOOST_CHECK_EQUAL(f3[i], -f1[i]);
  }

  // ---------------------------
  // Assignment of an expression
  // ---------------------------
  f3 = f1 + f2;

  for (Uint i = 0; i < M; ++i)
  {
    BOOST_CHECK_EQUAL(f3[i], 5. * i);
  }

  // -----------------------------
  // Accumulation of an expression
  // -----------------------------
  f3 += (f1 + f2);

  for (Uint i = 0; i < M; ++i)
  {
    BOOST_CHECK_EQUAL(f3[i], 10. * i);
  }

  // -----------------------------
  // Subtraction of an expression
  // -----------------------------
  f3 -= f2;

  for (Uint i = 0; i < M; ++i)
  {
    BOOST_CHECK_EQUAL(f3[i], 7. * i);
  }

  f3 -= (f1 - f2);

  for (Uint i = 0; i < M; ++i)
  {
    BOOST_CHECK_EQUAL(f3[i], 8. * i);
  }

  // ---------------------------------
  // Multiplication by a scalar factor
  // ---------------------------------
  f3 = 2. * (3. * f1 - 4. * f2);

  for (Uint i = 0; i < M; ++i)
  {
    BOOST_CHECK_EQUAL(f3[i], -12. * i);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(vector_mesh_function_utest)
{
  interpolation::VectorMeshFunction<Real> f1("space", "f1"), f2("space", "f2"), f3("space", "f3");
  typedef interpolation::VectorMeshFunction<Real>::entry_type entry_type;
  typedef interpolation::VectorMeshFunction<Real>::const_entry_type const_entry_type;

  const Uint nb_fields  = 3;
  const Uint nb_entries = 10;

  f1.resize(nb_fields, nb_entries);
  f2.resize(nb_fields, nb_entries);
  f3.resize(nb_fields, nb_entries);

  for (Uint i = 0; i < nb_entries; ++i)
  {
    entry_type node_value1 = f1.value(i);
    entry_type node_value2 = f2.value(i);
    for (Uint j = 0; j < nb_fields; ++j)
    {
      node_value1[j] = 2. * i;
      node_value2[j] = 3. * i;
    }
  }

  // --------------
  // Unary addition
  // --------------
  f3 = -1. * f1;

  for (Uint i = 0; i < nb_entries; ++i)
  {
    const const_entry_type node_value3 = f3.const_value(i);
    const const_entry_type node_value1 = f1.const_value(i);

    for (Uint j = 0; j < nb_fields; ++j)
    {
      BOOST_CHECK_EQUAL(node_value3[j], -1. * node_value1[j]);
    }
  }

  // ---------------------------
  // Assignment of an expression
  // ---------------------------
  f3 = f1 + f2;

  for (Uint i = 0; i < nb_entries; ++i)
  {
    const const_entry_type node_value1 = f1.const_value(i);
    const const_entry_type node_value2 = f2.const_value(i);
    const const_entry_type node_value3 = f3.const_value(i);

    for (Uint j = 0; j < nb_fields; ++j)
    {
      BOOST_CHECK_EQUAL(node_value3[j], node_value1[j] + node_value2[j]);
    }
  }

  // -----------------------------
  // Accumulation of an expression
  // -----------------------------
  f3 += (f1 + f2);

  for (Uint i = 0; i < nb_entries; ++i)
  {
    const const_entry_type node_value1 = f1.const_value(i);
    const const_entry_type node_value2 = f2.const_value(i);
    const const_entry_type node_value3 = f3.const_value(i);

    for (Uint j = 0; j < nb_fields; ++j)
    {
      BOOST_CHECK_EQUAL(node_value3[j], 2. * (node_value1[j] + node_value2[j]));
    }
  }

  // -----------------------------
  // Subtraction of an expression
  // -----------------------------
  f3 -= f2;

  for (Uint i = 0; i < nb_entries; ++i)
  {
    const const_entry_type node_value1 = f1.const_value(i);
    const const_entry_type node_value2 = f2.const_value(i);
    const const_entry_type node_value3 = f3.const_value(i);

    for (Uint j = 0; j < nb_fields; ++j)
    {
      BOOST_CHECK_EQUAL(node_value3[j], 2. * node_value1[j] + node_value2[j]);
    }
  }

  f3 -= (f1 - f2);

  for (Uint i = 0; i < nb_entries; ++i)
  {
    const const_entry_type node_value1 = f1.const_value(i);
    const const_entry_type node_value2 = f2.const_value(i);
    const const_entry_type node_value3 = f3.const_value(i);

    for (Uint j = 0; j < nb_fields; ++j)
    {
      BOOST_CHECK_EQUAL(node_value3[j], node_value1[j] + 2. * node_value2[j]);
    }
  }

  // ---------------------------------
  // Multiplication by a scalar factor
  // ---------------------------------
  f3 = 2. * (3. * f1 - 4. * f2);

  for (Uint i = 0; i < nb_entries; ++i)
  {
    const const_entry_type node_value1 = f1.const_value(i);
    const const_entry_type node_value2 = f2.const_value(i);
    const const_entry_type node_value3 = f3.const_value(i);

    for (Uint j = 0; j < nb_fields; ++j)
    {
      BOOST_CHECK_EQUAL(node_value3[j], 6. * node_value1[j] - 8. * node_value2[j]);
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(scalar_mesh_function_vector_mesh_function_utest)
{
  typedef interpolation::VectorMeshFunction<Real>::entry_type entry_type;
  typedef interpolation::VectorMeshFunction<Real>::const_entry_type const_entry_type;
  interpolation::VectorMeshFunction<Real> f1_vec("space", "f1_vec"), f2_vec("space", "f2_vec");
  interpolation::ScalarMeshFunction<Real> f1_scal("space", "f1_scal");

  const Uint nb_fields  = 3;
  const Uint nb_entries = 10;

  f1_vec.resize(nb_fields, nb_entries);
  f2_vec.resize(nb_fields, nb_entries);
  f1_scal.resize(nb_entries);

  for (Uint i = 0; i < nb_entries; ++i)
  {
    entry_type node_value1 = f1_vec.value(i);
    for (Uint j = 0; j < nb_fields; ++j)
    {
      node_value1[j] = 2. * i;
    }
  }

  f1_scal.fill(3.0);

  // --------------------------------------------
  // Element-wise multiplication between a scalar
  // and vector function
  // --------------------------------------------
  f2_vec = f1_scal * f1_vec;

  for (Uint i = 0; i < nb_entries; ++i)
  {
    const const_entry_type node_value1 = f1_vec.const_value(i);
    const const_entry_type node_value2 = f2_vec.const_value(i);

    for (Uint j = 0; j < nb_fields; ++j)
    {
      BOOST_CHECK_EQUAL(node_value2[j], f1_scal[i] * node_value1[j]);
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mesh_function_norms_utest)
{
  typedef interpolation::VectorMeshFunction<Real>::entry_type entry_type;

  interpolation::ScalarMeshFunction<Real> f_scal("space", "f_scal");
  interpolation::VectorMeshFunction<Real> f_vect("space", "f_vect");

  const Uint nb_entries = 10;
  const Uint nb_fields  = 2;

  f_scal.resize(nb_entries);
  f_vect.resize(nb_fields, nb_entries);

  Real norm_L1_ref = 0.0;
  Real norm_L2_ref = 0.0;

  for (Uint e = 0; e < nb_entries; ++e)
  {
    entry_type node_value = f_vect.value(e);

    f_scal[e] = 2. * e;

    norm_L1_ref += std::abs(2. * e);
    norm_L2_ref += (2. * e) * (2. * e);

    for (Uint f = 0; f < nb_fields; ++f)
    {
      node_value[f] = 2. * e;
    }
  }

  norm_L1_ref /= nb_entries;
  norm_L2_ref = std::sqrt(norm_L2_ref / nb_entries);

  math::DenseDVec<Real> norm;

  norm_L1(f_scal, norm);
  BOOST_CHECK_EQUAL(norm.size(), 1u);
  BOOST_CHECK_EQUAL(norm[0], norm_L1_ref);
  // std::cout << "Norm L1 (scalar) = " << norm << std::endl;

  norm_L2(f_scal, norm);
  BOOST_CHECK_EQUAL(norm.size(), 1u);
  BOOST_CHECK_EQUAL(norm[0], norm_L2_ref);
  // std::cout << "Norm L2 (scalar) = " << norm << std::endl;

  norm_L1(f_vect, norm);
  BOOST_CHECK_EQUAL(norm.size(), nb_fields);
  for (Uint f = 0; f < nb_fields; ++f)
  {
    BOOST_CHECK_EQUAL(norm[f], norm_L1_ref);
  }
  // std::cout << "Norm L1 (vector) = " << norm << std::endl;

  norm_L2(f_vect, norm);
  BOOST_CHECK_EQUAL(norm.size(), nb_fields);
  for (Uint f = 0; f < nb_fields; ++f)
  {
    BOOST_CHECK_EQUAL(norm[f], norm_L2_ref);
  }
  // std::cout << "Norm L2 (vector) = " << norm << std::endl;
}

// ----------------------------------------------------------------------------
