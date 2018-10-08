/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Carnevali_expansion_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iomanip>
#include <iostream>

/// PDEKIT headers
#include "common/Constants.hpp"
#include "math/DenseConstVecView.hpp"
#include "mesh/shape_function/CarnevaliExpansionTriag.hpp"
#include "mesh/shape_function/DubinerExpansionTetra.hpp"
#include "mesh/std_region/StdRegionFactory.hpp"
#include "test/utest-mesh/shape_function/fixture/ReferenceShapeFunctionsTetra.hpp"
#include "test/utest-mesh/shape_function/fixture/ReferenceShapeFunctionsTriag.hpp"

using namespace pdekit;
using namespace boost::unit_test;

struct CarnevaliBasis_Fixture
{
  /// common setup for each test case
  CarnevaliBasis_Fixture() : ref_topology_factory(StdRegionFactory::instance())
  {
  }

  /// common tear-down for each test case
  ~CarnevaliBasis_Fixture()
  {
  }

  // --------------------------------------------------------------------------

  template <typename ModalBasisType, typename NodalBasisType>
  void check_modal_basis(const std::string &ref_topology_name, const math::DenseDVec<Real> &point,
                         const Real eps, const bool verbose = false)
  {
    // set up some matrices to hold results
    ModalBasisType basis(NodalBasisType::poly_order);

    // point.resize(2);
    // point[KSI] = -0.7;
    // point[ETA] = -0.4;

    reference_sf_values.resize(basis.nb_modes());
    computed_sf_values.resize(basis.nb_modes());
    reference_sf_derivatives.resize(basis.nb_modes(), basis.topo_dim());
    computed_sf_derivatives.resize(basis.nb_modes(), basis.topo_dim());

    // get the reference topology
    const StdRegionFactory::instance_type::const_product_base_ptr rt =
        ref_topology_factory.create(mesh::PointSetTag::string_to_tag(ref_topology_name));

    rt->coordinates(rt_coords);

    basis.Vandermonde_matrix(rt_coords, V);
    V.transpose_in_place();

    // Inverse of Vandermonde matrix
    V_inv.resize(V.rows(), V.cols());
    V.inv(V_inv);

    UnitM         = V * V_inv;
    Real inf_norm = 0.0;
    for (Uint i = 0; i < UnitM.rows(); ++i)
    {
      for (Uint j = 0; j < UnitM.cols(); ++j)
      {
        inf_norm += std::abs(UnitM(i, j));
      }
    }

    /// Check that the off-diagonal entries of the product
    /// V * V_inv are almost zero
    BOOST_CHECK_LE(std::abs(inf_norm - UnitM.rows()), eps);

    // V_inv.transpose_in_place();

    // Check the result by hard-coded function
    NodalBasisType::eval(point, reference_sf_values);

    // Compute the shape function values
    basis.evaluate_in_one_point(point, mode_values);
    computed_sf_values = V_inv * mode_values;

    if (verbose)
    {
      std::cout << "[" << ref_topology_name
                << "] Reference shape function values = " << reference_sf_values << std::endl;
      std::cout << "[" << ref_topology_name
                << "] Computed shape function values  = " << computed_sf_values << std::endl;

      Real sum_values = 0.0;

      for (Uint i = 0; i < reference_sf_values.size(); ++i)
      {
        sum_values += reference_sf_values[i];
      }
      std::cout << "[" << ref_topology_name
                << "] Sum of shape function values (reference) : " << sum_values << std::endl;

      sum_values = 0.0;
      for (Uint i = 0; i < computed_sf_values.size(); ++i)
      {
        sum_values += computed_sf_values[i];
      }
      std::cout << "[" << ref_topology_name
                << "] Sum of shape function values (computed) : " << sum_values << std::endl
                << std::endl;
    }

    inf_norm = 0.0;
    for (Uint i = 0; i < reference_sf_values.size(); ++i)
    {
      inf_norm += std::abs(reference_sf_values[i] - computed_sf_values[i]);
    }
    /// Make the reference and computed values are the same
    BOOST_CHECK_LE(inf_norm, eps);

    /// Check the derivatives
    NodalBasisType::eval_deriv(point, reference_sf_derivatives);

    basis.evaluate_derivatives_in_one_point(point, mode_values, mode_derivatives);
    computed_sf_derivatives = V_inv * mode_derivatives;

    for (Uint d = 0; d < basis.topo_dim(); ++d)
    {
      inf_norm = 0.0;
      for (Uint i = 0; i < reference_sf_derivatives.rows(); ++i)
      {
        inf_norm += std::abs(reference_sf_derivatives(i, d) - computed_sf_derivatives(i, d));
      }
      /// Make the reference and computed values are the same
      BOOST_CHECK_LE(inf_norm, eps);
    }

    if (verbose)
    {
      std::string var_names[3] = {"xi", "eta", "zeta"};

      for (Uint d = 0; d < basis.topo_dim(); ++d)
      {
        std::cout << std::endl;
        std::cout << "[" << ref_topology_name << "] The derivatives with respect to "
                  << var_names[d] << " (reference) : " << std::endl;
        std::cout << reference_sf_derivatives.const_col(d) << std::endl;

        std::cout << "[" << ref_topology_name << "] The derivatives with respect to "
                  << var_names[d] << " (computed) : " << std::endl;
        // result = V_inv * mode_derivatives.col(0);
        std::cout << computed_sf_derivatives.const_col(d) << std::endl << std::endl;

        Real sum_deriv = 0.0;

        for (Uint i = 0; i < reference_sf_derivatives.rows(); ++i)
        {
          sum_deriv += reference_sf_derivatives(i, d);
        }
        std::cout << "[" << ref_topology_name << "] Sum of derivatives with respect to "
                  << var_names[d] << " (reference) : " << sum_deriv << std::endl;

        sum_deriv = 0.0;
        for (Uint i = 0; i < computed_sf_derivatives.rows(); ++i)
        {
          sum_deriv += computed_sf_derivatives(i, d);
        }
        std::cout << "[" << ref_topology_name << "] Sum of derivatives with respect to "
                  << var_names[d] << " (computed) : " << sum_deriv << std::endl
                  << std::endl;
      }
    }
  }

  // ----------------------------------------------------------------------------

  /// Member variables of fixture

  StdRegionFactory::instance_type &ref_topology_factory;
  math::DenseDMat<Real> rt_coords;

  /// Vandermonde matrix, its inverse and a matrix which is obtained
  /// as V * V_inv (should be unit matrix)
  math::DenseDMat<Real> V, V_inv, UnitM;

  /// Point in which the basis values and their derivatives are
  /// computed
  math::DenseDVec<Real> point;

  /// Values of modal basis computed in 'point'
  math::DenseDVec<Real> mode_values;

  /// Derivatives of modal basis in one point
  math::DenseDMat<Real> mode_derivatives;

  /// Reference values of shape functions or derivatives
  /// This should be the correct value against which we check
  /// data computed using orthogonal basis
  math::DenseDVec<Real> reference_sf_values;

  /// Result computed using Dubiner basis
  math::DenseDVec<Real> computed_sf_values;

  /// Reference values of sf derivatives
  math::DenseDMat<Real> reference_sf_derivatives;

  /// Values of sf derivatives computed with orthogonal basis
  math::DenseDMat<Real> computed_sf_derivatives;
};

BOOST_FIXTURE_TEST_SUITE(CarnevaliBasis_TestSuite, CarnevaliBasis_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Carnevali_modal_basis_triag_p1)
{
  point.resize(2);
  const Real xi0_test[3] = {-0.7, 0.7, -1.0};
  const Real xi1_test[3] = {-0.4, -0.4, 1.0};

  for (Uint i = 0; i < 3; ++i)
  {
    point[XI0] = xi0_test[i];
    point[XI1] = xi1_test[i];
    check_modal_basis<mesh::sf::CarnevaliExpansionTriag, mesh::utest_fixture::p1_triag_Lagrange_sf>(
        "Triag-P1-Equidist", point, 1e-14, true);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Carnevali_modal_basis_triag_p2)
{
  point.resize(2);
  const Real xi0_test[3] = {-0.7, 0.7, -1.0};
  const Real xi1_test[3] = {-0.4, -0.4, 1.0};

  for (Uint i = 0; i < 3; ++i)
  {
    point[XI0] = xi0_test[i];
    point[XI1] = xi1_test[i];
    check_modal_basis<mesh::sf::CarnevaliExpansionTriag, mesh::utest_fixture::p2_triag_Lagrange_sf>(
        "Triag-P2-Equidist", point, 1e-14, true);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Carnevali_modal_basis_triag_p3)
{
  point.resize(2);
  const Real xi0_test[3] = {-0.7, 0.7, -1.0};
  const Real xi1_test[3] = {-0.4, -0.4, 1.0};

  for (Uint i = 0; i < 3; ++i)
  {
    point[XI0] = xi0_test[i];
    point[XI1] = xi1_test[i];
    check_modal_basis<mesh::sf::CarnevaliExpansionTriag, mesh::utest_fixture::p3_triag_Lagrange_sf>(
        "Triag-P3-Equidist", point, 1e-14, true);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Carnevali_prime_basis_derivatives_triag_p2)
{
  math::DenseDMat<Real> coordinates(3, 2);
  coordinates(0, 0) = -1.;
  coordinates(0, 1) = -1.;
  coordinates(1, 0) = 0.0;
  coordinates(1, 1) = -1.0;
  coordinates(2, 0) = 0.0;
  coordinates(2, 1) = 0.0;

  mesh::sf::CarnevaliExpansionTriag det;
  det.set_polynomial_order(P2);

  std::vector<math::DenseDMat<Real>> dV;

  det.Vandermonde_matrix_derivatives(coordinates, dV);

  std::cout << "Derivatives with respect to xi0 = " << std::endl << dV[0] << std::endl << std::endl;
  std::cout << "Derivatives with respect to xi1 = " << std::endl << dV[1] << std::endl << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
