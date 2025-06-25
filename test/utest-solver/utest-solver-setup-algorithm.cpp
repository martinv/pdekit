/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE solver_setup_algorithm_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iomanip>
#include <iostream>

#include "mesh/io/MeshCreator.hpp"
#include "solver/SolverSetupAlgorithm.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

struct SolverSetupAlgorithm_Fixture
{
  /// common setup for each test case
  SolverSetupAlgorithm_Fixture()
  {
    /// arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~SolverSetupAlgorithm_Fixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;

  struct PolyRule
  {
    // Return: a tuple consisting of
    // 1) Shape function tag of source shape function (which we're
    // projecting FROM) 2) Shape function tag of target polynomial (i.e
    // space we're projecting TO)
    mesh::sf::SFTag sf_tag(const mesh::PointSetTag &geo_tag, const mesh::PointSetTag &sol_tag) const
    {
      const mesh::sf::SFTag poly_space_tag(sol_tag.elem_shape(), SFunc::Lagrange,
                                           sol_tag.poly_order(), ModalBasis::Modal);

      return poly_space_tag;
    }

    // Return: a tuple consisting of
    // 1) Quadrature order to be used in order to construct mass matrix
    // operator
    Uint quadrature_order(const mesh::PointSetTag &geo_tag, const mesh::PointSetTag &sol_tag) const
    {
      const Uint quad_order = std::max(2 * sol_tag.poly_order(), geo_tag.poly_order());
      return quad_order;
    }
  };
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(SolverSetupAlgorithm_TestSuite, SolverSetupAlgorithm_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(solve_setup_in_mesh_interior)
{
  using MeshType = Tria<Cart2D>;
  // using dof_range_t =
  // common::IteratorRange<MeshType::dof_storage_type::const_dof_iterator>;

  std::shared_ptr<MeshType> mesh2D(new MeshType("mesh2D"));

  // Create a quad domain [0,1] x [0,1] meshed with triangles
  MeshCreator::make_unit_quad(*mesh2D, "geo_dofs", 5, true);

  // Increase the polynomial order
  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");
  MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P2,
                                                  PointSetID::Warpblend);

  solver::SolverSetupAlgorithm setup_algorithm;

  std::vector<mesh::DiscreteElemKey> geo_fe_values;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> sol_fe_values;

  using geo_dof_iter_t = decltype((*geo_dofs).cbegin());
  using sol_dof_iter_t = decltype((*sol_dofs).cbegin());

  setup_algorithm.generate_fe_pairs<geo_dof_iter_t, sol_dof_iter_t>(
      (*geo_dofs).cbegin(), {common::make_iter_range((*sol_dofs).cbegin(), (*sol_dofs).cend())},
      PolyRule{}, geo_fe_values, sol_fe_values);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(solver_setup_on_mesh_boundary)
{
  using MeshType = Tria<Cart2D>;
  // using dof_range_t =
  // common::IteratorRange<MeshType::dof_storage_type::const_dof_iterator>;

  std::shared_ptr<MeshType> mesh2D(new MeshType("mesh2D"));

  // Create a quad domain [0,1] x [0,1] meshed with triangles
  MeshCreator::make_unit_quad(*mesh2D, "geo_dofs", 5, true);

  // Increase the polynomial order
  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");
  MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P2,
                                                  PointSetID::Warpblend);

  solver::SolverSetupAlgorithm setup_algorithm;

  std::vector<mesh::DiscreteElemKey> geo_fe_values;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> sol_fe_values;

  using boundary_iterator = mesh::BoundaryFacets<Cart2D, Cart2D::FACET_DIM>::const_dof_iterator;
  const boundary_iterator geo_bdry_begin =
      mesh2D->all_boundaries().domain("Edge0")->cbegin(*geo_dofs);
  // const boundary_iterator geo_bdry_end =
  // mesh2D->all_boundaries().domain("Edge0")->cend(*geo_dofs);
  const boundary_iterator sol_bdry_begin =
      mesh2D->all_boundaries().domain("Edge0")->cbegin(*sol_dofs);
  const boundary_iterator sol_bdry_end = mesh2D->all_boundaries().domain("Edge0")->cend(*sol_dofs);

  setup_algorithm.generate_fe_pairs<boundary_iterator, boundary_iterator>(
      geo_bdry_begin, {common::make_iter_range(sol_bdry_begin, sol_bdry_end)}, PolyRule{},
      geo_fe_values, sol_fe_values);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
