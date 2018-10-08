/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HDG_elliptic_solver_utest
#include <boost/test/unit_test.hpp>

#include <ctime>
#include <iostream>
#include <memory>

#include "math/MathConstants.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "solver/fe/HelmholtzSolverHDG.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig2D;
typedef Tria<MeshConfig2D> MeshType2D;
// typedef interpolation::FunctionSpace<MeshConfig>::scal_f scalar_function;
// typedef interpolation::FunctionSpace<MeshConfig>::vect_f vector_function;

// ----------------------------------------------------------------------------

struct HDGEllipticSolverUtestFixture
{
  /// common setup for each test case
  HDGEllipticSolverUtestFixture()
  {
    // uncomment if you want to use arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~HDGEllipticSolverUtestFixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(HDGEllipticSolverTestSuite, HDGEllipticSolverUtestFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(hdg_elliptic_solver_parenv_init)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(m_argc, m_argv);
  BOOST_CHECK_EQUAL(mpi_env.is_initialized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(HDG_elliptic_solver_2D_utest)
{
  MeshType2D::shared_ptr mesh = std::make_shared<MeshType2D>("mesh");

  MeshCreator::make_unit_quad(*mesh, "geo_dofs", 11, true);
  common::PtrHandle<MeshType2D::dof_storage_type> geo_dofs = mesh->dof_storage("geo_dofs");
  common::PtrHandle<MeshType2D::dof_storage_type> sol_dofs = mesh->create_dof_storage("sol_dofs");

  MeshType2D::dof_storage_type::clone_discontinuous(*mesh, *geo_dofs, *sol_dofs, P3,
                                                    PointSetID::Warpblend);

  fe::HelmholtzSolverHDG<MeshConfig2D> solver;

  interpolation::VectorMeshFunction<Real>::ptr solution(
      new interpolation::VectorMeshFunction<Real>("", "solution"));

  solver.setup(*mesh, *sol_dofs);
  solver.set_solution(solution);

  solver.solve(*mesh, *sol_dofs);

  gmsh::GmshWriter meshwriter;
  const std::string outfilename = "solution_HDG_elliptic_2D.msh";
  meshwriter.write_mesh_to_file(*mesh, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*mesh, "solution_HDG_elliptic_2D.msh", *solution,
                                           "solution");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(hdg_elliptic_solver_parenv_finalize)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance();
  mpi_env.finalize();
  BOOST_CHECK_EQUAL(mpi_env.is_finalized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
