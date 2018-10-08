/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RDSolver_test
#include <boost/test/unit_test.hpp>

#include <ctime>
#include <iostream>

#include "common/MPI/MPIEnv.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "math/MathConstants.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/MeshManipulator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/rdm/RDSolver.hpp"
#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::Euler2DCons phys_model;

// ----------------------------------------------------------------------------

struct RDSolverUtestGlobalFixture
{
  /// common setup for each test case
  RDSolverUtestGlobalFixture()
  {
    /// arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~RDSolverUtestGlobalFixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(RDSolverTestSuite, RDSolverUtestGlobalFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(RDSolverTestSuite_ParenvInit)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(m_argc, m_argv);
  BOOST_CHECK_EQUAL(mpi_env.is_initialized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(RD_solver_utest)
{
  using scheme_type = rdm::RDSolver<MeshConfig, phys_model>;
  scheme_type scheme;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(RDSolverTestSuite_ParenvFinalize)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance();
  mpi_env.finalize();
  BOOST_CHECK_EQUAL(mpi_env.is_finalized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
