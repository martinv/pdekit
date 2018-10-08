/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Poisson_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>
#include <map>
#include <memory>

/// PDEKIT headers
#include "common/MPI/MPIEnv.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "linear_system/LSTpetra.hpp"
#include "math/MathConstants.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

#include "solver/fe/HelmholtzSolverCG.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

#if PDEKIT_HAVE_TRILINOS
using namespace pdekit::ls;
#endif

// ----------------------------------------------------------------------------

struct MPIFixture
{
  MPIFixture()
  {
    int argc    = boost::unit_test::framework::master_test_suite().argc;
    char **argv = boost::unit_test::framework::master_test_suite().argv;

    common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);
    is_initialized                              = mpi_env.is_initialized();
    is_finalized                                = mpi_env.is_finalized();
  }

  ~MPIFixture()
  {
    common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance();
    mpi_env.finalize();

    std::cout << "MPI environment is finalized (should be true): " << mpi_env.is_finalized()
              << std::endl;
  }

  static bool is_initialized;
  static bool is_finalized;
};

bool MPIFixture::is_initialized = false;
bool MPIFixture::is_finalized   = false;

// ----------------------------------------------------------------------------

struct WeakDirBC_Fixture
{
  /// common setup for each test case
  WeakDirBC_Fixture()
  {
  }

  /// common tear-down for each test case
  ~WeakDirBC_Fixture()
  {
  }

  template <typename MeshConfig>
  void solve_ls(const mesh::Tria<MeshConfig> &in_mesh,
                std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
                std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs,
                interpolation::VectorMeshFunction<Real> &solution)
  {
    std::shared_ptr<ls::TpetraMultiVector<Real>> x(
        new ls::TpetraMultiVector<Real>(global_stiffness_matrix->map()));

    ls::LSTpetra<Real> lin_system;
    // lin_system.configure(global_stiffness_matrix, global_rhs, x, true,
    // false);
    lin_system.initialize_solver(global_stiffness_matrix, global_rhs, x, false);

    std::shared_ptr<ls::TrilinosPC<Real>> preconditioner = std::make_shared<ls::IfpackPC<Real>>();
    preconditioner->create("ILUT", global_stiffness_matrix);

    lin_system.connect_preconditioner(preconditioner);
    lin_system.update_after_mat_values_change(global_stiffness_matrix, global_rhs, x, false);
    lin_system.solve(550, 1.e-12);

    // Create a mesh function to copy the result

    typedef typename pdekit::result_of::dof_map_t<MeshConfig> cell_dofs_type;
    const cell_dofs_type &cell_dofs = *(in_mesh.dof_storage("geo_dofs"));
    const Uint nb_nodes_in_mesh     = cell_dofs.nb_nodes();

    solution.resize(1, nb_nodes_in_mesh);
    solution.fill(0.0);

    for (Uint n = 0; n < nb_nodes_in_mesh; ++n)
    {
      typename interpolation::VectorMeshFunction<Real>::entry_type nodal_value = solution.value(n);
      nodal_value[0]                                                           = (*x).value(n, 0);
    }
  }

  // ----------------------------------------------------------------------------

  // Dirichlet boundary condition for 2D Poisson problem
  struct DirichletBC2D
  {
    inline static Real eval(const math::DenseConstVecView<Real> &coord)
    {
      return 1.0 + std::sin(math::pi * coord[X0]) + std::sin(math::pi * coord[X1]);
    }
  };

  // ----------------------------------------------------------------------------

  // Dirichlet boundary condition for 3D Poisson problem
  struct DirichletBC3D
  {
    inline static Real eval(const math::DenseConstVecView<Real> &coord)
    {
      return 1.0 + std::sin(math::pi * coord[X0]) + std::sin(math::pi * coord[X1]) +
             +std::sin(math::pi * coord[X2]);
    }
  };

  // ----------------------------------------------------------------------------

  // RHS for 2D Poisson problem
  struct RHS2D
  {
    /*
    inline static const Real eval(const math::DenseConstVecView<Real>
    &coord)
    {
      return 2.0 * math::pi * math::pi * std::sin(math::pi * coord[X0]) *
             std::sin(math::pi * coord[X1]);
    }
    */
    inline static Real eval(const math::DenseConstVecView<Real> &coord)
    {
      return 0.0;
    }
  };

  // ----------------------------------------------------------------------------

  // RHS for 3D Poisson problem
  struct RHS3D
  {
    /*
    inline static const Real eval(const math::DenseConstVecView<Real>
    &coord)
    {
      return 3.0 * math::pi * math::pi * std::sin(math::pi * coord[X0]) *
             std::sin(math::pi * coord[X1]) * std::sin(math::pi *
    coord[X2]);
    }
    */
    inline static Real eval(const math::DenseConstVecView<Real> &coord)
    {
      return 0.0;
    }
  };
};

// ----------------------------------------------------------------------------

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_FIXTURE_TEST_SUITE(WeakDirBC_TestSuite, WeakDirBC_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Poisson_utest_2D)
{
  BOOST_CHECK_EQUAL(MPIFixture::is_initialized, true);
  BOOST_CHECK_EQUAL(MPIFixture::is_finalized, false);

  const std::string infilename  = "L_shape_2D_tri_p3.msh";
  const std::string outfilename = "Poisson_solution_tri.msh";

  const std::string assembly_method = "Nitsche"; // standard, Nitsche

  mesh::gmsh::GmshReader gmsh_reader;
  mesh::Tria<Cart2D>::shared_ptr input_mesh = std::make_shared<mesh::Tria<Cart2D>>("mesh");
  gmsh_reader.read_mesh_from_file(infilename, *input_mesh, "geo_dofs");

  std::shared_ptr<ls::TpetraCrsMatrix<Real>> global_stiffness_matrix;
  std::shared_ptr<ls::TpetraMultiVector<Real>> global_rhs;

  RHS2D rhs;
  DirichletBC2D bc;

  solver::fe::HelmholtzSolverCG<Cart2D> Helm_solver;
  const std::vector<std::string> boundary_names = {"Wall1", "Wall2", "Wall3",
                                                   "Wall4", "Wall5", "Wall6"};
  Helm_solver.assemble(assembly_method, *input_mesh, rhs, boundary_names, bc,
                       global_stiffness_matrix, global_rhs);

  interpolation::VectorMeshFunction<Real> solution("", "solution");
  solve_ls(*input_mesh, global_stiffness_matrix, global_rhs, solution);

  mesh::gmsh::GmshWriter gmsh_writer;
  gmsh_writer.write_mesh_to_file(*input_mesh, "geo_dofs", outfilename);
  gmsh_writer.append_nodal_function_to_file(*input_mesh, outfilename, solution, "solution");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Poisson_utest_3D)
{
  BOOST_CHECK_EQUAL(MPIFixture::is_initialized, true);
  BOOST_CHECK_EQUAL(MPIFixture::is_finalized, false);

  const std::string infilename  = "L_shape_3D_tet_p3.msh";
  const std::string outfilename = "Poisson_solution_tet.msh";

  const std::string assembly_method = "Nitsche"; // standard, Nitsche

  mesh::gmsh::GmshReader gmsh_reader;
  mesh::Tria<Cart3D>::shared_ptr input_mesh = std::make_shared<mesh::Tria<Cart3D>>("mesh");
  gmsh_reader.read_mesh_from_file(infilename, *input_mesh, "geo_dofs");

  std::shared_ptr<ls::TpetraCrsMatrix<Real>> global_stiffness_matrix;
  std::shared_ptr<ls::TpetraMultiVector<Real>> global_rhs;

  RHS3D rhs;
  DirichletBC3D bc;

  solver::fe::HelmholtzSolverCG<Cart3D> Helm_solver;
  const std::vector<std::string> boundary_names = {"Bottom", "Top", "Wall"};
  Helm_solver.assemble(assembly_method, *input_mesh, rhs, boundary_names, bc,
                       global_stiffness_matrix, global_rhs);

  interpolation::VectorMeshFunction<Real> solution("", "solution");
  solve_ls(*input_mesh, global_stiffness_matrix, global_rhs, solution);

  mesh::gmsh::GmshWriter gmsh_writer;
  gmsh_writer.write_mesh_to_file(*input_mesh, "geo_dofs", outfilename);
  gmsh_writer.append_nodal_function_to_file(*input_mesh, outfilename, solution, "solution");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
