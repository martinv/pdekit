#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>

// #include <sys/resource.h>

#include "common/PDEKit.hpp"
#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/ErrorEstimator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/scalar/RotationAdvection2D.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodImplicit.hpp"
#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/bc/WeakDirichletBC.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::RotationAdvection2D<physics::AroundOrigin2D> phys_model;

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // ------------------------------------------------------
  // INITIALIZE ENVIRONMENT
  // ------------------------------------------------------

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  // struct rusage r_usage;

  std::shared_ptr<common::Component> mesh_component = common::make_component<MeshType>("mesh2D");
  std::shared_ptr<MeshType> mesh2D                  = mesh_component->handle<MeshType>();

  const std::string infilename = "advection-tri.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                                  PointSetID::Warpblend);

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P4;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  vector_function::ptr solution = std::make_shared<vector_function>("", "solution");
  vector_function::ptr residual = std::make_shared<vector_function>("", "residual");

  // Choose scheme, configure the function spaces (solution, residual,
  // update coefficients

  const Uint nb_nodes = (*sol_dofs).nb_nodes();

  // const Uint nb_facets = sol_facets.nb_cells();

  solution->resize(phys_model::NEQ, nb_nodes);
  residual->resize(phys_model::NEQ, nb_nodes);

  //  interpolation::FunctionSpace::Function::Ptr f_ptr =
  // fs_solution->function("residual");
  //  std::cout << "Got function " << f_ptr->name() << std::endl;

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  typedef rdm::DRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA, rdm::FacetDG> scheme_type;
  scheme_type scheme;

  /*
  scheme.configure_spaces(fs_geometry_cells, fs_solution_cells,
  fs_geometry_cell_contours, fs_solution_cell_contours, fs_geometry_facets,
  fs_solution_facets);
  */

  scheme.configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme.initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);

  scheme.set_vec_function(SolverVecFn::solution, solution);
  scheme.set_vec_function(SolverVecFn::residuals, residual);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  /*
  // Initial condition
  solver::InitialCondition<MeshConfig> ic("hat_initial_condition", sol_mesh);

  ic.set_domain(_2D, "fluid");
  ic.set_expression(&CosineHat2D::value);
  ic.apply_function(*solution);
  */

  /*
  result_of::geometry<MeshConfig>::type const &sol_geometry =
  sol_mesh->geometry();

  for (Uint c = 0; c < sol_cells.nb_cells(); ++c)
  {
    const mesh::MeshEntity cell = sol_cells.cell(c);

    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      const math::ConstVectorBlock<Real> node_coords =
  sol_geometry.node(cell.vertex(v));
      // const Real init_value = std::sin(node_coords[X0] * node_coords[X1]);
      const Real init_value = node_coords[X0] * node_coords[X0] +
  node_coords[X1] * node_coords[X1];
      (*solution)(0, cell.vertex(v)) = init_value;
    }
  }
  */

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  std::shared_ptr<scheme_type::bc_base_type> inlet_bc =
      scheme.add_boundary_condition("StrongDirichlet", "inlet_bc", "inlet");
  inlet_bc->set_expression(&CosineHat2D::value);

  std::shared_ptr<scheme_type::bc_base_type> farfield_bc =
      scheme.add_boundary_condition("StrongDirichlet", "strong_farfield", "farfield");
  farfield_bc->set_expression(&Zero::value);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::cout.precision(15);

  math::DenseDVec<Real> res_L1_norm(1);

  start = std::chrono::system_clock::now();

  for (Uint iter = 0; iter < 50; ++iter)
  {
    /*
    nodal_residual->fill(0.0);
    time_update->reset_wave_speeds();
    */

    const Real CFL = iter < 5 ? 3.0 : std::min(25.0, 4.0 + std::pow(1.30, iter - 5));

    scheme.set_cfl(CFL);

    if (iter % 5 == 0)
    {
      scheme.assemble_lhs_and_rhs(CFL);
      scheme.solve({SolverOption::RecomputePreconditioner});
    }
    else
    {
      scheme.assemble_rhs();
      scheme.solve({});
    }

    scheme.compute_residual_norm(res_L1_norm);

    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L1_norm);
  }

  end                                 = std::chrono::system_clock::now();
  std::chrono::duration<Real> elapsed = end - start;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed.count() << " s" << std::endl;
  // getrusage(RUSAGE_SELF, &r_usage);
  // std::cout << "Memory usage: " << r_usage.ru_maxrss << " Kb" << std::endl;

  /// Compute the error in the solution
  Rotation2DSolution exact_solution;
  interpolation::ErrorEstimator<MeshConfig> error_estim(mesh2D);
  error_estim.compute_pointwise_L2_error(exact_solution, *sol_dofs, *solution);
  error_estim.compute_L1_error(exact_solution, *mesh2D, *sol_dofs, *solution);
  error_estim.compute_L2_error(exact_solution, *mesh2D, *sol_dofs, *solution);
  error_estim.compute_infty_error(exact_solution, *mesh2D, *sol_dofs, *solution);

  /// Write the output to file
  const std::string outfilename = "output_advection_p1_discontinuous_2D.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "solution");
  // meshwriter.append_nodal_function_to_file(*sol_mesh, outfilename,
  // *nodal_residual, "residual");
  // meshwriter.append_nodal_function_to_file(*sol_mesh, outfilename,
  // *update_coeff, "update_coeff");

  /*
  vtk::VtkWriter vtkwriter;
  vtkwriter.append_nodal_function_to_file(*mesh2D,
  "output_advection_p1_discontinuous_2D.vtu", *solution, "solution");
  */

  mpi_env.finalize();

  return 0;
}
