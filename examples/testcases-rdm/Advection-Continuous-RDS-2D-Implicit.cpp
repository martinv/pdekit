#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "common/MPI/MPIEnv.hpp"
#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/ErrorEstimator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/scalar/RotationAdvection2D.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/PGRDMethodImplicit.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"

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

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  const std::string infilename = "advection-tri.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_continuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                               PointSetID::Warpblend);

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P4;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  vector_function::ptr solution = std::make_shared<vector_function>("", "solution");
  vector_function::ptr residual = std::make_shared<vector_function>("", "residual");

  const Uint nb_nodes = (*sol_dofs).nb_nodes();

  solution->resize(phys_model::NEQ, nb_nodes);
  residual->resize(phys_model::NEQ, nb_nodes);

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;

  scheme_type scheme;

  scheme.configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme.initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme.set_vec_function(SolverVecFn::solution, solution);
  scheme.set_vec_function(SolverVecFn::residuals, residual);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);
  solver::InitialCondition<MeshConfig> ic("hat_initial_condition", mesh2D);

  ic.set_domain(_2D, "fluid");
  ic.set_expression(&CosineHat2D::value);
  ic.apply_function("sol_dofs", *solution);
  // ic.set_domain(_1D,"inlet");
  // ic.apply_function(hat,*solution);

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

  clock_t start, end;
  Real elapsed;
  std::cout.precision(15);

  math::DenseDVec<Real> res_L1_norm(1);

  start = clock();

  for (Uint iter = 0; iter < 25; ++iter)
  {
    // std::cout << "Iter = " << std::setw(5) << iter << ", ";

    const Real CFL = iter < 10 ? 3.0 : std::min(15.0, 4.0 + std::pow(1.30, iter - 10));
    // std::cout << "dt = " << std::setw(20) << CFL << " ";

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

    /*
    std::cout << ", res = " << std::setw(20) << res_L1_norm[0] << " , log(res) = " << std::setw(20)
              << std::log(res_L1_norm[0]) << std::endl;
    */

    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L1_norm);
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Compute the error in the solution
  Rotation2DSolution exact_solution;
  interpolation::ErrorEstimator<MeshConfig> error_estim(mesh2D);
  error_estim.compute_pointwise_L2_error(exact_solution, *sol_dofs, *solution);
  error_estim.compute_L1_error(exact_solution, *mesh2D, *sol_dofs, *solution);
  error_estim.compute_L2_error(exact_solution, *mesh2D, *sol_dofs, *solution);
  error_estim.compute_infty_error(exact_solution, *mesh2D, *sol_dofs, *solution);

  /// Write the output to file
  const std::string outfilename = "output_advection_p1_2D_implicit.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "solution");

  /*
  vtk::VtkWriter vtkwriter;
  vtkwriter.append_nodal_function_to_file(*sol_mesh, "output.vtu", *solution,
  "solution");
  */

  mpi_env.finalize();

  return 0;
}
