#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/ErrorEstimator.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/scalar/AdvDiffLinear2DFOS.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/time/ExplicitTimeStepper.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::AdvDiffLinear2DFOS phys_model;

// ----------------------------------------------------------------------------

int main()
{
  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  const std::string infilename = "unit_square_tri_p1.msh";

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

  typedef rdm::PGRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme = std::make_shared<scheme_type>();
  // scheme_type scheme;

  scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residual);

  time::ExplicitEuler time_stepper;
  time_stepper.setup(phys_model::NEQ, nb_nodes);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);
  solver::InitialCondition<MeshConfig> ic("diffusion_initial_condition", mesh2D);

  ic.set_domain(_2D, "interior");
  ic.set_expression(&Zero::value);
  ic.apply_function("sol_dofs", *solution);
  // ic.set_domain(_1D,"inlet");
  // ic.apply_function(hat,*solution);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  std::shared_ptr<scheme_type::bc_base_type> bottom_bc =
      scheme->add_boundary_condition("WeakDirichlet", "bottom_bc", "bottom");
  bottom_bc->set_expression(&DiffusionBC2D::value);

  std::shared_ptr<scheme_type::bc_base_type> top_bc =
      scheme->add_boundary_condition("WeakDirichlet", "top_bc", "top");
  top_bc->set_expression(&DiffusionBC2D::value);

  std::shared_ptr<scheme_type::bc_base_type> left_bc =
      scheme->add_boundary_condition("WeakDirichlet", "left_bc", "left");
  left_bc->set_expression(&DiffusionBC2D::value);

  std::shared_ptr<scheme_type::bc_base_type> right_bc =
      scheme->add_boundary_condition("WeakDirichlet", "right_bc", "right");
  right_bc->set_expression(&DiffusionBC2D::value);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;
  std::cout.precision(15);

  math::DenseDVec<Real> res_L1_norm(1);

  start = clock();

  for (Uint iter = 0; iter < 150; ++iter)
  {
    const Real CFL = 0.5;

    /*
    residual->fill(0.0);
    time_update->reset_wave_speeds();

    // weak_inlet_bc.apply();
    // strong_inlet_bc.apply();
    // farfield.apply();

    scheme->iterate_by_std_region_type();
    scheme->apply_boundary_conditions();

    time_update->compute_local_time_step(0.5);
    const scalar_function &dt = time_update->time_update_function();
    (*solution) -= dt * (*residual);
    */

    // Blending coefficient for nonlinear schemes
    /*
    blending_coeff->calculate_blending_coeff(topology.dof_storage(),
    sol_topology.dof_storage(), *solution);
    */

    // time_stepper.advance_in_time(scheme, 0.5);

    scheme->assemble_lhs_and_rhs(CFL);
    scheme->solve({});

    norm_L1(*residual, res_L1_norm);
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
  const std::string outfilename = "output_diffusion_p1_2D.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "solution");
  meshwriter.append_nodal_function_to_file(
      *mesh2D, outfilename, *scheme->scal_function(SolverScalFn::time_step), "time_step");

  vtk::VtkWriter vtkwriter;
  vtkwriter.append_nodal_function_to_file(*mesh2D, *sol_dofs, "output_diffusion_p1_2D.vtu",
                                          *solution, "solution");

  return 0;
}
