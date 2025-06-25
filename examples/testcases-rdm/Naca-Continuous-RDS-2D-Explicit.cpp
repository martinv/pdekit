#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/time/ExplicitTimeStepper.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::interpolation;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::Euler2DCons phys_model;

// ----------------------------------------------------------------------------

int main()
{
  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  const std::string infilename = "naca0012_p1_tri.msh";
  // const std::string infilename = "naca0012_p2_mixed_elem.msh";
  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_continuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                               PointSetID::Equidist);

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P2;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  vector_function::ptr solution = std::make_shared<vector_function>("", "solution");
  vector_function::ptr residual = std::make_shared<vector_function>("", "residual");

  const Uint nb_nodes = (*sol_dofs).nb_nodes();

  solution->resize(phys_model::NEQ, nb_nodes);
  residual->resize(phys_model::NEQ, nb_nodes);

  //  interpolation::FunctionSpace::Function::Ptr f_ptr =
  // fs_solution->function("residual");
  //  std::cout << "Got function " << f_ptr->name() << std::endl;

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  typedef rdm::PGRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA> lda_scheme_type;
  lda_scheme_type lda_scheme;

  typedef rdm::PGRDMethodExplicit<MeshConfig, phys_model, rdm::PGB> b_scheme_type;
  b_scheme_type b_scheme;

  lda_scheme.configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  lda_scheme.initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);

  lda_scheme.set_vec_function(SolverVecFn::solution, solution);
  lda_scheme.set_vec_function(SolverVecFn::residuals, residual);

  b_scheme.configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  b_scheme.initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);

  b_scheme.set_vec_function(SolverVecFn::solution, solution);
  b_scheme.set_vec_function(SolverVecFn::residuals, residual);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh2D);
  ic.set_domain(_2D, "InnerCells");
  math::DenseDVec<Real> ic_state(phys_model::NEQ);

  ic_state[0] = 1.22503;
  ic_state[1] = 208.30559;
  ic_state[2] = 7.27419;
  ic_state[3] = 271044.38;

  ic.apply_values("sol_dofs", ic_state,
                  *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  lda_scheme.add_boundary_condition("WeakWall", "wall_bc", "Wall");

  // std::shared_ptr<bc_base_type> farfield_bc =
  auto farfield_bc = lda_scheme.add_boundary_condition("WeakFarfield", "farfield_bc", "Farfield");
  farfield_bc->set_reference_state(ic_state);
  farfield_bc->print_bc_parameters();

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  math::DenseDVec<Real> residual_norm(phys_model::NEQ);

  clock_t start, end;
  Real elapsed;

  start = clock();

  time::RK3TVD time_stepper;
  time_stepper.setup(phys_model::NEQ, nb_nodes);

  for (Uint iter = 0; iter <= 300; ++iter)
  {
    const Real CFL = 0.8;

    /*
    lda_scheme.assemble_lhs_and_rhs(0.5);
    lda_scheme.solve({});
    */

    time_stepper.advance_in_time(lda_scheme, CFL);

    if (iter % 10 == 0)
    {
      interpolation::norm_L2(*residual, residual_norm);
      solver::SolverIO::print_iter_and_res_norm(iter, CFL, residual_norm);
    }
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  // ------------------------------------------------------
  // WRITE THE OUTPUT TO FILE
  // ------------------------------------------------------

  const std::string outfilename = "output_naca_continuous.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "u");

  // Write Mach number values to file
  scalar_function::ptr Mach_number = std::make_shared<scalar_function>("", "Ma");
  Mach_number->resize(nb_nodes);

  physics::Euler2DCons physical_model;
  physics::Euler2DProperties properties;
  physics::Euler2DProperties::SolGradM gradient_matrix;

  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution, [](const phys_model::Properties &props) { return props.Ma; },
      *Mach_number);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *Mach_number, "Ma");
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *residual, "Res");

  return 0;
}
