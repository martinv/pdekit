#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/ErrorEstimator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/scalar/RotationAdvection2D.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodImplicit.hpp"
#include "solver/rdm/PGRDMethodImplicit.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::RotationAdvection2D<physics::ConstAdvection2D> phys_model;

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // ------------------------------------------------------
  // INITIALIZE ENVIRONMENT
  // ------------------------------------------------------

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // SELECT SCHEME
  // ------------------------------------------------------

  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;
  // typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::LF>
  // scheme_type; typedef rdm::DRDMethodImplicit<MeshConfig, phys_model,
  // rdm::PGLDA, rdm::FacetDG> scheme_type;

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme = std::make_shared<scheme_type>();

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  const std::string infilename = "unit_square_tri_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  if (scheme->is_continuous())
  {
    MeshType::dof_storage_type::clone_continuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                                 PointSetID::Warpblend);
  }
  else
  {
    MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                                    PointSetID::Warpblend);
  }

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
  // BLENDING COEFFICIENT FOR NONLINEAR SCHEMES
  // ------------------------------------------------------

  /*
  rdm::RDBlendingCoeff<MeshConfig> rd_blending_coeff_operator;

  // Types of viscosity coefficient: "LF_Blend", "ViscosityIndicator"
  const std::string blend_coeff_type = "LF_Blend";
  rd_blending_coeff_operator.setup(geo_mesh->topology().dof_storage(),
                                   sol_mesh->topology().dof_storage(),
  blend_coeff_type);

  scalar_function::ptr blending_coeff =
      std::make_shared<scalar_function>("", "blending_coefficient");
  blending_coeff->resize(geo_mesh->topology().dof_storage().nb_active_cells());
  */

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);

  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residual);
  // scheme->set_blending_coeff(blending_coeff);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);
  solver::InitialCondition<MeshConfig> ic("zero_ic", mesh2D);

  ic.set_domain(_2D, "interior");
  ic.set_expression(&Zero::value);
  ic.apply_function("sol_dofs", *solution);
  // ic.set_domain(_1D,"inlet");
  // ic.apply_function(hat,*solution);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  std::shared_ptr<scheme_type::bc_base_type> left_bc =
      scheme->add_boundary_condition("StrongDirichlet", "left_bc", "left");
  left_bc->set_expression(&Zero::value);

  std::shared_ptr<scheme_type::bc_base_type> right_bc =
      scheme->add_boundary_condition("StrongDirichlet", "right_bc", "right");
  right_bc->set_expression(&Zero::value);

  std::shared_ptr<scheme_type::bc_base_type> inlet_bc =
      scheme->add_boundary_condition("StrongDirichlet", "inlet_bc", "bottom");
  inlet_bc->set_expression(&SineWave2D::value);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;
  std::cout.precision(15);

  math::DenseDVec<Real> res_L2_norm(1);
  math::DenseDVec<Real> res_L2_norm0(1);

  const Real CFL0   = 1.1;
  const Real CFLmax = 100.0;
  const Real alpha  = 1.05;

  start = clock();

  for (Uint iter = 0; iter < 70; ++iter)
  {
    // Blending coefficient for nonlinear schemes
    /*
    rd_blending_coeff_operator.calculate(geo_mesh->topology().dof_storage(),
                                         sol_mesh->topology().dof_storage(),
    *solution, *blending_coeff);
    */

    // const Real CFL = iter < 30 ? 3.0 : std::min(50.0, 4.0 +
    // std::pow(1.30, iter - 30));
    const Real CFL =
        iter < 1 ? CFL0
                 : std::max(10.0, std::min(CFL0 * std::pow(res_L2_norm0[0] / res_L2_norm[0], alpha),
                                           CFLmax));

    if (iter % 5 == 0)
    {
      scheme->set_cfl(CFL);
      scheme->assemble_lhs_and_rhs(CFL);
      scheme->solve({SolverOption::RecomputePreconditioner});
    }
    else
    {
      scheme->assemble_rhs();
      scheme->solve({});
    }

    scheme->compute_residual_norm(res_L2_norm);

    // Store the initial norm
    if (iter == 0)
    {
      res_L2_norm0 = res_L2_norm;
    }

    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Compute the error in the solution
  SineWave2DSolution exact_solution;
  interpolation::ErrorEstimator<MeshConfig> error_estim(mesh2D);
  error_estim.compute_pointwise_L2_error(exact_solution, *sol_dofs, *solution);
  error_estim.compute_L1_error(exact_solution, *mesh2D, *sol_dofs, *solution);
  error_estim.compute_L2_error(exact_solution, *mesh2D, *sol_dofs, *solution);
  error_estim.compute_infty_error(exact_solution, *mesh2D, *sol_dofs, *solution);

  /// Write the output to file
  const std::string outfilename = "output_lin_advection_2D_implicit.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "solution");

  /*
  meshwriter.append_cell_function_to_file(*sol_mesh, outfilename,
  *blending_coeff, "blending_coefficient");
  */

  /*
  vtk::VtkWriter vtkwriter;
  vtkwriter.append_nodal_function_to_file(*sol_mesh,
  "output_lin_advection_2D_implicit.vtu", *solution, "solution");

  */

  mpi_env.finalize();

  return 0;
}
