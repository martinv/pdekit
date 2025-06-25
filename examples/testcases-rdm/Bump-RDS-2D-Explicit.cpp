#include <ctime>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodExplicit.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/bc/WeakSubInlet.hpp"
#include "solver/rdm/bc/WeakSubOutlet.hpp"
#include "solver/rdm/bc/WeakWall.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"

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

  typedef rdm::PGRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;
  // typedef rdm::DRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA,
  // rdm::FacetDG> scheme_type;
  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme = std::make_shared<scheme_type>();

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  const std::string infilename = "bump_p1_tri.msh";
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

  interpolation::FunctionSpace<MeshConfig>::ptr fs_solution =
      std::make_shared<interpolation::FunctionSpace<MeshConfig>>();

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

  scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residual);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  // Material values
  const Real gamma = 1.4;
  const Real R     = 287.046;

  // Inlet values
  const Real Ttot = 307.6488978;
  const Real Ptot = 120195.4453;
  const Real M_in = 0.5;

  // Outlet values
  const Real Pout = 101325.0;

  const Real T_in = Ttot / (1. + 0.5 * (gamma - 1) * M_in * M_in);
  const Real p_in = Ptot / std::pow(1. + 0.5 * (gamma - 1) * M_in * M_in, gamma / (gamma - 1));

  // SinusBumpInit  init_func;

  /// Initial condition

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh2D);
  ic.set_domain(_2D, "domain");
  math::DenseDVec<Real> ic_state(physics::Euler2DCons::NEQ);

  const Real rho_in = p_in / (R * T_in);
  const Real u_in   = rho_in * M_in * sqrt(gamma * R * T_in);
  const Real v_in   = 0.0;

  ic_state[0] = rho_in;
  ic_state[1] = rho_in * u_in;
  ic_state[2] = rho_in * v_in;
  ic_state[3] = p_in / (gamma - 1.) + 0.5 * rho_in * (u_in * u_in + v_in * v_in);

  ic.apply_values("sol_dofs", ic_state,
                  *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  scheme->add_boundary_condition("WeakWall", "bottom_wall_bc", "bottom_wall");

  scheme->add_boundary_condition("WeakWall", "top_wall_bc", "top_wall");

  std::shared_ptr<scheme_type::bc_base_type> inlet_bc =
      scheme->add_boundary_condition("WeakSubInlet", "inlet_bc", "inlet");
  inlet_bc->set_parameter("total_temperature", Ttot);
  inlet_bc->set_parameter("total_pressure", Ptot);
  inlet_bc->print_bc_parameters();

  std::shared_ptr<scheme_type::bc_base_type> outlet_bc =
      scheme->add_boundary_condition("WeakSubOutlet", "outlet_bc", "outlet");
  outlet_bc->set_parameter("total_pressure", Pout);

#if 0

  mesh::MeshBoundarySet<MeshConfig> const &boundaries = mesh->topology().all_boundaries();

// ---

  BoundaryFacets<MeshConfig> const &bottom_wall_bdry = *(boundaries.domain("bottom_wall"));
  FunctionSpace<MeshConfig, _1D>::ptr bottom_wall_space =
      std::make_shared<FunctionSpace<MeshConfig, _1D>>();
  bottom_wall_space->set_reference_fe_values(bottom_wall_bdry, mesh->topology().dof_storage(),
                                             SFunc::Lagrange, quadrature_order, QuadType::Gauss);

  rdm::WeakWall<MeshConfig, physics::Euler2DCons> bottom_wall_bc("bottom_wall_bc");
  bottom_wall_bc.configure_mesh_data(
      topology.cells_ptr(), topology.all_boundaries_ptr(), mesh->geometry_ptr(),
      topology.dof_storage_ptr(), mesh->geometry_ptr(), topology.dof_storage_ptr(), "bottom_wall");
  bottom_wall_bc.configure_spaces(bottom_wall_space, bottom_wall_space, solution);
  bottom_wall_bc.set_residuals(residual);
  bottom_wall_bc.set_update_coeffs(update_coeff);

  // ---

  BoundaryFacets<MeshConfig> const &top_wall_bdry = *(boundaries.domain("top_wall"));
  FunctionSpace<MeshConfig, _1D>::ptr top_wall_space =
      std::make_shared<FunctionSpace<MeshConfig, _1D>>();
  top_wall_space->set_reference_fe_values(top_wall_bdry, mesh->topology().dof_storage(),
                                          SFunc::Lagrange, quadrature_order, QuadType::Gauss);

  rdm::WeakWall<MeshConfig, physics::Euler2DCons> top_wall_bc("top_wall_bc");
  top_wall_bc.configure_mesh_data(topology.cells_ptr(), topology.all_boundaries_ptr(),
                                  mesh->geometry_ptr(), topology.dof_storage_ptr(),
                                  mesh->geometry_ptr(), topology.dof_storage_ptr(), "top_wall");
  top_wall_bc.configure_spaces(top_wall_space, top_wall_space, solution);
  top_wall_bc.set_residuals(residual);
  top_wall_bc.set_update_coeffs(update_coeff);

  // ---

  BoundaryFacets<MeshConfig> const &inlet_bdry = *(boundaries.domain("inlet"));
  FunctionSpace<MeshConfig, _1D>::ptr inlet_bdry_space =
      std::make_shared<FunctionSpace<MeshConfig, _1D>>();
  inlet_bdry_space->set_reference_fe_values(inlet_bdry, mesh->topology().dof_storage(),
                                            SFunc::Lagrange, quadrature_order, QuadType::Gauss);

  rdm::WeakSubInlet<MeshConfig, physics::Euler2DCons> inlet_bc("inlet_bc");
  inlet_bc.configure_mesh_data(topology.cells_ptr(), topology.all_boundaries_ptr(),
                               mesh->geometry_ptr(), topology.dof_storage_ptr(),
                               mesh->geometry_ptr(), topology.dof_storage_ptr(), "inlet");
  inlet_bc.configure_spaces(inlet_bdry_space, inlet_bdry_space, solution);
  inlet_bc.set_residuals(residual);
  inlet_bc.set_update_coeffs(update_coeff);
  inlet_bc.set_parameter("total_temperature", Ttot);
  inlet_bc.set_parameter("total_pressure", Ptot);
  inlet_bc.print_bc_parameters();

  // ---

  BoundaryFacets<MeshConfig> const &outlet_bdry = *(boundaries.domain("outlet"));
  FunctionSpace<MeshConfig, _1D>::ptr outlet_bdry_space =
      std::make_shared<FunctionSpace<MeshConfig, _1D>>();
  outlet_bdry_space->set_reference_fe_values(outlet_bdry, mesh->topology().dof_storage(),
                                             SFunc::Lagrange, quadrature_order, QuadType::Gauss);

  rdm::WeakSubOutlet<MeshConfig, physics::Euler2DCons> outlet_bc("outlet_bc");
  outlet_bc.configure_mesh_data(topology.cells_ptr(), topology.all_boundaries_ptr(),
                                mesh->geometry_ptr(), topology.dof_storage_ptr(),
                                mesh->geometry_ptr(), topology.dof_storage_ptr(), "outlet");
  outlet_bc.configure_spaces(outlet_bdry_space, outlet_bdry_space, solution);
  outlet_bc.set_residuals(residual);
  outlet_bc.set_update_coeffs(update_coeff);
  outlet_bc.set_parameter("total_pressure", Pout);
#endif

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;

  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  start = clock();
  for (Uint iter = 0; iter < 1000; ++iter)
  {
    const Real CFL = 0.5;

    scheme->assemble_lhs_and_rhs(CFL);
    scheme->solve({});
    scheme->compute_residual_norm(res_L2_norm);

    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Write the output to file

  const std::string outfilename = "output_bump_rds_2D_explicit.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "u");

  /// Write Mach number values to file

  std::shared_ptr<interpolation::ScalarMeshFunction<Real>> Mach_number =
      std::make_shared<interpolation::ScalarMeshFunction<Real>>("", "Ma");

  Mach_number->resize(nb_nodes);

  solver::PostprocessingUtils<MeshConfig>::compute_Euler_Ma<phys_model>(*sol_dofs, *solution,
                                                                        *Mach_number);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *Mach_number, "Ma");

  return 0;
}
