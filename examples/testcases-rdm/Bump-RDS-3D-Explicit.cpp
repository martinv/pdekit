#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler3DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/bc/WeakSubInlet.hpp"
#include "solver/rdm/bc/WeakSubOutlet.hpp"
#include "solver/rdm/bc/WeakWall.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::interpolation;
using namespace pdekit::solver;

typedef Cart3D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::Euler3DCons phys_model;

// ----------------------------------------------------------------------------

int main()
{
  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh3D = std::make_shared<MeshType>("mesh3D");

  const std::string infilename = "bump_p1_tet.msh";
  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh3D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh3D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh3D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_continuous(*mesh3D, *geo_dofs, *sol_dofs, P1,
                                               PointSetID::Warpblend);

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

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  typedef rdm::PGRDMethodExplicit<MeshConfig, physics::Euler3DCons, rdm::PGLDA> scheme_type;
  scheme_type scheme;

  scheme.configure_mesh_data(mesh3D, geo_dofs, sol_dofs);
  scheme.initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme.set_vec_function(SolverVecFn::solution, solution);
  scheme.set_vec_function(SolverVecFn::residuals, residual);

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

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh3D);
  ic.set_domain(_3D, "Fluid");
  math::DenseDVec<Real> ic_state(physics::Euler3DCons::NEQ);

  const Real rho_in = p_in / (R * T_in);
  const Real u_in   = rho_in * M_in * sqrt(gamma * R * T_in);
  const Real v_in   = 0.0;
  const Real w_in   = 0.0;

  ic_state[0] = rho_in;
  ic_state[1] = rho_in * u_in;
  ic_state[2] = rho_in * v_in;
  ic_state[3] = rho_in * w_in;
  ic_state[4] = p_in / (gamma - 1.) + 0.5 * rho_in * (u_in * u_in + v_in * v_in + w_in * w_in);

  ic.apply_values("sol_dofs", ic_state,
                  *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  scheme.add_boundary_condition("WeakWall", "bottom_wall_bc", "Bottom");
  scheme.add_boundary_condition("WeakWall", "top_wall_bc", "Top");
  scheme.add_boundary_condition("WeakSymmetry", "left_wall_bc", "SymmetryLeft");
  scheme.add_boundary_condition("WeakSymmetry", "right_wall_bc", "SymmetryRight");

  std::shared_ptr<scheme_type::bc_base_type> inlet_bc =
      scheme.add_boundary_condition("WeakSubInlet", "inlet_bc", "Inlet");
  inlet_bc->set_parameter("total_temperature", Ttot);
  inlet_bc->set_parameter("total_pressure", Ptot);
  inlet_bc->print_bc_parameters();

  std::shared_ptr<scheme_type::bc_base_type> outlet_bc =
      scheme.add_boundary_condition("WeakSubOutlet", "outlet_bc", "Outlet");
  outlet_bc->set_parameter("total_pressure", Pout);

  /*
  // Boundary conditions for the multi-threaded case
  // Add boundary conditions to thread 0
  scheme_type &worker0 = work_stream.worker(0);

  worker0.add_boundary_condition("WeakWall", "bottom_wall_bc", "Bottom");
  worker0.add_boundary_condition("WeakWall", "top_wall_bc", "Top");
  worker0.add_boundary_condition("WeakWall", "left_wall_bc", "SymmetryLeft");
  worker0.add_boundary_condition("WeakWall", "right_wall_bc",
  "SymmetryRight");

  std::shared_ptr<scheme_type::bc_base_type> inlet_bc =
      worker0.add_boundary_condition("WeakSubInlet", "inlet_bc", "Inlet");
  inlet_bc->set_parameter("total_temperature", Ttot);
  inlet_bc->set_parameter("total_pressure", Ptot);
  inlet_bc->print_bc_parameters();

  std::shared_ptr<scheme_type::bc_base_type> outlet_bc =
      worker0.add_boundary_condition("WeakSubOutlet", "outlet_bc", "Outlet");
  outlet_bc->set_parameter("total_pressure", Pout);
  */

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;

  math::DenseDVec<Real> residual_norm(phys_model::NEQ);

  start = clock();

  for (Uint iter = 0; iter < 100; ++iter)
  {
    const Real CFL = 0.5;
    scheme.assemble_lhs_and_rhs(CFL);
    scheme.solve({});

    if (iter % 5 == 0)
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

  /// Write the output to file

  const std::string outfilename = "output_bump_tet_p1.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh3D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh3D, outfilename, *solution, "u");

  /// Write Mach number values to file
  scalar_function::ptr Mach_number = std::make_shared<scalar_function>("", "Ma");
  Mach_number->resize(nb_nodes);

  physics::Euler3DCons physical_model;
  physics::Euler3DProperties properties;
  physics::Euler3DProperties::SolGradM gradient_matrix;

  /*
  for (Uint c = 0; c < (*sol_dofs).nb_active_cells(); ++c)
  {
    const mesh::MeshEntity active_cell = (*sol_dofs).active_cell(mesh::ActiveIdx(c));
    const mesh::DofCoordinates<MeshConfig::GDIM> active_coords =
        (*sol_dofs).active_cell_coords(active_cell);

    for (Uint n = 0; n < active_cell.nb_vert(); ++n)
    {
      physical_model.compute_properties(active_coords.c(n),
                                        solution->const_value(active_cell.vertex(n)),
                                        gradient_matrix, properties);
      (*Mach_number)[active_cell.vertex(n)] = properties.Ma;
    }
  }
  */

  PostprocessingUtils<MeshConfig>::compute_Euler_Ma<phys_model>(*sol_dofs, *solution, *Mach_number);

  meshwriter.append_nodal_function_to_file(*mesh3D, outfilename, *Mach_number, "Ma");

  meshwriter.save_mesh_boundary(*mesh3D, *sol_dofs, "bump3D_Ma_wall.msh", *Mach_number, "Ma",
                                {"Bottom"});

  // --------------------------------------------------------------------------
  // Compute entropy in the whole domain
  // --------------------------------------------------------------------------

  std::cout << "Evaluating entropy ... " << std::endl;
  interpolation::ScalarMeshFunction<Real> entropy_dev("", "entropy_dev");
  entropy_dev.resize((*sol_dofs).nb_nodes());

  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution,
      [rho_in, p_in](const phys_model::Properties &props) {
        return (props.P / p_in) / (std::pow(props.rho / rho_in, 1.4)) - 1.;
      },
      entropy_dev);

  meshwriter.append_nodal_function_to_file(*mesh3D, outfilename, entropy_dev, "entropy");

  interpolation::FunctionSpace<MeshConfig> geo_cell_space;
  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto quad_generator = [=](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order, PointSetID::Gauss);
  };

  geo_cell_space.set_reference_fe_values((*geo_dofs).as_range(), sf_generator, quad_generator);
  interpolation::FunctionSpace<MeshConfig> sol_cell_space;
  sol_cell_space.set_reference_fe_values((*sol_dofs).as_range(), sf_generator, quad_generator);

  const Real entropy_norm = PostprocessingUtils<MeshConfig>::compute_function_norm(
      *geo_dofs, *sol_dofs, geo_cell_space, sol_cell_space, entropy_dev);

  std::cout.precision(15);
  std::cout.setf(std::ios::fixed);
  std::cout << "Entropy L2 norm = " << entropy_norm << std::endl;

  meshwriter.save_mesh_boundary(*mesh3D, *sol_dofs, "bump3D_entropy_wall.msh", entropy_dev, "ds",
                                {"Bottom"});

  return 0;
}
