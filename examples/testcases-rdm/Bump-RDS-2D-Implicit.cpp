#include <chrono>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodImplicit.hpp"
#include "solver/rdm/PGRDMethodImplicit.hpp"
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

int main(int argc, char *argv[])
{
  // ------------------------------------------------------
  // INITIALIZE ENVIRONMENT
  // ------------------------------------------------------

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // DEFINE METHOD TYPE
  // ------------------------------------------------------

  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGN> n_scheme_type;
  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA> lda_scheme_type;
  // typedef rdm::DRDMethodImplicit<MeshConfig, phys_model, rdm::PGN,
  // rdm::FacetDG> n_scheme_type; typedef rdm::DRDMethodImplicit<MeshConfig,
  // phys_model, rdm::PGLDA, rdm::FacetDG> lda_scheme_type;

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> n_scheme =
      std::make_shared<n_scheme_type>();
  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> lda_scheme =
      std::make_shared<lda_scheme_type>();

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("geo_mesh");

  const std::string infilename = "bump_p1_tri.msh";
  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  if (lda_scheme->is_continuous())
  {
    MeshType::dof_storage_type::clone_continuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                                 PointSetID::Warpblend);
  }
  else
  {
    MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                                    PointSetID::Warpblend);
  }

  /*
  /// TEST - REORDER NODES
  {
    result_of::dof_map_t<MeshConfig> &sol_dofs =
  sol_mesh->topology().dof_storage(); std::vector<Uint>
  new_dof_ids(sol_dofs.nb_nodes()); for (Uint i = 0; i < new_dof_ids.size();
  ++i)
    {
      new_dof_ids[i] = i;
    }
    std::random_shuffle(new_dof_ids.begin(), new_dof_ids.end());
    sol_dofs.renumber_dofs(new_dof_ids);
  }
  */

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

  n_scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  n_scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  n_scheme->set_vec_function(SolverVecFn::solution, solution);
  n_scheme->set_vec_function(SolverVecFn::residuals, residual);

  {
    /*
    interpolation::ScalarMeshFunction<Uint> cell_id("", "cell_id");
    cell_id.resize((*sol_dofs).nb_active_cells());
    for (Uint ac = 0; ac < (*sol_dofs).nb_active_cells(); ++ac)
    {
      const mesh::MeshEntity active_cell = (*sol_dofs).active_cell(ac);
      cell_id[ac] = active_cell.idx() + 1;
    }

    gmsh::GmshWriter writer;
    writer.write_mesh_to_file(*mesh2D, "sol_dofs", "before_reordering.msh");
    writer.append_cell_function_to_file(*mesh2D, "before_reordering.msh",
    cell_id, "cell_id_before_reordering");
    */

    std::vector<Int> reordering;
    if (lda_scheme->is_continuous())
    {
      lda_scheme->compute_node_reordering(*mesh2D, *sol_dofs, reordering);
      (*sol_dofs).renumber_dofs(reordering);
    }
    else
    {
      lda_scheme->compute_cell_reordering(*mesh2D, *sol_dofs, reordering);
      (*sol_dofs).renumber_dofs_blockwise(reordering);
    }
  }

  lda_scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  lda_scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  lda_scheme->set_vec_function(SolverVecFn::solution, solution);
  lda_scheme->set_vec_function(SolverVecFn::residuals, residual);

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

  // Boundary conditions for N scheme
  {
    n_scheme->add_boundary_condition("WeakWall", "bottom_wall_bc", "bottom_wall");
    n_scheme->add_boundary_condition("WeakWall", "top_wall_bc", "top_wall");

    std::shared_ptr<n_scheme_type::bc_base_type> inlet_bc =
        n_scheme->add_boundary_condition("WeakSubInlet", "inlet_bc", "inlet");
    inlet_bc->set_parameter("total_temperature", Ttot);
    inlet_bc->set_parameter("total_pressure", Ptot);
    // inlet_bc->print_reference_values();

    std::shared_ptr<n_scheme_type::bc_base_type> outlet_bc =
        n_scheme->add_boundary_condition("WeakSubOutlet", "outlet_bc", "outlet");
    outlet_bc->set_parameter("total_pressure", Pout);
  }

  lda_scheme->add_boundary_condition("WeakWall", "bottom_wall_bc", "bottom_wall");
  lda_scheme->add_boundary_condition("WeakWall", "top_wall_bc", "top_wall");

  std::shared_ptr<lda_scheme_type::bc_base_type> inlet_bc =
      lda_scheme->add_boundary_condition("WeakSubInlet", "inlet_bc", "inlet");
  inlet_bc->set_parameter("total_temperature", Ttot);
  inlet_bc->set_parameter("total_pressure", Ptot);
  // inlet_bc->print_reference_values();

  std::shared_ptr<lda_scheme_type::bc_base_type> outlet_bc =
      lda_scheme->add_boundary_condition("WeakSubOutlet", "outlet_bc", "outlet");
  outlet_bc->set_parameter("total_pressure", Pout);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  std::cout.precision(15);
  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  const Real CFLMax = lda_scheme->is_continuous() ? 50.0 : 20.0;

  // Pointer to active scheme - swich between N and LDA as needed
  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> active_scheme = lda_scheme;

  const std::clock_t cpu_start = std::clock();
  auto wall_start              = std::chrono::high_resolution_clock::now();

  for (Uint iter = 0; iter < 70; ++iter)
  {
    const Real CFL = iter < 10 ? 15.0 : std::min(CFLMax, 20.0 + std::pow(1.50, iter - 10));

    active_scheme->set_cfl(CFL);

    if (iter % 5 == 0)
    {
      active_scheme->assemble_lhs_and_rhs(CFL);
      active_scheme->solve({SolverOption::RecomputePreconditioner});
    }
    else
    {
      active_scheme->assemble_rhs();
      active_scheme->solve({});
    }
    active_scheme->compute_residual_norm(res_L2_norm);

    const std::clock_t cpu_end = std::clock();
    const double cpu_duration  = (cpu_end - cpu_start) / (double)CLOCKS_PER_SEC;
    solver::SolverIO::print_iter_and_res_norm_w_timing(iter, CFL, res_L2_norm, cpu_duration);

    // print_res_norm(iter, CFL, res_L2_norm);
  }

  const std::clock_t cpu_end = std::clock();
  auto wall_end              = std::chrono::high_resolution_clock::now();

  const double cpu_duration = (cpu_end - cpu_start) / (double)CLOCKS_PER_SEC;
  auto wclock_duration = std::chrono::duration<double, std::milli>(wall_end - wall_start).count();

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop)        = " << cpu_duration << " s" << std::endl;
  std::cout << "Wall clock time (main iteration loop) = " << wclock_duration * 1.e-3 << " s"
            << std::endl;

  /// Write the output to file

  const std::string outfilename = "output_bump_rds_2D_implicit.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "u");

  /// Write Mach number values to file
  scalar_function::ptr Mach_number = std::make_shared<scalar_function>("", "Ma");
  Mach_number->resize(nb_nodes);

  physics::Euler2DCons physical_model;
  physics::Euler2DProperties properties;
  physics::Euler2DProperties::SolGradM gradient_matrix;

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

  meshwriter.append_nodal_function_to_file((*mesh2D), outfilename, *Mach_number, "Ma");
  meshwriter.append_nodal_function_to_file(
      (*mesh2D), outfilename, *lda_scheme->vec_function(SolverVecFn::residuals), "residuals");

  // --------------------------------------------------------------------------
  // Solution and entropy on bottom wall
  // --------------------------------------------------------------------------

  interpolation::VectorMeshFunction<Real> bottom_wall_coords("", "bottom_coords");
  interpolation::VectorMeshFunction<Real> bottom_wall_solution_values("", "bottom_solution");

  solver::extract_boundary_data(*mesh2D, *sol_dofs,
                                *active_scheme->vec_function(SolverVecFn::solution), "bottom_wall",
                                bottom_wall_coords, bottom_wall_solution_values);

  std::ofstream outstream;
  outstream.open("bump_bottom_wall_data.dat");
  outstream.precision(15);
  outstream.setf(std::ios::fixed);

  for (Uint i = 0; i < bottom_wall_coords.nb_entries(); ++i)
  {
    typedef interpolation::VectorMeshFunction<Real>::const_entry_type entry_type;

    const entry_type coords = bottom_wall_coords.const_value(i);
    const entry_type values = bottom_wall_solution_values.const_value(i);

    physical_model.compute_properties(coords, values, gradient_matrix, properties);

    const Real ds = (properties.P / p_in) / (std::pow(properties.rho / rho_in, 1.4)) - 1.;

    outstream << std::setw(15) << coords[X0] << " " << std::setw(15) << coords[X1] << " "
              << std::setw(15) << values[0] << " " << std::setw(15) << values[1] << " "
              << std::setw(15) << values[2] << " " << std::setw(15) << values[3] << " "
              << std::setw(15) << ds << std::endl;
  }
  outstream.close();

  // --------------------------------------------------------------------------
  // Compute entropy in the whole domain
  // --------------------------------------------------------------------------

  std::cout << "Evaluating entropy ... " << std::endl;
  interpolation::ScalarMeshFunction<Real> entropy_dev("", "entropy_dev");
  entropy_dev.resize((*sol_dofs).nb_nodes());

  /*
  PostprocessingUtils<MeshConfig>::compute_Euler_entropy<phys_model>(*sol_dofs, *solution, rho_in,
                                                                     p_in, entropy_dev);
  */

  PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution,
      [rho_in, p_in](const phys_model::Properties &props) {
        return (props.P / p_in) / (std::pow(props.rho / rho_in, 1.4)) - 1.;
      },
      entropy_dev);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, entropy_dev, "entropy");

  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto quad_generator = [=](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order, PointSetID::Gauss);
  };

  interpolation::FunctionSpace<MeshConfig> geo_cell_space;
  geo_cell_space.set_reference_fe_values((*geo_dofs).as_range(), sf_generator, quad_generator);
  interpolation::FunctionSpace<MeshConfig> sol_cell_space;
  sol_cell_space.set_reference_fe_values((*sol_dofs).as_range(), sf_generator, quad_generator);

  const Real entropy_norm = PostprocessingUtils<MeshConfig>::compute_function_norm(
      *geo_dofs, *sol_dofs, geo_cell_space, sol_cell_space, entropy_dev);

  std::cout.precision(15);
  std::cout.setf(std::ios::fixed);
  std::cout << "Entropy L2 norm = " << entropy_norm << std::endl;

  mpi_env.finalize();

  return 0;
}
