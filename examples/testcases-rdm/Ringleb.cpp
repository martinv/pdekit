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

  const std::string infilename = "ringleb_p1.msh";
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

  mesh::adapt::LocalInterpolator loc_interpolator;
  std::vector<Real> tmp(phys_model::NEQ);

  for (mesh::DofMap<MeshConfig>::const_dof_iterator it = (*sol_dofs).cbegin();
       it != (*sol_dofs).cend(); ++it)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = it->tcell();
    const mesh::MeshEntity cell                         = it->mesh_entity();

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < cell_coords.rows(); ++n)
    {
      interpolation::VectorMeshFunction<Real>::entry_type u_node =
          (*solution).value(cell.vertex(n));

      u_node[0] = 1.2; // rho
      u_node[1] = 0.0;
      u_node[2] = -0.1 * std::abs(cell_coords.row(n)[X1]);
      u_node[3] = 253318.0;
    }
  }

  /*
  // Initialize the field
  for (Uint i = 0; i < solution->nb_entries(); ++i)
  {
    interpolation::VectorMeshFunction<Real>::entry_type node_value =
  solution->value(i); node_value[0] = 1.2;    // rho node_value[1] = 0.01; //
  rho * u node_value[2] = 0.01;   //-1.2; // rho * v node_value[3] = 253318;
  //
  e
  }
  */

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  // Boundary conditions for N scheme
  {
    n_scheme->add_boundary_condition("WeakWall", "inner_wall_bc", "inner_wall");
    n_scheme->add_boundary_condition("WeakWall", "outer_wall_bc", "outer_wall");

    std::shared_ptr<n_scheme_type::bc_base_type> inlet_bc =
        n_scheme->add_boundary_condition("WeakDirichlet", "inlet_bc", "inlet");
    inlet_bc->set_expression(Ringleb2DSolution::value);

    std::shared_ptr<n_scheme_type::bc_base_type> outlet_bc =
        n_scheme->add_boundary_condition("WeakDirichlet", "outlet_bc", "outlet");
    outlet_bc->set_expression(Ringleb2DSolution::value);
  }

  lda_scheme->add_boundary_condition("WeakWall", "inner_wall_bc", "inner_wall");
  lda_scheme->add_boundary_condition("WeakWall", "outer_wall_bc", "outer_wall");

  std::shared_ptr<lda_scheme_type::bc_base_type> inlet_bc =
      lda_scheme->add_boundary_condition("WeakDirichlet", "inlet_bc", "inlet");
  inlet_bc->set_expression(Ringleb2DSolution::value);

  std::shared_ptr<lda_scheme_type::bc_base_type> outlet_bc =
      lda_scheme->add_boundary_condition("WeakDirichlet", "outlet_bc", "outlet");
  outlet_bc->set_expression(Ringleb2DSolution::value);

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

  const std::string outfilename = "output_ringleb_2D.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "u");

  /// Write Mach number values to file
  scalar_function::ptr Mach_number = std::make_shared<scalar_function>("", "Ma");
  Mach_number->resize(nb_nodes);

  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution, [](const phys_model::Properties &props) { return props.Ma; },
      *Mach_number);

  meshwriter.append_nodal_function_to_file((*mesh2D), outfilename, *Mach_number, "Ma");
  meshwriter.append_nodal_function_to_file(
      (*mesh2D), outfilename, *lda_scheme->vec_function(SolverVecFn::residuals), "residuals");

  // --------------------------------------------------------------------------
  // Compute entropy in the whole domain
  // --------------------------------------------------------------------------

  /*
  std::cout << "Evaluating entropy ... " << std::endl;
  interpolation::VectorMeshFunction<Real> entropy_dev("", "entropy_dev");
  entropy_dev.resize(1, (*sol_dofs).nb_nodes());

  for (Uint c = 0; c < (*sol_dofs).nb_active_cells(); ++c)
  {
    const mesh::MeshEntity active_cell =
  (*sol_dofs).active_cell(mesh::ActiveIdx(c)); const
  mesh::DofCoordinates<MeshConfig::GDIM> active_coords =
        (*sol_dofs).active_cell_coords(active_cell);

    for (Uint n = 0; n < active_cell.nb_vert(); ++n)
    {
      physical_model.compute_properties(active_coords.c(n),
                                        solution->const_value(active_cell.vertex(n)),
                                        gradient_matrix, properties);
      interpolation::VectorMeshFunction<Real>::entry_type node_value =
          entropy_dev.value(active_cell.vertex(n));

      node_value[0] = (properties.P / p_in) / (std::pow(properties.rho /
  rho_in, 1.4)) - 1.;
    }
  }

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, entropy_dev,
  "entropy");


  interpolation::FunctionSpace<MeshConfig> geo_cell_space;
  geo_cell_space.set_reference_fe_values(*geo_dofs, SFunc::Lagrange,
  quadrature_order, PointSetID::Gauss);
  interpolation::FunctionSpace<MeshConfig> sol_cell_space;
  sol_cell_space.set_reference_fe_values(*sol_dofs, SFunc::Lagrange,
  quadrature_order, PointSetID::Gauss);

  const Real entropy_norm =
      SolverUtils<MeshConfig>::compute_function_norm(*geo_dofs, *sol_dofs, geo_cell_space,
  sol_cell_space, entropy_dev);

  std::cout.precision(15);
  std::cout.setf(std::ios::fixed);
  std::cout << "Entropy L2 norm = " << entropy_norm << std::endl;
  */

  mpi_env.finalize();

  return 0;
}
