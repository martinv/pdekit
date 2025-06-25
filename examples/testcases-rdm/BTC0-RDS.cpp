#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/euler/Euler3DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodImplicit.hpp"
#include "solver/rdm/PGRDMethodImplicit.hpp"
#include "solver/rdm/bc/WeakFarfield.hpp"
#include "solver/rdm/bc/WeakSymmetry.hpp"
#include "solver/rdm/bc/WeakWall.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"

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

int main(int argc, char *argv[])
{
  // ------------------------------------------------------
  // INITIALIZE ENVIRONMENT
  // ------------------------------------------------------

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // DEFINE METHOD TYPE
  // ------------------------------------------------------

  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGN> NSchemeType;
  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA> LDASchemeType;
  // typedef rdm::DRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA,
  // rdm::FacetDG> LDASchemeType;
  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::SUPG> SUPGSchemeType;
  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGB> BSchemeType;

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme(new LDASchemeType());

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh3D = std::make_shared<MeshType>("mesh3D");

  const std::string infilename = "btc0-mesh-p2-coarse.msh";
  gmsh::GmshReader meshreader;
  gmsh::GmshWriter mshwriter;
  // vtk::VtkWriter vtuwriter;
  meshreader.read_mesh_from_file(infilename, *mesh3D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh3D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh3D->create_dof_storage("sol_dofs");

  if (scheme->is_continuous())
  {
    /*
    MeshType::dof_storage_type::clone_continuous(*mesh3D, *geo_dofs,
    *sol_dofs, P2, PointSetID::Equidist);
    */

    MeshType::dof_storage_type::make_identical_copy(*geo_dofs, *sol_dofs);
  }
  else
  {
    MeshType::dof_storage_type::clone_discontinuous(*mesh3D, *geo_dofs, *sol_dofs, P2,
                                                    PointSetID::Equidist);
  }

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P4;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  const Uint nb_nodes = (*sol_dofs).nb_nodes();

  vector_function::ptr solution = std::make_shared<vector_function>("", "solution");
  vector_function::ptr residual = std::make_shared<vector_function>("", "residual");

  solution->resize(phys_model::NEQ, nb_nodes);
  residual->resize(phys_model::NEQ, nb_nodes);

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  {
    // scheme->reorder_dofs(*mesh3D, *sol_dofs);

    std::vector<Int> reordering;
    scheme->compute_node_reordering(*mesh3D, *sol_dofs, reordering);
    (*sol_dofs).renumber_dofs(reordering);
  }

  scheme->configure_mesh_data(mesh3D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residual);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh3D);
  ic.set_domain(_3D, "Fluid");
  math::DenseDVec<Real> ic_state(phys_model::NEQ);

  const Real gamma    = 1.4;
  const Real alpha_in = 1.;
  const Real rho_in   = 1.22503;   // 1.0;
  const Real rhou_in  = 208.30559; // 1.0;
  const Real rhov_in  = 0.0;
  const Real rhow_in  = rhou_in * std::tan(alpha_in * math::pi / 180);

  const Real v2_in =
      (rhou_in * rhou_in + rhov_in * rhov_in + rhow_in * rhow_in) / (rho_in * rho_in);

  const Real M_in = 0.5;
  const Real p_in = rho_in / gamma * v2_in / (M_in * M_in);
  const Real e_in = p_in / (gamma - 1.) + 0.5 * rho_in * v2_in;

  ic_state[0] = rho_in;
  ic_state[1] = rhou_in;
  ic_state[2] = rhov_in;
  ic_state[3] = rhow_in;
  ic_state[4] = e_in;

  std::cout << "Initial state: " << ic_state << std::endl;

  // Initialize by inlet state
  ic.apply_values("sol_dofs", ic_state,
                  *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  std::shared_ptr<LDASchemeType::bc_base_type> wall1_bc =
      scheme->add_boundary_condition("WeakWall", "wall1_bc", "WallBottom");

  // This is boundary condition at the fuselage end - wall bc does not work
  // here We prescribe symmetry instead
  std::shared_ptr<LDASchemeType::bc_base_type> wall2_bc =
      scheme->add_boundary_condition("WeakWall", "wall2_bc", "WallTop");

  std::shared_ptr<LDASchemeType::bc_base_type> farfield_bc =
      scheme->add_boundary_condition("WeakFarfield", "farfield_bc", "Farfield");
  farfield_bc->set_reference_state(ic_state);
  farfield_bc->print_bc_parameters();

  // Global residual for convergence monitoring
  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  // Variables for writing Mach number values to file
  scalar_function Mach_number("", "Ma");
  Mach_number.resize(nb_nodes);

  // Variables for writing pressure distribution to file
  scalar_function Cp("", "Cp");
  Cp.resize(nb_nodes);

  phys_model physical_model;
  physics::Euler3DProperties properties;
  physics::Euler3DProperties::SolGradM gradient_matrix;

  // Variables for timing
  clock_t start, end;
  Real elapsed;

  std::cout.precision(15);

  auto t_wall_start = std::chrono::high_resolution_clock::now();

  /// Main iteration loop
  start = clock();

  for (Uint iter = 0; iter < 5; ++iter)
  {
    const Real CFL = (iter < 100) ? 10.0 : 20.0;

    if (iter % 5 == 0)
    {
      scheme->assemble_lhs_and_rhs(CFL);
      scheme->solve({SolverOption::RecomputePreconditioner});
    }
    else
    {
      scheme->assemble_rhs();
      scheme->solve({});
    }

    scheme->compute_residual_norm(res_L2_norm);
    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);

    // Write the output to file
    if (iter % 100 == 0) // 5000
    {
      // Compute Mach number
      solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
          *sol_dofs, *solution, [](const phys_model::Properties &props) { return props.Ma; },
          Mach_number);

      // Compute pressure coefficient
      solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
          *sol_dofs, *solution,
          [rho_in, p_in, v2_in](const phys_model::Properties &props) {
            return (props.P - p_in) / (0.5 * rho_in * v2_in);
          },
          Cp);

      std::stringstream ss;
      ss.str(std::string());
      ss << std::setfill('0') << std::setw(5) << iter;

      // Output writing

      const std::string msh_outfilename = "output_btc0_rds_" + ss.str() + "_iter.msh";
      mshwriter.write_mesh_to_file(*mesh3D, "sol_dofs", msh_outfilename);
      mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, *solution, "u");
      mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, Mach_number, "Ma");
      mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, Cp, "Cp");

      mshwriter.save_mesh_boundary(*mesh3D, *sol_dofs, "btc0_wall.msh", *solution, "u",
                                   {"WallBottom", "WallTop"});

      /*
      const std::string vtk_outfilename = "output_onera_m6_tet_p2_" +
      ss.str() +
      "_iter.vtu";
      vtuwriter.write_data_to_file(*mesh,vtk_outfilename,*solution,"solution");
      */
    }

  } // Loop for iterations

  end = clock();

  auto t_wall_end = std::chrono::high_resolution_clock::now();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;
  std::cout << "Wall clock time (main iteration loop) = "
            << std::chrono::duration_cast<std::chrono::seconds>(t_wall_end - t_wall_start).count()
            << " s" << std::endl;

  const std::string msh_outfilename = "output_btc0_rds_final.msh";
  mshwriter.write_mesh_to_file(*mesh3D, "sol_dofs", msh_outfilename);
  mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, *solution, "u");

  // Compute Mach number
  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution, [](const phys_model::Properties &props) { return props.Ma; },
      Mach_number);

  // Compute pressure coefficient
  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution,
      [rho_in, p_in, v2_in](const phys_model::Properties &props) {
        return (props.P - p_in) / (0.5 * rho_in * v2_in);
      },
      Cp);

  mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, Mach_number, "Ma");
  mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, Cp, "Cp");

  mshwriter.save_mesh_boundary(*mesh3D, *sol_dofs, "btc0_u_wall.msh", *solution, "u",
                               {"WallBottom", "WallTop"});
  mshwriter.save_mesh_boundary(*mesh3D, *sol_dofs, "btc0_cp_wall.msh", Cp, "Cp",
                               {"WallBottom", "WallTop"});

  // --------------------------------------------------------------------------
  // Compute entropy in the whole domain
  // --------------------------------------------------------------------------

  std::cout << "Evaluating entropy ... " << std::endl;
  interpolation::ScalarMeshFunction<Real> entropy_dev("", "entropy_dev");
  entropy_dev.resize((*sol_dofs).nb_nodes());

  // Compute Entropy
  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution,
      [rho_in, p_in](const phys_model::Properties &props) {
        return (props.P / p_in) / (std::pow(props.rho / rho_in, 1.4)) - 1.;
      },
      entropy_dev);

  mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, entropy_dev, "entropy");

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

  mshwriter.save_mesh_boundary(*mesh3D, *sol_dofs, "btc0_entropy_wall.msh", entropy_dev, "ds",
                               {"WallBottom", "WallTop"});

  mpi_env.finalize();

  return 0;
}
