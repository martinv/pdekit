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
#include "solver/rdm/PGRDMethodImplicit.hpp"
#include "solver/rdm/bc/WeakFarfield.hpp"
#include "solver/rdm/bc/WeakSymmetry.hpp"
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

int main(int argc, char *argv[])
{
  // Settings to make the simulation run:
  // Use "M6_unstr_euler_c2_p2_rescaled.msh"
  // The quadrature order should be P4
  // SUPG scheme works, LDA is not stable enough
  // The interpolation order for solution and flux should be the same (to be
  // verified)

  // ------------------------------------------------------
  // INITIALIZE ENVIRONMENT
  // ------------------------------------------------------

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh3D = std::make_shared<MeshType>("mesh3D");

  const std::string infilename = "M6_unstr_euler_c2_p2_rescaled.msh";
  gmsh::GmshReader meshreader;
  gmsh::GmshWriter mshwriter;
  // vtk::VtkWriter vtuwriter;
  meshreader.read_mesh_from_file(infilename, *mesh3D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh3D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh3D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_continuous(*mesh3D, *geo_dofs, *sol_dofs, P2,
                                               PointSetID::Warpblend);

  // const result_of::dof_handler<MeshConfig>::type &geo_dofs =
  // mesh->topology().dof_storage();

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

  //  interpolation::FunctionSpace::Function::Ptr f_ptr =
  // fs_solution->function("residual");
  //  std::cout << "Got function " << f_ptr->name() << std::endl;

  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGN> NSchemeType;
  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA> LDASchemeType;
  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::SUPG> SUPGSchemeType;
  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGB> BSchemeType;

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme =
      std::make_shared<SUPGSchemeType>();

  scheme->configure_mesh_data(mesh3D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residual);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh3D);
  ic.set_domain(_3D, "Fluid");
  math::DenseDVec<Real> ic_state(phys_model::NEQ);

  // This corresponds to angle of attack = 3.06 and Mach number on inlet Ma_in
  // = 0.84
  const Real alpha = 3.06;
  // const Real tan_alpha = 0.05345791105765804881; // tan(3.06*PI/180)
  const Real Ma_in = 0.8395;
  const Real kappa = 1.4;

  const Real rho    = 1.4;
  const Real u      = 1.0; // 1.7432/rho;
  const Real v      = 0.0;
  const Real w      = u * std::tan(PI * alpha / 180.0); // u*tan_alpha;
  const Real speed2 = u * u + v * v + w * w;

  const Real e = speed2 / (Ma_in * Ma_in) * rho / (kappa * (kappa - 1.)) + 0.5 * rho * speed2;

  ic_state[0] = rho;
  ic_state[1] = rho * u;
  ic_state[2] = rho * v;
  ic_state[3] = rho * w;
  ic_state[4] = e;

  phys_model::Properties props;

  math::DenseSVec<Real, _3D> coord;
  math::DenseSMat<Real, 1, 5> sol;
  sol.insert_row(0, ic_state);
  phys_model::SolGradM grad_vars;
  grad_vars.fill(0.0);

  phys_model::compute_properties(coord, sol.const_row_transp(0), grad_vars, props);

  std::cout << "Inlet mach number = " << props.Ma << std::endl;
  std::cout << "Inlet density     = " << props.rho << std::endl;
  std::cout << "angle = " << std::atan(w / u) * 180.0 / PI << " deg" << std::endl;
  std::cout << "ic_state = " << ic_state << std::endl;

  /*
  // This was working
  ic_state[0] = 1.40000;
  ic_state[1] = 1.17432;
  ic_state[2] = 0.0;
  ic_state[3] = 0.06278;
  ic_state[4] = 2.99392;
  */

  // solution->fill(0.0);
  // Initialize by inlet state
  ic.apply_values("sol_dofs", ic_state,
                  *solution); // The inlet function is also used as bc

  // Initialize from file
  // meshreader.read_nodal_function_from_file(*mesh,
  // "onera_m6_tet_p2_9000_iter_N.msh",
  //                                         solution->data(),
  // {"\"u_0\"","\"u_\"1","\"u_2\"","\"u_3\"","\"u_4\""});

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  typedef rdm::RDMBCBase<MeshConfig, phys_model, phys_model::DIM - 1> bc_base_type;

  scheme->add_boundary_condition("WeakWall", "wing_bc", "Wing");

  /*
  std::shared_ptr<bc_base_type> symmetry_bc =
      scheme->add_boundary_condition("WeakSymmetry", "symmetry_bc",
  "Symmetry");
  */

  scheme->add_boundary_condition("WeakSymmetry", "symmetry_bc", "Symmetry");

  std::shared_ptr<bc_base_type> farfield_bc =
      scheme->add_boundary_condition("WeakFarfield", "farfield_bc", "Farfield");
  farfield_bc->set_reference_state(ic_state);
  farfield_bc->print_bc_parameters();

  // ------------------------

  // Global residual for convergence monitoring
  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  // Variables for writing Mach number values to file
  scalar_function Mach_number("", "Ma");
  Mach_number.resize(nb_nodes);

  phys_model physical_model;
  physics::Euler3DProperties properties;
  physics::Euler3DProperties::SolGradM gradient_matrix;

  // Variables for timing
  clock_t start, end;
  Real elapsed;

  std::cout.precision(15);

  auto t_wall_start = std::chrono::high_resolution_clock::now();

  /// Main iteration loop
  typedef vector_function::entry_type node_value_type;
  typedef vector_function::const_entry_type const_node_value_type;

  start = clock();

  for (Uint iter = 0; iter <= 10; ++iter)
  {
    const Real CFL = 10.0;

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

    if (iter % 10 == 0)
    {
      res_L2_norm.fill(0.0);
      norm_L2(*residual, res_L2_norm);
      solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);
    }

    // Write the output to file
    if (iter % 5000 == 0) // 5000
    {
      solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
          *sol_dofs, *solution, [](const phys_model::Properties &props) { return props.Ma; },
          Mach_number);

      std::stringstream ss;
      ss.str(std::string());
      ss << std::setfill('0') << std::setw(5) << iter;

      // Output writing

      const std::string msh_outfilename = "output_onera_m6_tet_p2_" + ss.str() + "_iter.msh";
      mshwriter.write_mesh_to_file(*mesh3D, "sol_dofs", msh_outfilename);
      mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, *solution, "u");
      mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, Mach_number, "Ma");

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

  const std::string msh_outfilename = "output_onera_m6_tet_p2_final.msh";
  mshwriter.write_mesh_to_file(*mesh3D, "sol_dofs", msh_outfilename);
  mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, *solution, "u");
  mshwriter.append_nodal_function_to_file(*mesh3D, msh_outfilename, Mach_number, "Ma");

  // Save the aircraft surface and symmetry surface to a separate file
  mshwriter.save_mesh_boundary(*mesh3D, *sol_dofs, "onera_m6_wall.msh", *solution, "u",
                               {"Wing", "Symmetry"});

  /*
  const std::string vtk_outfilename = "output_onera_m6_tet_p2_final.vtu";
  vtuwriter.write_data_to_file(*mesh,vtk_outfilename,*solution,"solution");
  */

  mpi_env.finalize();

  return 0;
}
