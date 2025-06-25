#include <ctime>
#include <iostream>
#include <memory>

#include "common/PDEKit.hpp"
#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/euler/Euler3DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodImplicit.hpp"
#include "solver/rdm/PGRDMethodImplicit.hpp"
#include "solver/rdm/bc/WeakSubInlet.hpp"
#include "solver/rdm/bc/WeakSubOutlet.hpp"
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

class BoundaryData
{
  public:
  BoundaryData()
  {
  }

  ~BoundaryData()
  {
  }

  static Real inlet_density()
  {
    return p_in / (R * T_in);
  }

  static Real inlet_pressure()
  {
    return p_in;
  }

  static constexpr Real total_pressure()
  {
    return Ptot;
  }

  static constexpr Real inlet_temp_tot()
  {
    return Ttot;
  }

  static constexpr Real outlet_press_tot()
  {
    return Pout;
  }

  static void inlet_state(math::DenseDVec<Real> &in_state)
  {
    const Real rho_in = p_in / (R * T_in);
    const Real u_in   = rho_in * M_in * sqrt(gamma * R * T_in);
    const Real v_in   = 0.0;
    const Real w_in   = 0.0;

    in_state.resize(phys_model::NEQ);

    in_state[0] = rho_in;
    in_state[1] = rho_in * u_in;
    in_state[2] = rho_in * v_in;
    in_state[3] = rho_in * w_in;
    in_state[4] = p_in / (gamma - 1.) + 0.5 * rho_in * (u_in * u_in + v_in * v_in + w_in * w_in);
  }

  private:
  // Material values
  static constexpr Real gamma = 1.4;
  static constexpr Real R     = 287.046;

  // Inlet values
  static constexpr Real Ttot = 307.6488978;
  static constexpr Real Ptot = 120195.4453;
  static constexpr Real M_in = 0.5;

  // Outlet values
  static constexpr Real Pout = 101325.0;

  static constexpr Real T_in = Ttot / (1. + 0.5 * (gamma - 1) * M_in * M_in);
  static const Real p_in;
};

const Real BoundaryData::p_in =
    BoundaryData::Ptot /
    std::pow(1. + 0.5 * (BoundaryData::gamma - 1) * BoundaryData::M_in * BoundaryData::M_in,
             BoundaryData::gamma / (BoundaryData::gamma - 1));

// ----------------------------------------------------------------------------

std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> build_solver(
    const std::shared_ptr<MeshType const> &mesh,
    common::PtrHandle<MeshType::dof_storage_type> geo_dofs,
    common::PtrHandle<MeshType::dof_storage_type> sol_dofs, const Uint quad_order,
    vector_function::ptr solution, vector_function::ptr residuals, const std::string &scheme_type,
    const bool reorder_dofs)
{
  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  // typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA>
  // scheme_type;

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme;

  if (scheme_type == "LDA")
  {
    scheme = std::make_shared<rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA>>();
    // scheme = std::make_shared<rdm::DRDMethodImplicit<MeshConfig,
    // phys_model, rdm::PGLDA, rdm::FacetDG>>();
  }
  else
  {
    scheme = std::make_shared<rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGN>>();
    // scheme = std::make_shared<rdm::DRDMethodImplicit<MeshConfig,
    // phys_model, rdm::PGN, FacetDG>>();
  }

  if (reorder_dofs)
  {
    std::vector<Int> reordering;
    scheme->compute_node_reordering(*mesh, *sol_dofs, reordering);
    (*sol_dofs).renumber_dofs(reordering);
  }

  scheme->configure_mesh_data(mesh, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quad_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residuals);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  scheme->add_boundary_condition("WeakWall", "bottom_wall_bc", "Bottom");
  scheme->add_boundary_condition("WeakWall", "top_wall_bc", "Top");
  scheme->add_boundary_condition("WeakSymmetry", "left_symmetry_bc", "SymmetryLeft");
  scheme->add_boundary_condition("WeakSymmetry", "right_symmetry_bc", "SymmetryRight");

  std::shared_ptr<rdm::RDMBCBase<MeshConfig, phys_model>> inlet_bc =
      scheme->add_boundary_condition("WeakSubInlet", "inlet_bc", "Inlet");
  inlet_bc->set_parameter("total_temperature", BoundaryData::inlet_temp_tot());
  inlet_bc->set_parameter("total_pressure", BoundaryData::total_pressure());
  // inlet_bc->print_reference_values();

  std::shared_ptr<rdm::RDMBCBase<MeshConfig, phys_model>> outlet_bc =
      scheme->add_boundary_condition("WeakSubOutlet", "outlet_bc", "Outlet");
  outlet_bc->set_parameter("total_pressure", BoundaryData::outlet_press_tot());

  return scheme;
}

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

  MeshType::shared_ptr mesh3D = std::make_shared<MeshType>("mesh3D");

  const std::string infilename = "bump_p1_tet.msh";
  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh3D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh3D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh3D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_continuous(*mesh3D, *geo_dofs, *sol_dofs, P1,
                                               PointSetID::Warpblend);
  /*
  MeshType::dof_storage_type::clone_discontinuous(*mesh3D, *geo_dofs,
  *sol_dofs, P1, PointSetID::Warpblend);
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
  // BUILD SOLVER
  // ------------------------------------------------------

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> solver =
      build_solver(mesh3D, geo_dofs, sol_dofs, quadrature_order, solution, residual, "LDA", true);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  /// Initial condition

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh3D);
  ic.set_domain(_3D, "Fluid");
  math::DenseDVec<Real> ic_state(phys_model::NEQ);

  BoundaryData::inlet_state(ic_state);

  ic.apply_values("sol_dofs", ic_state,
                  *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  const std::clock_t cpu_start = std::clock();
  auto wall_start              = std::chrono::high_resolution_clock::now();

  for (Uint iter = 0; iter < 70; ++iter)
  {
    const Real CFL = iter < 10 ? 15.0 : std::min(1.e2, 20.0 + std::pow(1.50, iter - 10));

    solver->set_cfl(CFL);

    if (iter % 5 == 0)
    {
      solver->assemble_lhs_and_rhs(CFL);
      solver->solve({SolverOption::RecomputePreconditioner});
    }
    else
    {
      solver->assemble_rhs();
      solver->solve({});
    }

    solver->compute_residual_norm(res_L2_norm);

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

  const std::string outfilename = "output_bump_rds_3D_implicit.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh3D, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*mesh3D, outfilename, *solution, "u");

  /// Write Mach number values to file
  scalar_function::ptr Mach_number = std::make_shared<scalar_function>("", "Ma");
  Mach_number->resize(nb_nodes);

  phys_model physical_model;
  physics::Euler3DProperties properties;
  physics::Euler3DProperties::SolGradM gradient_matrix;

  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution, [](const phys_model::Properties &props) { return props.Ma; },
      *Mach_number);

  meshwriter.append_nodal_function_to_file(*mesh3D, outfilename, *Mach_number, "Ma");

  meshwriter.save_mesh_boundary(*mesh3D, *sol_dofs, "bump3D_Ma_wall.msh", *Mach_number, "Ma",
                                {"Bottom"});
  /*
  vtk::VtkWriter vtkwriter;
  vtkwriter.append_nodal_function_to_file(
      *mesh3D, *sol_dofs, "output_bump_continuous_rds_3D_implicit.vtu",
  *solution, "u");
  */

  // --------------------------------------------------------------------------
  // Compute entropy in the whole domain
  // --------------------------------------------------------------------------

  std::cout << "Evaluating entropy ... " << std::endl;
  interpolation::ScalarMeshFunction<Real> entropy_dev("", "entropy_dev");
  entropy_dev.resize((*sol_dofs).nb_nodes());

  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution,
      [](const phys_model::Properties &props) {
        return (props.P / BoundaryData::inlet_pressure()) /
                   (std::pow(props.rho / BoundaryData::inlet_density(), 1.4)) -
               1.;
      },
      entropy_dev);

  meshwriter.append_nodal_function_to_file(*mesh3D, outfilename, entropy_dev, "entropy");

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

  meshwriter.save_mesh_boundary(*mesh3D, *sol_dofs, "bump3D_entropy_wall.msh", entropy_dev, "ds",
                                {"Bottom"});

  // --------------------------------------------------------------------------
  mpi_env.finalize();
  // --------------------------------------------------------------------------

  return 0;
}
