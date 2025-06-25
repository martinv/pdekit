#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/L2ProjectionGlobal.hpp"
#include "interpolation/L2ProjectionLocal.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/art_visc/ArtificialViscosity.hpp"
#include "solver/rdm/DRDMethodImplicit.hpp"
#include "solver/rdm/PGRDMethodImplicit.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"
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

#define TEST_CASE_IS_TRANSONIC 0

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

  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;
  // typedef rdm::DRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA,
  // rdm::FacetDG> scheme_type;

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme = std::make_shared<scheme_type>();

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
  // DISCONTINUITY SENSOR SETUP - FOR TRANSONIC FLOW
  // ------------------------------------------------------

#if TEST_CASE_IS_TRANSONIC

  solver::ArtificialViscosity<MeshConfig> art_visc_operator;
  art_visc_operator.prepare_indicator_operators(*geo_dofs, *sol_dofs);

  scalar_function::ptr indicator = std::make_shared<scalar_function>("", "smoothness_indicator");
  indicator->resize((*sol_dofs).nb_active_cells());

  scalar_function::ptr art_visc = std::make_shared<scalar_function>("", "artificial_viscosity");
  art_visc->resize((*sol_dofs).nb_nodes());

  // ------------------------------------------------------
  // BLENDING COEFFICIENT FOR NONLINEAR SCHEMES
  // ------------------------------------------------------

  rdm::RDBlendingCoeff<MeshConfig> rd_blending_coeff_operator;

  // Types of viscosity coefficient: "LF_Blend", "ViscosityIndicator",
  // "WongJansen", "Bonanni2D"
  const std::string blend_coeff_type = "Bonanni2D";
  rd_blending_coeff_operator.setup(*geo_dofs, *sol_dofs, blend_coeff_type);

  // This is reference value for Wong-Jansen blending coefficient
  rd_blending_coeff_operator.set_param("u_inf", 1.0);

  scalar_function::ptr blending_coeff =
      std::make_shared<scalar_function>("", "blending_coefficient");
  blending_coeff->resize((*sol_dofs).nb_nodes());

#endif

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);

  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residual);
#if TEST_CASE_IS_TRANSONIC
  scheme->set_artificial_viscosity(art_visc);
// scheme->set_blending_coeff(blending_coeff);
#endif

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh2D);
  ic.set_domain(_2D, "InnerCells");
  math::DenseDVec<Real> ic_state(phys_model::NEQ);

  const Real gamma = 1.4;

#if TEST_CASE_IS_TRANSONIC
  const Real alpha_in = 1.25;
  const Real M_in     = 0.8;
#else
  const Real alpha_in = 2.0;
  const Real M_in     = 0.5;
#endif

  const Real rho_in  = 1.22503;   // 1.0;
  const Real rhou_in = 208.30559; // 1.0;
  const Real rhov_in = rhou_in * std::tan(alpha_in * math::pi / 180);

  const Real v2_in = (rhou_in * rhou_in + rhov_in * rhov_in) / (rho_in * rho_in);

  const Real p_in = rho_in / gamma * v2_in / (M_in * M_in);
  const Real e_in = p_in / (gamma - 1.) + 0.5 * rho_in * v2_in;

  ic_state[0] = rho_in;
  ic_state[1] = rhou_in;
  ic_state[2] = rhov_in;
  ic_state[3] = e_in;

  /*
  ic_state[0] = 1.22503;
  ic_state[1] = 208.30559;
  ic_state[2] = 7.27419;
  ic_state[3] = 271044.38;
  */

  // std::cout << "MA = " << std::sqrt(gamma * p_in / rho_in) /
  // std::sqrt(v2_in)
  // << std::endl; std::cout << "IC STATE = " << ic_state << std::endl;

  ic.apply_values("sol_dofs", ic_state,
                  *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  scheme->add_boundary_condition("WeakWall", "wall_bc", "Wall");

  // std::shared_ptr<bc_base_type> farfield_bc =
  auto farfield_bc = scheme->add_boundary_condition("WeakFarfield", "farfield_bc", "Farfield");
  farfield_bc->set_reference_state(ic_state);
  farfield_bc->print_bc_parameters();

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  clock_t start, end;
  Real elapsed;

  const Real CFLMax = 30.0;

  std::cout.precision(15);

  start = clock();

  for (Uint iter = 0; iter <= 70; ++iter)
  {
    const Real CFL = iter < 10 ? 2.0 : std::min(CFLMax, 2.0 + std::pow(1.50, iter - 10));

    scheme->set_cfl(CFL);

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

#if TEST_CASE_IS_TRANSONIC
    if ((iter > 0) && ((iter % 5) == 0))
    {
      /*
      art_visc_operator.compute_indicator(*mesh2D, *geo_dofs, *sol_dofs,
      *solution, *indicator);
      */
      if (iter < 140)
      {
        art_visc_operator.compute_artificial_viscosity_nodal(*mesh2D, *geo_dofs, *sol_dofs,
                                                             *solution, *art_visc);
      }

      rd_blending_coeff_operator.calculate(*mesh2D, *geo_dofs, *sol_dofs, *solution,
                                           *blending_coeff);
    }
#endif

    scheme->compute_residual_norm(res_L2_norm);
    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);
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

  Mach_number->resize((*sol_dofs).nb_nodes());

  physics::Euler2DCons physical_model;
  physics::Euler2DProperties properties;
  physics::Euler2DProperties::SolGradM gradient_matrix;

  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution, [](const phys_model::Properties &props) { return props.Ma; },
      *Mach_number);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *Mach_number, "Ma");
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *residual, "Res");

#if TEST_CASE_IS_TRANSONIC
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *art_visc, "art_visc");
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *blending_coeff, "theta");
#endif

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

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, entropy_dev, "entropy_dev");

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
  std::cout << "Entropy deviation L2 norm = " << entropy_norm << std::endl;

  // --------------------------------------------------------------------------
  // WRITE DATA ON THE WALL
  // --------------------------------------------------------------------------

  std::cout << "Writing 'naca2d_wall_data.dat' - data on profile boundary ..." << std::endl;
  interpolation::VectorMeshFunction<Real> wall_coords("", "wall_coords");
  interpolation::VectorMeshFunction<Real> wall_solution_values("", "wall_solution");

  solver::extract_boundary_data(*mesh2D, *sol_dofs, *scheme->vec_function(SolverVecFn::solution),
                                "Wall", wall_coords, wall_solution_values);

  std::ofstream outstream;
  outstream.open("naca2d_wall_data.dat");
  outstream << "#    X coord          Y coord               rho              rhou  "
               "  "
               "         "
            << "rhov                    E                     Cp                "
            << "Ma           entropy" << std::endl;
  outstream.precision(15);
  outstream.setf(std::ios::fixed);

  for (Uint i = 0; i < wall_coords.nb_entries(); ++i)
  {
    typedef interpolation::VectorMeshFunction<Real>::const_entry_type entry_type;

    const entry_type coords = wall_coords.const_value(i);
    const entry_type values = wall_solution_values.const_value(i);

    physical_model.compute_properties(coords, values, gradient_matrix, properties);
    const Real cp      = (properties.P - p_in) / (0.5 * rho_in * v2_in);
    const Real Ma      = properties.Ma;
    const Real entropy = (properties.P / p_in) / (std::pow(values[0] / rho_in, 1.4)) - 1.;

    outstream << std::setw(12) << coords[X0] << " " << std::setw(12) << coords[X1] << " "
              << std::setw(15) << values[0] << " " << std::setw(15) << values[1] << " "
              << std::setw(15) << values[2] << " " << std::setw(15) << values[3] << " "
              << std::setw(15) << cp << " " << std::setw(15) << Ma << " " << std::setw(15)
              << entropy << std::endl;
  }
  outstream.close();

  mpi_env.finalize();

  return 0;
}
