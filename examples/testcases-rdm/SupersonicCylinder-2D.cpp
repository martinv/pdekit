#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/MeshManipulator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/art_visc//ArtificialViscosity.hpp"
#include "solver/rdm/DRDMethodExplicit.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/bc/WeakFarfield.hpp"
#include "solver/rdm/bc/WeakSuperInlet.hpp"
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

void evaluate_inlet_state(const Real Ma_in, math::DenseDVec<Real> &inlet_state)
{
  // Material constants
  const Real gamma = 1.4;

  // Inlet values
  const Real rho_in = 1.0;
  const Real u_in   = 1.0;
  const Real v_in   = 0.0;

  const Real p_in = (u_in * u_in + v_in * v_in) / (Ma_in * Ma_in) * rho_in / gamma;
  // std::cout << "Pressure at inlet = " << p_in << std::endl;

  inlet_state.resize(phys_model::NEQ);

  inlet_state[0] = rho_in;
  inlet_state[1] = rho_in * u_in;
  inlet_state[2] = rho_in * v_in;
  inlet_state[3] = p_in / (gamma - 1) + 0.5 * rho_in * (u_in * u_in + v_in * v_in);
}

// ----------------------------------------------------------------------------

void print_res_norm(const math::DenseDVec<Real> &norm)
{
  std::cout.precision(15);
  std::cout << "res =";
  for (Uint i = 0; i < norm.size(); ++i)
  {
    std::cout << " " << std::setw(20) << norm[i];
  }

  std::cout << " , log(res) =";
  for (Uint i = 0; i < norm.size(); ++i)
  {
    std::cout << " " << std::setw(20) << std::log(norm[i]);
  }
  std::cout << std::endl;
}

// ----------------------------------------------------------------------------

class CylinderInlet
{
  public:
  static void set_Ma_inlet(const Real Ma)
  {
    M_in = Ma;
  }

  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component)
  {
    const Real gamma = 1.4;

    // Inlet values
    const Real rho_in = 1.0;
    const Real u_in   = 1.0;
    const Real v_in   = 0.0;
    const Real p_in   = (u_in * u_in + v_in * v_in) / (M_in * M_in) * rho_in / gamma;

    Real result = 0.0;

    switch (component)
    {
      case 0:
        result = rho_in;
        break;
      case 1:
        result = rho_in * u_in;
        break;
      case 2:
        result = rho_in * v_in;
        break;
      case 3:
        const Real e_in = p_in / (gamma - 1) + 0.5 * rho_in * (u_in * u_in + v_in * v_in);
        result          = e_in;
        break;
    };

    return result;
  }

  private:
  static Real M_in;
};

Real CylinderInlet::M_in = 2.0;

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // ------------------------------------------------------
  // DEFINE METHOD TYPE
  // ------------------------------------------------------

  typedef rdm::PGRDMethodExplicit<MeshConfig, phys_model, rdm::PGB> scheme_type;
  // typedef rdm::DRDMethodExplicit<MeshConfig, phys_model, rdm::PGB,
  // rdm::FacetDG> scheme_type;

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme = std::make_shared<scheme_type>();

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  const std::string infilename = "supersonic_cylinder_tri_p1.msh";
  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  // Scale the coordinates - in case cylinder4.CFmesh is used

  // MeshManipulator::scale_mesh(*mesh, 1.e-6);

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  if (scheme->is_continuous())
  {
    MeshType::dof_storage_type::clone_continuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                                 PointSetID::Equidist);
  }
  else
  {
    MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P3,
                                                    PointSetID::Equidist);
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

  // vector_function::ptr cell_lc =
  // fs_solution->create_vector_function("cell_lc");
  // solver::ComputeUpwindLc::rdm_cellwise(*mesh, *cell_lc);

  //  interpolation::FunctionSpace::Function::Ptr f_ptr =
  // fs_solution->function("residual");
  //  std::cout << "Got function " << f_ptr->name() << std::endl;

  // ------------------------------------------------------
  // DISCONTINUITY SENSOR SETUP
  // ------------------------------------------------------

  solver::ArtificialViscosity<MeshConfig> art_visc_operator;
  art_visc_operator.prepare_indicator_operators(*sol_dofs);

  scalar_function::ptr indicator = std::make_shared<scalar_function>("", "smoothness_indicator");
  indicator->resize((*sol_dofs).nb_active_cells());

  scalar_function::ptr art_visc = std::make_shared<scalar_function>("", "artificial_viscosity");
  // art_visc->resize((*sol_dofs).nb_active_cells());
  art_visc->resize((*sol_dofs).nb_nodes());

  // ------------------------------------------------------
  // BLENDING COEFFICIENT FOR NONLINEAR SCHEMES
  // ------------------------------------------------------

  rdm::RDBlendingCoeff<MeshConfig> rd_blending_coeff_operator;

  // Types of viscosity coefficient: "LF_Blend", "ViscosityIndicator",
  // "WongJansen", "Bonanni2D" const std::string blend_coeff_type =
  // "WongJansen";
  const std::string blend_coeff_type = "Bonanni2D";
  rd_blending_coeff_operator.setup(*sol_dofs, blend_coeff_type);

  // This is inlet density - parameter needed to configure the Wong-Jansen
  // discontinuity indicator
  rd_blending_coeff_operator.set_param("u_inf", 1.0);

  scalar_function::ptr blending_coeff =
      std::make_shared<scalar_function>("", "blending_coefficient");
  blending_coeff->resize((*sol_dofs).nb_nodes());
  // Initially, the blending coefficient is set to 1 so that the full
  // dissipation is activated
  blending_coeff->fill(0.0);

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residual);
  // scheme->set_artificial_viscosity(art_visc);
  scheme->set_blending_coeff(blending_coeff);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  // Initial condition (also supersonic inlet condition)

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh2D);
  ic.set_domain(_2D, "InnerCells");
  math::DenseDVec<Real> ic_state(phys_model::NEQ);

  // Initial much number is 2.0
  const Real Ma_init_cond = 2.0;
  evaluate_inlet_state(Ma_init_cond, ic_state);

  ic.apply_values("sol_dofs", ic_state,
                  *solution); // This state is used as both initial
                              // and boundary condition

  CylinderInlet::set_Ma_inlet(Ma_init_cond);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  typedef rdm::RDMBCBase<MeshConfig, phys_model, phys_model::DIM - 1> bc_base_type;

  scheme->add_boundary_condition("WeakWall", "wall_bc", "Wall");

  std::shared_ptr<bc_base_type> inlet_bc =
      scheme->add_boundary_condition("WeakDirichlet", "inlet_bc", "Inlet");
  inlet_bc->set_expression(&CylinderInlet::value);
  inlet_bc->print_bc_parameters();

  /*
  std::shared_ptr<bc_base_type> inlet_bc =
      scheme->add_boundary_condition("WeakSuperInlet", "inlet_bc", "Inlet");
  inlet_bc->set_expression(&CylinderInlet::value);
  */

  std::shared_ptr<bc_base_type> outlet_top_bc =
      scheme->add_boundary_condition("WeakFarfield", "top_bc", "OutletTop");
  outlet_top_bc->set_reference_state(ic_state);
  outlet_top_bc->print_bc_parameters();

  std::shared_ptr<bc_base_type> outlet_bottom_bc =
      scheme->add_boundary_condition("WeakFarfield", "bottom_bc", "OutletBottom");
  outlet_bottom_bc->set_reference_state(ic_state);
  outlet_bottom_bc->print_bc_parameters();

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  // typedef vector_function::entry_type node_value_type;
  // typedef vector_function::const_entry_type const_node_value_type;

  clock_t start, end;
  Real elapsed;

  std::cout.precision(15);
  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  start = clock();
  for (Uint iter = 0; iter < 5000; ++iter)
  {
    /*
    if ((iter > 0) && ((iter % 5) == 0))
    {
      art_visc_operator.compute_artificial_viscosity(geo_mesh->topology().dof_storage(),
                                                     sol_mesh->topology().dof_storage(),
    *solution, *art_visc);
    }
    */

    const Real CFL = 0.5;

    scheme->set_cfl(CFL);
    scheme->assemble_lhs_and_rhs(CFL);
    scheme->solve({});

    // After 2500 iterations, start updating the blending coefficient
    // if ((iter > 2500) && (iter < 4000) &&  (iter % 10 == 0))
    // if ((iter > 2500) && (iter % 10 == 0))
    if (iter % 10 == 0)
    {
      rd_blending_coeff_operator.calculate(*mesh2D, *sol_dofs, *solution, *blending_coeff);
    }

    // At 3000 iterations, increase the Mach number

    if (iter == 3000)
    {
      CylinderInlet::set_Ma_inlet(4.0);

      auto inlet_bc = scheme->boundary_condition("inlet_bc");
      inlet_bc->set_expression(&CylinderInlet::value);

      evaluate_inlet_state(4.0, ic_state);

      auto outlet_top_bc = scheme->boundary_condition("top_bc");
      outlet_top_bc->set_reference_state(ic_state);

      auto outlet_bottom_bc = scheme->boundary_condition("bottom_bc");
      outlet_bottom_bc->set_reference_state(ic_state);
    }

    interpolation::norm_L2(*residual, res_L2_norm);
    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  art_visc_operator.compute_indicator_normalized(*sol_dofs, *solution, *indicator);
  art_visc_operator.compute_artificial_viscosity_nodal(*mesh2D, *sol_dofs, *solution, *art_visc);

  // Compute the blending coefficient and store its values for visualization
  rd_blending_coeff_operator.calculate(*mesh2D, *sol_dofs, *solution, *blending_coeff);

  // --------------------------------------------------------------------------
  // Write the output to file
  // --------------------------------------------------------------------------

  const std::string outfilename = "output_supersonic_cylinder_tri_p1.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "u");
  meshwriter.append_cell_function_to_file(*mesh2D, outfilename, *indicator, "smoothness_indicator");
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *art_visc, "artificial_viscosity");
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *blending_coeff,
                                           "blending_coefficient");

  // --------------------------------------------------------------------------
  // Write Mach number values to file
  // --------------------------------------------------------------------------
  scalar_function::ptr Mach_number = std::make_shared<scalar_function>("", "Ma");
  Mach_number->resize(nb_nodes);

  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution, [](const phys_model::Properties &props) { return props.Ma; },
      *Mach_number);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *Mach_number, "Ma");

  return 0;
}
