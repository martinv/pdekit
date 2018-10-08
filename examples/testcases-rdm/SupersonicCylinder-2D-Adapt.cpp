#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/mesh_function/MeshFunctionSnapshot.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/MeshManipulator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/art_visc/ArtificialViscosity.hpp"
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

template <typename MeshConfig, typename Physics>
class Simulation
{
  public:
  Simulation();

  ~Simulation();

  void read_mesh(const std::string &mesh_name, const PolyOrderID poly_order);

  void setup_artifical_viscosity();

  void setup_blending_coefficient(const std::string &blend_coeff_type);

  void setup_viscosity_indicator();

  void prepare_solver(const bool reset_solution);

  void iterate(const Uint nb_iter);

  void h_adapt();

  void p_adapt();

  void write_output(const std::string &outfilename);

  private:
  typedef Tria<MeshConfig> MeshType;
  typedef typename Tria<MeshConfig>::dof_storage_type DofStorageType;
  typedef typename interpolation::ScalarMeshFunction<Real> scalar_function;
  typedef typename interpolation::VectorMeshFunction<Real> vector_function;

  std::shared_ptr<rdm::RDMethod<MeshConfig, Physics>> build_solver(
      const std::shared_ptr<MeshType const> &mesh, common::PtrHandle<DofStorageType> geo_dofs,
      common::PtrHandle<DofStorageType> sol_dofs, const Uint quad_order,
      typename vector_function::ptr solution, typename vector_function::ptr residuals,
      const std::string &scheme_type, const bool reset_solution, const bool reorder_dofs);

  // Mesh
  typename MeshType::shared_ptr m_mesh;

  // Degrees of freedom for geometry and solution
  common::PtrHandle<DofStorageType> m_geo_dofs;
  common::PtrHandle<DofStorageType> m_sol_dofs;

  // Numerical solver
  std::shared_ptr<rdm::RDMethod<MeshConfig, Physics>> m_solver;

  // Mesh functions for solution and residuals
  typename vector_function::ptr m_solution;
  typename vector_function::ptr m_residual;

  // Operator for artifical viscosity
  solver::ArtificialViscosity<MeshConfig> m_art_visc_operator;
  // Artifical viscosity function
  typename scalar_function::ptr m_art_visc;
  // Indicator based on artifical viscosity
  typename scalar_function::ptr m_indicator;

  // Blending coefficient operator
  rdm::RDBlendingCoeff<MeshConfig> rd_blending_coeff_operator;
  // Blending coefficient function
  typename scalar_function::ptr m_blending_coeff;

  // Types of viscosity coefficient: "LF_Blend", "ViscosityIndicator",
  // "WongJansen"
  std::string m_blend_coeff_type;

  // -------------------
  // Data for adaptation
  // -------------------

  MeshAdaptSequence<Cart2D> adapt_schedule;
  std::vector<CellTransform> cell_adapt_ops;

  std::unique_ptr<interpolation::ScalarMeshFunction<Uint>>
      m_adapt_cell_marker; //("", "adapt_cell_marker");
  interpolation::MeshFunctionSnapshot<Real> mesh_function_snapshot;

  // Iteration counter
  Uint m_iter;

  // Quadrature order
  PolyOrderID m_quadrature_order;
};

template <typename MeshConfig, typename Physics>
Simulation<MeshConfig, Physics>::Simulation() : m_iter(0u), m_quadrature_order(P0)
{
  m_blend_coeff_type = "LF_Blend";

  m_adapt_cell_marker = std::move(std::unique_ptr<interpolation::ScalarMeshFunction<Uint>>(
      new interpolation::ScalarMeshFunction<Uint>("", "adapt_cell_marker")));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
Simulation<MeshConfig, Physics>::~Simulation()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::read_mesh(const std::string &mesh_name,
                                                const PolyOrderID poly_order)
{
  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  m_mesh = typename MeshType::shared_ptr(new MeshType("mesh2D"));

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(mesh_name, *m_mesh, "geo_dofs");

  m_geo_dofs = m_mesh->dof_storage("geo_dofs");
  m_sol_dofs = m_mesh->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_discontinuous(*m_mesh, *m_geo_dofs, *m_sol_dofs, poly_order,
                                                  PointSetID::Equidist);

  /*
  MeshType::dof_storage_type::clone_continuous(m_mesh->topology(),
  *m_geo_dofs, *m_sol_dofs, poly_order, PointSetID::Equidist);
  */

  m_quadrature_order = static_cast<PolyOrderID>(2 * poly_order);

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  m_solution = std::make_shared<vector_function>("", "solution");
  m_residual = std::make_shared<vector_function>("", "residual");

  const Uint nb_nodes = (*m_sol_dofs).nb_nodes();

  m_solution->resize(Physics::NEQ, nb_nodes);
  m_residual->resize(Physics::NEQ, nb_nodes);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::setup_artifical_viscosity()
{
  m_art_visc_operator.prepare_indicator_operators(*m_sol_dofs);

  m_art_visc = std::make_shared<scalar_function>("", "artificial_viscosity");
  m_art_visc->resize((*m_sol_dofs).nb_nodes());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::setup_blending_coefficient(
    const std::string &blend_coeff_type)
{
  m_blend_coeff_type = blend_coeff_type;

  rd_blending_coeff_operator.setup(*m_sol_dofs, m_blend_coeff_type);

  if (m_blend_coeff_type == "WongJansen")
  {
    // This is reference value for Wong-Jansen blending coefficient
    rd_blending_coeff_operator.set_param("u_inf", 1.0);
  }

  m_blending_coeff = std::make_shared<scalar_function>("", "blending_coefficient");
  m_blending_coeff->resize((*m_sol_dofs).nb_nodes());
  m_blending_coeff->fill(1.0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::setup_viscosity_indicator()
{
  m_indicator = std::make_shared<scalar_function>("", "smoothness_indicator");
  m_indicator->resize((*m_geo_dofs).nb_active_cells());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::prepare_solver(const bool reset_solution)
{
  m_solver = build_solver(m_mesh, m_geo_dofs, m_sol_dofs, m_quadrature_order, m_solution,
                          m_residual, "B", reset_solution, true);

  // m_solver->set_artificial_viscosity(m_art_visc);
  m_solver->set_blending_coeff(m_blending_coeff);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::iterate(const Uint nb_iter)
{
  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  // Initial condition (also supersonic inlet condition)

  /* THIS IS NOT NEEDED
  solver::InitialCondition<MeshConfig> ic("initial_cond", m_mesh);
  ic.set_domain(_2D, "InnerCells");
  */
  math::DenseDVec<Real> ic_state(phys_model::NEQ);

  // Initial much number is 2.0
  const Real Ma_init_cond = 2.0;
  evaluate_inlet_state(Ma_init_cond, ic_state);

  for (Uint iter = 0; iter < nb_iter; ++iter)
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

    m_solver->set_cfl(CFL);
    m_solver->assemble_lhs_and_rhs(CFL);
    m_solver->solve({});

    // After 2500 iterations, start updating the blending coefficient
    // if ((iter > 2500) && (iter < 4000) &&  (iter % 10 == 0))
    if ((iter > 2500) && (iter % 10 == 0))
    {
      rd_blending_coeff_operator.calculate(*m_mesh, *m_sol_dofs, *m_solution, *m_blending_coeff);
    }

    /*
    // At 3000 iterations, increase the Mach number
    if (iter == 3000)
    {
      evaluate_inlet_state(4.0, ic_state);

      auto outlet_top_bc = m_solver->boundary_condition("top_bc");
      outlet_top_bc->set_reference_state(ic_state);

      auto outlet_bottom_bc = m_solver->boundary_condition("bottom_bc");
      outlet_bottom_bc->set_reference_state(ic_state);

      CylinderInlet::set_Ma_inlet(4.0);
    }
    */

    interpolation::norm_L2(*m_residual, res_L2_norm);
    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::h_adapt()
{
  const Uint old_nb_active_cells = (*m_sol_dofs).nb_active_cells();

  m_adapt_cell_marker->resize(old_nb_active_cells);
  cell_adapt_ops.resize(old_nb_active_cells);

  m_indicator->resize(old_nb_active_cells);
  m_art_visc_operator.compute_indicator(*m_sol_dofs, *m_solution, *m_indicator, 5.0);

  for (Uint c = 0; c < old_nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > 5.e-5)
    {
      (*m_adapt_cell_marker)[c] = 1;
      cell_adapt_ops[c]         = CellTransform::UNIFORM_REFINE;
    }
    else
    {
      (*m_indicator)[c] = 0;
      cell_adapt_ops[c] = CellTransform::NO_TRANS;
    }
  }

  adapt_schedule.define_h_adapt_ops(*m_mesh, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);

  mesh_function_snapshot.create<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);

  m_mesh->adapt(adapt_schedule);
  (*m_geo_dofs).adapt(adapt_schedule);
  (*m_sol_dofs).adapt(adapt_schedule);

  /*
  tools::MeshInspector<MeshConfig> mesh_inspector;
  const Real volume =
      mesh_inspector.check_mesh_consistency(*m_mesh, "geo_dofs",
  m_quadrature_order); std::cout << "Volume computed by mesh inspector = " <<
  volume << std::endl;
  */

  mesh_function_snapshot.restore_function<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);

  const Uint nb_nodes        = (*m_sol_dofs).nb_nodes();
  const Uint nb_active_cells = (*m_sol_dofs).nb_active_cells();

  // SOLUTION IS NOT RESIZED: this was done by mesh function snapshot above!
  // solution->resize(phys_model::NEQ, nb_nodes);

  m_residual->resize(Physics::NEQ, nb_nodes);
  m_indicator->resize(nb_active_cells);
  m_art_visc->resize(nb_nodes);
  m_blending_coeff->resize(nb_nodes);

  rd_blending_coeff_operator.clear();
  rd_blending_coeff_operator.setup(*m_sol_dofs, m_blend_coeff_type);
  // This is reference value for Wong-Jansen blending coefficient
  rd_blending_coeff_operator.set_param("u_inf", 1.0);

  /*
  m_solver = build_solver(m_mesh, m_geo_dofs, m_sol_dofs, m_quadrature_order,
  m_solution, m_residual, "B", false);
  // m_solver->set_artificial_viscosity(m_art_visc);

  m_solver->set_blending_coeff(m_blending_coeff);
  */

  const bool reset_solution = false;
  prepare_solver(reset_solution);

  const Real CFL = 0.8;

  m_solver->assemble_lhs_and_rhs(CFL);
  m_solver->solve({SolverOption::RecomputePreconditioner});
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::p_adapt()
{
  const Uint nb_active_cells = (*m_sol_dofs).nb_active_cells();

  std::vector<Uint> old_cell_p_order(nb_active_cells);
  std::vector<Uint> adapt_cell_p_order(nb_active_cells);
  std::vector<Real> adapt_quality(nb_active_cells);

  m_indicator->resize(nb_active_cells);
  m_art_visc_operator.compute_indicator(*m_geo_dofs, *m_sol_dofs, *m_solution, *m_indicator, 5.0);

  Uint poly_order = 0;

  for (Uint c = 0; c < nb_active_cells; ++c)
  {
    const mesh::MeshEntity active_cell = (*m_sol_dofs).active_cell(mesh::ActiveIdx(c));
    old_cell_p_order[c]                = active_cell.pt_set_id().poly_order();

    adapt_quality[c] = (*m_indicator)[c];

    if ((*m_indicator)[c] > 5.e-5)
    {
      adapt_cell_p_order[c] = old_cell_p_order[c];
    }
    else
    {
      adapt_cell_p_order[c] = active_cell.pt_set_id().poly_order() + P1;
      poly_order            = std::max(poly_order, adapt_cell_p_order[c]);
    }
  }

  adapt_schedule.define_p_adapt_ops(adapt_cell_p_order);
  // adapt_schedule.define_p_adapt_ops(old_cell_p_order, adapt_cell_p_order,
  // adapt_quality, 0.3);

  // Make sure that the quadrature order to set up the solver is sufficiently
  // high:
  m_quadrature_order = static_cast<PolyOrderID>(poly_order);

  mesh_function_snapshot.create<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);

  // m_mesh->adapt(adapt_schedule);
  // (*m_geo_dofs).adapt(adapt_schedule);
  (*m_sol_dofs).adapt(adapt_schedule);

  mesh_function_snapshot.restore_function<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);

  const Uint nb_nodes = (*m_sol_dofs).nb_nodes();

  m_art_visc_operator.clear();
  m_art_visc_operator.prepare_indicator_operators(*m_geo_dofs, *m_sol_dofs);

  // SOLUTION IS NOT RESIZED: this was done by mesh function snapshot above!
  // solution->resize(phys_model::NEQ, nb_nodes);

  m_residual->resize(Physics::NEQ, nb_nodes);
  m_indicator->resize(nb_active_cells);
  m_art_visc->resize((*m_sol_dofs).nb_nodes());
  m_blending_coeff->resize(nb_nodes);

  rd_blending_coeff_operator.clear();
  rd_blending_coeff_operator.setup(*m_geo_dofs, *m_sol_dofs, m_blend_coeff_type);
  // This is reference value for Wong-Jansen blending coefficient
  rd_blending_coeff_operator.set_param("u_inf", 1.0);

  /*
  m_solver = build_solver(m_mesh, m_geo_dofs, m_sol_dofs, m_quadrature_order,
  m_solution, m_residual, "B", false);
  // m_solver->set_artificial_viscosity(m_art_visc);

  m_solver->set_blending_coeff(m_blending_coeff);
  */

  // ######################################################################################

  // gmsh::GmshWriter writer;
  // writer.write_mesh_to_file(*m_mesh, "sol_dofs", "debug.msh");
  // writer.append_nodal_function_to_file(*m_mesh, "debug.msh", *m_solution,
  // "solution");

  // ######################################################################################

  const bool reset_solution = false;
  prepare_solver(reset_solution);

  const Real CFL = 0.8;

  m_solver->assemble_lhs_and_rhs(CFL);
  m_solver->solve({SolverOption::RecomputePreconditioner});
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::write_output(const std::string &outfilename)
{
  /*
  // Compute the smoothness indicator and artificial viscosity once more to
  store for
  // visualization
  m_art_visc_operator.compute_indicator(*m_geo_dofs, *m_sol_dofs, *m_solution,
  *m_indicator, 7.0);
  m_art_visc_operator.compute_artificial_viscosity(*m_geo_dofs, *m_sol_dofs,
  *m_solution, *m_art_visc, 7.0);
  */

  // Compute the blending coefficient and store its values for visualization
  rd_blending_coeff_operator.calculate(*m_mesh, *m_sol_dofs, *m_solution, *m_blending_coeff);

  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*m_mesh, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *m_solution, "solution");

  /*
  meshwriter.append_cell_function_to_file(*m_mesh, outfilename, *m_indicator,
                                          "smoothness_indicator");
  meshwriter.append_cell_function_to_file(*m_mesh, outfilename, *m_art_visc,
                                          "artificial_viscosity");
                                          */
  meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *m_blending_coeff,
                                           "blending_coefficient");

  // --------------------------------------------------------------------------
  // Write Mach number values to file
  // --------------------------------------------------------------------------
  typename scalar_function::ptr Mach_number = std::make_shared<scalar_function>("", "Ma");
  Mach_number->resize((*m_sol_dofs).nb_nodes());

  solver::PostprocessingUtils<MeshConfig>::template compute_Euler_quantity<Physics>(
      *m_sol_dofs, *m_solution, [](const typename Physics::Properties &props) { return props.Ma; },
      *Mach_number);

  meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *Mach_number, "Ma");

  /*
  meshwriter.append_cell_function_to_file(*m_mesh, outfilename,
  *m_adapt_cell_marker, "adapt_cell_marker");

  Uint p_order_min = 10;
  Uint p_order_max = 0;

  const Uint nb_active_cells = (*m_sol_dofs).nb_active_cells();

  for (Uint c = 0; c < (*m_sol_dofs).nb_active_cells(); ++c)
  {
    const mesh::MeshEntity active_cell =
  (*m_sol_dofs).active_cell(mesh::ActiveIdx(c)); const Uint p_active =
  active_cell.std_region_id().poly_order();

    p_order_min = std::min(p_order_min, p_active);
    p_order_max = std::max(p_order_max, p_active);
  }

  if (p_order_min != p_order_max)
  {
    std::vector<Uint> cell_p_order(nb_active_cells);
    cell_p_order.assign(nb_active_cells, p_order_max);
    adapt_schedule.define_p_adapt_ops(cell_p_order);

    mesh_function_snapshot.create<MeshConfig>(m_mesh->topology(), *m_sol_dofs,
  *m_solution);
    (*m_sol_dofs).adapt(adapt_schedule);
    mesh_function_snapshot.restore_function<MeshConfig>(m_mesh->topology(),
  *m_sol_dofs, *m_solution);

    meshwriter.write_mesh_to_file(*m_mesh, "sol_dofs", "viz_" + outfilename);
    meshwriter.append_nodal_function_to_file(*m_mesh, "viz_" + outfilename,
  *m_solution, "solution");
  }
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
std::shared_ptr<rdm::RDMethod<MeshConfig, Physics>> Simulation<MeshConfig, Physics>::build_solver(
    const std::shared_ptr<MeshType const> &mesh, common::PtrHandle<DofStorageType> geo_dofs,
    common::PtrHandle<DofStorageType> sol_dofs, const Uint quad_order,
    typename vector_function::ptr solution, typename vector_function::ptr residuals,
    const std::string &scheme_type, const bool reset_solution, const bool reorder_dofs)
{
  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  std::shared_ptr<rdm::RDMethod<MeshConfig, Physics>> scheme;

  if (scheme_type == "B")
  {
    scheme =
        std::make_shared<rdm::DRDMethodExplicit<MeshConfig, Physics, rdm::PGB, rdm::FacetDG>>();
    // std::make_shared<rdm::PGRDMethodExplicit<MeshConfig, Physics,
    // rdm::PGB>>();
  }
  else
  {
    scheme =
        std::make_shared<rdm::DRDMethodExplicit<MeshConfig, Physics, rdm::PGLDA, rdm::FacetDG>>();
  }

  /*
  if (reorder_dofs)
  {
    scheme->reorder_dofs(mesh.topology(), *sol_dofs);
  }
  */

  scheme->configure_mesh_data(mesh, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quad_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residuals);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  // Initial condition (also supersonic inlet condition)

  solver::InitialCondition<MeshConfig> ic("initial_cond", m_mesh);
  ic.set_domain(_2D, "InnerCells");
  math::DenseDVec<Real> ic_state(phys_model::NEQ);

  // Initial Mach number is 2.0
  const Real Ma_init_cond = 2.0;
  evaluate_inlet_state(Ma_init_cond, ic_state);

  if (reset_solution)
  {
    //  solution->fill(0.0);
    ic.apply_values("sol_dofs", ic_state,
                    *solution); // This state is used as both initial
                                // and boundary condition
  }

  CylinderInlet::set_Ma_inlet(Ma_init_cond);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  typedef rdm::RDMBCBase<MeshConfig, phys_model, phys_model::DIM - 1> bc_base_type;

  scheme->add_boundary_condition("WeakWall", "wall_bc", "Wall");

  /*
  std::shared_ptr<bc_base_type> inlet_bc =
      scheme->add_boundary_condition("WeakDirichlet", "inlet_bc", "Inlet");
  inlet_bc->set_expression(&CylinderInlet::value);
  inlet_bc->print_bc_parameters();
  */

  std::shared_ptr<bc_base_type> inlet_bc =
      scheme->add_boundary_condition("WeakSuperInlet", "inlet_bc", "Inlet");

  std::shared_ptr<bc_base_type> outlet_top_bc =
      scheme->add_boundary_condition("WeakFarfield", "top_bc", "OutletTop");
  outlet_top_bc->set_reference_state(ic_state);
  outlet_top_bc->print_bc_parameters();

  std::shared_ptr<bc_base_type> outlet_bottom_bc =
      scheme->add_boundary_condition("WeakFarfield", "bottom_bc", "OutletBottom");
  outlet_bottom_bc->set_reference_state(ic_state);
  outlet_bottom_bc->print_bc_parameters();

  return scheme;
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Simulation<Cart2D, physics::Euler2DCons> simulation;
  simulation.read_mesh("supersonic_cylinder_tri_p2.msh", P2);

  simulation.setup_artifical_viscosity();
  simulation.setup_viscosity_indicator();

  // Types of viscosity coefficient: "LF_Blend", "ViscosityIndicator",
  // "WongJansen", "Bonanni2D"
  const std::string blend_coeff_type = "WongJansen";
  simulation.setup_blending_coefficient(blend_coeff_type);

  const bool reset_solution = true;
  simulation.prepare_solver(reset_solution);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;

  start = clock();

  simulation.iterate(1500);

  for (Uint adapt_pass = 0; adapt_pass < 1; ++adapt_pass)
  {
    simulation.h_adapt();
    // simulation.p_adapt();
    simulation.iterate(700);
  }

  /*
  simulation.iterate(600);
  */

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Write the output to file
  simulation.write_output("output_supersonic_cylinder_tri_adapt.msh");

  return 0;
}
