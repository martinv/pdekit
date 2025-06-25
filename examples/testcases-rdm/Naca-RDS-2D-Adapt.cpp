#include <ctime>
#include <iomanip>
#include <iostream>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/OperatorResidualCellwise.hpp"
#include "interpolation/mesh_function/MeshFunctionSnapshot.hpp"
#include "interpolation/mesh_function/MeshFunctionTools.hpp"
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
// typedef physics::Euler2DCons phys_model;

#define TEST_CASE_IS_TRANSONIC 1

// ----------------------------------------------------------------------------

struct SimulationOptions
{
  Uint nb_iter;
  Uint freeze_artif_visc_at_iter;
  Real CFLMin;
  Real CFLMax;
  bool use_artificial_viscosity;
  bool use_blending;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
class Simulation
{
  public:
  Simulation();

  ~Simulation();

  void read_mesh(const std::string &mesh_name, const PolyOrderID poly_order);

  void setup_artifical_viscosity_operator();

  void setup_artifical_viscosity_field();

  void setup_blending_coefficient(const std::string &blend_coeff_type);

  void setup_viscosity_indicator_field();

  void prepare_solver(const bool reset_solution, const bool reorder_dofs,
                      const bool use_artificial_viscosity, const bool use_blending);

  void iterate(const SimulationOptions &options);

  void h_adapt(const h_AdaptStrategy adapt_strategy, const SimulationOptions &options);

  void p_adapt(const Uint p_max, const SimulationOptions &options);

  void hp_adapt(const h_AdaptStrategy adapt_strategy, const SimulationOptions &options);

  void write_output(const std::string &outfilename, const SimulationOptions &options);

  private:
  typedef Tria<MeshConfig> MeshType;
  typedef typename Tria<MeshConfig>::dof_storage_type DofStorageType;
  typedef typename interpolation::ScalarMeshFunction<Real> scalar_function;
  typedef typename interpolation::VectorMeshFunction<Real> vector_function;

  /// METHODS

  std::shared_ptr<rdm::RDMethod<MeshConfig, Physics>> build_solver(
      const std::shared_ptr<MeshType const> &mesh, common::PtrHandle<DofStorageType> geo_dofs,
      common::PtrHandle<DofStorageType> sol_dofs, const Uint quad_order,
      typename vector_function::ptr solution, typename vector_function::ptr residuals,
      const std::string &scheme_type, const bool reset_solution, const bool reorder_dofs);

  void set_initial_condition(const std::shared_ptr<MeshType const> &mesh,
                             vector_function &solution);

  void print_res_norm(const Uint iter, const Real CFL, const math::DenseDVec<Real> &norm) const;

  /// DATA
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

  std::unique_ptr<interpolation::ScalarMeshFunction<Uint>> m_cell_dist;

  // Vector of new node ids after reordering
  std::vector<Int> m_reordering;

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

  m_cell_dist = std::move(std::unique_ptr<interpolation::ScalarMeshFunction<Uint>>(
      new interpolation::ScalarMeshFunction<Uint>("", "cell_distance")));
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
void Simulation<MeshConfig, Physics>::setup_artifical_viscosity_operator()
{
  m_art_visc_operator.prepare_indicator_operators(*m_sol_dofs);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::setup_artifical_viscosity_field()
{
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
void Simulation<MeshConfig, Physics>::setup_viscosity_indicator_field()
{
  m_indicator = std::make_shared<scalar_function>("", "smoothness_indicator");
  m_indicator->resize((*m_geo_dofs).nb_active_cells());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::prepare_solver(const bool reset_solution,
                                                     const bool reorder_dofs,
                                                     const bool use_artificial_viscosity,
                                                     const bool use_blending)
{
  const std::string rdm_scheme_type = use_artificial_viscosity ? "LDA" : "B";
  m_solver = build_solver(m_mesh, m_geo_dofs, m_sol_dofs, m_quadrature_order, m_solution,
                          m_residual, rdm_scheme_type, reset_solution, reorder_dofs);

  if (use_artificial_viscosity)
  {
    m_solver->set_artificial_viscosity(m_art_visc);
  }

  if (use_blending)
  {
    m_solver->set_blending_coeff(m_blending_coeff);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::iterate(const SimulationOptions &options)
{
  std::cout.precision(15);
  math::DenseDVec<Real> res_L2_norm(Physics::NEQ);

  for (Uint iter = 0; iter <= options.nb_iter; ++iter)
  {
    const Real CFL =
        iter < 10 ? options.CFLMin : std::min(options.CFLMax, 2.0 + std::pow(1.50, iter - 10));

    m_solver->set_cfl(CFL);

    if (iter % 20 == 0)
    {
      m_solver->assemble_lhs_and_rhs(CFL);
      m_solver->solve({SolverOption::RecomputePreconditioner});
    }
    else
    {
      m_solver->assemble_rhs();
      m_solver->solve({});
    }

    m_solver->compute_residual_norm(res_L2_norm);

    print_res_norm(iter, CFL, res_L2_norm);

#if TEST_CASE_IS_TRANSONIC
    if ((iter < options.freeze_artif_visc_at_iter) && (iter > 0) && ((iter % 5) == 0))
    {
      /*
      art_visc_operator.compute_indicator(*mesh2D, *geo_dofs, *sol_dofs,
      *solution, *indicator);
      */

      if (options.use_artificial_viscosity)
      {
        m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_sol_dofs, *m_solution,
                                                               *m_art_visc);
      }

      if (options.use_blending)
      {
        rd_blending_coeff_operator.calculate(*m_mesh, *m_sol_dofs, *m_solution, *m_blending_coeff);
      }
    }
#endif
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::h_adapt(const h_AdaptStrategy adapt_strategy,
                                              const SimulationOptions &options)
{
  const Uint old_nb_active_cells = (*m_sol_dofs).nb_active_cells();

  m_adapt_cell_marker->resize(old_nb_active_cells);
  cell_adapt_ops.resize(old_nb_active_cells);

  m_indicator->resize(old_nb_active_cells);
  m_art_visc_operator.compute_indicator(*m_sol_dofs, *m_solution, *m_indicator, 5.0);

  /*
  // Adapt based on artificial viscosity operator
  m_art_visc_operator.compute_indicator(*m_geo_dofs, *m_sol_dofs, *m_solution,
  *m_indicator, 5.0);

  for (Uint c = 0; c < old_nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > 5.e-5)
    {
      (*m_adapt_cell_marker)[c] = 1;
      cell_adapt_ops[c] = CellTransform::UNIFORM_REFINE;
    }
    else
    {
      (*m_indicator)[c] = 0;
      cell_adapt_ops[c] = CellTransform::DO_NOTHING;
    }
  }
  */

  // Adapt based on blending operator
  if (options.use_blending)
  {
    rd_blending_coeff_operator.calculate_in_cells(*m_mesh, *m_sol_dofs, *m_solution, *m_indicator);
  }

#define USE_CELL_DISTANCE 0
#if USE_CELL_DISTANCE

  const Uint max_cell_dist = 2;

  std::vector<mesh::ActiveIdx> level_zero_cells;

  for (Uint c = 0; c < old_nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > 5.e-2)
    {
      level_zero_cells.push_back(mesh::ActiveIdx(c));
    }
  }

  interpolation::ScalarMeshFunction<Uint> cell_dist("", "cell_distance");
  const common::ArrayView<const mesh::ActiveIdx, Uint> level_zero_view(level_zero_cells.data(),
                                                                       level_zero_cells.size());
  interpolation::MeshFunctionTools::compute_cell_distance(*m_mesh, level_zero_view, cell_dist);

  for (Uint c = 0; c < cell_dist.nb_entries(); ++c)
  {
    if (cell_dist[c] < max_cell_dist)
    {
      (*m_adapt_cell_marker)[c] = 1;
      cell_adapt_ops[c]         = CellTransform::UNIFORM_REFINE;
    }
    else
    {
      (*m_indicator)[c] = 0;
      cell_adapt_ops[c] = CellTransform::DO_NOTHING;
    }
  }

#else

  for (Uint c = 0; c < old_nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > 5.e-2)
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
#endif

  // ######################################################################################

  if (options.use_blending)
  {
    rd_blending_coeff_operator.calculate(*m_mesh, *m_sol_dofs, *m_solution, *m_blending_coeff);
  }

  if (options.use_artificial_viscosity)
  {
    m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_sol_dofs, *m_solution,
                                                           *m_art_visc);
  }

  write_output("debug_before_h_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_before_h_adapt.msh");

  // ######################################################################################

  adapt_schedule.define_h_adapt_ops(*m_mesh, cell_adapt_ops, adapt_strategy);

  if (adapt_strategy == h_AdaptStrategy::w_hanging_nodes)
  {
    mesh_function_snapshot.create<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);
  }

  m_mesh->adapt(adapt_schedule);
  (*m_geo_dofs).adapt_update(*m_mesh);
  (*m_sol_dofs).adapt_update(*m_mesh);

  /*
  tools::MeshInspector<MeshConfig> mesh_inspector;
  const Real volume =
      mesh_inspector.check_mesh_consistency(*m_mesh, "geo_dofs",
  m_quadrature_order); std::cout << "Volume computed by mesh inspector = " <<
  volume << std::endl;
  */

  if (adapt_strategy == h_AdaptStrategy::w_hanging_nodes)
  {
    mesh_function_snapshot.restore_function<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);
  }
  else if (adapt_strategy == h_AdaptStrategy::red_green)
  {
    const Uint sol_dofs_nb_nodes = (*m_sol_dofs).nb_nodes();
    m_solution->resize(Physics::NEQ, sol_dofs_nb_nodes);
    set_initial_condition(m_mesh, *m_solution);
  }

  const Uint nb_nodes        = (*m_sol_dofs).nb_nodes();
  const Uint nb_active_cells = (*m_sol_dofs).nb_active_cells();

  // SOLUTION IS NOT RESIZED: this was done by mesh function snapshot above!
  // solution->resize(phys_model::NEQ, nb_nodes);

  m_residual->resize(Physics::NEQ, nb_nodes);
  m_indicator->resize(nb_active_cells);

  if (options.use_blending)
  {
    m_blending_coeff->resize(nb_nodes);
    rd_blending_coeff_operator.clear();
    rd_blending_coeff_operator.setup(*m_sol_dofs, m_blend_coeff_type);
    // This is reference value for Wong-Jansen blending coefficient
    rd_blending_coeff_operator.set_param("u_inf", 1.0);

    rd_blending_coeff_operator.calculate(*m_mesh, *m_sol_dofs, *m_solution, *m_blending_coeff);
  }

  if (options.use_artificial_viscosity)
  {
    m_art_visc->resize(nb_nodes);
    m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_sol_dofs, *m_solution,
                                                           *m_art_visc);
  }

  write_output("debug_after_h_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_after_h_adapt.msh");

  // ######################################################################################

  const bool reset_solution = false;
  const bool reorder_dofs   = true;
  prepare_solver(reset_solution, reorder_dofs, options.use_artificial_viscosity,
                 options.use_blending);

  if (reorder_dofs)
  {
    if (options.use_blending)
    {
      m_blending_coeff->apply_reordering(m_reordering);
    }
    if (options.use_artificial_viscosity)
    {
      m_art_visc->apply_reordering(m_reordering);
    }
  }

  const Real CFL = 0.8;

  m_solver->assemble_lhs_and_rhs(CFL);
  m_solver->solve({SolverOption::RecomputePreconditioner});
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::p_adapt(const Uint p_max, const SimulationOptions &options)
{
  // --------------------------------------------

  const Uint nb_active_cells = (*m_sol_dofs).nb_active_cells();
  m_indicator->resize(nb_active_cells);

  m_art_visc_operator.clear();
  m_art_visc_operator.prepare_indicator_operators(*m_geo_dofs, *m_sol_dofs);

  if (options.use_blending)
  {
    rd_blending_coeff_operator.calculate_in_cells(*m_mesh, *m_geo_dofs, *m_sol_dofs, *m_solution,
                                                  *m_indicator);
  }

  if (options.use_artificial_viscosity)
  {
    m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_geo_dofs, *m_sol_dofs,
                                                           *m_solution, *m_art_visc);
    m_art_visc_operator.compute_indicator(*m_geo_dofs, *m_sol_dofs, *m_solution, *m_indicator, 5.0);
  }

  // --------------------------------------------
  write_output("debug_before_p_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_before_p_adapt.msh");
  // --------------------------------------------

  std::vector<Uint> adapt_cell_p_order(nb_active_cells);

  Uint poly_order             = 0;
  const Real refine_threshold = (options.use_artificial_viscosity) ? 1.e-3 : 5.e-2;

  for (Uint c = 0; c < nb_active_cells; ++c)
  {
    const mesh::MeshEntity active_cell = (*m_sol_dofs).active_cell(mesh::ActiveIdx(c));
    const Uint old_cell_p_order        = active_cell.pt_set_id().poly_order();

    if ((*m_indicator)[c] > refine_threshold)
    {
      adapt_cell_p_order[c] = old_cell_p_order;
    }
    else
    {
      adapt_cell_p_order[c] = old_cell_p_order + P1;
      poly_order            = std::max(poly_order, adapt_cell_p_order[c]);
    }
  }

  // ---------

  /*
  std::vector<mesh::ActiveIdx> level_zero_cells;

  for (Uint c = 0; c < nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > 1.e-3)
    {
      level_zero_cells.push_back(mesh::ActiveIdx(c));
    }
  }

  const common::ArrayView<const mesh::ActiveIdx, _1D, Uint>
  level_zero_view(level_zero_cells.data(), level_zero_cells.size());
  interpolation::MeshFunctionTools::compute_cell_distance(*m_mesh,
  level_zero_view, *m_cell_dist);

  for (Uint c = 0; c < m_cell_dist->nb_entries(); ++c)
  {
    adapt_cell_p_order[c] = std::min(p_max, (*m_cell_dist)[c] + P1);
    poly_order = std::max(poly_order, adapt_cell_p_order[c]);
  }
  */

  // ---------

  adapt_schedule.define_p_adapt_ops(adapt_cell_p_order);
  // adapt_schedule.define_p_adapt_ops(old_cell_p_order, adapt_cell_p_order,
  // adapt_quality, 0.3);

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

  // Make sure that the quadrature order to set up the solver is sufficiently
  // high:
  m_quadrature_order = static_cast<PolyOrderID>(2 * poly_order);

  if (options.use_blending)
  {
    m_blending_coeff->resize(nb_nodes);
    rd_blending_coeff_operator.clear();
    rd_blending_coeff_operator.setup(*m_geo_dofs, *m_sol_dofs, m_blend_coeff_type);
    // This is reference value for Wong-Jansen blending coefficient
    rd_blending_coeff_operator.set_param("u_inf", 1.0);

    /*
    rd_blending_coeff_operator.calculate(*m_mesh, *m_geo_dofs, *m_sol_dofs,
    *m_solution, *m_blending_coeff);
    */
    m_blending_coeff->fill(0.0);
  }

  if (options.use_artificial_viscosity)
  {
    m_art_visc->resize((*m_sol_dofs).nb_nodes());
    m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_geo_dofs, *m_sol_dofs,
                                                           *m_solution, *m_art_visc);
  }

  const bool reset_solution = false;
  const bool reorder_dofs   = true;
  prepare_solver(reset_solution, reorder_dofs, options.use_artificial_viscosity,
                 options.use_blending);

  if (reorder_dofs)
  {
    if (options.use_blending)
    {
      m_blending_coeff->apply_reordering(m_reordering);
    }
    if (options.use_artificial_viscosity)
    {
      m_art_visc->apply_reordering(m_reordering);
    }
  }

  // --------------------------------------------
  write_output("debug_after_p_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_after_p_adapt.msh");
  // --------------------------------------------

  m_solver->assemble_lhs_and_rhs(options.CFLMin);
  m_solver->solve({SolverOption::RecomputePreconditioner});
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::hp_adapt(const h_AdaptStrategy adapt_strategy,
                                               const SimulationOptions &options)
{
  const Uint old_nb_active_cells = (*m_sol_dofs).nb_active_cells();

  interpolation::ScalarMeshFunction<Uint> cell_p_order("", "cell_p_order");
  cell_p_order.resize(old_nb_active_cells);

  m_adapt_cell_marker->resize(old_nb_active_cells);
  cell_adapt_ops.resize(old_nb_active_cells);

  m_indicator->resize(old_nb_active_cells);

  // If we're running with Blended scheme, the artificial viscosity operator
  // is NOT automatically updated between adaptation loops to handle new
  // element types added to the mesh - we have to do it
  // here manually
  if (options.use_blending)
  {
    m_art_visc_operator.clear();
    m_art_visc_operator.prepare_indicator_operators(*m_geo_dofs, *m_sol_dofs);
  }

  m_art_visc_operator.compute_indicator(*m_geo_dofs, *m_sol_dofs, *m_solution, *m_indicator, 5.0);

  /*
  // Adapt based on artificial viscosity operator
  m_art_visc_operator.compute_indicator(*m_geo_dofs, *m_sol_dofs, *m_solution,
  *m_indicator, 5.0);

  for (Uint c = 0; c < old_nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > 5.e-5)
    {
      (*m_adapt_cell_marker)[c] = 1;
      cell_adapt_ops[c] = CellTransform::UNIFORM_REFINE;
    }
    else
    {
      (*m_indicator)[c] = 0;
      cell_adapt_ops[c] = CellTransform::DO_NOTHING;
    }
  }
  */

  // Adapt based on blending operator
  if (options.use_blending)
  {
    rd_blending_coeff_operator.calculate_in_cells(*m_mesh, *m_geo_dofs, *m_sol_dofs, *m_solution,
                                                  *m_indicator);
  }

  Uint poly_order = 0;

#define USE_CELL_DISTANCE 0
#if USE_CELL_DISTANCE

  const Uint max_cell_dist = 2;

  std::vector<mesh::ActiveIdx> level_zero_cells;

  for (Uint c = 0; c < old_nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > 5.e-2)
    {
      level_zero_cells.push_back(mesh::ActiveIdx(c));
    }
  }

  interpolation::ScalarMeshFunction<Uint> cell_dist("", "cell_distance");
  const common::ArrayView<const mesh::ActiveIdx, Uint> level_zero_view(level_zero_cells.data(),
                                                                       level_zero_cells.size());
  interpolation::MeshFunctionTools::compute_cell_distance(*m_mesh, level_zero_view, cell_dist);

  for (Uint c = 0; c < cell_dist.nb_entries(); ++c)
  {
    if (cell_dist[c] < max_cell_dist)
    {
      (*m_adapt_cell_marker)[c] = 1;
      cell_adapt_ops[c]         = CellTransform::UNIFORM_REFINE;
    }
    else
    {
      (*m_indicator)[c] = 0;
      cell_adapt_ops[c] = CellTransform::DO_NOTHING;
    }
  }

#else

  std::vector<Real> indicator_vector;
  indicator_vector.resize(m_indicator->nb_entries());

  // const Real refine_threshold = (options.use_artificial_viscosity) ? 5.e-2
  // : 0.1;
  const Real refine_threshold = (options.use_artificial_viscosity) ? 1.e-3 : 5.e-2;

  for (Uint c = 0; c < old_nb_active_cells; ++c)
  {
    const mesh::MeshEntity active_cell = (*m_sol_dofs).active_cell(mesh::ActiveIdx(c));
    cell_p_order[c]                    = active_cell.pt_set_id().poly_order();

    if ((*m_indicator)[c] > refine_threshold)
    {
      (*m_adapt_cell_marker)[c] = 1;
      cell_adapt_ops[c]         = CellTransform::UNIFORM_REFINE;
    }
    else
    {
      (*m_indicator)[c] = 0;
      cell_adapt_ops[c] = CellTransform::NO_TRANS;

      cell_p_order[c] = cell_p_order[c] + P1;
      poly_order      = std::max(poly_order, cell_p_order[c]);
    }
    indicator_vector[c] = (*m_indicator)[c];
  }

  m_quadrature_order  = static_cast<PolyOrderID>(2 * poly_order);
#endif

  // ######################################################################################

  if (options.use_blending)
  {
    rd_blending_coeff_operator.calculate(*m_mesh, *m_geo_dofs, *m_sol_dofs, *m_solution,
                                         *m_blending_coeff);
  }

  if (options.use_artificial_viscosity)
  {
    m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_geo_dofs, *m_sol_dofs,
                                                           *m_solution, *m_art_visc);
  }

  write_output("debug_before_hp_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_before_hp_adapt.msh");

  // ######################################################################################

  adapt_schedule.define_h_adapt_ops(*m_mesh, cell_adapt_ops, adapt_strategy);

  mesh_function_snapshot.create_cellwise<MeshConfig>(*m_mesh, *m_sol_dofs, cell_p_order);

  if (adapt_strategy == h_AdaptStrategy::w_hanging_nodes)
  {
    mesh_function_snapshot.create<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);
  }

  m_mesh->adapt(adapt_schedule);

  mesh_function_snapshot.restore_function_cellwise<MeshConfig>(*m_mesh, *m_sol_dofs, cell_p_order);

  const common::ArrayView<const Uint, _1D, Uint> cell_p_order_view(&cell_p_order[0],
                                                                   cell_p_order.nb_entries());

  (*m_geo_dofs).adapt_update(*m_mesh, cell_p_order_view);
  (*m_sol_dofs).adapt_update(*m_mesh, cell_p_order_view);

  /*
  tools::MeshInspector<MeshConfig> mesh_inspector;
  const Real volume =
      mesh_inspector.check_mesh_consistency(*m_mesh, "geo_dofs",
  m_quadrature_order); std::cout << "Volume computed by mesh inspector = " <<
  volume << std::endl;
  */

  if (adapt_strategy == h_AdaptStrategy::w_hanging_nodes)
  {
    mesh_function_snapshot.restore_function<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);
  }
  else if (adapt_strategy == h_AdaptStrategy::red_green)
  {
    const Uint sol_dofs_nb_nodes = (*m_sol_dofs).nb_nodes();
    m_solution->resize(Physics::NEQ, sol_dofs_nb_nodes);
    set_initial_condition(m_mesh, *m_solution);
  }

  const Uint nb_nodes        = (*m_sol_dofs).nb_nodes();
  const Uint nb_active_cells = (*m_sol_dofs).nb_active_cells();

  // SOLUTION IS NOT RESIZED: this was done by mesh function snapshot above!
  // solution->resize(phys_model::NEQ, nb_nodes);

  m_residual->resize(Physics::NEQ, nb_nodes);
  m_indicator->resize(nb_active_cells);

  if (options.use_blending)
  {
    m_blending_coeff->resize(nb_nodes);
    rd_blending_coeff_operator.clear();
    rd_blending_coeff_operator.setup(*m_geo_dofs, *m_sol_dofs, m_blend_coeff_type);
    // This is reference value for Wong-Jansen blending coefficient
    rd_blending_coeff_operator.set_param("u_inf", 1.0);

    rd_blending_coeff_operator.calculate(*m_mesh, *m_geo_dofs, *m_sol_dofs, *m_solution,
                                         *m_blending_coeff);
  }

  if (options.use_artificial_viscosity)
  {
    m_art_visc->resize(nb_nodes);
    m_art_visc_operator.clear();
    m_art_visc_operator.prepare_indicator_operators(*m_geo_dofs, *m_sol_dofs);
    m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_geo_dofs, *m_sol_dofs,
                                                           *m_solution, *m_art_visc);
  }

  write_output("debug_after_hp_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_after_hp_adapt.msh");

  // ######################################################################################

  const bool reset_solution = false;
  const bool reorder_dofs   = true;
  prepare_solver(reset_solution, reorder_dofs, options.use_artificial_viscosity,
                 options.use_blending);

  if (reorder_dofs)
  {
    if (options.use_blending)
    {
      m_blending_coeff->apply_reordering(m_reordering);
    }
    if (options.use_artificial_viscosity)
    {
      m_art_visc->apply_reordering(m_reordering);
    }
  }

  m_solver->assemble_lhs_and_rhs(options.CFLMin);
  m_solver->solve({SolverOption::RecomputePreconditioner});
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::write_output(const std::string &outfilename,
                                                   const SimulationOptions &options)
{
  // ------------------------------------------------------
  // FREESTREAM STATE NEEDED FOR WRITING THE OUTPUT TO FILE
  // ------------------------------------------------------

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

  // ------------------------------------------------------
  // WRITE THE OUTPUT TO FILE
  // ------------------------------------------------------

  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*m_mesh, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *m_solution, "u");

  // Write Mach number values to file
  scalar_function::ptr Mach_number = std::make_shared<scalar_function>("", "Ma");

  Mach_number->resize((*m_sol_dofs).nb_nodes());

  using post_utils = solver::PostprocessingUtils<MeshConfig>;

  post_utils::template compute_Euler_quantity<Physics>(
      *m_sol_dofs, *m_solution, [](const typename Physics::Properties &props) { return props.Ma; },
      *Mach_number);

  meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *Mach_number, "Ma");
  meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *m_residual, "Res");

#if TEST_CASE_IS_TRANSONIC
  if (options.use_artificial_viscosity)
  {
    meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *m_art_visc, "art_visc");
  }

  if (options.use_blending)
  {
    meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *m_blending_coeff, "theta");
  }
#endif

  std::vector<Uint> p_order(m_mesh->nb_active_cells());
  for (Uint ac = 0; ac < m_mesh->nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity active_cell = (*m_sol_dofs).active_cell(mesh::ActiveIdx(ac));
    p_order[ac]                        = active_cell.pt_set_id().poly_order();
  }

  const common::ArrayView<const Uint, _1D, Uint> poly_orders_view(p_order.data(), p_order.size());
  meshwriter.append_cell_function_to_file(*m_mesh, outfilename, poly_orders_view, "poly_order");

  /*
  interpolation::VectorMeshFunction<Real> op_residual_function("",
  "operator_residual"); interpolation::OperatorResidualCellwise<MeshConfig,
  Physics> operator_res; operator_res.setup(*m_geo_dofs, *m_sol_dofs);
  operator_res.evaluate(*m_geo_dofs, *m_sol_dofs, *m_solution,
  op_residual_function); meshwriter.append_cell_function_to_file(*m_mesh,
  outfilename, op_residual_function, "residuals");
  */

  // --------------------------------------------------------------------------

  std::cout << "Evaluating entropy ... " << std::endl;
  interpolation::ScalarMeshFunction<Real> entropy_dev("", "entropy_dev");

  entropy_dev.resize((*m_sol_dofs).nb_nodes());

  post_utils::template compute_Euler_quantity<Physics>(
      *m_sol_dofs, *m_solution,
      [rho_in, p_in](const typename Physics::Properties &props) {
        return (props.P / p_in) / (std::pow(props.rho / rho_in, 1.4)) - 1.;
      },
      entropy_dev);

  meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, entropy_dev, "entropy_dev");

  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto quad_generator = [this](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, this->m_quadrature_order, PointSetID::Gauss);
  };

  interpolation::FunctionSpace<MeshConfig> geo_cell_space;
  geo_cell_space.set_reference_fe_values((*m_geo_dofs).as_range(), sf_generator, quad_generator);
  interpolation::FunctionSpace<MeshConfig> sol_cell_space;
  sol_cell_space.set_reference_fe_values((*m_sol_dofs).as_range(), sf_generator, quad_generator);

  const Real entropy_norm = PostprocessingUtils<MeshConfig>::compute_function_norm(
      *m_geo_dofs, *m_sol_dofs, geo_cell_space, sol_cell_space, entropy_dev);

  std::cout.precision(15);
  std::cout.setf(std::ios::fixed);
  std::cout << "Entropy deviation L2 norm = " << entropy_norm << std::endl;

  // --------------------------------------------------------------------------
  // WRITE DATA ON THE WALL
  // --------------------------------------------------------------------------

  Physics physical_model;
  typename Physics::Properties properties;
  typename Physics::Properties::SolGradM gradient_matrix;

  std::cout << "Writing 'naca2d_adapt_wall_data.dat' - data on profile boundary ..." << std::endl;
  interpolation::VectorMeshFunction<Real> wall_coords("", "wall_coords");
  interpolation::VectorMeshFunction<Real> wall_solution_values("", "wall_solution");

  solver::extract_boundary_data(*m_mesh, *m_sol_dofs, *m_solution, "Wall", wall_coords,
                                wall_solution_values);

  std::ofstream outstream;
  outstream.open("naca2d_adapt_wall_data.dat");
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

  if (scheme_type == "LDA")
  {
    scheme =
        std::make_shared<rdm::DRDMethodImplicit<MeshConfig, Physics, rdm::PGLDA, rdm::FacetDG>>();
  }
  else if (scheme_type == "B")
  {
    scheme =
        std::make_shared<rdm::DRDMethodImplicit<MeshConfig, Physics, rdm::PGB, rdm::FacetDG>>();
  }
  else
  {
    // Set the pointer to zero and let the simulation crash ...
    scheme.reset();
  }

  if (reorder_dofs)
  {
    scheme->compute_node_reordering(*mesh, *sol_dofs, m_reordering);
    (*sol_dofs).renumber_dofs(m_reordering);

    solution->apply_reordering(m_reordering);
    residuals->apply_reordering(m_reordering);
  }

  scheme->configure_mesh_data(mesh, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quad_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residuals);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------
  solver::InitialCondition<MeshConfig> ic("initial_cond", m_mesh);
  ic.set_domain(_2D, "InnerCells");
  math::DenseDVec<Real> ic_state(Physics::NEQ);

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

  if (reset_solution)
  {
    ic.apply_values("sol_dofs", ic_state,
                    *m_solution); // The inlet function is also used as bc
  }

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  scheme->add_boundary_condition("WeakWall", "wall_bc", "Wall");

  // std::shared_ptr<bc_base_type> farfield_bc =
  auto farfield_bc = scheme->add_boundary_condition("WeakFarfield", "farfield_bc", "Farfield");
  farfield_bc->set_reference_state(ic_state);
  farfield_bc->print_bc_parameters();

  return scheme;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::set_initial_condition(
    const std::shared_ptr<MeshType const> &mesh, vector_function &solution)
{
  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------
  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh);
  ic.set_domain(_2D, "InnerCells");
  math::DenseDVec<Real> ic_state(Physics::NEQ);

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

  ic.apply_values("sol_dofs", ic_state, solution);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::print_res_norm(const Uint iter, const Real CFL,
                                                     const math::DenseDVec<Real> &norm) const
{
  solver::SolverIO::print_iter_and_res_norm(iter, CFL, norm);
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // ------------------------------------------------------
  // INITIALIZE ENVIRONMENT
  // ------------------------------------------------------

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // SET INITIAL OPTION VALUES
  // ------------------------------------------------------

  SimulationOptions options;
  options.nb_iter                   = 250;
  options.freeze_artif_visc_at_iter = 10000;
  options.CFLMin                    = 2.0;
  options.CFLMax                    = 20.0;
  options.use_artificial_viscosity  = true;
  options.use_blending              = false;

  const bool reset_solution = true;
  const bool reorder_dofs   = true;

  // ------------------------------------------------------
  // READ MESH AND PREPARE SOLVER
  // ------------------------------------------------------

  Simulation<Cart2D, physics::Euler2DCons> simulation;
  simulation.read_mesh("naca0012_p2_tri.msh", P2);

#if TEST_CASE_IS_TRANSONIC

  simulation.setup_artifical_viscosity_operator();

  if (options.use_artificial_viscosity)
  {
    simulation.setup_artifical_viscosity_field();
    simulation.setup_viscosity_indicator_field();
  }
  else
  {
    // Types of viscosity coefficient: "LF_Blend", "ViscosityIndicator",
    // "WongJansen", "Bonanni2D"
    const std::string blend_coeff_type = "Bonanni2D";
    simulation.setup_blending_coefficient(blend_coeff_type);
    simulation.setup_viscosity_indicator_field();
  }

#endif

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;

  simulation.prepare_solver(reset_solution, reorder_dofs, options.use_artificial_viscosity,
                            options.use_blending);

  start = clock();

  simulation.iterate(options);

  options.nb_iter                   = 100;
  options.freeze_artif_visc_at_iter = 3000;

  options.CFLMin = 0.7;
  options.CFLMax = 0.8;

  for (Uint adapt_pass = 0; adapt_pass < 1; ++adapt_pass)
  {
    simulation.h_adapt(h_AdaptStrategy::w_hanging_nodes, options);
    // simulation.p_adapt(P3, options);
    // simulation.hp_adapt(h_AdaptStrategy::w_hanging_nodes, options);
    simulation.iterate(options);
  }

  // simulation.iterate(6000, 20.0);

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  // ------------------------------------------------------
  // WRITE THE OUTPUT TO FILE
  // ------------------------------------------------------
  simulation.write_output("output_naca_rds_adapt.msh", options);

  mpi_env.finalize();

  return 0;
}
