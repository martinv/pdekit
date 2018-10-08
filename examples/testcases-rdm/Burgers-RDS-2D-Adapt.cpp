#include <ctime>
#include <iomanip>
#include <iostream>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/mesh_function/MeshFunctionSnapshot.hpp"
#include "interpolation/mesh_function/MeshFunctionTools.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/scalar/Burgers2D.hpp"
#include "solver/SolverIO.hpp"
#include "solver/art_visc/ArtificialViscosity.hpp"
#include "solver/rdm/DRDMethodImplicit.hpp"
#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/blending_coeff/RDBlendingCoeff.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"
#include "tools/MeshInspector.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

// ----------------------------------------------------------------------------

struct SimulationOptions
{
  SimulationOptions();

  Uint nb_iter;
  Uint freeze_artif_visc_at_iter;
  Real CFLMin;
  Real CFLMax;
  Uint p_max;
  bool use_artificial_viscosity;
  bool use_blending;
};

// Constructor to set default values
SimulationOptions::SimulationOptions()
    : nb_iter(1500), freeze_artif_visc_at_iter(10000), CFLMin(0.7), CFLMax(0.8), p_max(5),
      use_artificial_viscosity(true), use_blending(false)
{
}

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

  void prepare_solver(const bool reset_solution, const bool reorder_dofs,
                      const bool use_artificial_viscosity, const bool use_blending);

  void iterate(const SimulationOptions &options);

  void h_adapt(const h_AdaptStrategy adapt_strategy, const SimulationOptions &options);

  void p_adapt(const Uint p_max, const SimulationOptions &options);

  void hp_adapt(const h_AdaptStrategy adapt_strategy, const SimulationOptions &options);

  void write_output(const std::string &outfilename, const SimulationOptions &options);

  void write_final_viz_output(const std::string &outfilename);

  void check_mesh() const;

  private:
  typedef Tria<MeshConfig> MeshType;
  typedef typename Tria<MeshConfig>::dof_storage_type DofStorageType;
  typedef typename interpolation::ScalarMeshFunction<Real> scalar_function;
  typedef typename interpolation::VectorMeshFunction<Real> vector_function;

  static Real left_bc_expr(
      const math::DenseConstVecView<Real> &point_coord,
      const typename interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
      const Uint component);

  static Real right_bc_expr(
      const math::DenseConstVecView<Real> &point_coord,
      const typename interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
      const Uint component);

  static Real top_bc_expr(
      const math::DenseConstVecView<Real> &point_coord,
      const typename interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
      const Uint component);

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

  std::unique_ptr<interpolation::ScalarMeshFunction<Uint>> m_cell_dist;

  interpolation::MeshFunctionSnapshot<Real> mesh_function_snapshot;

  // Vector of new node ids after reordering
  std::vector<Int> m_cell_reordering;
  std::vector<Int> m_node_reordering;

  // Iteration counter
  Uint m_iter;

  // Quadrature order
  PolyOrderID m_quadrature_order;
};

// ----------------------------------------------------------------------------

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

  // const std::string infilename = "square_quad_p3.msh";
  const std::string infilename = "square_tri_p2.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(mesh_name, *m_mesh, "geo_dofs");

  m_geo_dofs = m_mesh->dof_storage("geo_dofs");
  m_sol_dofs = m_mesh->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_discontinuous(*m_mesh, *m_geo_dofs, *m_sol_dofs, poly_order,
                                                  PointSetID::Equidist);

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
  // m_art_visc->resize((*m_geo_dofs).nb_active_cells());
  m_art_visc->resize((*m_sol_dofs).nb_nodes());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::setup_blending_coefficient(
    const std::string &blend_coeff_type)
{
  m_blend_coeff_type = blend_coeff_type;

  rd_blending_coeff_operator.setup(*m_sol_dofs, m_blend_coeff_type);

  // This is reference value for Wong-Jansen blending coefficient
  rd_blending_coeff_operator.set_param("u_inf", 1.0);

  m_blending_coeff = std::make_shared<scalar_function>("", "blending_coefficient");
  m_blending_coeff->resize((*m_sol_dofs).nb_nodes());
  m_blending_coeff->fill(0.0);
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
  math::DenseDVec<Real> res_L1_norm(1);

  for (Uint iter = 0; iter < options.nb_iter; ++iter)
  {
    // const Real CFL = iter < 100 ? 3.0 : std::min(15.0, 4.0 +
    // std::pow(1.30, iter - 100));
    const Real CFL = iter < 100
                         ? options.CFLMin
                         : std::min(options.CFLMax, options.CFLMin + std::pow(1.30, iter - 100));

    m_solver->set_cfl(CFL);

    if (iter % 5 == 0)
    {
      m_solver->assemble_lhs_and_rhs(CFL);
      m_solver->solve({SolverOption::RecomputePreconditioner});
    }
    else
    {
      m_solver->assemble_rhs();
      m_solver->solve({});
    }

    /*
    rd_blending_coeff_operator.calculate(geo_mesh->topology().dof_storage(),
                                         sol_mesh->topology().dof_storage(),
    *solution, *blending_coeff);
    */

    if ((iter < options.freeze_artif_visc_at_iter) && (iter > 0) && ((iter % 5) == 0))
    {
      if (options.use_artificial_viscosity)
      {
        m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_sol_dofs, *m_solution,
                                                               *m_art_visc, 7.0);
      }

      if (options.use_blending)
      {
        rd_blending_coeff_operator.calculate(*m_mesh, *m_sol_dofs, *m_solution, *m_blending_coeff);
      }
    }

    if (iter % 10 == 0)
    {
      m_solver->compute_residual_norm(res_L1_norm);
      solver::SolverIO::print_iter_and_res_norm(m_iter, CFL, res_L1_norm);
    }

    m_iter++;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::h_adapt(const h_AdaptStrategy adapt_strategy,
                                              const SimulationOptions &options)
{
  // --------------------------------------------
  write_output("debug_before_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_before_adapt.msh");
  // --------------------------------------------

  const Uint old_nb_active_cells = (*m_sol_dofs).nb_active_cells();

  m_adapt_cell_marker->resize(old_nb_active_cells);
  cell_adapt_ops.resize(old_nb_active_cells);

  m_indicator->resize(old_nb_active_cells);
  m_art_visc_operator.compute_indicator(*m_sol_dofs, *m_solution, *m_indicator, 5.0);

  // Adapt based on blending operator
  if (options.use_blending)
  {
    rd_blending_coeff_operator.calculate_in_cells(*m_mesh, *m_sol_dofs, *m_solution, *m_indicator);
  }

#define USE_CELL_DISTANCE 0
#if USE_CELL_DISTANCE

  std::vector<mesh::ActiveIdx> level_zero_cells;

  for (Uint c = 0; c < old_nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > 5.e-5)
    {
      level_zero_cells.push_back(mesh::ActiveIdx(c));
    }
  }

  interpolation::ScalarMeshFunction<Uint> cell_dist("", "cell_distance");
  const common::ArrayView<const mesh::ActiveIdx, Uint> level_zero_view(level_zero_cells.data(),
                                                                       level_zero_cells.size());
  interpolation::MeshFunctionTools::compute_cell_distance(*m_mesh, level_zero_view, cell_dist);

  std::vector<Real> indicator_vector;
  indicator_vector.resize(m_indicator->nb_entries());

  for (Uint c = 0; c < cell_dist.nb_entries(); ++c)
  {
    if (cell_dist[c] < 3)
    {
      (*m_adapt_cell_marker)[c] = 1;
      cell_adapt_ops[c]         = CellTransform::UNIFORM_REFINE;
    }
    else
    {
      (*m_indicator)[c] = 0;
      cell_adapt_ops[c] = CellTransform::DO_NOTHING;
    }
    indicator_vector[c] = (*m_indicator)[c];
  }

#else

  std::vector<Real> indicator_vector;
  indicator_vector.resize(m_indicator->nb_entries());

  const Real refine_threshold = (options.use_artificial_viscosity) ? 5.e-5 : 0.1;

  for (Uint c = 0; c < old_nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > refine_threshold)
    {
      (*m_adapt_cell_marker)[c] = 1;
      cell_adapt_ops[c]         = CellTransform::UNIFORM_REFINE;
    }
    else
    {
      (*m_indicator)[c] = 0;
      cell_adapt_ops[c] = CellTransform::NO_TRANS;
    }
    indicator_vector[c] = (*m_indicator)[c];
  }
#endif

  // --------------------------------------------

  if (options.use_blending)
  {
    rd_blending_coeff_operator.calculate(*m_mesh, *m_sol_dofs, *m_solution, *m_blending_coeff);
  }

  if (options.use_artificial_viscosity)
  {
    m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_sol_dofs, *m_solution,
                                                           *m_art_visc, 7.0);
  }

  // --------------------------------------------

  // adapt_schedule.define_h_adapt_ops(*m_mesh, cell_adapt_ops);

  const Real refinement_ratio = 0.1;
  adapt_schedule.define_h_adapt_ops(*m_mesh, cell_adapt_ops, indicator_vector, refinement_ratio,
                                    adapt_strategy);

  if (adapt_strategy == h_AdaptStrategy::w_hanging_nodes)
  {
    mesh_function_snapshot.create<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);
  }

  m_mesh->adapt(adapt_schedule);

  /*
  (*m_geo_dofs).adapt(adapt_schedule);
  (*m_sol_dofs).adapt(adapt_schedule);
  */

  (*m_geo_dofs).adapt_update(*m_mesh);
  (*m_sol_dofs).adapt_update(*m_mesh);

  tools::MeshInspector<MeshConfig> mesh_inspector;
  const Real volume =
      mesh_inspector.check_mesh_consistency(*m_mesh, "geo_dofs", m_quadrature_order);
  std::cout << "Volume computed by mesh inspector = " << volume << std::endl;

  if (adapt_strategy == h_AdaptStrategy::w_hanging_nodes)
  {
    mesh_function_snapshot.restore_function<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);
  }
  else if (adapt_strategy == h_AdaptStrategy::red_green)
  {
    const Uint sol_dofs_nb_nodes = (*m_sol_dofs).nb_nodes();
    m_solution->resize(Physics::NEQ, sol_dofs_nb_nodes);
    m_solution->fill(0.0);
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
                                                           *m_art_visc, 7.0);
  }

  const bool reset_solution = false;
  const bool reorder_dofs   = true;
  prepare_solver(reset_solution, reorder_dofs, options.use_artificial_viscosity,
                 options.use_blending);

  if (reorder_dofs)
  {
    if (options.use_blending)
    {
      m_blending_coeff->apply_reordering(m_node_reordering);
    }
    if (options.use_artificial_viscosity)
    {
      m_art_visc->apply_reordering(m_node_reordering);
    }
  }

  // --------------------------------------------
  // m_mesh->print_complete_skeleton(_1D);
  write_output("debug_after_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_after_adapt.msh");
  // --------------------------------------------

  m_solver->assemble_lhs_and_rhs(options.CFLMin);
  m_solver->solve({SolverOption::RecomputePreconditioner});
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::p_adapt(const Uint p_max, const SimulationOptions &options)
{
  // --------------------------------------------

  write_output("debug_before_p_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_before_p_adapt.msh");

  // --------------------------------------------

  const Uint nb_active_cells = (*m_sol_dofs).nb_active_cells();

  std::vector<Uint> old_cell_p_order(nb_active_cells);
  std::vector<Uint> adapt_cell_p_order(nb_active_cells);
  std::vector<Real> adapt_quality(nb_active_cells);

  m_indicator->resize(nb_active_cells);

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

  // Adapt based on blending operator
  if (options.use_blending)
  {
    rd_blending_coeff_operator.calculate_in_cells(*m_mesh, *m_geo_dofs, *m_sol_dofs, *m_solution,
                                                  *m_indicator);
  }

  Uint poly_order = 0;

  // ---------

  std::vector<mesh::ActiveIdx> level_zero_cells;

  const Uint p1_cell_region_width = 0;
  const Real refine_threshold     = (options.use_artificial_viscosity) ? 5.e-5 : 0.1;

  for (Uint c = 0; c < nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > refine_threshold)
    {
      level_zero_cells.push_back(mesh::ActiveIdx(c));
    }
  }

  const common::ArrayView<const mesh::ActiveIdx, _1D, Uint> level_zero_view(
      level_zero_cells.data(), level_zero_cells.size());
  interpolation::MeshFunctionTools::compute_cell_distance(*m_mesh, level_zero_view, *m_cell_dist);

  /*
  std::vector<Real> indicator_vector;
  indicator_vector.resize(m_indicator->nb_entries());
  */

  for (Uint c = 0; c < m_cell_dist->nb_entries(); ++c)
  {
    const mesh::MeshEntity active_cell = (*m_sol_dofs).active_cell(mesh::ActiveIdx(c));
    adapt_cell_p_order[c]              = active_cell.pt_set_id().poly_order();

    if ((*m_cell_dist)[c] > p1_cell_region_width)
    {
      adapt_cell_p_order[c] = std::min(p_max, 2 + (*m_cell_dist)[c] - p1_cell_region_width);

      /*
      adapt_cell_p_order[c] =
          std::min(p_max, adapt_cell_p_order[c] + (*m_cell_dist)[c] -
      p1_cell_region_width + 1);
      */
    }

    poly_order = std::max(poly_order, adapt_cell_p_order[c]);
  }

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
                                                           *m_solution, *m_art_visc, 7.0);
  }

  // Make sure that the quadrature order to set up the solver is sufficiently
  // high:
  m_quadrature_order        = static_cast<PolyOrderID>(2 * poly_order);
  const bool reset_solution = false;
  const bool reorder_dofs   = true;

  prepare_solver(reset_solution, reorder_dofs, options.use_artificial_viscosity,
                 options.use_blending);

  if (reorder_dofs)
  {
    if (options.use_blending)
    {
      m_blending_coeff->apply_reordering(m_node_reordering);
    }
    if (options.use_artificial_viscosity)
    {
      m_art_visc->apply_reordering(m_node_reordering);
    }
  }

  // --------------------------------------------

  // m_mesh->print_complete_skeleton(_1D);
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

  // Adapt based on blending operator
  if (options.use_blending)
  {
    rd_blending_coeff_operator.calculate_in_cells(*m_mesh, *m_geo_dofs, *m_sol_dofs, *m_solution,
                                                  *m_indicator);
  }

  Uint poly_order = 0;

#define USE_CELL_DISTANCE_HP 0
#if USE_CELL_DISTANCE_HP
  const Uint p1_cell_region_width = 1;

  std::vector<mesh::ActiveIdx> level_zero_cells;
  const Real refine_threshold = (options.use_artificial_viscosity) ? 5.e-5 : 0.1;

  for (Uint c = 0; c < old_nb_active_cells; ++c)
  {
    if ((*m_indicator)[c] > refine_threshold)
    {
      level_zero_cells.push_back(mesh::ActiveIdx(c));
    }
  }

  interpolation::ScalarMeshFunction<Uint> cell_dist("", "cell_distance");
  const common::ArrayView<const mesh::ActiveIdx, _1D, Uint> level_zero_view(
      level_zero_cells.data(), level_zero_cells.size());
  interpolation::MeshFunctionTools::compute_cell_distance(*m_mesh, level_zero_view, cell_dist);

  std::vector<Real> indicator_vector;
  indicator_vector.resize(m_indicator->nb_entries());

  for (Uint c = 0; c < cell_dist.nb_entries(); ++c)
  {
    const mesh::MeshEntity active_cell = (*m_sol_dofs).active_cell(mesh::ActiveIdx(c));
    cell_p_order[c]                    = active_cell.std_region_id().poly_order();

    if (cell_dist[c] < p1_cell_region_width)
    {
      (*m_adapt_cell_marker)[c] = 1;
      cell_adapt_ops[c]         = CellTransform::UNIFORM_REFINE;
    }
    else
    {
      (*m_indicator)[c] = 0;
      cell_adapt_ops[c] = CellTransform::DO_NOTHING;
      cell_p_order[c]   = std::min(options.p_max, cell_dist[c] - p1_cell_region_width + 2);
    }
    indicator_vector[c] = (*m_indicator)[c];
    poly_order          = std::max(poly_order, cell_p_order[c]);
  }

#else

  std::vector<Real> indicator_vector;
  indicator_vector.resize(m_indicator->nb_entries());

  const Real refine_threshold = (options.use_artificial_viscosity) ? 5.e-5 : 0.1;

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

      cell_p_order[c] = active_cell.pt_set_id().poly_order() + P1;
      poly_order      = std::max(poly_order, cell_p_order[c]);
    }
    indicator_vector[c] = (*m_indicator)[c];
  }
#endif

  m_quadrature_order = static_cast<PolyOrderID>(2 * poly_order);

  // --------------------------------------------

  if (options.use_blending)
  {
    rd_blending_coeff_operator.calculate(*m_mesh, *m_geo_dofs, *m_sol_dofs, *m_solution,
                                         *m_blending_coeff);
  }

  if (options.use_artificial_viscosity)
  {
    m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_geo_dofs, *m_sol_dofs,
                                                           *m_solution, *m_art_visc, 7.0);
  }

  write_output("debug_before_hp_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_before_hp_adapt.msh");

  // --------------------------------------------

  // adapt_schedule.define_h_adapt_ops(*m_mesh, cell_adapt_ops);

  const Real refinement_ratio = 0.1;
  adapt_schedule.define_h_adapt_ops(*m_mesh, cell_adapt_ops, indicator_vector, refinement_ratio,
                                    adapt_strategy);

  mesh_function_snapshot.create_cellwise<MeshConfig>(*m_mesh, *m_sol_dofs, cell_p_order);

  if (adapt_strategy == h_AdaptStrategy::w_hanging_nodes)
  {
    mesh_function_snapshot.create<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);
  }

  m_mesh->adapt(adapt_schedule);

  /*
  (*m_geo_dofs).adapt(adapt_schedule);
  (*m_sol_dofs).adapt(adapt_schedule);
  */

  mesh_function_snapshot.restore_function_cellwise<MeshConfig>(*m_mesh, *m_sol_dofs, cell_p_order);

  const common::ArrayView<const Uint, _1D, Uint> cell_p_order_view(&cell_p_order[0],
                                                                   cell_p_order.nb_entries());

  (*m_geo_dofs).adapt_update(*m_mesh, cell_p_order_view);
  (*m_sol_dofs).adapt_update(*m_mesh, cell_p_order_view);

  tools::MeshInspector<MeshConfig> mesh_inspector;
  const Real volume =
      mesh_inspector.check_mesh_consistency(*m_mesh, "geo_dofs", m_quadrature_order);
  std::cout << "Volume computed by mesh inspector = " << volume << std::endl;

  if (adapt_strategy == h_AdaptStrategy::w_hanging_nodes)
  {
    mesh_function_snapshot.restore_function<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);
  }
  else if (adapt_strategy == h_AdaptStrategy::red_green)
  {
    const Uint sol_dofs_nb_nodes = (*m_sol_dofs).nb_nodes();
    m_solution->resize(Physics::NEQ, sol_dofs_nb_nodes);
    m_solution->fill(0.0);
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
                                                           *m_solution, *m_art_visc, 7.0);
  }

  const bool reset_solution = false;
  const bool reorder_dofs   = true;

  prepare_solver(reset_solution, reorder_dofs, options.use_artificial_viscosity,
                 options.use_blending);

  if (reorder_dofs)
  {
    if (options.use_blending)
    {
      m_blending_coeff->apply_reordering(m_node_reordering);
    }
    if (options.use_artificial_viscosity)
    {
      m_art_visc->apply_reordering(m_node_reordering);
    }
  }

  // --------------------------------------------

  // m_mesh->print_complete_skeleton(_1D);
  write_output("debug_after_hp_adapt.msh", options);
  m_mesh->write_dual_graph_to_gmsh("sol_dofs", "skeleton_after_hp_adapt.msh");

  // --------------------------------------------

  m_solver->assemble_lhs_and_rhs(options.CFLMin);
  m_solver->solve({SolverOption::RecomputePreconditioner});
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::write_output(const std::string &outfilename,
                                                   const SimulationOptions &options)
{
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*m_mesh, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *m_solution, "solution");

  if (options.use_artificial_viscosity)
  {
    // Compute the smoothness indicator and artificial viscosity once more
    // to store for visualization
    m_art_visc_operator.compute_indicator_normalized(*m_sol_dofs, *m_solution, *m_indicator, 7.0);
    m_art_visc_operator.compute_artificial_viscosity_nodal(*m_mesh, *m_sol_dofs, *m_solution,
                                                           *m_art_visc, 7.0);

    interpolation::ScalarMeshFunction<Real> cellwise_art_visc("", "artificial_viscosity_cellwise");
    cellwise_art_visc.resize((*m_sol_dofs).nb_active_cells());
    m_art_visc_operator.compute_artificial_viscosity(*m_sol_dofs, *m_solution, cellwise_art_visc,
                                                     7.0);

    meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *m_art_visc,
                                             "artificial_viscosity");
    meshwriter.append_cell_function_to_file(*m_mesh, outfilename, cellwise_art_visc,
                                            "artificial_viscosity_cellwise");
    meshwriter.append_cell_function_to_file(*m_mesh, outfilename, *m_indicator,
                                            "smoothness_indicator");
    meshwriter.append_cell_function_to_file(*m_mesh, outfilename, *m_adapt_cell_marker,
                                            "adapt_cell_marker");
  }

  if (options.use_blending)
  {
    // Compute the blending coefficient and store its values for
    // visualization
    rd_blending_coeff_operator.calculate(*m_mesh, *m_sol_dofs, *m_solution, *m_blending_coeff);
    meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *m_blending_coeff,
                                             "blending_coefficient");
  }

  std::vector<Uint> p_order(m_mesh->nb_active_cells());
  for (Uint ac = 0; ac < m_mesh->nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity active_cell = (*m_sol_dofs).active_cell(mesh::ActiveIdx(ac));
    p_order[ac]                        = active_cell.pt_set_id().poly_order();
  }

  const common::ArrayView<const Uint, _1D, Uint> poly_orders_view(p_order.data(), p_order.size());
  meshwriter.append_cell_function_to_file(*m_mesh, outfilename, poly_orders_view, "poly_order");
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::write_final_viz_output(const std::string &outfilename)
{
  Uint p_order_min = 10;
  Uint p_order_max = 0;

  const Uint nb_active_cells = (*m_sol_dofs).nb_active_cells();

  for (Uint c = 0; c < (*m_sol_dofs).nb_active_cells(); ++c)
  {
    const mesh::MeshEntity active_cell = (*m_sol_dofs).active_cell(mesh::ActiveIdx(c));
    const Uint p_active                = active_cell.pt_set_id().poly_order();

    p_order_min = std::min(p_order_min, p_active);
    p_order_max = std::max(p_order_max, p_active);
  }

  if (p_order_min != p_order_max)
  {
    std::vector<Uint> cell_p_order(nb_active_cells);
    cell_p_order.assign(nb_active_cells, p_order_max);
    adapt_schedule.define_p_adapt_ops(cell_p_order);

    mesh_function_snapshot.create<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);
    (*m_sol_dofs).adapt(adapt_schedule);
    mesh_function_snapshot.restore_function<MeshConfig>(*m_mesh, *m_sol_dofs, *m_solution);
  }

  gmsh::GmshWriter meshwriter;

  meshwriter.write_mesh_to_file(*m_mesh, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*m_mesh, outfilename, *m_solution, "solution");
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void Simulation<MeshConfig, Physics>::check_mesh() const
{
  tools::MeshInspector<MeshConfig>::check_mesh_consistency(*m_mesh, "sol_dofs", P5);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
Real Simulation<MeshConfig, Physics>::left_bc_expr(
    const math::DenseConstVecView<Real> &point_coord,
    const typename interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
    const Uint component)
{
  return 1.5;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
Real Simulation<MeshConfig, Physics>::right_bc_expr(
    const math::DenseConstVecView<Real> &point_coord,
    const typename interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
    const Uint component)
{
  return -0.5;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
Real Simulation<MeshConfig, Physics>::top_bc_expr(
    const math::DenseConstVecView<Real> &point_coord,
    const typename interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
    const Uint component)
{
  return 1.5;
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
    std::cout << "Building LDA scheme ..." << std::endl;
    scheme =
        std::make_shared<rdm::DRDMethodImplicit<MeshConfig, Physics, rdm::PGLDA, rdm::FacetDG>>();
  }
  else if (scheme_type == "B")
  {
    std::cout << "Building B scheme ..." << std::endl;
    scheme =
        std::make_shared<rdm::DRDMethodImplicit<MeshConfig, Physics, rdm::PGB, rdm::FacetDG>>();
  }
  else
  {
    scheme.reset();

    /*
    scheme =
        std::make_shared<rdm::DRDMethodImplicit<MeshConfig, Physics,
    rdm::PGN, rdm::FacetDG>>();
    */
  }

  if (reorder_dofs)
  {
    scheme->compute_cell_reordering(*mesh, *sol_dofs, m_cell_reordering);
    (*sol_dofs).renumber_dofs_blockwise(m_cell_reordering, m_node_reordering);

    solution->apply_reordering(m_node_reordering);
    residuals->apply_reordering(m_node_reordering);
  }

  scheme->configure_mesh_data(mesh, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quad_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residuals);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  /*
  std::shared_ptr<rdm::RDMBCBase<MeshConfig, Physics>> bottom =
      scheme->add_boundary_condition("WeakDirichlet", "bottom_bc", "bottom");
  bottom->set_expression(&BurgersInlet2D::value);
  */

  std::shared_ptr<rdm::RDMBCBase<MeshConfig, Physics>> bottom =
      scheme->add_boundary_condition("StrongDirichlet", "bottom_bc", "bottom");
  bottom->set_expression(&BurgersInlet2D::value);

  /*
  std::shared_ptr<rdm::RDMBCBase<MeshConfig, phys_model>> top =
      scheme->add_boundary_condition("StrongDirichlet", "top_bc", "top");
  top->set_expression(&top_bc_expr);
  */

  std::shared_ptr<rdm::RDMBCBase<MeshConfig, Physics>> left =
      scheme->add_boundary_condition("WeakDirichlet", "left_bc", "left");
  left->set_expression(&left_bc_expr);

  /*
  std::shared_ptr<rdm::RDMBCBase<MeshConfig, Physics>> right =
      scheme->add_boundary_condition("WeakDirichlet", "right_bc", "right");
  */

  std::shared_ptr<rdm::RDMBCBase<MeshConfig, Physics>> right =
      scheme->add_boundary_condition("StrongDirichlet", "right_bc", "right");

  right->set_expression(&right_bc_expr);

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
  // SET INITIAL OPTION VALUES
  // ------------------------------------------------------

  SimulationOptions options;
  options.nb_iter                   = 1500;
  options.freeze_artif_visc_at_iter = 10000;
  options.CFLMin                    = 0.7;
  options.CFLMax                    = 0.8;
  options.use_artificial_viscosity  = true;
  options.use_blending              = false;

  const bool reset_solution = true;
  const bool reorder_dofs   = true;

  Simulation<Cart2D, physics::Burgers2D> simulation;
  simulation.read_mesh("square_tri_p2.msh", P2);

  simulation.setup_artifical_viscosity();
  simulation.setup_viscosity_indicator();
  simulation.setup_blending_coefficient("WongJansen");

  simulation.prepare_solver(reset_solution, reorder_dofs, options.use_artificial_viscosity,
                            options.use_blending);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;

  start = clock();

  simulation.iterate(options);

  // ------------------------
  // h- and hp-adaptation
  // ------------------------

  options.CFLMax  = 0.8;
  options.nb_iter = 700;

  for (Uint adapt_pass = 0; adapt_pass < 4; ++adapt_pass)
  {
    simulation.h_adapt(h_AdaptStrategy::w_hanging_nodes, options);
    // simulation.hp_adapt(h_AdaptStrategy::w_hanging_nodes, options);
    simulation.iterate(options);
  }

  // ------------------------
  // p-adaptation
  // ------------------------

  /*
  options.CFLMin = 0.1;
  options.CFLMax = 0.7;
  options.nb_iter = 700;

  for (Uint adapt_pass = 0; adapt_pass < 2; ++adapt_pass)
  {
    simulation.p_adapt(P5, options);
    simulation.iterate(options);
  }
  */

  options.CFLMin  = options.CFLMax;
  options.CFLMax  = 2.0;
  options.nb_iter = 700;
  simulation.iterate(options);

  // ------------------------

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Write the output to file
  simulation.write_output("output_Burgers_2D_adapt.msh", options);
  simulation.write_final_viz_output("viz_output_Burgers_2D_adapt.msh");

  /*
  vtk::VtkWriter vtkwriter;
  vtkwriter.append_nodal_function_to_file(*geo_mesh, "output.vtu", *solution,
  "solution");
  */

  mpi_env.finalize();

  return 0;
}
