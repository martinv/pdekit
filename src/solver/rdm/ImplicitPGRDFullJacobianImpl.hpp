#ifndef PDEKIT_Implicit_PGRD_Full_Jacobian_Implementation_hpp
#define PDEKIT_Implicit_PGRD_Full_Jacobian_Implementation_hpp

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <thread>

#include "graph/GraphReordering.hpp"
#include "linear_system/LSAssembler.hpp"
#include "linear_system/LSTpetra.hpp"
#include "solver/rdm/assembly/PGRDMFinDiffCellAssemblyWorker.hpp"
#include "solver/rdm/assembly/PGRDMPicardJacCellAssemblyWorker.hpp"
#include "solver/rdm/bc/RDMBCManager.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

namespace detail
{

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
class ImplicitPGRDFullJacobianImpl
{
  public:
  using f_space_cells = interpolation::FunctionSpace<MeshConfig>;
  using bc_base_type  = RDMBCBase<MeshConfig, Physics, Physics::DIM - 1>;

  /// TYPEDEFS
  using sol_cache_type = interpolation::SolutionCache;
  /// TYPES
  struct BCScratchDataType
  {
    sol_cache_type residual_cache;
    sol_cache_type residual_cache_perturbed;
  };

  /// Default constructor
  ImplicitPGRDFullJacobianImpl();

  /// Default destructor
  ~ImplicitPGRDFullJacobianImpl();

  void configure_cell_spaces(typename result_of::tria_t<MeshConfig> const &tria,
                             typename result_of::dof_map_t<MeshConfig> const &sol_cells,
                             ls::TpetraCrsMatrix<Real> &mat, ls::TpetraMultiVector<Real> &rhs,
                             ls::TpetraMultiVector<Real> &du, const Uint nb_threads,
                             const Uint nb_blocks, const SFunc sf_type, const PointSetID quad_type,
                             const Uint quad_order);

  /// Reorder degrees of freedom to improve conditioning of system matrix
  void compute_node_reordering_impl(const mesh::Tria<MeshConfig> &tria,
                                    typename result_of::dof_map_t<MeshConfig> const &dofs,
                                    std::vector<Int> &reordering);

  /// Reorder degrees of freedom to improve conditioning of system matrix
  void compute_cell_reordering_impl(const mesh::Tria<MeshConfig> &tria,
                                    typename result_of::dof_map_t<MeshConfig> const &dofs,
                                    std::vector<Int> &reordering);

  void set_sources(const interpolation::VectorMeshFunction<Real>::ptr &sources);

  void set_blending_coeff(const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff);

  void set_artificial_viscosity(
      const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity);

  void assemble_lhs_and_rhs_impl(
      const mesh::Tria<MeshConfig> &mesh,
      typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
      std::vector<bool> const &is_Dirichlet_node, RDMBCManager<MeshConfig, Physics> &bc_manager,
      std::map<std::string, BCScratchDataType *> &weak_bc_data,
      interpolation::VectorMeshFunction<Real> &solution, ls::TpetraCrsMatrix<Real> &mat,
      ls::TpetraMultiVector<Real> &rhs, RDTimeUpdate &time_update, const Real CFL);

  void assemble_rhs_impl(const mesh::Tria<MeshConfig> &mesh,
                         typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
                         std::vector<bool> const &is_Dirichlet_node,
                         RDMBCManager<MeshConfig, Physics> &bc_manager,
                         std::map<std::string, BCScratchDataType *> &weak_bc_data,
                         interpolation::VectorMeshFunction<Real> &solution,
                         ls::TpetraMultiVector<Real> &rhs, RDTimeUpdate &time_update);

  void solve(typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
             std::shared_ptr<ls::TpetraCrsMatrix<Real>> &mat,
             std::shared_ptr<ls::TpetraMultiVector<Real>> &rhs,
             std::shared_ptr<ls::TpetraMultiVector<Real>> &du, ls::LSTpetra<Real> &linear_system,
             interpolation::VectorMeshFunction<Real> &solution,
             const bool recompute_preconditioner);

  /// This method transfers the accumulated residuals from the RHS of the
  /// linear system to a vector mesh function
  void update_residuals(ls::TpetraMultiVector<Real> const &rhs,
                        interpolation::VectorMeshFunction<Real> &residuals);

  void compute_residual_norm(typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
                             ls::TpetraMultiVector<Real> const &rhs,
                             math::DenseDVec<Real> &norm_L2) const;

  private:
  using assembly_worker_type = PGRDMFinDiffCellAssemblyWorker<MeshConfig, Physics, SchemeTraits>;
  // using assembly_worker_type = PGRDMPicardJacCellAssemblyWorker<MeshConfig, Physics,
  // SchemeTraits>;

  /// METHODS
  void create_system_matrix_structure(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                                      ls::TpetraCrsMatrix<Real> &mat,
                                      ls::TpetraMultiVector<Real> &rhs,
                                      ls::TpetraMultiVector<Real> &du);

  /// DATA
  std::vector<std::thread> m_threads;
  std::vector<std::unique_ptr<assembly_worker_type>> m_assembly_workers;

  std::vector<std::tuple<Uint, Uint, Real>> m_mat_buffer;
  std::vector<std::tuple<Uint, Real>> m_rhs_buffer;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::ImplicitPGRDFullJacobianImpl()
{
  m_threads.resize(0);
  m_assembly_workers.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::~ImplicitPGRDFullJacobianImpl()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::configure_cell_spaces(
    typename result_of::tria_t<MeshConfig> const &tria,
    typename result_of::dof_map_t<MeshConfig> const &sol_cells, ls::TpetraCrsMatrix<Real> &mat,
    ls::TpetraMultiVector<Real> &rhs, ls::TpetraMultiVector<Real> &du, const Uint nb_threads,
    const Uint nb_blocks, const SFunc sf_type, const PointSetID quad_type, const Uint quad_order)
{
  m_threads.resize(nb_threads);
  m_assembly_workers.clear();
  m_assembly_workers.resize(0);

  const Uint total_nb_cells      = sol_cells.nb_active_cells();
  const Uint nb_cells_per_thread = (sol_cells.nb_active_cells() % nb_threads == 0)
                                       ? sol_cells.nb_active_cells() / nb_threads
                                       : (sol_cells.nb_active_cells() + nb_threads) / nb_threads;

  /*
  std::cout << "Going to use " << m_workers.size() << " threads" << std::endl;
  std::cout << "Total number of cells = " << total_nb_cells << std::endl;
  std::cout << "Number of cells per thread = " << nb_cells_per_thread <<
  std::endl;
  */

  for (Uint t = 0; t < nb_threads; ++t)
  {
    std::unique_ptr<assembly_worker_type> worker(new assembly_worker_type(t));

    const Uint first_cell = t * nb_cells_per_thread;
    const Uint last_cell  = std::min((t + 1) * nb_cells_per_thread - 1, total_nb_cells - 1);

    // std::cout << "Configuring worker to process cells [" << first_cell <<
    // ","
    // << last_cell << "]"
    //           << std::endl;
    worker->configure_cell_spaces(tria, sol_cells, first_cell, last_cell, nb_blocks, sf_type,
                                  quad_type, quad_order);

    m_assembly_workers.push_back(std::move(worker));
  }

  // Prepare the linear system
  create_system_matrix_structure(sol_cells, mat, rhs, du);

  return;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::compute_node_reordering_impl(
    const mesh::Tria<MeshConfig> &tria, typename result_of::dof_map_t<MeshConfig> const &dofs,
    std::vector<Int> &reordering)
{
  ls::LSAssembler assembler;

  graph::Graph<Int> mesh_dof_pattern(dofs.nb_nodes());
  assembler.build_dof_sparsity_pattern_continuous<MeshConfig>(dofs, mesh_dof_pattern);

  // Save the mesh graph before reordering
  // mesh_dof_pattern.write_to_file("mesh_sparsity_graph_before.dat");

  graph::GraphReordering::compute_reverse_cuthill_mckee(mesh_dof_pattern, reordering);

  // Reorder dofs
  // dofs.renumber_dofs(reordering);

  // For postprocessing purposes, apply reordering to the mesh graph as well
  mesh_dof_pattern.apply_reordering(reordering);
  // mesh_dof_pattern.write_to_file("mesh_sparsity_graph_after.dat");
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::compute_cell_reordering_impl(
    const mesh::Tria<MeshConfig> &tria, typename result_of::dof_map_t<MeshConfig> const &dofs,
    std::vector<Int> &reordering)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::set_sources(
    const interpolation::VectorMeshFunction<Real>::ptr &sources)
{
  /*
  m_source_cache.allocate(elem_type_map, m_nb_blocks, Physics::NEQ);
  m_source_metric.allocate_buffer(elem_type_map, m_nb_blocks, Physics::NEQ);
  */

  for (Uint t = 0; t < m_assembly_workers.size(); ++t)
  {
    m_assembly_workers[t]->set_sources(sources);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::set_blending_coeff(
    const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff)
{
  for (Uint t = 0; t < m_assembly_workers.size(); ++t)
  {
    m_assembly_workers[t]->set_blending_coeff(blending_coeff);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::set_artificial_viscosity(
    const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity)
{
  for (Uint t = 0; t < m_assembly_workers.size(); ++t)
  {
    m_assembly_workers[t]->set_artificial_viscosity(artificial_viscosity);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::assemble_lhs_and_rhs_impl(
    const mesh::Tria<MeshConfig> &mesh,
    typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
    std::vector<bool> const &is_Dirichlet_node, RDMBCManager<MeshConfig, Physics> &bc_manager,
    std::map<std::string, BCScratchDataType *> &weak_bc_data,
    interpolation::VectorMeshFunction<Real> &solution, ls::TpetraCrsMatrix<Real> &mat,
    ls::TpetraMultiVector<Real> &rhs, RDTimeUpdate &time_update, const Real CFL)
{

  // using cell_dofs_type = typename result_of::dof_map_t<MeshConfig>;
  using bc_manager_type = RDMBCManager<MeshConfig, Physics>;

  // **************************************************************************
  // APPLY STRONG BCs
  // **************************************************************************
  for (typename bc_manager_type::const_iterator it = bc_manager.cbegin(); it != bc_manager.cend();
       ++it)
  {
    if (it->second->bc_type() == BC_TYPE_STRONG)
    {
      // Apply the boundary condition
      it->second->apply(time_update);
    }
  }

  // Prepare the system matrix and RHS for assembly
  mat.unlock();
  mat.fill(0.0);
  rhs.fill(0.0);

  // m_residuals->fill(0.0);
  time_update.reset_wave_speeds();

  // **************************************************************************
  // PROCESS TERMS COMING INTERIOR CELLS
  // **************************************************************************

  /*
  const std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  */

  for (Uint t = 0; t < m_assembly_workers.size(); ++t)
  {
    m_threads[t] = std::thread(
        &assembly_worker_type::assemble_mat_and_rhs_part, &(*m_assembly_workers[t]), std::ref(mesh),
        std::ref(sol_cell_connectivity), std::ref(solution), std::ref(time_update),
        std::ref(is_Dirichlet_node), std::ref(mat), std::ref(rhs));
  }

  std::for_each(m_threads.begin(), m_threads.end(), std::mem_fn(&std::thread::join));

  /*
  const std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::ofstream outfile;
  outfile.open("cell_workers_timing.dat", std::ios::app);
  const double cpu_duration = (c_end - c_start) / (double)CLOCKS_PER_SEC;
  auto wclock_duration = 0.001 * std::chrono::duration<double,
  std::milli>(t_end
  - t_start).count(); outfile << "full cpu/wall " << std::setw(10) <<
  std::setprecision(10) << cpu_duration << " "
          << std::setw(10) << std::setprecision(10) << wclock_duration <<
  std::endl; outfile.close();
  */

  std::vector<Real> residual_finite_diff(1); //(Physics::NEQ);
  std::vector<Int> global_col_indexes(1);    //(Physics::NEQ);

  // **************************************************************************
  // PROCESS TERMS COMING FROM WEAK BOUNDARY CONDITIONS
  // **************************************************************************

  for (typename std::map<std::string, BCScratchDataType *>::const_iterator it =
           weak_bc_data.cbegin();
       it != weak_bc_data.cend(); ++it)
  {
    const std::string &bc_condition_name = it->first;

    std::shared_ptr<bc_base_type> weak_bc    = bc_manager.get_bc(bc_condition_name);
    BCScratchDataType &bc_data               = *it->second;
    sol_cache_type &residual_cache           = bc_data.residual_cache;
    sol_cache_type &residual_cache_perturbed = bc_data.residual_cache_perturbed;

    m_mat_buffer.resize(0);
    m_rhs_buffer.resize(0);

    weak_bc->global_lsys_fin_diff_jacobian_entries(
        time_update, residual_cache, residual_cache_perturbed, m_mat_buffer, m_rhs_buffer);

    mat.add_values(m_mat_buffer);
    rhs.add_values(0, m_rhs_buffer);
  }

  // **************************************************************************
  // 6) ADD DIAGONAL TERMS TO THE SYSTEM JACOBIAN
  // **************************************************************************

  const interpolation::ScalarMeshFunction<Real> &nodal_dual_volume =
      time_update.nodal_dual_volume();

  // using node_value_type = typename
  // interpolation::VectorMeshFunction<Real>::entry_type;

  const Int tot_nb_nodes = sol_cell_connectivity.nb_nodes();

  mat.unlock();

#if 1

  for (Int n = 0; n < tot_nb_nodes; ++n)
  {
    // Note that dt = |S|_i / (sum k+) = nodal_dual_volume[n] /
    // time_update.wave_speed(n) Hence time_update.wave_speed(n) /
    // nodal_dual_volume[n] = 1/ dt
    residual_finite_diff[0] = time_update.wave_speed(n) / (nodal_dual_volume[n] * CFL);

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      const Int global_row_idx = n * Physics::NEQ + eq;
      global_col_indexes[0]    = global_row_idx;

      mat.add_values_to_row(global_row_idx, residual_finite_diff, global_col_indexes);
    }
  }

#else

  Real one_over_min_time_step = time_update.wave_speed(0) / (nodal_dual_volume[0] * CFL);

  for (Int n = 0; n < tot_nb_nodes; ++n)
  {
    const Real time_step   = time_update.wave_speed(n) / (nodal_dual_volume[n] * CFL);
    one_over_min_time_step = std::max(one_over_min_time_step, time_step);
  }

  residual_finite_diff[0] = one_over_min_time_step;

  for (Int n = 0; n < tot_nb_nodes; ++n)
  {
    // Note that dt = |S|_i / (sum k+) = nodal_dual_volume[n] /
    // time_update.wave_speed(n) Hence time_update.wave_speed(n) /
    // nodal_dual_volume[n] = 1/ dt

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      const Int global_row_idx = n * Physics::NEQ + eq;
      global_col_indexes[0]    = global_row_idx;

      mat.add_values_to_row(global_row_idx, residual_finite_diff, global_col_indexes);
    }
  }

#endif

  mat.lock();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::assemble_rhs_impl(
    const mesh::Tria<MeshConfig> &mesh,
    typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
    std::vector<bool> const &is_Dirichlet_node, RDMBCManager<MeshConfig, Physics> &bc_manager,
    std::map<std::string, BCScratchDataType *> &weak_bc_data,
    interpolation::VectorMeshFunction<Real> &solution, ls::TpetraMultiVector<Real> &rhs,
    RDTimeUpdate &time_update)
{
  // using cell_dofs_type = typename result_of::dof_map_t<MeshConfig>;
  using bc_manager_type = RDMBCManager<MeshConfig, Physics>;

  // **************************************************************************
  // APPLY STRONG BCs
  // **************************************************************************

  for (typename bc_manager_type::const_iterator it = bc_manager.cbegin(); it != bc_manager.cend();
       ++it)
  {
    if (it->second->bc_type() == BC_TYPE_STRONG)
    {
      // Apply the boundary condition
      it->second->apply(time_update);
    }
  }

  rhs.fill(0.0);

  // m_residuals->fill(0.0);
  time_update.reset_wave_speeds();

  // **************************************************************************
  // PROCESS CELL VOLUME TERMS
  // **************************************************************************

  const std::clock_t c_start = std::clock();
  auto t_start               = std::chrono::high_resolution_clock::now();

  for (Uint t = 0; t < m_assembly_workers.size(); ++t)
  {
    m_threads[t] = std::thread(&assembly_worker_type::assemble_rhs_part, &(*m_assembly_workers[t]),
                               std::ref(mesh), std::ref(sol_cell_connectivity), std::ref(solution),
                               std::ref(time_update), std::ref(is_Dirichlet_node), std::ref(rhs));
  }

  std::for_each(m_threads.begin(), m_threads.end(), std::mem_fn(&std::thread::join));

  const std::clock_t c_end = std::clock();
  auto t_end               = std::chrono::high_resolution_clock::now();

  std::ofstream outfile;
  outfile.open("cell_workers_timing.dat", std::ios::app);
  const double cpu_duration = (c_end - c_start) / (double)CLOCKS_PER_SEC;
  auto wclock_duration = 0.001 * std::chrono::duration<double, std::milli>(t_end - t_start).count();
  outfile << "rhs cpu/wall  " << std::setw(10) << std::setprecision(10) << cpu_duration << " "
          << std::setw(10) << std::setprecision(10) << wclock_duration << std::endl;
  outfile.close();

  // **************************************************************************
  // PROCESS TERMS COMING FROM WEAK BOUNDARY CONDITIONS
  // **************************************************************************

  for (typename std::map<std::string, BCScratchDataType *>::const_iterator it =
           weak_bc_data.cbegin();
       it != weak_bc_data.cend(); ++it)
  {
    const std::string &bc_condition_name = it->first;

    std::shared_ptr<bc_base_type> weak_bc = bc_manager.get_bc(bc_condition_name);
    BCScratchDataType &bc_data            = *it->second;
    sol_cache_type &residual_cache        = bc_data.residual_cache;

    m_rhs_buffer.resize(0);
    weak_bc->global_lsys_rhs_entries(time_update, residual_cache, m_rhs_buffer);
    rhs.add_values(0, m_rhs_buffer);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::solve(
    typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
    std::shared_ptr<ls::TpetraCrsMatrix<Real>> &mat,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &rhs,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &du, ls::LSTpetra<Real> &linear_system,
    interpolation::VectorMeshFunction<Real> &solution, const bool recompute_preconditioner)
{
  const Real omega_relax = 1.0;

  // **************************************************************************
  // SOLVE THE LINEAR SYSTEM
  // **************************************************************************

  // m_lin_system.configure(m_mat, m_rhs, m_du, recompute_preconditioner,
  // false);
  if (recompute_preconditioner)
  {
    linear_system.update_after_mat_values_change(mat, rhs, du, false);
  }
  linear_system.solve();

  using node_value_type = typename interpolation::VectorMeshFunction<Real>::entry_type;

  for (Uint n = 0; n < sol_cell_connectivity.nb_nodes(); ++n)
  {
    node_value_type u_node = solution.value(n);

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      u_node[eq] += omega_relax * (*du).value(n * Physics::NEQ + eq, 0);
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::update_residuals(
    ls::TpetraMultiVector<Real> const &rhs, interpolation::VectorMeshFunction<Real> &residuals)
{
  residuals.fill(0.0);

  using entry_type     = interpolation::VectorMeshFunction<Real>::entry_type;
  const Uint nb_fields = residuals.nb_fields(); // This should be equal to Physics::NEQ

  for (Uint n = 0; n < residuals.nb_entries(); ++n)
  {
    entry_type node_value = residuals.value(n);
    for (Uint f = 0; f < nb_fields; ++f)
    {
      node_value[f] = -rhs.value(n * nb_fields + f);
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::compute_residual_norm(
    typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
    ls::TpetraMultiVector<Real> const &rhs, math::DenseDVec<Real> &norm_L2) const
{
  const Uint nb_nodes = sol_cell_connectivity.nb_nodes();

  if (norm_L2.size() != Physics::NEQ)
  {
    norm_L2.resize(Physics::NEQ);
  }

  norm_L2.fill(0.0);

  for (Uint n = 0; n < nb_nodes; ++n)
  {
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      norm_L2[eq] += rhs.value(n * Physics::NEQ + eq, 0) * rhs.value(n * Physics::NEQ + eq, 0);
    }
  }

  for (Uint eq = 0; eq < Physics::NEQ; ++eq)
  {
    norm_L2[eq] = std::sqrt(norm_L2[eq] / nb_nodes);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits>::
    create_system_matrix_structure(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                                   ls::TpetraCrsMatrix<Real> &mat, ls::TpetraMultiVector<Real> &rhs,
                                   ls::TpetraMultiVector<Real> &du)
{

  ls::LSAssembler assembler;
  graph::Graph<Int> ls_dofs(sol_dofs.nb_nodes() * Physics::NEQ);
  assembler.build_matrix_sparsity_pattern_continuous<MeshConfig>(Physics::NEQ, sol_dofs, ls_dofs);

  /*
  ls_dofs.write_to_svg("matrix_graph_pg_rd_implicit.svg");
  ls_dofs.write_to_file("matrix_graph_pg_rd_implicit.dat");

  std::unique_ptr<common::BlockArray<Int, Int>> nonzero_coords(new
  common::BlockArray<Int, Int>()); ls_dofs.compress_to_crs(*nonzero_coords);

  math::MatrixSparsityPattern<Int> spattern;
  spattern.build_sparsity(std::move(nonzero_coords));
  spattern.print_vtu("matrix_graph_pg_rd_implicit.vtu");
  */

  /*
  graph::Graph<Int> mesh_dof_pattern(sol_dofs.nb_nodes());
  assembler.template
  build_dof_sparsity_pattern_continuous<MeshConfig>(sol_dofs,
  mesh_dof_pattern);
  mesh_dof_pattern.write_to_file("mesh_sparsity_graph_before.dat");
  std::vector<Int> reordering;
  graph::GraphReordering::compute_reverse_cuthill_mckee(mesh_dof_pattern,
  reordering); mesh_dof_pattern.apply_reordering(reordering);
  mesh_dof_pattern.write_to_file("mesh_sparsity_graph_after.dat");
  */

  /*
  std::vector<Int> reordering;
  graph::GraphReordering::compute_reverse_cuthill_mckee(ls_dofs, reordering);
  ls_dofs.apply_reordering(reordering);
  ls_dofs.write_to_file("matrix_graph_after.dat");
  */

  // using cell_dofs_type = typename result_of::dof_map_t<MeshConfig>;

#if 0
  const Uint ls_size = sol_dofs.nb_nodes() * Physics::NEQ;

  Uint max_row_len = 0;
  for (Uint i = 0; i < ls_dofs.size(); ++i)
  {
    max_row_len = std::max(max_row_len, static_cast<Uint>(ls_dofs[i].size()));
  }
#endif

  std::cout << "  Creating matrix structure in continuous implicit solver" << std::endl;

  // m_mat->init(ls_size, max_row_len);
  mat.init(ls_dofs);
  du.init(mat.map(), 1);
  rhs.init(mat.map(), 1);

  // m_mat->unlock();
  // m_mat->fill(0.0);
  rhs.fill(0.0);

  // m_lin_system.initialize_solver(m_mat, m_rhs, m_du, false);

  std::cout << "  Matrix structure created" << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
