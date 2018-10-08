#ifndef PDEKIT_Implicit_D_RD_Full_Jacobian_Implementation_hpp
#define PDEKIT_Implicit_D_RD_Full_Jacobian_Implementation_hpp

#include <thread>

#include "graph/GraphReordering.hpp"
#include "interpolation/CoarseScaleCorrectionOpBuilder.hpp"
#include "linear_system/CoarseScaleCorrection.hpp"
#include "linear_system/LSAssembler.hpp"
#include "solver/rdm/assembly/PGRDMFinDiffCellAssemblyWorker.hpp"
#include "solver/rdm/assembly/PGRDMFullJacFacetAssemblyWorker.hpp"
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

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
class ImplicitDRDFullJacobianImpl
{
  public:
  using f_space_cells  = interpolation::FunctionSpace<MeshConfig>;
  using f_space_facets = interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1>;
  using bc_base_type   = RDMBCBase<MeshConfig, Physics, Physics::DIM - 1>;

  /// TYPEDEFS
  using sol_cache_type = interpolation::SolutionCache;

  /*
  enum
  {
    needs_volume_derivatives_on_trace =
  FacetScheme::needs_volume_derivatives_on_trace
  };
  */

  /// TYPES
  struct BCScratchDataType
  {
    sol_cache_type residual_cache;
    sol_cache_type residual_cache_perturbed;
  };

  /// Default constructor
  ImplicitDRDFullJacobianImpl();

  /// Default destructor
  ~ImplicitDRDFullJacobianImpl();

  void configure_cell_spaces(const mesh::Tria<MeshConfig> &tria,
                             typename result_of::dof_map_t<MeshConfig> const &sol_cells,
                             ls::TpetraCrsMatrix<Real> &mat, ls::TpetraMultiVector<Real> &rhs,
                             ls::TpetraMultiVector<Real> &du, const Uint nb_threads,
                             const Uint nb_blocks, const SFunc sf_type, const PointSetID quad_type,
                             const Uint quad_order);

  void configure_facet_spaces(const mesh::Tria<MeshConfig> &tria,
                              typename result_of::dof_map_t<MeshConfig> const &sol_cells,
                              const SFunc sf_type, const PointSetID quad_type,
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
      const mesh::Tria<MeshConfig> &tria,
      typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
      std::vector<bool> const &is_Dirichlet_node, RDMBCManager<MeshConfig, Physics> &bc_manager,
      std::map<std::string, BCScratchDataType *> &weak_bc_data,
      interpolation::VectorMeshFunction<Real> &solution, ls::TpetraCrsMatrix<Real> &mat,
      ls::TpetraMultiVector<Real> &rhs, interpolation::VectorMeshFunction<Real> &nodal_residuals,
      RDTimeUpdate &time_update, const Real CFL);

  void assemble_rhs_impl(const mesh::Tria<MeshConfig> &tria,
                         typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
                         std::vector<bool> const &is_Dirichlet_node,
                         RDMBCManager<MeshConfig, Physics> &bc_manager,
                         std::map<std::string, BCScratchDataType *> &weak_bc_data,
                         interpolation::VectorMeshFunction<Real> &solution,
                         ls::TpetraMultiVector<Real> &rhs,
                         interpolation::VectorMeshFunction<Real> &nodal_residuals,
                         RDTimeUpdate &time_update);

  void solve(typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
             std::shared_ptr<ls::TpetraCrsMatrix<Real>> &mat,
             std::shared_ptr<ls::TpetraMultiVector<Real>> &rhs,
             std::shared_ptr<ls::TpetraMultiVector<Real>> &du, ls::LSTpetra<Real> &linear_system,
             interpolation::VectorMeshFunction<Real> &solution,
             const bool recompute_preconditioner);

  void compute_residual_norm(typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
                             ls::TpetraMultiVector<Real> const &rhs,
                             math::DenseDVec<Real> &norm_L2) const;

  private:
  /// TYPEDEFS
  using cell_assembly_worker_type = PGRDMFinDiffCellAssemblyWorker<MeshConfig, Physics, CellScheme>;
  // using cell_assembly_worker_type = PGRDMPicardJacCellAssemblyWorker<MeshConfig, Physics,
  // CellScheme>;

  using facet_assembly_worker_type =
      PGRDMFullJacFacetAssemblyWorker<MeshConfig, Physics, CellScheme, FacetScheme>;

  /// METHODS
  void create_system_matrix_structure(const mesh::Tria<MeshConfig> &tria,
                                      typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                                      ls::TpetraCrsMatrix<Real> &mat,
                                      ls::TpetraMultiVector<Real> &rhs,
                                      ls::TpetraMultiVector<Real> &du);

  /// DATA
  std::vector<std::thread> m_threads;
  std::vector<std::unique_ptr<cell_assembly_worker_type>> m_cell_assembly_workers;
  std::vector<std::unique_ptr<facet_assembly_worker_type>> m_facet_assembly_workers;

  std::vector<std::tuple<Uint, Uint, Real>> m_mat_buffer;
  std::vector<std::tuple<Uint, Real>> m_rhs_buffer;

  /// Flag to determine whether coarse scale correction should be used
  bool m_use_coarse_scale_corr;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme,
                            FacetScheme>::ImplicitDRDFullJacobianImpl()
{
  m_threads.resize(0);
  m_cell_assembly_workers.resize(0);
  m_facet_assembly_workers.resize(0);
  m_use_coarse_scale_corr = false;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme,
                            FacetScheme>::~ImplicitDRDFullJacobianImpl()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::
    configure_cell_spaces(const mesh::Tria<MeshConfig> &tria,
                          typename result_of::dof_map_t<MeshConfig> const &sol_cells,
                          ls::TpetraCrsMatrix<Real> &mat, ls::TpetraMultiVector<Real> &rhs,
                          ls::TpetraMultiVector<Real> &du, const Uint nb_threads,
                          const Uint nb_blocks, const SFunc sf_type, const PointSetID quad_type,
                          const Uint quad_order)
{
  m_threads.resize(nb_threads);
  m_cell_assembly_workers.clear();
  m_cell_assembly_workers.resize(0);
  m_facet_assembly_workers.clear();
  m_facet_assembly_workers.resize(0);

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
    std::unique_ptr<cell_assembly_worker_type> cell_worker(new cell_assembly_worker_type(t));
    std::unique_ptr<facet_assembly_worker_type> facet_worker(new facet_assembly_worker_type(t));

    const Uint first_cell = t * nb_cells_per_thread;
    const Uint last_cell  = std::min((t + 1) * nb_cells_per_thread - 1, total_nb_cells - 1);

    // std::cout << "Configuring worker to process cells [" << first_cell <<
    // ","
    // << last_cell << "]"
    //           << std::endl;
    cell_worker->configure_cell_spaces(tria, sol_cells, first_cell, last_cell, nb_blocks, sf_type,
                                       quad_type, quad_order);
    facet_worker->configure_cell_spaces(sol_cells, first_cell, last_cell, nb_blocks, sf_type,
                                        quad_type, quad_order);

    m_cell_assembly_workers.push_back(std::move(cell_worker));
    m_facet_assembly_workers.push_back(std::move(facet_worker));
  }

  // Prepare the linear system
  create_system_matrix_structure(tria, sol_cells, mat, rhs, du);

  return;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::
    configure_facet_spaces(const mesh::Tria<MeshConfig> &tria,
                           typename result_of::dof_map_t<MeshConfig> const &sol_cells,
                           const SFunc sf_type, const PointSetID quad_type, const Uint quad_order)
{
  const Uint total_nb_facets = tria.active_skeleton_size(MeshConfig::TDIM - 1);
  const Uint nb_threads      = m_facet_assembly_workers.size();

  const Uint nb_facets_per_thread = (total_nb_facets % nb_threads == 0)
                                        ? total_nb_facets / nb_threads
                                        : (total_nb_facets + nb_threads) / nb_threads;

  for (Uint t = 0; t < nb_threads; ++t)
  {
    const Uint first_facet = t * nb_facets_per_thread;
    const Uint last_facet  = std::min((t + 1) * nb_facets_per_thread - 1, total_nb_facets - 1);

    m_facet_assembly_workers[t]->configure_facet_spaces(tria, sol_cells, first_facet, last_facet,
                                                        sf_type, quad_type, quad_order);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::
    compute_node_reordering_impl(const mesh::Tria<MeshConfig> &tria,
                                 typename result_of::dof_map_t<MeshConfig> const &dofs,
                                 std::vector<Int> &reordering)
{
  ls::LSAssembler assembler;

  if (!m_use_coarse_scale_corr)
  {
    graph::Graph<Int> mesh_dof_pattern(dofs.nb_nodes());
    assembler.build_dof_sparsity_pattern_discontinuous<MeshConfig>(tria, dofs, mesh_dof_pattern);

    // Save the mesh graph before reordering
    // mesh_dof_pattern.write_to_file("mesh_sparsity_graph_before.dat");

    graph::GraphReordering::compute_reverse_cuthill_mckee(mesh_dof_pattern, reordering);
    // Reorder dofs
    // dofs.renumber_dofs(reordering);

    // For postprocessing purposes, apply reordering to the mesh graph as
    // well
    mesh_dof_pattern.apply_reordering(reordering);
    // mesh_dof_pattern.write_to_file("mesh_sparsity_graph_after.dat");
  }
  else
  {
    common::BlockArray<Uint, Uint> mesh_dual_graph_crs;

    tria.build_dual_graph_undirected(mesh_dual_graph_crs);

    // For each old active cell, 'cell_reordering_old_to_new' stores the new
    // active index on position i Note that we're actually NOT going to
    // renumber the cells, only change the order in which the cells are
    // assigned their dof numbers
    std::vector<Uint> cell_reordering_old_to_new;

    // Inverse map to the above: for each new active index on position i,
    // store what the old active index was
    std::vector<Int> cell_reordering_new_to_old;
    graph::GraphReordering::compute_reverse_cuthill_mckee(mesh_dual_graph_crs,
                                                          cell_reordering_old_to_new);

    // Generate the inverse map
    cell_reordering_new_to_old.resize(cell_reordering_old_to_new.size());

    for (Uint i = 0; i < cell_reordering_old_to_new.size(); ++i)
    {
      cell_reordering_new_to_old[cell_reordering_old_to_new[i]] = i;
    }

    reordering.resize(dofs.nb_nodes());

    Uint dof_idx = 0;
    for (Uint ac = 0; ac < cell_reordering_new_to_old.size(); ++ac)
    {
      const Uint old_ac_id        = cell_reordering_new_to_old[ac];
      const mesh::MeshEntity cell = dofs.active_cell(mesh::ActiveIdx(old_ac_id));
      for (Uint v = 0; v < cell.nb_vert(); ++v)
      {
        reordering[cell.vertex(v)] = dof_idx++;
      }
    }

    /*
    std::cout << "Cell reordering: " << std::endl;
    for (Uint i = 0; i < cell_reordering_old_to_new.size(); ++i)
    {
      std::cout << cell_reordering_old_to_new[i] << " ";
    }
    std::cout << std::endl;
    */

  } // if m_use_coarse_scale_corr
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::
    compute_cell_reordering_impl(const mesh::Tria<MeshConfig> &tria,
                                 typename result_of::dof_map_t<MeshConfig> const &dofs,
                                 std::vector<Int> &reordering)
{
  common::BlockArray<Int, Uint> mesh_dual_graph_crs;
  tria.build_dual_graph_undirected(mesh_dual_graph_crs);
  graph::GraphReordering::compute_reverse_cuthill_mckee(mesh_dual_graph_crs, reordering);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::set_sources(
    const interpolation::VectorMeshFunction<Real>::ptr &sources)
{
  /*
  m_source_cache.allocate(elem_type_map, m_nb_blocks, Physics::NEQ);
  m_source_metric.allocate_buffer(elem_type_map, m_nb_blocks, Physics::NEQ);
  */

  for (Uint t = 0; t < m_cell_assembly_workers.size(); ++t)
  {
    m_cell_assembly_workers[t]->set_sources(sources);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::set_blending_coeff(
    const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff)
{
  for (Uint t = 0; t < m_cell_assembly_workers.size(); ++t)
  {
    m_cell_assembly_workers[t]->set_blending_coeff(blending_coeff);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::
    set_artificial_viscosity(
        const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity)
{
  for (Uint t = 0; t < m_cell_assembly_workers.size(); ++t)
  {
    m_cell_assembly_workers[t]->set_artificial_viscosity(artificial_viscosity);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::
    assemble_lhs_and_rhs_impl(
        const mesh::Tria<MeshConfig> &tria,
        typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
        std::vector<bool> const &is_Dirichlet_node, RDMBCManager<MeshConfig, Physics> &bc_manager,
        std::map<std::string, BCScratchDataType *> &weak_bc_data,

        interpolation::VectorMeshFunction<Real> &solution, ls::TpetraCrsMatrix<Real> &mat,
        ls::TpetraMultiVector<Real> &rhs, interpolation::VectorMeshFunction<Real> &nodal_residuals,
        RDTimeUpdate &time_update, const Real CFL)
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
  // PROCESS TERMS COMING INTERIOR FACETS
  // **************************************************************************

  for (Uint t = 0; t < m_facet_assembly_workers.size(); ++t)
  {
    /*
    m_threads[t] = std::thread(
        &assembly_worker_type::assemble_mat_and_rhs_facets_part,
    std::ref(*m_assembly_workers[t]), std::ref(tria),
    std::ref(geo_cell_connectivity), std::ref(sol_cell_connectivity),
        std::ref(solution), std::ref(time_update),
        std::ref(is_Dirichlet_node), std::ref(*m_mat), std::ref(*m_rhs));
    */

    m_threads[t] =
        std::thread(&facet_assembly_worker_type::assemble_mat_and_rhs_facets_part,
                    &(*m_facet_assembly_workers[t]), std::ref(tria),
                    std::ref(sol_cell_connectivity), std::ref(solution), std::ref(time_update),
                    std::ref(is_Dirichlet_node), std::ref(mat), std::ref(rhs));
  }

  std::for_each(m_threads.begin(), m_threads.end(), std::mem_fn(&std::thread::join));

  // **************************************************************************
  // PROCESS TERMS COMING INTERIOR CELLS
  // **************************************************************************

  for (Uint t = 0; t < m_cell_assembly_workers.size(); ++t)
  {
    /*
    m_threads[t] =
    std::thread(&assembly_worker_type::assemble_mat_and_rhs_cells_part,
                               std::ref(*m_assembly_workers[t]),
    std::ref(geo_cell_connectivity), std::ref(sol_cell_connectivity),
    std::ref(solution), std::ref(time_update), std::ref(is_Dirichlet_node),
    std::ref(*m_mat),
                               std::ref(*m_rhs));
    */

    m_threads[t] = std::thread(
        &cell_assembly_worker_type::assemble_mat_and_rhs_part, &(*m_cell_assembly_workers[t]),
        std::ref(tria), std::ref(sol_cell_connectivity), std::ref(solution), std::ref(time_update),
        std::ref(is_Dirichlet_node), std::ref(mat), std::ref(rhs));
  }

  std::for_each(m_threads.begin(), m_threads.end(), std::mem_fn(&std::thread::join));

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

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::assemble_rhs_impl(
    const mesh::Tria<MeshConfig> &tria,
    typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
    std::vector<bool> const &is_Dirichlet_node, RDMBCManager<MeshConfig, Physics> &bc_manager,
    std::map<std::string, BCScratchDataType *> &weak_bc_data,
    interpolation::VectorMeshFunction<Real> &solution, ls::TpetraMultiVector<Real> &rhs,
    interpolation::VectorMeshFunction<Real> &nodal_residuals, RDTimeUpdate &time_update)
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
  // PROCESS TERMS COMING INTERIOR FACETS
  // **************************************************************************

  for (Uint t = 0; t < m_facet_assembly_workers.size(); ++t)
  {
    /*
    m_threads[t] =
    std::thread(&assembly_worker_type::assemble_rhs_facets_part,
                               std::ref(*m_assembly_workers[t]),
    std::ref(tria), std::ref(geo_cell_connectivity),
    std::ref(sol_cell_connectivity), std::ref(solution),
    std::ref(time_update), std::ref(is_Dirichlet_node), std::ref(*m_rhs));
    */

    m_threads[t] = std::thread(&facet_assembly_worker_type::assemble_rhs_facets_part,
                               &(*m_facet_assembly_workers[t]), std::ref(tria),
                               std::ref(sol_cell_connectivity), std::ref(solution),
                               std::ref(time_update), std::ref(is_Dirichlet_node), std::ref(rhs));
  }

  std::for_each(m_threads.begin(), m_threads.end(), std::mem_fn(&std::thread::join));

  // **************************************************************************
  // PROCESS TERMS COMING INTERIOR CELLS
  // **************************************************************************

  for (Uint t = 0; t < m_cell_assembly_workers.size(); ++t)
  {
    /*
    m_threads[t] = std::thread(
        &assembly_worker_type::assemble_rhs_cells_part,
    std::ref(*m_assembly_workers[t]), std::ref(geo_cell_connectivity),
    std::ref(sol_cell_connectivity), std::ref(solution),
    std::ref(time_update), std::ref(is_Dirichlet_node), std::ref(*m_rhs));
    */

    m_threads[t] =
        std::thread(&cell_assembly_worker_type::assemble_rhs_part, &(*m_cell_assembly_workers[t]),
                    std::ref(tria), std::ref(sol_cell_connectivity), std::ref(solution),
                    std::ref(time_update), std::ref(is_Dirichlet_node), std::ref(rhs));
  }

  std::for_each(m_threads.begin(), m_threads.end(), std::mem_fn(&std::thread::join));

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

  // **************************************************************************
  // SOLVE THE LINEAR SYSTEM
  // **************************************************************************

  /*
  m_lin_system.configure(m_mat, m_rhs, m_du);
  m_lin_system.solve();

  typedef typename interpolation::VectorMeshFunction<Real>::entry_type
  node_value_type;

  for (Uint n = 0; n < sol_cell_connectivity.nb_nodes(); ++n)
  {
    node_value_type u_node = solution.value(n);

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      u_node[eq] += (*m_du)(n *Physics::NEQ + eq, 0);
    }
  }
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::solve(
    typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
    std::shared_ptr<ls::TpetraCrsMatrix<Real>> &mat,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &rhs,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &du, ls::LSTpetra<Real> &linear_system,
    interpolation::VectorMeshFunction<Real> &solution, const bool recompute_preconditioner)
{
  // --------------------------------------------
  // SOLVE THE LINEAR SYSTEM
  // --------------------------------------------

  const Real omega_relax = 1.0;

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

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::
    compute_residual_norm(typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity,
                          ls::TpetraMultiVector<Real> const &rhs,
                          math::DenseDVec<Real> &norm_L2) const
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

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme>::
    create_system_matrix_structure(const mesh::Tria<MeshConfig> &tria,
                                   typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                                   ls::TpetraCrsMatrix<Real> &mat, ls::TpetraMultiVector<Real> &rhs,
                                   ls::TpetraMultiVector<Real> &du)
{
  ls::LSAssembler assembler;
  graph::Graph<Int> ls_dofs(sol_dofs.nb_nodes() * Physics::NEQ);
  assembler.build_matrix_sparsity_pattern_discontinuous<MeshConfig>(Physics::NEQ, tria, sol_dofs,
                                                                    ls_dofs);

  // ls_dofs.write_to_file("matrix_graph_drd_implicit.dat");

  // using cell_dofs_type = typename result_of::dof_map_t<MeshConfig>;

  std::cout << "  Creating matrix structure in discontinuous implicit solver" << std::endl;

  // m_mat->init(sol_dofs.nb_nodes() * Physics::NEQ);
  mat.init(ls_dofs);
  du.init(mat.map(), 1);
  rhs.init(mat.map(), 1);

  // m_mat->unlock();
  // m_mat->fill(0.0);
  rhs.fill(0.0);

  std::cout << "  Matrix structure created" << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
