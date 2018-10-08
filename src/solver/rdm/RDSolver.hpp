#ifndef PDEKIT_PG_RD_Solver_hpp
#define PDEKIT_PG_RD_Solver_hpp

#include <memory>
#include <thread>

#include "common/OptionMap.hpp"
#include "interpolation/FluxSpaceMetric.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "solver/SolverBase.hpp"
#include "solver/rdm/RDTimeUpdate.hpp"
#include "solver/rdm/assembly/PGRDMExplicitCellWorker.hpp"
#include "solver/rdm/bc/RDMBCManager.hpp"
#include "solver/rdm/bc/RDMBCMetricData.hpp"
#include "solver/rdm/cellsplitters/PGLDA.hpp"

namespace pdekit
{

namespace mesh
{
template <typename MeshConfig>
class Mesh;
}

namespace solver
{

namespace rdm
{

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
class RDMBCBase;

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
class RDSolver : public SolverBase
{
  public:
  /// TYPEDEFS
  using mesh_type     = typename mesh::Tria<MeshConfig>;
  using mesh_boundary = typename result_of::mesh_boundary_set_t<MeshConfig>;
  using cell_dofs     = typename result_of::dof_map_t<MeshConfig>;

  using f_space_cells = interpolation::FunctionSpace<MeshConfig>;
  using f_space_bdry  = interpolation::FunctionSpace<MeshConfig, Physics::DIM - 1>;
  using bc_base_type  = RDMBCBase<MeshConfig, Physics, Physics::DIM - 1>;

  /// Default constructor
  RDSolver();

  /// Destructor
  ~RDSolver() override;

  /// Set the data:
  /// a) Topology:
  ///   - cell topology
  ///   - mesh boundaries
  /// b) Geometrical support:
  ///   - geometry dofs
  /// c) Solution support:
  ///   - solution dofs
  void configure_mesh_data(const std::shared_ptr<const mesh_type> &mesh,
                           const common::PtrHandle<const cell_dofs> geo_dofs,
                           const common::PtrHandle<const cell_dofs> sol_dofs);

  /// Configure underlying function space for geometry and solution
  void initialize_work_data(const SFunc sf_type, const PointSetID quad_type,
                            const Uint quadrature_order, const Uint nb_threads,
                            const Uint nb_blocks);

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  void compute_node_reordering(const mesh_type &mesh, const cell_dofs &dofs,
                               std::vector<Int> &reordering);

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  void compute_cell_reordering(const mesh_type &mesh, const cell_dofs &dofs,
                               std::vector<Int> &reordering);

  /// Return pointer to the underlying geometry space
  typename f_space_cells::ptr geo_space() const;

  /// Return pointer to the underlying solution space
  typename f_space_cells::ptr sol_space() const;

  /// Configure a vector function field
  void set_vec_function(const SolverVecFn vec_fn_id,
                        const interpolation::VectorMeshFunction<Real>::ptr &vec_fn) override;

  /// Set the computation of blending coefficient
  void set_blending_coeff(const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff);

  /// Set the smoothness/discontinuity sensor
  void set_artificial_viscosity(
      const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity);

  /// Add a new boundary condition to all levels
  /// @param bc_type_name   ... type of this boundary condition (e.g.
  /// WeakSubInlet)
  /// @param condition_name ... (user defined) name of the boundary condition
  /// @param boundary_name  ... name of the mesh boundary segment where the bc
  /// should be applied
  std::shared_ptr<bc_base_type> add_boundary_condition(const std::string &bc_type_name,
                                                       const std::string &condition_name,
                                                       const std::string &boundary_name);

  /// Return existing boundary condition
  /// @param condition_name ... (user defined) name of the boundary condition
  std::shared_ptr<bc_base_type> boundary_condition(const std::string &condition_name) const;

  /// Remove a boundary condition
  void remove_boundary_condition(const std::string &condition_name);

  /// Remove all boundary conditions
  void remove_all_boundary_conditions();

  /// Compute L1 norm of the residuals (to monitor convergence)
  void compute_residual_norm(math::DenseDVec<Real> &norm_L2) const;

  /// Flag saying whether solver is continuous or not
  bool is_continuous() const override;

  /// Flag saying whether solver is explicit
  bool is_explicit() const override;

  /// Assemble both lhs and rhs simultaneously
  void assemble_lhs_and_rhs(const Real CFL) override;

  /// Assemble only rhs (assuming that left-hand side) is already
  /// assembled
  void assemble_rhs() override;

  /// Solve the underlying system
  void solve(std::vector<SolverOption> const &solver_options = {}) override;

  private:
  /// METHODS
  /// Apply all boundary conditions
  void apply_boundary_conditions();

  /// TYPEDEFS

  using geo_cache_type = interpolation::GeometryCache<MeshConfig::GDIM>;
  using geo_metric_type =
      interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, MeshConfig::GDIM>;

  using sol_cache_type  = interpolation::SolutionCache;
  using sol_metric_type = interpolation::SolutionSpaceMetric<MeshConfig, MeshConfig::GDIM>;

  using flux_metric_type = interpolation::FluxSpaceMetric<MeshConfig, Physics, MeshConfig::GDIM>;

  using bc_metric_data_type = RDMBCMetricData<MeshConfig, Physics, MeshConfig::TDIM - 1>;

  using worker_type = typename detail::PGRDMExplicitCellWorker<MeshConfig, Physics, PGLDA>;

  /// DATA
  /// Hold all solver options
  std::unique_ptr<common::OptionMap> m_options;

  /// Handles to mesh data: support and solution
  std::shared_ptr<const mesh_type> m_tria;
  common::PtrHandle<const cell_dofs> m_geo_dofs;
  common::PtrHandle<const cell_dofs> m_sol_dofs;

  /// Handles to approximation spaces
  typename f_space_cells::ptr m_geo_cell_space;
  typename f_space_cells::ptr m_sol_cell_space;

  /// Type of shape functions
  SFunc m_sf_type;

  /// Quadrature type
  PointSetID m_quad_type;

  /// Quadrature order
  Uint m_quad_order;

  /// Manager to create boundary conditions
  RDMBCManager<MeshConfig, Physics> m_bc_manager;

  std::vector<std::thread> m_threads;
  std::vector<std::unique_ptr<worker_type>> m_workers;

  /// Pointers to the functions holding update coefficient,
  /// blending coefficient and artificial viscosity
  RDTimeUpdate m_rd_time_update;
  interpolation::ScalarMeshFunction<Real>::ptr m_blending_coeff;
  interpolation::ScalarMeshFunction<Real>::ptr m_art_viscosity;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
RDSolver<MeshConfig, Physics>::RDSolver()
    : SolverBase(), m_options(new common::OptionMap()), m_tria(nullptr), m_geo_dofs(nullptr),
      m_sol_dofs(nullptr), m_geo_cell_space(nullptr), m_sol_cell_space(nullptr),
      m_sf_type(SFunc::Undefined), m_quad_type(PointSetID::Undefined), m_quad_order(0)
{
  m_options->create<std::string>("MeshFilePath");
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
RDSolver<MeshConfig, Physics>::~RDSolver()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::configure_mesh_data(
    const std::shared_ptr<const mesh_type> &mesh, const common::PtrHandle<const cell_dofs> geo_dofs,
    const common::PtrHandle<const cell_dofs> sol_dofs)
{
  m_tria     = mesh;
  m_geo_dofs = geo_dofs;
  m_sol_dofs = sol_dofs;

  // Setup the time update object
  m_rd_time_update.setup<MeshConfig>(*mesh, *sol_dofs);

  auto time_step = scal_function(SolverScalFn::time_step);
  time_step->resize((*sol_dofs).nb_nodes());
  time_step->fill(0.0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::initialize_work_data(const SFunc sf_type,
                                                         const PointSetID quad_type,
                                                         const Uint quadrature_order,
                                                         const Uint nb_threads,
                                                         const Uint nb_blocks)
{
  m_sf_type    = sf_type;
  m_quad_type  = quad_type;
  m_quad_order = quadrature_order;

  auto sf_generator_geo = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto sf_generator_sol = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, sf_type, order, ModalBasis::Modal);
  };

  auto cell_quad_generator = [quadrature_order, quad_type](const ElemShape shape,
                                                           const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order, quad_type);
  };

  m_geo_cell_space = std::make_shared<f_space_cells>();
  m_geo_cell_space->set_reference_fe_values((*m_tria).as_active_cell_range(), sf_generator_geo,
                                            cell_quad_generator);
  m_sol_cell_space = std::make_shared<f_space_cells>();
  m_sol_cell_space->set_reference_fe_values((*m_sol_dofs).as_range(), sf_generator_sol,
                                            cell_quad_generator);

  const typename result_of::dof_map_t<MeshConfig> &sol_cells = *m_sol_dofs;

  m_threads.resize(nb_threads);
  m_workers.clear();
  m_workers.resize(0);

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
    std::unique_ptr<worker_type> worker(new worker_type());

    const Uint first_cell = t * nb_cells_per_thread;
    const Uint last_cell  = std::min((t + 1) * nb_cells_per_thread - 1, total_nb_cells - 1);

    // std::cout << "Configuring worker to process cells [" << first_cell <<
    // ","
    // << last_cell << "]"
    //           << std::endl;
    worker->configure_cell_spaces(*m_tria, sol_cells, first_cell, last_cell, nb_blocks, m_sf_type,
                                  m_quad_type, m_quad_order);

    m_workers.push_back(std::move(worker));
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::compute_node_reordering(const mesh_type &mesh,
                                                            const cell_dofs &dofs,
                                                            std::vector<Int> &reordering)
{
  // For the moment, do nothing
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::compute_cell_reordering(const mesh_type &mesh,
                                                            const cell_dofs &dofs,
                                                            std::vector<Int> &reordering)
{
  // For the moment, do nothing
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
typename RDSolver<MeshConfig, Physics>::f_space_cells::ptr RDSolver<MeshConfig,
                                                                    Physics>::geo_space() const
{
  return m_geo_cell_space;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
typename RDSolver<MeshConfig, Physics>::f_space_cells::ptr RDSolver<MeshConfig,
                                                                    Physics>::sol_space() const
{
  return m_sol_cell_space;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::set_vec_function(
    const SolverVecFn vec_fn_id, const interpolation::VectorMeshFunction<Real>::ptr &vec_fn)
{
  SolverBase::set_vec_function(vec_fn_id, vec_fn);

  if (SolverVecFn::solution == vec_fn_id)
  {
    // Make sure that the boundary conditions are updated as well regarding
    // the active solution field
    using bc_iterator = typename RDMBCManager<MeshConfig, Physics>::const_iterator;
    for (bc_iterator bc_it = m_bc_manager.cbegin(); bc_it != m_bc_manager.cend(); ++bc_it)
    {
      bc_it->second->set_solution(vec_fn);
    }
  }
  else if (SolverVecFn::residuals == vec_fn_id)
  {
    // Make sure that the boundary conditions are updated as well regarding
    // the active residual field
    using bc_iterator = typename RDMBCManager<MeshConfig, Physics>::const_iterator;
    for (bc_iterator bc_it = m_bc_manager.cbegin(); bc_it != m_bc_manager.cend(); ++bc_it)
    {
      bc_it->second->set_residuals(vec_fn);
    }
  }
  else if (SolverVecFn::sources == vec_fn_id)
  {
    for (Uint t = 0; t < m_workers.size(); ++t)
    {
      m_workers[t]->set_sources(vec_fn);
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::set_blending_coeff(
    const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff)
{
  m_blending_coeff = blending_coeff;
  for (Uint t = 0; t < m_workers.size(); ++t)
  {
    m_workers[t]->set_blending_coeff(blending_coeff);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::set_artificial_viscosity(
    const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity)
{
  m_art_viscosity = artificial_viscosity;
  for (Uint t = 0; t < m_workers.size(); ++t)
  {
    m_workers[t]->set_artificial_viscosity(artificial_viscosity);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
std::shared_ptr<typename RDSolver<MeshConfig, Physics>::bc_base_type> RDSolver<
    MeshConfig, Physics>::add_boundary_condition(const std::string &bc_type_name,
                                                 const std::string &condition_name,
                                                 const std::string &boundary_name)
{
  std::shared_ptr<bc_base_type> new_bc = m_bc_manager.create_bc(bc_type_name, condition_name);

  new_bc->configure_mesh_data(m_tria, m_geo_dofs, m_sol_dofs, boundary_name);

  new_bc->configure_spaces(m_sf_type, m_quad_type, m_quad_order);
  new_bc->set_solution(vec_function(SolverVecFn::solution));
  new_bc->set_residuals(vec_function(SolverVecFn::residuals));

  return new_bc;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
std::shared_ptr<typename RDSolver<MeshConfig, Physics>::bc_base_type> RDSolver<
    MeshConfig, Physics>::boundary_condition(const std::string &condition_name) const
{
  return m_bc_manager.get_bc(condition_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::remove_boundary_condition(const std::string &condition_name)
{
  m_bc_manager.remove_bc(condition_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::remove_all_boundary_conditions()
{
  m_bc_manager.clear_all_bcs();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::compute_residual_norm(math::DenseDVec<Real> &norm_L2) const
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
bool RDSolver<MeshConfig, Physics>::is_continuous() const
{
  return true;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
bool RDSolver<MeshConfig, Physics>::is_explicit() const
{
  return true;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::assemble_lhs_and_rhs(const Real CFL)
{
  assemble_rhs();

  RDTimeUpdate &time_update = m_rd_time_update;

  time_update.compute_local_time_step(CFL, *scal_function(SolverScalFn::time_step));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::assemble_rhs()
{
  /// @TODO: check that the ptr_handles/shared_ptrs on the RHS of the
  /// assignment are valid
  const mesh::Tria<MeshConfig> &tria = *m_tria;

  const typename result_of::dof_map_t<MeshConfig> &sol_cells = *m_sol_dofs;

  const interpolation::VectorMeshFunction<Real> &solution =
      *SolverBase::vec_function(SolverVecFn::solution);

  interpolation::VectorMeshFunction<Real> &nodal_residuals =
      *SolverBase::vec_function(SolverVecFn::residuals);
  nodal_residuals.fill(0.0);

  RDTimeUpdate &time_update = m_rd_time_update;
  time_update.reset_wave_speeds();

  for (Uint t = 0; t < m_workers.size(); ++t)
  {
    /*
    m_threads[t] =
        std::thread(&worker_type::assemble_lhs_and_rhs_impl,
    std::ref(*m_workers[t]), std::ref(cell_topology), std::ref(sol_cells),
                    std::ref(solution), std::ref(nodal_residuals),
    std::ref(time_update));
    */

    m_threads[t] = std::thread(&worker_type::assemble_lhs_and_rhs_impl, &(*m_workers[t]),
                               std::ref(tria), std::ref(sol_cells), std::ref(solution),
                               std::ref(nodal_residuals), std::ref(time_update));
  }

  std::for_each(m_threads.begin(), m_threads.end(), std::mem_fn(&std::thread::join));

  apply_boundary_conditions();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::solve(std::vector<SolverOption> const &solver_options)
{
  interpolation::VectorMeshFunction<Real> &u         = *vec_function(SolverVecFn::solution);
  const interpolation::VectorMeshFunction<Real> &res = *vec_function(SolverVecFn::residuals);

  u -= ((*scal_function(SolverScalFn::time_step)) * res);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDSolver<MeshConfig, Physics>::apply_boundary_conditions()
{
  using bc_iterator = typename RDMBCManager<MeshConfig, Physics>::const_iterator;
  for (bc_iterator bc_it = m_bc_manager.cbegin(); bc_it != m_bc_manager.cend(); ++bc_it)
  {
    // std::cout << "Applying boundary condition \'" << bc_it->first << "\'"
    // << std::endl;
    RDTimeUpdate &time_update = m_rd_time_update;
    bc_it->second->apply(time_update);
  }
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
