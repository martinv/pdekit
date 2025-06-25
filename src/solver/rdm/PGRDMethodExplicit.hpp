#ifndef PDEKIT_PG_RD_Method_Explicit_hpp
#define PDEKIT_PG_RD_Method_Explicit_hpp

#include <thread>

#include "solver/ExplicitSolver.hpp"
#include "solver/rdm/RDMethod.hpp"
#include "solver/rdm/assembly/PGRDMExplicitCellWorker.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, typename SchemeTraits>
class PGRDMethodExplicit : public RDMethod<MeshConfig, Physics>, protected ExplicitSolver
{

  public:
  /// TYPEDEFS

  using tria_t         = result_of::tria_t<MeshConfig>;
  using dof_map_t      = result_of::dof_map_t<MeshConfig>;
  using f_space_cells  = interpolation::FunctionSpace<MeshConfig>;
  using f_space_facets = interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1>;
  using bc_base_type   = RDMBCBase<MeshConfig, Physics, Physics::DIM - 1>;

  /// Constructor
  PGRDMethodExplicit();

  /// Destructor
  ~PGRDMethodExplicit() override;

  /// Flag saying whether solver is continuous or not
  bool is_continuous() const override;

  /// Flag saying whether solver is explicit
  bool is_explicit() const override;

  /// Prepare the method - allocate data etc.
  void initialize_work_data(const SFunc sf_type, const PointSetID quad_type,
                            const Uint quadrature_order, const Uint nb_threads,
                            const Uint nb_blocks) override;

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  void compute_node_reordering(typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
                               typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs,
                               std::vector<Int> &reordering) override;

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  void compute_cell_reordering(typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
                               typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs,
                               std::vector<Int> &reordering) override;

  /// Configure a vector function field
  void set_vec_function(const SolverVecFn vec_fn_id,
                        const interpolation::VectorMeshFunction<Real>::ptr &vec_fn) override;

  void set_blending_coeff(
      const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff) override;

  void set_artificial_viscosity(
      const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity) override;

  /// Return existing boundary condition
  /// @param condition_name ... (user defined) name of the boundary condition
  std::shared_ptr<bc_base_type> boundary_condition(
      const std::string &condition_name) const override;

  /// Remove a boundary condition with given name
  void remove_boundary_condition(const std::string &condition_name) override;

  /// Remove all boundary conditions
  void remove_all_boundary_conditions() override;

  /// Assemble both lhs and rhs simultaneously
  void assemble_lhs_and_rhs(const Real CFL) override;

  /// Assemble only rhs (assuming that left-hand side) is already
  /// assembled
  void assemble_rhs() override;

  /// Solve the underlying system
  void solve(std::vector<SolverOption> const &solver_options = {}) override;

  /// Compute L2 norm of the residuals (to monitor convergence)
  void compute_residual_norm(math::DenseDVec<Real> &norm_L2) const override;

  private:
  /// TYPES
  /*
  class SFGenerator
  {
  public:
    constexpr SFGenerator(const SFunc sf_type) : m_sf_type(sf_type) {}

    constexpr mesh::sf::SFTag operator()(const ElemShape shape, const Uint
  order) const
    {
      return mesh::sf::SFTag(shape, m_sf_type, order, ModalBasis::Modal);
    }

  private:
    const SFunc m_sf_type;
  };

  class QuadGenerator
  {
  public:
    constexpr QuadGenerator(const PointSetID quad_type, const Uint quad_order)
      : m_quad_type(quad_type), m_quad_order(quad_order)
    {
    }

    constexpr mesh::PointSetTag operator()(const ElemShape shape, const Uint
  elem_order)
    {
      return mesh::PointSetTag(shape, m_quad_order, m_quad_type);
    }

  private:
    const PointSetID m_quad_type;
    const Uint m_quad_order;
  };
  */

  class SFGenerator
  {
public:
    constexpr SFGenerator(const SFunc sf_type) : m_sf_type(sf_type)
    {
    }

    mesh::sf::SFTag operator()(const ElemShape shape, const Uint order) const
    {
      return mesh::sf::SFTag(shape, m_sf_type, order, ModalBasis::Modal);
    }

private:
    const SFunc m_sf_type;
  };

  class QuadGenerator
  {
public:
    constexpr QuadGenerator(const PointSetID quad_type, const Uint quad_order)
        : m_quad_type(quad_type), m_quad_order(quad_order)
    {
    }

    mesh::PointSetTag operator()(const ElemShape shape, const Uint elem_order) const
    {
      return mesh::PointSetTag(shape, m_quad_order, m_quad_type);
    }

private:
    const PointSetID m_quad_type;
    const Uint m_quad_order;
  };

  /// METHODS

  /// Make one iteration, loop over groups of cells of the same type
  void iterate_by_std_region_type(const Uint first_cell_idx = 0, const Uint last_cell_idx = 1e9);

  /// Apply all boundary conditions
  void apply_boundary_conditions();

  /// TYPEDEFS
  using rd_method_base = RDMethod<MeshConfig, Physics>;

  using geo_cache_type  = interpolation::GeometryCache<MeshConfig::GDIM>;
  using geo_metric_type = interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>;

  using sol_cache_type  = interpolation::SolutionCache;
  using sol_metric_type = interpolation::SolutionSpaceMetric<MeshConfig, MeshConfig::GDIM>;

  using flux_metric_type = interpolation::FluxSpaceMetric<MeshConfig, Physics, MeshConfig::GDIM>;

  using bc_metric_data_type = RDMBCMetricData<MeshConfig, Physics, MeshConfig::TDIM - 1>;

  using worker_type = typename detail::PGRDMExplicitCellWorker<MeshConfig, Physics, SchemeTraits>;

  /// DATA
  std::vector<std::thread> m_threads;
  std::vector<std::unique_ptr<worker_type>> m_workers;
};

// ----------------------------------------------------------------------------
// // METHODS OF CONTINUOUS (PG) RD METHOD - EXPLICIT
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::PGRDMethodExplicit()
    : RDMethod<MeshConfig, Physics>()
{
  m_threads.resize(0);
  m_workers.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::~PGRDMethodExplicit()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
bool PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::is_continuous() const
{
  return true;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
bool PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::is_explicit() const
{
  return true;
}
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::initialize_work_data(
    const SFunc sf_type, const PointSetID quad_type, const Uint quadrature_order,
    const Uint nb_threads, const Uint nb_blocks)
{
  rd_method_base::initialize_work_data(sf_type, quad_type, quadrature_order, nb_threads, nb_blocks);

  const tria_t &tria         = *(rd_method_base::m_tria);
  const dof_map_t &sol_cells = *(rd_method_base::m_sol_dofs);

#if 0
  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [quadrature_order, quad_type](const ElemShape shape,
                                                      const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order, quad_type);
  };

  /*
  const SFGenerator sf_generator{sf_type};
  const QuadGenerator quad_generator{quad_type, quadrature_order};
  */

  rd_method_base::m_geo_cell_space = std::make_shared<f_space_cells>();

  rd_method_base::m_geo_cell_space->set_reference_fe_values(
    (*rd_method_base::m_geo_dofs).as_range(), sf_generator, quad_generator);
  rd_method_base::m_sol_cell_space = std::make_shared<f_space_cells>();

  rd_method_base::m_sol_cell_space->set_reference_fe_values(
    (*rd_method_base::m_sol_dofs).as_range(), sf_generator, quad_generator);
#endif

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
    worker->configure_cell_spaces(tria, sol_cells, first_cell, last_cell, nb_blocks,
                                  rd_method_base::m_sf_type, rd_method_base::m_quad_type,
                                  rd_method_base::m_quad_order);

    m_workers.push_back(std::move(worker));
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::compute_node_reordering(
    typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
    typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs, std::vector<Int> &reordering)
{
  // For the moment, do nothing
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::compute_cell_reordering(
    typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
    typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs, std::vector<Int> &reordering)
{
  // For the moment, do nothing
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::set_vec_function(
    const SolverVecFn vec_fn_id, const interpolation::VectorMeshFunction<Real>::ptr &vec_fn)
{
  rd_method_base::set_vec_function(vec_fn_id, vec_fn);

  if (SolverVecFn::solution == vec_fn_id)
  {
    // Make sure that the boundary conditions are updated as well regarding
    // the active solution field
    using bc_iterator = typename RDMBCManager<MeshConfig, Physics>::const_iterator;
    for (bc_iterator bc_it = rd_method_base::m_bc_manager.cbegin();
         bc_it != rd_method_base::m_bc_manager.cend(); ++bc_it)
    {
      bc_it->second->set_solution(vec_fn);
    }
  }
  else if (SolverVecFn::residuals == vec_fn_id)
  {
    // Make sure that the boundary conditions are updated as well regarding
    // the active residual field
    using bc_iterator = typename RDMBCManager<MeshConfig, Physics>::const_iterator;
    for (bc_iterator bc_it = rd_method_base::m_bc_manager.cbegin();
         bc_it != rd_method_base::m_bc_manager.cend(); ++bc_it)
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
template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::set_blending_coeff(
    const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff)
{
  rd_method_base::set_blending_coeff(blending_coeff);

  for (Uint t = 0; t < m_workers.size(); ++t)
  {
    m_workers[t]->set_blending_coeff(blending_coeff);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::set_artificial_viscosity(
    const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity)
{
  rd_method_base::set_artificial_viscosity(artificial_viscosity);

  for (Uint t = 0; t < m_workers.size(); ++t)
  {
    m_workers[t]->set_artificial_viscosity(artificial_viscosity);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
std::shared_ptr<typename PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::bc_base_type>
PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::boundary_condition(
    const std::string &condition_name) const
{
  return rd_method_base::m_bc_manager.get_bc(condition_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::remove_boundary_condition(
    const std::string &condition_name)
{
  rd_method_base::m_bc_manager.remove_bc(condition_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::remove_all_boundary_conditions()
{
  rd_method_base::m_bc_manager.clear_all_bcs();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::assemble_lhs_and_rhs(const Real CFL)
{
  assemble_rhs();

  RDTimeUpdate &time_update = rd_method_base::m_rd_time_update;
  time_update.compute_local_time_step(CFL, *rd_method_base::scal_function(SolverScalFn::time_step));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::assemble_rhs()
{
  /// @TODO: check that the ptr_handles/shared_ptrs on the RHS of the
  /// assignment are valid
  const mesh::Tria<MeshConfig> &tria = *(rd_method_base::m_tria);
  const dof_map_t &sol_dofs          = *(rd_method_base::m_sol_dofs);

  const interpolation::VectorMeshFunction<Real> &solution =
      (*rd_method_base::vec_function(SolverVecFn::solution));

  interpolation::VectorMeshFunction<Real> &nodal_residuals =
      (*rd_method_base::vec_function(SolverVecFn::residuals));
  nodal_residuals.fill(0.0);

  RDTimeUpdate &time_update = rd_method_base::m_rd_time_update;
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
                               std::ref(tria), std::ref(sol_dofs), std::ref(solution),
                               std::ref(nodal_residuals), std::ref(time_update));
  }

  std::for_each(m_threads.begin(), m_threads.end(), std::mem_fn(&std::thread::join));

  apply_boundary_conditions();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::solve(
    std::vector<SolverOption> const &solver_options)
{
  interpolation::VectorMeshFunction<Real> &u =
      (*rd_method_base::vec_function(SolverVecFn::solution));
  const interpolation::VectorMeshFunction<Real> &res =
      (*rd_method_base::vec_function(SolverVecFn::residuals));

  u -= ((*rd_method_base::scal_function(SolverScalFn::time_step)) * res);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::compute_residual_norm(
    math::DenseDVec<Real> &norm_L2) const
{
  const interpolation::VectorMeshFunction<Real> &res =
      (*rd_method_base::vec_function(SolverVecFn::residuals));

  const interpolation::ScalarMeshFunction<Real> &dt =
      (*rd_method_base::scal_function(SolverScalFn::time_step));

  const Uint nb_nodes = res.nb_entries();

  if (norm_L2.size() != Physics::NEQ)
  {
    norm_L2.resize(Physics::NEQ);
  }

  norm_L2.fill(0.0);

  for (Uint n = 0; n < nb_nodes; ++n)
  {
    const interpolation::VectorMeshFunction<Real>::const_entry_type res_in_node =
        res.const_value(n);
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      norm_L2[eq] += (res_in_node[eq] * dt[n]) * (res_in_node[eq] * dt[n]);
    }
  }

  for (Uint eq = 0; eq < Physics::NEQ; ++eq)
  {
    norm_L2[eq] = std::sqrt(norm_L2[eq] / nb_nodes);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::iterate_by_std_region_type(
    const Uint first_cell_idx, const Uint last_cell_idx)
{
  /// @TODO: check that the ptr_handles/shared_ptrs on the RHS of the
  /// assignment are valid
  const mesh::Tria<MeshConfig> &cell_topology = *(rd_method_base::m_tria);
  const dof_map_t &sol_dofs                   = *(rd_method_base::m_sol_dofs);

  const typename f_space_cells::vect_f &solution  = (*rd_method_base::m_solution);
  typename f_space_cells::vect_f &nodal_residuals = (*rd_method_base::m_residuals);

  // typename f_space::scal_f &update_coeff =
  // (*rd_method_base::m_rd_time_update);
  RDTimeUpdate &time_update = (*rd_method_base::m_rd_time_update);

  for (Uint t = 0; t < m_workers.size(); ++t)
  {
    m_workers[t]->iterate_by_std_region_type(cell_topology, sol_dofs, solution, nodal_residuals,
                                             time_update, first_cell_idx, last_cell_idx);
  }

  apply_boundary_conditions();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>::apply_boundary_conditions()
{
  using bc_iterator = typename RDMBCManager<MeshConfig, Physics>::const_iterator;
  for (bc_iterator bc_it = rd_method_base::m_bc_manager.cbegin();
       bc_it != rd_method_base::m_bc_manager.cend(); ++bc_it)
  {
    // std::cout << "Applying boundary condition \'" << bc_it->first << "\'"
    // << std::endl;
    bc_it->second->apply(rd_method_base::m_rd_time_update);
  }
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
