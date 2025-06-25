#ifndef PDEKIT_D_RD_Method_Explicit_hpp
#define PDEKIT_D_RD_Method_Explicit_hpp

#include "solver/ExplicitSolver.hpp"
#include "solver/rdm/DRDMethodLFImpl.hpp"
#include "solver/rdm/RDMethod.hpp"
#include "solver/rdm/RDMethodScratchData.hpp"
#include "solver/rdm/assembly/PGRDMExplicitCellWorker.hpp"
#include "solver/rdm/assembly/PGRDMExplicitFacetWorker.hpp"
#include "solver/rdm/cellsplitters/CellSchemeSelector.hpp"
#include "solver/rdm/facetsplitters/FacetSchemeSelector.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
class DRDMethodExplicit : public RDMethod<MeshConfig, Physics>, protected ExplicitSolver
{
  private:
  using cell_scheme_type  = typename internal::CellSchemeSelector<Physics, CellScheme>::type;
  using facet_scheme_type = typename internal::FacetSchemeSelector<Physics, FacetScheme>::type;

  public:
  /// TYPEDEFS

  using tria_t         = result_of::tria_t<MeshConfig>;
  using dof_map_t      = result_of::dof_map_t<MeshConfig>;
  using f_space_cells  = interpolation::FunctionSpace<MeshConfig>;
  using f_space_facets = interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1>;
  using bc_base_type   = RDMBCBase<MeshConfig, Physics, Physics::DIM - 1>;

  enum
  {
    is_var_beta_type_rdm = cell_scheme_type::is_var_beta_type_rdm
  };
  enum
  {
    needs_volume_derivatives_on_trace = facet_scheme_type::needs_volume_derivatives_on_trace
  };

  /// Constructor
  DRDMethodExplicit();

  /// Destructor
  ~DRDMethodExplicit() override;

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

  /// Compute L1 norm of the residuals (to monitor convergence)
  void compute_residual_norm(math::DenseDVec<Real> &norm_L2) const override;

  private:
  /// METHODS

  /// Apply all boundary conditions
  void apply_boundary_conditions();

  enum
  {
    FACET_DIM = MeshConfig::TDIM - _1D,
    EDGE_DIM  = _1D,
    FACET_DATA_TOPO_DIM =
        needs_volume_derivatives_on_trace ? MeshConfig::TDIM : MeshConfig::TDIM - _1D
  };

  /// TYPEDEFS
  using rd_method_base = RDMethod<MeshConfig, Physics>;

  // Select scheme for discretization in cell interiors
  using cell_loop_implementation_type = typename common::SelectType<
      is_var_beta_type_rdm,
      typename detail::PGRDMExplicitCellWorker<MeshConfig, Physics, CellScheme>,
      typename detail::DRDMethodLFImpl<MeshConfig, Physics, CellScheme, FacetScheme>>::type;

  // Select scheme for discretization on facets
  using facet_loop_implementation_type = typename common::SelectType<
      is_var_beta_type_rdm,
      typename detail::PGRDMExplicitFacetWorker<MeshConfig, Physics, FacetScheme>,
      typename detail::DRDMethodLFImpl<MeshConfig, Physics, CellScheme, FacetScheme>>::type;

  using scal_function_ptr = interpolation::ScalarMeshFunction<Real>::ptr;
  using vect_function_ptr = interpolation::VectorMeshFunction<Real>::ptr;

  /// DATA
  cell_loop_implementation_type m_cell_worker;
  facet_loop_implementation_type m_facet_worker;
};

// ----------------------------------------------------------------------------
// // METHODS OF DISCONTINUOUS RD METHOD - EXPLICIT
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::DRDMethodExplicit()
    : RDMethod<MeshConfig, Physics>()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::~DRDMethodExplicit()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
bool DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::is_continuous() const
{
  return false;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
bool DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::is_explicit() const
{
  return true;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::initialize_work_data(
    const SFunc sf_type, const PointSetID quad_type, const Uint quadrature_order,
    const Uint nb_threads, const Uint nb_blocks)
{
  rd_method_base::initialize_work_data(sf_type, quad_type, quadrature_order, nb_threads, nb_blocks);

  const tria_t &tria        = *(rd_method_base::m_tria);
  const dof_map_t &sol_dofs = *(rd_method_base::m_sol_dofs);

  const Uint total_nb_cells = sol_dofs.nb_active_cells();

  const Uint first_cell = 0;
  const Uint last_cell  = total_nb_cells - 1;

  // Configure the actual discretization
  m_cell_worker.configure_cell_spaces(tria, sol_dofs, first_cell, last_cell, nb_blocks,
                                      rd_method_base::m_sf_type, rd_method_base::m_quad_type,
                                      rd_method_base::m_quad_order);

  m_facet_worker.configure_facet_spaces(*rd_method_base::m_tria, sol_dofs,
                                        rd_method_base::m_sf_type, rd_method_base::m_quad_type,
                                        rd_method_base::m_quad_order);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::compute_node_reordering(
    typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
    typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs, std::vector<Int> &reordering)
{
  // This is explicit solver - do nothing
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::compute_cell_reordering(
    typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
    typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs, std::vector<Int> &reordering)
{
  // This is explicit solver - do nothing
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::set_vec_function(
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
    m_cell_worker.set_sources(vec_fn);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::set_blending_coeff(
    const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff)
{
  rd_method_base::set_blending_coeff(blending_coeff);
  m_cell_worker.set_blending_coeff(blending_coeff);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::set_artificial_viscosity(
    const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity)
{
  rd_method_base::set_artificial_viscosity(artificial_viscosity);
  m_cell_worker.set_artificial_viscosity(artificial_viscosity);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
std::shared_ptr<
    typename DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::bc_base_type>
DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::boundary_condition(
    const std::string &condition_name) const
{
  return rd_method_base::m_bc_manager.get_bc(condition_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::remove_boundary_condition(
    const std::string &condition_name)
{
  rd_method_base::m_bc_manager.remove_bc(condition_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme,
                       FacetScheme>::remove_all_boundary_conditions()
{
  rd_method_base::m_bc_manager.clear_all_bcs();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::assemble_lhs_and_rhs(
    const Real CFL)
{
  assemble_rhs();

  RDTimeUpdate &time_update = rd_method_base::m_rd_time_update;
  time_update.compute_local_time_step(CFL, *rd_method_base::scal_function(SolverScalFn::time_step));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::assemble_rhs()
{
  const mesh::Tria<MeshConfig> &tria = *(rd_method_base::m_tria);
  const dof_map_t &sol_dofs          = *(rd_method_base::m_sol_dofs);

  const interpolation::VectorMeshFunction<Real> &solution =
      (*rd_method_base::vec_function(SolverVecFn::solution));
  interpolation::VectorMeshFunction<Real> &nodal_residuals =
      (*rd_method_base::vec_function(SolverVecFn::residuals));
  nodal_residuals.fill(0.0);

  RDTimeUpdate &time_update = rd_method_base::m_rd_time_update;
  time_update.reset_wave_speeds();

  // For the Lax-Friedrchis scheme, it is important that 'iterate_over_facets'
  // is called before 'iterate_over_cells'

  const Uint first_facet_idx = 0;
  const Uint last_facet_idx  = tria.active_skeleton_size(MeshConfig::TDIM - 1) - 1;
  m_facet_worker.iterate_over_facets(tria, sol_dofs, solution, nodal_residuals, time_update,
                                     first_facet_idx, last_facet_idx);

  m_cell_worker.assemble_lhs_and_rhs_impl(tria, sol_dofs, solution, nodal_residuals, time_update);

  apply_boundary_conditions();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::solve(
    std::vector<SolverOption> const &solver_options)
{
  interpolation::VectorMeshFunction<Real> &u =
      (*rd_method_base::vec_function(SolverVecFn::solution));
  const interpolation::VectorMeshFunction<Real> &res =
      (*rd_method_base::vec_function(SolverVecFn::residuals));

  u -= ((*rd_method_base::scal_function(SolverScalFn::time_step)) * res);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::compute_residual_norm(
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

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodExplicit<MeshConfig, Physics, CellScheme, FacetScheme>::apply_boundary_conditions()
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
