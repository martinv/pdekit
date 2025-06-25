#ifndef PDEKIT_DG_Hyperbolic_Solver_hpp
#define PDEKIT_DG_Hyperbolic_Solver_hpp

#include <array>

#include "interpolation/FluxSpaceMetric.hpp"
#include "mesh/Tria.hpp"
#include "mesh/containers/DofMap.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "solver/SolverBase.hpp"
#include "solver/fe/DGSolver.hpp"
#include "solver/fe/assembly/DGExplicitCellWorker.hpp"
#include "solver/fe/assembly/DGExplicitFacetWorker.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
class DGHyperbolicSolver : public DGSolver<MeshConfig, Physics>
{
  public:
  /// Constructor
  DGHyperbolicSolver();

  /// Destructor
  ~DGHyperbolicSolver();

  /// Configure underlying function space for geometry and solution
  void build_approximation_spaces(const SFunc sf_type, const PointSetID quad_type,
                                  const Uint quadrature_order) override;

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  void compute_node_reordering(typename DGSolver<MeshConfig, Physics>::mesh_type const &topology,
                               typename DGSolver<MeshConfig, Physics>::cell_dofs const &dofs,
                               std::vector<Int> &reordering) override;

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  void compute_cell_reordering(typename DGSolver<MeshConfig, Physics>::mesh_type const &topology,
                               typename DGSolver<MeshConfig, Physics>::cell_dofs const &dofs,
                               std::vector<Int> &reordering) override;

  /// Prepare the method - allocate data etc.
  void initialize_work_data(const Uint nb_threads, const Uint nb_blocks) override;

  /// Assemble both lhs and rhs simultaneously
  void assemble_lhs_and_rhs(const Real CFL);

  /// Assemble only rhs (assuming that left-hand side) is already
  /// assembled
  void assemble_rhs();

  /// Solve the underlying system
  void solve(std::vector<SolverOption> const &solver_options = {});

  private:
  /// TYPEDEFS
  using dg_method_base = DGSolver<MeshConfig, Physics>;

  /// Scheme for discretization in cell interiors
  using cell_loop_implementation_type =
      detail::DGExplicitCellWorker<MeshConfig, Physics, SchemeTraits>;

  /// Scheme for discretization on facets
  using facet_loop_implementation_type =
      detail::DGExplicitFacetWorker<MeshConfig, Physics, SchemeTraits>;

  using f_space_facets = interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1>;

  typename f_space_facets::ptr m_geo_facet_space;
  typename f_space_facets::ptr m_sol_facet_space;

  cell_loop_implementation_type m_cell_loop_impl;
  facet_loop_implementation_type m_facet_loop_impl;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
DGHyperbolicSolver<MeshConfig, Physics, SchemeTraits>::DGHyperbolicSolver()
    : DGSolver<MeshConfig, Physics>(), m_geo_facet_space(nullptr), m_sol_facet_space(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
DGHyperbolicSolver<MeshConfig, Physics, SchemeTraits>::~DGHyperbolicSolver()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void DGHyperbolicSolver<MeshConfig, Physics, SchemeTraits>::build_approximation_spaces(
    const SFunc sf_type, const PointSetID quad_type, const Uint quadrature_order)
{
  dg_method_base::m_sf_type    = sf_type;
  dg_method_base::m_quad_type  = quad_type;
  dg_method_base::m_quad_order = quadrature_order;

  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [quad_type, quadrature_order](const ElemShape shape,
                                                      const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order, quad_type);
  };

  dg_method_base::m_geo_cell_space = std::make_shared<interpolation::FunctionSpace<MeshConfig>>();

  dg_method_base::m_geo_cell_space->set_reference_fe_values(
      (*dg_method_base::m_geo_dofs).as_range(), sf_generator, quad_generator);

  m_geo_facet_space = std::make_shared<f_space_facets>();
  m_geo_facet_space->set_reference_fe_values_on_skeleton(
      *dg_method_base::m_mesh, *dg_method_base::m_geo_dofs, sf_generator, quad_generator);

  dg_method_base::m_sol_cell_space = std::make_shared<interpolation::FunctionSpace<MeshConfig>>();

  dg_method_base::m_sol_cell_space->set_reference_fe_values(
      (*dg_method_base::m_sol_dofs).as_range(), sf_generator, quad_generator);

  m_sol_facet_space = std::make_shared<f_space_facets>();
  m_sol_facet_space->set_reference_fe_values_on_skeleton(
      *dg_method_base::m_mesh, *dg_method_base::m_sol_dofs, sf_generator, quad_generator);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void DGHyperbolicSolver<MeshConfig, Physics, SchemeTraits>::compute_node_reordering(
    typename DGSolver<MeshConfig, Physics>::mesh_type const &topology,
    typename DGSolver<MeshConfig, Physics>::cell_dofs const &dofs, std::vector<Int> &reordering)
{
  // For the moment, do nothing
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void DGHyperbolicSolver<MeshConfig, Physics, SchemeTraits>::compute_cell_reordering(
    typename DGSolver<MeshConfig, Physics>::mesh_type const &topology,
    typename DGSolver<MeshConfig, Physics>::cell_dofs const &dofs, std::vector<Int> &reordering)
{
  // For the moment, do nothing
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void DGHyperbolicSolver<MeshConfig, Physics, SchemeTraits>::initialize_work_data(
    const Uint nb_threads, const Uint nb_blocks)
{
  const typename result_of::dof_map_t<MeshConfig> &geo_cells = *(dg_method_base::m_geo_dofs);

  const Uint total_nb_cells = geo_cells.nb_active_cells();

  const Uint first_cell = 0;
  const Uint last_cell  = total_nb_cells - 1;

  // Configure the actual discretization
  m_cell_loop_impl.configure_cell_spaces(dg_method_base::m_geo_cell_space,
                                         dg_method_base::m_sol_cell_space, first_cell, last_cell,
                                         nb_blocks, dg_method_base::m_sf_type,
                                         dg_method_base::m_quad_type, dg_method_base::m_quad_order);

  m_facet_loop_impl.configure_facet_spaces(m_geo_facet_space, m_sol_facet_space,
                                           dg_method_base::m_sf_type, dg_method_base::m_quad_type,
                                           dg_method_base::m_quad_order);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void DGHyperbolicSolver<MeshConfig, Physics, SchemeTraits>::assemble_lhs_and_rhs(const Real CFL)
{
  assemble_rhs();

  DGTimeUpdate &time_update = dg_method_base::m_dg_time_update;
  time_update.compute_local_time_step(CFL, dg_method_base::m_timestep);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void DGHyperbolicSolver<MeshConfig, Physics, SchemeTraits>::assemble_rhs()
{
  const mesh::Tria<MeshConfig> &cell_topology               = *(dg_method_base::m_mesh);
  const typename result_of::dof_map_t<MeshConfig> &geo_dofs = *(dg_method_base::m_geo_dofs);
  const typename result_of::dof_map_t<MeshConfig> &sol_dofs = *(dg_method_base::m_sol_dofs);

  const interpolation::VectorMeshFunction<Real> &solution  = (*dg_method_base::m_solution);
  interpolation::VectorMeshFunction<Real> &nodal_residuals = (*dg_method_base::m_residuals);
  nodal_residuals.fill(0.0);

  DGTimeUpdate &time_update = dg_method_base::m_dg_time_update;
  time_update.reset_wave_speeds();

  // For the Lax-Friedrchis scheme, it is important that 'iterate_over_facets'
  // is called before 'iterate_over_cells'

  const Uint first_facet_idx = 0;
  const Uint last_facet_idx  = cell_topology.active_skeleton_size(MeshConfig::TDIM - 1) - 1;
  m_facet_loop_impl.iterate_over_facets(cell_topology, sol_dofs, solution, nodal_residuals,
                                        time_update, first_facet_idx, last_facet_idx);

  m_cell_loop_impl.assemble_lhs_and_rhs_impl(cell_topology, geo_dofs, sol_dofs, solution,
                                             nodal_residuals, time_update);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void DGHyperbolicSolver<MeshConfig, Physics, SchemeTraits>::solve(
    std::vector<SolverOption> const &solver_options)
{
  interpolation::VectorMeshFunction<Real> &u         = (*dg_method_base::m_solution);
  const interpolation::VectorMeshFunction<Real> &res = (*dg_method_base::m_residuals);

  u -= (dg_method_base::m_timestep * res);
}

// ----------------------------------------------------------------------------

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
