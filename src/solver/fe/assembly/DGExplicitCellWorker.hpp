#ifndef PDEKIT_Solver_DG_Explicit_Cell_Worker_hpp
#define PDEKIT_Solver_DG_Explicit_Cell_Worker_hpp

#include <memory>

#include "interpolation/FluxSpaceMetric.hpp"
#include "solver/fe/DGMethodConstData.hpp"
#include "solver/fe/DGMethodScratchData.hpp"
#include "solver/fe/DGTimeUpdate.hpp"
#include "solver/fe/cell_terms/DGCellTerm.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

namespace detail
{

template <typename MeshConfig, typename Physics, typename SchemeTraits>
class DGExplicitCellWorker
{
  public:
  using method_data_type = DGMethodScratchData<Physics>;
  using f_space_cells    = interpolation::FunctionSpace<MeshConfig>;

  /// Default constructor
  DGExplicitCellWorker();

  /// Default destructor
  ~DGExplicitCellWorker();

  /// Set the function space (reference elements)
  void configure_cell_spaces(const typename f_space_cells::ptr &geo_cell_space,
                             const typename f_space_cells::ptr &sol_cell_space,
                             const Uint first_cell_idx, const Uint last_cell_idx,
                             const Uint nb_blocks, const SFunc sf_type, const PointSetID quad_type,
                             const Uint quadrature_order);

  void assemble_lhs_and_rhs_impl(const mesh::Tria<MeshConfig> &tria,
                                 const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
                                 const typename result_of::dof_map_t<MeshConfig> &sol_dofs,
                                 const interpolation::VectorMeshFunction<Real> &solution,
                                 interpolation::VectorMeshFunction<Real> &nodal_residuals,
                                 DGTimeUpdate &time_update);

  private:
  enum
  {
    NEQ = Physics::NEQ
  };

  using geo_cache_type  = interpolation::GeometryCache<MeshConfig::GDIM>;
  using geo_metric_type = interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>;

  using sol_cache_type  = interpolation::SolutionCache;
  using sol_metric_type = interpolation::SolutionSpaceMetric<MeshConfig, MeshConfig::GDIM>;

  using flux_metric_type = interpolation::FluxSpaceMetric<MeshConfig, Physics, MeshConfig::GDIM>;

  /// DATA

  /// Pointer to the geometry space
  typename f_space_cells::ptr m_geometry_space;

  /// Pointer to the solution space
  typename f_space_cells::ptr m_solution_space;

  /// Geometry cache
  geo_cache_type m_geo_cache;

  /// Geometry metric
  geo_metric_type m_geo_metric;

  /// Solution cache
  sol_cache_type m_sol_cache;

  /// Interpolated values and derivatives of the solution u_h
  sol_metric_type m_sol_metric;

  /// Metric for the fluxes
  flux_metric_type m_flux_metric;

  /// Index of first cell on which we should iterate
  Uint m_first_cell_idx;

  /// Index of last cell on which we should iterate
  Uint m_last_cell_idx;

  /// Number of blocks to fill in the metric
  Uint m_nb_blocks;

  /// Type of shape functions
  SFunc m_sf_type;

  /// Quadrature type
  PointSetID m_quad_type;

  /// Quadrature order
  Uint m_quad_order;

  /// Constant method data
  DGMethodConstData<typename geo_metric_type::cellwise_metric,
                    typename sol_metric_type::cellwise_metric,
                    typename flux_metric_type::cellwise_metric>
      m_const_method_data;

  /// A map which associates to each element type in solution space a
  /// corresponding
  /// reference element for fluxes
  common::DataMap<mesh::PointSetTagExt, method_data_type> m_rdm_method_data;

  /// Object that will be used to discretize element-interior part of weak
  /// form
  typename internal::DGCellTerm<Physics> m_scheme;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
DGExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::DGExplicitCellWorker()
    : m_first_cell_idx(1), m_last_cell_idx(0), m_nb_blocks(1), m_sf_type(SFunc::Undefined),
      m_quad_type(PointSetID::Undefined), m_quad_order(0)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
DGExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::~DGExplicitCellWorker()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void DGExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::configure_cell_spaces(
    const typename f_space_cells::ptr &geo_space, const typename f_space_cells::ptr &sol_space,
    const Uint first_cell_idx, const Uint last_cell_idx, const Uint nb_blocks, const SFunc sf_type,
    const PointSetID quad_type, const Uint quadrature_order)
{
  // ---------------------
  // GEOMETRY SPACE CONFIG
  // ---------------------

  // typename result_of::dof_map_t<MeshConfig> const &geo_cells =
  // *rd_method_base::m_geo_dofs;

  m_geometry_space = geo_space;
  m_first_cell_idx = first_cell_idx;
  m_last_cell_idx  = last_cell_idx;
  m_nb_blocks      = nb_blocks;
  m_sf_type        = sf_type;
  m_quad_type      = quad_type;
  m_quad_order     = quadrature_order;

  m_geo_cache.allocate(m_geometry_space->discrete_elements().cbegin(),
                       m_geometry_space->discrete_elements().cend(), m_nb_blocks);
  m_geo_metric.allocate_buffer(m_geometry_space->discrete_elements().cbegin(),
                               m_geometry_space->discrete_elements().cend(), m_nb_blocks);

  // ---------------------
  // SOLUTION SPACE CONFIG
  // ---------------------
  m_solution_space = sol_space;

  // std::set<mesh::PointSetTag> elem_types;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> const &elem_type_map =
      m_solution_space->reference_elements();

  // typename result_of::cells<MeshConfig>::type const &sol_cells =
  // m_sol_mesh->topology().cells(); CHECK THAT m_nb_blocks ==
  // sol_cells.nb_active_cells()

  m_sol_cache.allocate(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, NEQ);
  m_sol_metric.allocate_buffer(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, NEQ);

  const common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> &geo_ref_elems =
      m_geometry_space->reference_elements();

  const typename common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>::const_iterator
      first_ref_elem = geo_ref_elems.cbegin();

  const mesh::PointSetTag first_geo_cell_tag = first_ref_elem.key_value().std_region_tag();
  const Uint first_geo_cell_order            = first_geo_cell_tag.poly_order();

  /*
  const mesh::MeshEntity first_geo_cell = (*geo_cells.begin());
  const mesh::PointSetTag first_geo_cell_tag = first_geo_cell.std_region_id();
  const Uint first_geo_cell_order = first_geo_cell_tag.poly_order();
  */

  // ---------------------
  // FLUX SPACE CONFIG
  // ---------------------
  m_flux_metric.allocate_buffer(SFunc::Lagrange, first_geo_cell_order, elem_type_map.cbegin(),
                                elem_type_map.cend(), m_nb_blocks);

  m_rdm_method_data.clear();

  for (common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>::const_iterator it =
           elem_type_map.cbegin();
       it != elem_type_map.cend(); ++it)
  {
    common::PtrHandle<method_data_type> method_data = m_rdm_method_data.create(it.key_value());

    (*method_data).resize_variables(*it.data_ptr());
  }

  return;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void DGExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::assemble_lhs_and_rhs_impl(
    const mesh::Tria<MeshConfig> &tria, const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &solution,
    interpolation::VectorMeshFunction<Real> &nodal_residuals, DGTimeUpdate &time_update)
{
  if (m_last_cell_idx < m_first_cell_idx)
  {
    return;
  }

  const Uint nb_cells_to_process = m_last_cell_idx - m_first_cell_idx + 1;

  const Uint nb_buffer_blocks = (nb_cells_to_process % m_nb_blocks == 0)
                                    ? nb_cells_to_process / m_nb_blocks
                                    : nb_cells_to_process / m_nb_blocks + 1;

  const interpolation::ScalarMeshFunction<Real> &nodal_dual_volume =
      time_update.nodal_dual_volume();

  /*
  auto geo_sf_generator = [this](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, m_sf_type, order, ModalBasis::Modal);
  };
  */
  auto geo_sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto quad_generator = [this](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, this->m_quad_order, this->m_quad_type);
  };

  /*
  std::cout << "Number of cells to process = " << nb_cells_to_process <<
  std::endl; std::cout << "First cell index (global) = " << first_cell_idx <<
  std::endl; std::cout << "Last cell index (global) = " << last_cell_idx <<
  std::endl; std::cout << "Nb. of blocks = " << nb_buffer_blocks << std::endl;
  */

  for (Uint bb = 0; bb < nb_buffer_blocks; ++bb)
  {
    const Uint first_cell_in_block = m_first_cell_idx + bb * m_nb_blocks;
    const Uint last_cell_in_block =
        std::min(m_first_cell_idx + (bb + 1) * m_nb_blocks - 1, m_last_cell_idx);

    /*
    std::cout << "  First cell in block: " << first_cell_in_block <<
    std::endl; std::cout << "  Last cell in block: " << last_cell_in_block
    << std::endl;
    */

    m_geo_cache.flush();
    m_sol_cache.flush();

    for (Uint c = first_cell_in_block; c <= last_cell_in_block; ++c)
    {
      const mesh::CellTopologyView<MeshConfig> tcell = tria.active_cell(mesh::ActiveIdx(c));
      const mesh::MeshEntity active_sol_cell         = sol_dofs.active_cell(mesh::ActiveIdx(c));

      const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell.coordinates();

      const mesh::PointSetTagExt active_geo_cell_tag_ext(tcell.pt_set_id(), P0,
                                                         mesh::CellTransform::NO_TRANS, 0u);

      const ElemShape elem_shape              = tcell.pt_set_id().elem_shape();
      const Uint geo_poly_order               = tcell.pt_set_id().poly_order();
      const mesh::sf::SFTag geo_basis         = geo_sf_generator(elem_shape, geo_poly_order);
      const mesh::PointSetTag geo_eval_pt_set = quad_generator(elem_shape, geo_poly_order);
      const mesh::PointSetTagExt geo_eval_pt_set_ext(geo_eval_pt_set, P0,
                                                     mesh::CellTransform::NO_TRANS, 0u);
      const mesh::DiscreteElemKey geo_key(active_geo_cell_tag_ext, geo_basis, geo_eval_pt_set_ext);

      const mesh::PointSetTagExt active_sol_cell_tag_ext(active_sol_cell.pt_set_id(), P0,
                                                         mesh::CellTransform::NO_TRANS, 0u);

      m_geo_cache.push_back_to_buffer(geo_cell_coords, geo_key);
      m_sol_cache.push_back_to_buffer(active_sol_cell, solution, active_sol_cell_tag_ext);
    }

    m_geo_metric.empty_buffer();
    m_sol_metric.empty_buffer();
    m_flux_metric.empty_buffer();

    m_geo_metric.evaluate(m_geo_cache, interpolation::RebuildMetricIndex{true});
    m_sol_metric.evaluate(m_geo_metric, m_sol_cache, interpolation::ComputeMetricDerivs{true},
                          interpolation::RebuildMetricIndex{true});
    m_flux_metric.evaluate(m_geo_cache, m_geo_metric, m_sol_cache, m_sol_metric,
                           interpolation::RebuildMetricIndex{true});

    using node_value_type = typename interpolation::VectorMeshFunction<Real>::entry_type;
    // using const_node_value_type = typename
    // interpolation::VectorMeshFunction<Real>::const_entry_type;

    Uint cell_idx_in_metric = 0;

    for (Uint c = first_cell_in_block; c <= last_cell_in_block; ++c)
    {

      const mesh::MeshEntity active_sol_cell = sol_dofs.active_cell(mesh::ActiveIdx(c));

      common::PtrHandle<method_data_type> method_data = m_rdm_method_data.std_region_data(
          mesh::PointSetTagExt(active_sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

      const math::DenseConstMatView<Real> sol_nodal_values =
          m_sol_cache.cell_values(cell_idx_in_metric);

      /*
      const typename geo_metric_type::cellwise_metric geo_cell_met =
          m_geo_metric.cellwise_values(cell_idx_in_metric);

      const typename sol_metric_type::cellwise_metric sol_cell_met =
          m_sol_metric.cellwise_values(cell_idx_in_metric);

      const typename flux_metric_type::cellwise_metric flux_cell_met =
          m_flux_metric.cellwise_values(cell_idx_in_metric);
      */

      m_const_method_data.CGM = m_geo_metric.cellwise_values(cell_idx_in_metric);
      m_const_method_data.CSM = m_sol_metric.cellwise_values(cell_idx_in_metric);
      m_const_method_data.CFM = m_flux_metric.cellwise_values(cell_idx_in_metric);

      m_scheme.compute_adv_residuals(sol_nodal_values, m_const_method_data, *method_data);

      math::DenseDVec<Real> const &elem_res          = (*method_data).m_elem_node_res;
      math::DenseDVec<Real> const &elem_update_coeff = (*method_data).m_elem_wave_speed;

      for (Uint n = 0; n < (*method_data).m_nb_nodes; ++n)
      {
        const Uint sol_vertex_id   = active_sol_cell.vertex(n);
        node_value_type node_res   = nodal_residuals.value(sol_vertex_id);
        const Real inv_node_volume = 1. / nodal_dual_volume[sol_vertex_id];

        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          node_res[eq] += inv_node_volume * elem_res[n * NEQ + eq];
        }

        time_update.accumulate_nodal_wave_speed(sol_vertex_id, elem_update_coeff[n]);

        /*
        const Uint sol_vertex_id = active_sol_cell.vertex(n);
        time_update.accumulate_nodal_wave_speed(sol_vertex_id,
        elem_update_coeff[n]);
        */
      }

      cell_idx_in_metric++;

    } // Loop over cells in given range

  } // Loop over buffer blocks
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
