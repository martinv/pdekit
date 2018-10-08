#ifndef PDEKIT_PG_RDM_Explicit_Cell_Worker_hpp
#define PDEKIT_PG_RDM_Explicit_Cell_Worker_hpp

#include <memory>

#include "interpolation/FluxSpaceMetric.hpp"
#include "solver/rdm/KeyGenerators.hpp"
#include "solver/rdm/RDMethodConstData.hpp"
#include "solver/rdm/RDTimeUpdate.hpp"
#include "solver/rdm/blending_coeff/RDBlendingCoeff.hpp"
#include "solver/rdm/cellsplitters/CellSchemeSelector.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

namespace detail
{

template <typename MeshConfig, typename Physics, typename SchemeTraits>
class PGRDMExplicitCellWorker
{
  public:
  using tria_t    = typename result_of::tria_t<MeshConfig>;
  using dof_map_t = typename result_of::dof_map_t<MeshConfig>;
  using method_data_type =
      typename internal::CellSchemeSelector<Physics, SchemeTraits>::type::method_data;
  using f_space_cells = interpolation::FunctionSpace<MeshConfig>;

  /// Default constructor
  PGRDMExplicitCellWorker();

  /// Default destructor
  ~PGRDMExplicitCellWorker();

  /// Set the function space (reference elements)
  void configure_cell_spaces(const tria_t &tria, const dof_map_t &sol_dofs,
                             const Uint first_cell_idx, const Uint last_cell_idx,
                             const Uint nb_blocks, const SFunc sf_type, const PointSetID quad_type,
                             const Uint quadrature_order);

  void set_sources(const interpolation::VectorMeshFunction<Real>::ptr &sources);

  void set_blending_coeff(const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff);

  void set_artificial_viscosity(
      const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity);

  void assemble_lhs_and_rhs_impl(const mesh::Tria<MeshConfig> &tria, const dof_map_t &sol_dofs,
                                 const interpolation::VectorMeshFunction<Real> &solution,
                                 interpolation::VectorMeshFunction<Real> &nodal_residuals,
                                 RDTimeUpdate &time_update);

  void iterate_by_std_region_type(const mesh::Tria<MeshConfig> &tria, const dof_map_t &sol_cells,
                                  const interpolation::VectorMeshFunction<Real> &solution,
                                  interpolation::VectorMeshFunction<Real> &nodal_residuals,
                                  RDTimeUpdate &time_update, const Uint first_cell_idx,
                                  const Uint last_cell_idx);

  private:
  enum
  {
    NEQ = Physics::NEQ
  };

  using geo_cache_type = interpolation::GeometryCache<MeshConfig::GDIM>;
  using geo_metric_type =
      interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, MeshConfig::GDIM>;

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

  /// Cache for source terms
  sol_cache_type m_source_cache;

  /// Interpolated values and derivatives of the source field
  sol_metric_type m_source_metric;

  /// Cache for stabilization coefficient (artificial viscosity, blending ...)
  sol_cache_type m_stab_coeff_cache;

  /// Metric for stabilization coefficient
  sol_metric_type m_stab_coeff_metric;

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
  RDMethodConstData<typename geo_metric_type::cellwise_metric,
                    typename sol_metric_type::cellwise_metric,
                    typename flux_metric_type::cellwise_metric>
      m_const_method_data;

  /// A map which associates to each element type in solution space a
  /// corresponding
  /// reference element for fluxes
  common::DataMap<mesh::PointSetTagExt, method_data_type> m_rdm_method_data;

  /// RD scheme (splitter) that will be used to compute nodal residuals
  /// in each cell
  typename internal::CellSchemeSelector<Physics, SchemeTraits>::type m_scheme;

  /// Pointer to the function that determines source terms
  interpolation::VectorMeshFunction<Real>::ptr m_sources;

  /// Blending coefficient for nonlinear schemes
  interpolation::ScalarMeshFunction<Real>::ptr m_blending_coeff;

  /// Pointer to the function that determines discontinuity sensor
  interpolation::ScalarMeshFunction<Real>::ptr m_art_viscosity;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
PGRDMExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::PGRDMExplicitCellWorker()
    : m_first_cell_idx(1), m_last_cell_idx(0), m_nb_blocks(1), m_sf_type(SFunc::Undefined),
      m_quad_type(PointSetID::Undefined), m_quad_order(0), m_sources(nullptr),
      m_blending_coeff(nullptr), m_art_viscosity(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
PGRDMExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::~PGRDMExplicitCellWorker()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::configure_cell_spaces(
    const tria_t &tria, const dof_map_t &sol_dofs, const Uint first_cell_idx,
    const Uint last_cell_idx, const Uint nb_blocks, const SFunc sf_type, const PointSetID quad_type,
    const Uint quadrature_order)
{
  // ---------------------
  // GEOMETRY SPACE CONFIG
  // ---------------------

  m_first_cell_idx = first_cell_idx;
  m_last_cell_idx  = last_cell_idx;
  m_nb_blocks      = nb_blocks;
  m_sf_type        = sf_type;
  m_quad_type      = quad_type;
  m_quad_order     = quadrature_order;

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

  m_geometry_space = std::make_shared<f_space_cells>();

  m_geometry_space->set_reference_fe_values(tria.as_active_cell_range(), sf_generator,
                                            quad_generator);

  m_geo_cache.allocate(m_geometry_space->discrete_elements().cbegin(),
                       m_geometry_space->discrete_elements().cend(), m_nb_blocks);
  m_geo_cache.print_types();

  // NEW
#if 0
  /*
  using geo_iter_t = typename dof_map_t::const_dof_iterator;
  const auto geo_range = common::make_iter_range(geo_dofs.cbegin(), geo_dofs.cend());
  DiscreteElemKeyGenerator<geo_iter_t> geo_key_generator{geo_range, sf_type, quad_type,
  quadrature_order};
  */

  const common::Range1D<Uint> cell_range{m_first_cell_idx, m_last_cell_idx};
  DiscreteElemKeyGenerator2<MeshConfig> geo_key_generator{geo_dofs, cell_range, sf_type, quad_type,
                                                         quadrature_order};

  m_geo_cache.clear();
  m_geo_cache.record_insertion_pattern(geo_key_generator, m_nb_blocks);
  m_geo_cache.print_types();
#endif
  // NEW

  m_geo_metric.allocate_buffer(m_geometry_space->discrete_elements().cbegin(),
                               m_geometry_space->discrete_elements().cend(), m_nb_blocks);

  // ---------------------
  // SOLUTION SPACE CONFIG
  // ---------------------
  m_solution_space = std::make_shared<f_space_cells>();

  m_solution_space->set_reference_fe_values(sol_dofs.as_range(), sf_generator, quad_generator);

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
void PGRDMExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::set_sources(
    const interpolation::VectorMeshFunction<Real>::ptr &sources)
{
  m_sources = sources;

  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> const &elem_type_map =
      m_solution_space->reference_elements();

  // typename result_of::cells<MeshConfig>::type const &sol_cells =
  // m_sol_mesh->topology().cells(); CHECK THAT m_nb_blocks ==
  // sol_cells.nb_active_cells()

  m_source_cache.allocate(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, NEQ);
  m_source_metric.allocate_buffer(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, NEQ);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::set_blending_coeff(
    const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff)
{
  m_blending_coeff = blending_coeff;

  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> const &elem_type_map =
      m_solution_space->reference_elements();

  m_stab_coeff_cache.allocate(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, 1u);
  m_stab_coeff_metric.allocate_buffer(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks,
                                      1u);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::set_artificial_viscosity(
    const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity)
{
  m_art_viscosity = artificial_viscosity;

  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> const &elem_type_map =
      m_solution_space->reference_elements();

  m_stab_coeff_cache.allocate(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, 1u);
  m_stab_coeff_metric.allocate_buffer(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks,
                                      1u);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::assemble_lhs_and_rhs_impl(
    const mesh::Tria<MeshConfig> &tria, const dof_map_t &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &solution,
    interpolation::VectorMeshFunction<Real> &nodal_residuals, RDTimeUpdate &time_update)
{
  if (m_last_cell_idx < m_first_cell_idx)
  {
    return;
  }

  const bool has_sources         = m_sources ? true : false;
  const bool has_blending_coeff  = m_blending_coeff ? true : false;
  const bool has_artif_viscosity = m_art_viscosity ? true : false;

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
    if (has_sources)
    {
      m_source_cache.flush();
    }
    if (has_artif_viscosity)
    {
      m_stab_coeff_cache.flush();
    }

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

      if (has_sources)
      {
        m_source_cache.push_back_to_buffer(active_sol_cell, *m_sources, active_sol_cell_tag_ext);
      }
      if (has_artif_viscosity)
      {
        m_stab_coeff_cache.push_back_to_buffer(active_sol_cell, *m_art_viscosity,
                                               active_sol_cell_tag_ext);
      }
    }

    m_geo_metric.empty_buffer();
    m_sol_metric.empty_buffer();
    m_flux_metric.empty_buffer();

    if (has_sources)
    {
      m_source_metric.empty_buffer();
    }

    if (has_artif_viscosity)
    {
      m_stab_coeff_metric.empty_buffer();
    }

    m_geo_metric.evaluate(m_geo_cache, interpolation::RebuildMetricIndex{true});
    m_sol_metric.evaluate(m_geo_metric, m_sol_cache, interpolation::ComputeMetricDerivs{true},
                          interpolation::RebuildMetricIndex{true});
    m_flux_metric.evaluate(m_geo_cache, m_geo_metric, m_sol_cache, m_sol_metric,
                           interpolation::RebuildMetricIndex{true});

    if (has_sources)
    {
      m_source_metric.evaluate(m_geo_metric, m_source_cache,
                               interpolation::ComputeMetricDerivs{false},
                               interpolation::RebuildMetricIndex{true});
    }

    if (has_artif_viscosity)
    {
      m_stab_coeff_metric.evaluate(m_geo_metric, m_stab_coeff_cache,
                                   interpolation::ComputeMetricDerivs{false},
                                   interpolation::RebuildMetricIndex{true});
    }

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

      (*method_data).m_use_external_theta = false;
      if (has_blending_coeff)
      {
        (*method_data).m_use_external_theta = true;
        // (*method_data).m_blending_coeff.fill((*m_blending_coeff)[active_sol_cell.idx()]);

        for (Uint v = 0; v < active_sol_cell.nb_vert(); ++v)
        {
          (*method_data).m_blending_coeff[v] = (*m_blending_coeff)[active_sol_cell.vertex(v)];
        }
      }

      if (has_artif_viscosity)
      {
        /*
        const Real art_visc = (*m_art_viscosity)[active_sol_cell.idx()];
        (*method_data).m_art_visc.fill(art_visc);
        */

        const typename sol_metric_type::cellwise_metric stab_coeff_cell_met =
            m_stab_coeff_metric.cellwise_values(cell_idx_in_metric);

        const math::DenseConstMatView<Real> stab_coeff_at_qd_pt =
            stab_coeff_cell_met.field_values();

        (*method_data).m_art_visc = stab_coeff_at_qd_pt.col(0);
      }

      if (!has_sources)
      {
        m_scheme.compute_adv_residuals(sol_nodal_values, m_const_method_data, *method_data);
      }
      else
      {
        m_const_method_data.CSrcM = m_source_metric.cellwise_values(cell_idx_in_metric);
        m_scheme.compute_adv_reaction_residuals(sol_nodal_values, m_const_method_data,
                                                *method_data);
      }

      math::DenseVecView<Real> const elem_res        = (*method_data).elem_node_res();
      math::DenseDVec<Real> const &elem_update_coeff = (*method_data).m_elem_wave_speed;

      for (Uint n = 0; n < (*method_data).nb_nodes(); ++n)
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

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMExplicitCellWorker<MeshConfig, Physics, SchemeTraits>::iterate_by_std_region_type(
    const mesh::Tria<MeshConfig> &tria, const dof_map_t &sol_cells,
    const interpolation::VectorMeshFunction<Real> &solution,
    interpolation::VectorMeshFunction<Real> &nodal_residuals, RDTimeUpdate &time_update,
    const Uint first_cell_idx, const Uint last_cell_idx)
{
  const interpolation::ScalarMeshFunction<Real> &nodal_dual_volume =
      time_update.nodal_dual_volume();

  m_geo_cache.flush();
  m_sol_cache.flush();
  if (m_sources)
  {
    m_source_cache.flush();
  }

  for (const typename dof_map_t::const_dof_range_typed &dof_group : sol_cells.all_dof_groups())
  {
    for (typename dof_map_t::const_dof_iterator_typed sol_dof_iter = dof_group.begin();
         sol_dof_iter != dof_group.end(); ++sol_dof_iter)
    {
      const mesh::MeshEntity sol_cell                = *sol_dof_iter;
      const mesh::CellTopologyView<MeshConfig> tcell = tria.active_cell(sol_dof_iter->active_idx());

      const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell.coordinates();

      m_geo_cache.push_back_to_buffer(
          geo_cell_coords,
          mesh::PointSetTagExt(tcell.cell_type(), P0, mesh::CellTransform::NO_TRANS, 0u));

      // std::cout << "Solution cell [" << cell_iter->idx() << "] " <<
      // *cell_iter << std::endl;
      m_sol_cache.push_back_to_buffer(
          sol_cell, solution,
          mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));
    }
    if (m_sources)
    {
      for (typename dof_map_t::const_dof_iterator_typed sol_dof_iter = dof_group.begin();
           sol_dof_iter != dof_group.end(); ++sol_dof_iter)
      {
        const mesh::MeshEntity sol_cell = *sol_dof_iter;
        m_source_cache.push_back_to_buffer(
            sol_cell, *m_sources,
            mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));
      }
    }
  }

  m_geo_metric.empty_buffer();
  m_sol_metric.empty_buffer();
  m_flux_metric.empty_buffer();
  if (m_sources)
  {
    m_source_metric.empty_buffer();
  }

  m_geo_metric.evaluate(m_geo_cache, true);
  m_sol_metric.evaluate(m_geo_metric, m_sol_cache, true, true);
  m_flux_metric.evaluate(m_geo_cache, m_geo_metric, m_sol_cache, m_sol_metric, true);

  if (m_sources)
  {
    m_source_metric.evaluate(m_geo_metric, m_source_cache, false, true);
  }

  using node_value_type = typename interpolation::VectorMeshFunction<Real>::entry_type;
  // using const_node_value_type = typename
  // interpolation::VectorMeshFunction<Real>::const_entry_type;

  Uint cell_idx_in_metric = 0;

  for (const typename dof_map_t::const_dof_range_typed &cell_group : sol_cells.all_dof_groups())
  {

    const mesh::PointSetTag cell_type_id = cell_group.begin()->std_region_id();

    common::PtrHandle<method_data_type> method_data = m_rdm_method_data.std_region_data(
        mesh::PointSetTagExt(cell_type_id, P0, mesh::CellTransform::NO_TRANS, 0u));

    for (typename dof_map_t::const_dof_iterator_typed cell_iter = cell_group.begin();
         cell_iter != cell_group.end(); ++cell_iter)
    {

      const math::DenseConstMatView<Real> sol_nodal_values =
          m_sol_cache.cell_values(cell_idx_in_metric);

      const typename geo_metric_type::cellwise_metric geo_cell_met =
          m_geo_metric.cellwise_values(cell_idx_in_metric);

      const typename sol_metric_type::cellwise_metric sol_cell_met =
          m_sol_metric.cellwise_values(cell_idx_in_metric);

      const typename flux_metric_type::cellwise_metric flux_cell_met =
          m_flux_metric.cellwise_values(cell_idx_in_metric);

      if (m_art_viscosity)
      {
        const Real art_visc       = (*m_art_viscosity)[cell_iter->idx()];
        (*method_data).m_art_visc = art_visc;
      }

      if (!m_sources)
      {
        m_scheme.compute_adv_residuals(sol_nodal_values, geo_cell_met, sol_cell_met, flux_cell_met,
                                       *method_data);
      }
      else
      {
        const typename sol_metric_type::cellwise_metric source_cell_met =
            m_source_metric.cellwise_values(cell_idx_in_metric);
        m_scheme.compute_adv_reaction_residuals(sol_nodal_values, geo_cell_met, sol_cell_met,
                                                flux_cell_met, source_cell_met, *method_data);
      }

      std::vector<typename Physics::FluxV> const &elem_res = (*method_data).m_elem_node_res;

      math::DenseDVec<Real> const &elem_update_coeff = (*method_data).m_elem_wave_speed;

      const mesh::MeshEntity solution_elem = *cell_iter;

      for (Uint n = 0; n < (*method_data).m_nb_nodes; ++n)
      {
        const Uint sol_vertex_id   = solution_elem.vertex(n);
        node_value_type node_res   = nodal_residuals.value(sol_vertex_id);
        const Real inv_node_volume = 1. / nodal_dual_volume[sol_vertex_id];

        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          node_res[eq] += inv_node_volume * elem_res[n][eq];
        }

        time_update.accumulate_nodal_wave_speed(solution_elem.vertex(n), elem_update_coeff[n]);

        /*
        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          (*m_residuals)(eq, solution_elem.vertex(n)) +=
        elem_res[n][eq];
        }
        (*m_update_coeff)(0, solution_elem.vertex(n)) +=
        elem_update_coeff[n];
        */
      }

      cell_idx_in_metric++;

    } // Loop over all cells of one group
  }   // Loop over all cell groups
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
