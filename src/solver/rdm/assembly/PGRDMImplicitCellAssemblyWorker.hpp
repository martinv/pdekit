#ifndef PDEKIT_PG_RDM_Implicit_Cell_Assembly_Worker_hpp
#define PDEKIT_PG_RDM_Implicit_Cell_Assembly_Worker_hpp

#include "interpolation/FluxSpaceMetric.hpp"
#include "linear_system/TpetraFwd.hpp"
#include "mesh/KeyCache.hpp"
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

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
class PGRDMImplicitCellAssemblyWorker
{
  public:
  /// TYPEDEFS
  using tria_t    = typename result_of::tria_t<MeshConfig>;
  using dof_map_t = typename result_of::dof_map_t<MeshConfig>;

  using f_space_cells = interpolation::FunctionSpace<MeshConfig>;

  PGRDMImplicitCellAssemblyWorker(const Uint worker_id);

  virtual ~PGRDMImplicitCellAssemblyWorker();

  virtual void configure_cell_spaces(tria_t const &tria, dof_map_t const &sol_dofs,
                                     const Uint first_cell_idx, const Uint last_cell_idx,
                                     const Uint nb_blocks, const SFunc sf_type,
                                     const PointSetID quad_type, const Uint quad_order);

  void set_sources(const interpolation::VectorMeshFunction<Real>::ptr &sources);

  void set_blending_coeff(const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff);

  void set_artificial_viscosity(
      const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity);

  virtual void assemble_mat_and_rhs_part(const tria_t &tria, const dof_map_t &sol_dofs,
                                         const interpolation::VectorMeshFunction<Real> &solution,
                                         RDTimeUpdate &time_update,
                                         const std::vector<bool> &is_Dirichlet_node,
                                         ls::TpetraCrsMatrix<Real> &mat,
                                         ls::TpetraMultiVector<Real> &rhs) = 0;

  virtual void assemble_rhs_part(const tria_t &tria, const dof_map_t &sol_dofs,
                                 const interpolation::VectorMeshFunction<Real> &solution,
                                 RDTimeUpdate &time_update,
                                 const std::vector<bool> &is_Dirichlet_node,
                                 ls::TpetraMultiVector<Real> &rhs) = 0;

  protected:
  enum
  {
    NEQ = Physics::NEQ
  };

  /// TYPEDEFS
  using geo_cache_type  = interpolation::GeometryCache<MeshConfig::GDIM>;
  using geo_metric_type = interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>;

  using sol_cache_type  = interpolation::SolutionCache;
  using sol_metric_type = interpolation::SolutionSpaceMetric<MeshConfig, MeshConfig::GDIM>;

  using flux_metric_type = interpolation::FluxSpaceMetric<MeshConfig, Physics, MeshConfig::GDIM>;

  using scal_function_ptr = interpolation::ScalarMeshFunction<Real>::ptr;
  using vect_function_ptr = interpolation::VectorMeshFunction<Real>::ptr;

  /// Method to fill discrete element key cache
  void fill_key_cache(const tria_t &tria, const dof_map_t &sol_dofs);

  /// Identifier of this worker
  const Uint m_worker_id;

  /// Space for geometry cells
  f_space_cells m_geo_cell_space;

  /// Space for solution cells
  f_space_cells m_sol_cell_space;

  /// Geometry cache
  geo_cache_type m_geo_cache;

  /// Geometry metric
  geo_metric_type m_geo_metric;

  /// Solution cache
  sol_cache_type m_sol_cache;

  /// Cache for computed nodal residuals
  sol_cache_type m_res_cache;

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

  /// Pointer to source field (if there is any)
  vect_function_ptr m_sources;

  /// Blending coefficient for nonlinear schemes
  scal_function_ptr m_blending_coeff;

  /// Pointer to the function that determines discontinuity sensor
  scal_function_ptr m_art_viscosity;

  /// Index of first cell on which we should iterate
  Uint m_first_cell_idx;

  /// Index of last cell on which we should iterate
  Uint m_last_cell_idx;

  /// Number of blocks to fill in the metric
  Uint m_nb_blocks;

  /// Shape function type
  SFunc m_sf_type;

  /// Quadrature type
  PointSetID m_quad_type;

  /// Quadrature oder
  Uint m_quad_order;

  /// Storage for discrete element keys
  mesh::KeyCache<mesh::DiscreteElemKey> m_geo_key_cache;

  /// Buffer for values that should be accumulated into system matrix
  std::vector<std::tuple<Uint, Uint, Real>> m_mat_buffer;

  /// Buffer for values that should be accumulated into RHS vector
  std::vector<std::tuple<Uint, Real>> m_rhs_buffer;

  /// Buffer for values that should be accumulated into 'time update' object
  std::vector<std::tuple<Uint, Real>> m_time_update_buffer;

  /// Constant method data
  RDMethodConstData<typename geo_metric_type::cellwise_metric,
                    typename sol_metric_type::cellwise_metric,
                    typename flux_metric_type::cellwise_metric>
      m_const_method_data;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>::PGRDMImplicitCellAssemblyWorker(
    const Uint worker_id)
    : m_worker_id(worker_id), m_sources(nullptr), m_blending_coeff(nullptr),
      m_art_viscosity(nullptr), m_first_cell_idx(1), m_last_cell_idx(0), m_nb_blocks(1),
      m_sf_type(SFunc::Undefined), m_quad_type(PointSetID::Undefined), m_quad_order(0)
{
  m_mat_buffer.resize(0);
  m_rhs_buffer.resize(0);
  m_time_update_buffer.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>::~PGRDMImplicitCellAssemblyWorker()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>::configure_cell_spaces(
    tria_t const &tria, dof_map_t const &sol_dofs, const Uint first_cell_idx,
    const Uint last_cell_idx, const Uint nb_blocks, const SFunc sf_type, const PointSetID quad_type,
    const Uint quad_order)
{

  m_first_cell_idx = first_cell_idx;
  m_last_cell_idx  = last_cell_idx;
  m_nb_blocks      = nb_blocks;
  m_sf_type        = sf_type;
  m_quad_type      = quad_type;
  m_quad_order     = quad_order;

  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [quad_order, quad_type](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, quad_order, quad_type);
  };

  m_geo_cell_space.set_reference_fe_values(tria.as_active_cell_range(), sf_generator,
                                           quad_generator);
  m_sol_cell_space.set_reference_fe_values(sol_dofs.as_range(), sf_generator, quad_generator);

  m_geo_cache.allocate(m_geo_cell_space.discrete_elements().cbegin(),
                       m_geo_cell_space.discrete_elements().cend(), m_nb_blocks);
  m_geo_metric.allocate_buffer(m_geo_cell_space.discrete_elements().cbegin(),
                               m_geo_cell_space.discrete_elements().cend(), m_nb_blocks);

  // std::set<mesh::PointSetTag> elem_types;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> const &elem_type_map =
      m_sol_cell_space.reference_elements();

  // CHECK THAT m_nb_blocks == sol_cells.nb_cells()

  m_sol_cache.allocate(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, NEQ);
  m_res_cache.allocate(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, NEQ);

  m_sol_metric.allocate_buffer(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, NEQ);

  const mesh::PointSetTag first_geo_cell_tag = tria.begin_cells_active()->pt_set_id();
  const Uint first_geo_cell_order            = first_geo_cell_tag.poly_order();

  m_flux_metric.allocate_buffer(SFunc::Lagrange, first_geo_cell_order, elem_type_map.cbegin(),
                                elem_type_map.cend(), m_nb_blocks);

  fill_key_cache(tria, sol_dofs);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>::set_sources(
    const interpolation::VectorMeshFunction<Real>::ptr &sources)
{
  /*
  m_source_cache.allocate(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, NEQ);
  m_source_metric.allocate_buffer(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, NEQ);
  */

  const typename f_space_cells::ref_elem_map &ref_elems = m_sol_cell_space.reference_elements();
  m_source_cache.allocate(ref_elems.cbegin(), ref_elems.cend(), m_nb_blocks, NEQ);
  m_source_metric.allocate_buffer(ref_elems.cbegin(), ref_elems.cend(), m_nb_blocks, NEQ);

  m_sources = sources;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>::set_blending_coeff(
    const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff)
{
  m_blending_coeff = blending_coeff;

  /*
  m_stab_coeff_cache.allocate(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, 1u);
  m_stab_coeff_metric.allocate_buffer(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks,
                                      1u);
  */

  const typename f_space_cells::ref_elem_map &ref_elems = m_sol_cell_space.reference_elements();
  m_stab_coeff_cache.allocate(ref_elems.cbegin(), ref_elems.cend(), m_nb_blocks, 1u);
  m_stab_coeff_metric.allocate_buffer(ref_elems.cbegin(), ref_elems.cend(), m_nb_blocks, NEQ);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>::set_artificial_viscosity(
    const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity)
{
  m_art_viscosity = artificial_viscosity;

  /*
  m_stab_coeff_cache.allocate(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks, 1u);
  m_stab_coeff_metric.allocate_buffer(elem_type_map.cbegin(), elem_type_map.cend(), m_nb_blocks,
                                      1u);
  */

  const typename f_space_cells::ref_elem_map &ref_elems = m_sol_cell_space.reference_elements();
  m_stab_coeff_cache.allocate(ref_elems.cbegin(), ref_elems.cend(), m_nb_blocks, 1u);
  m_stab_coeff_metric.allocate_buffer(ref_elems.cbegin(), ref_elems.cend(), m_nb_blocks, 1u);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>::assemble_mat_and_rhs_part(
    const tria_t &tria, const dof_map_t &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &solution, RDTimeUpdate &time_update,
    const std::vector<bool> &is_Dirichlet_node, ls::TpetraCrsMatrix<Real> &mat,
    ls::TpetraMultiVector<Real> &rhs)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>::assemble_rhs_part(
    const tria_t &tria, const dof_map_t &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &solution, RDTimeUpdate &time_update,
    const std::vector<bool> &is_Dirichlet_node, ls::TpetraMultiVector<Real> &rhs)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>::fill_key_cache(const tria_t &tria,
                                                                          const dof_map_t &sol_dofs)
{
  auto sf_generator = [this](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, m_sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [this](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, m_quad_order, m_quad_type);
  };

  m_geo_key_cache.clear();

  Uint cell_idx_in_metric = 0;

  for (const typename dof_map_t::const_dof_range_typed &dof_group :
       sol_dofs.all_active_dof_groups())
  {
    typename dof_map_t::const_dof_iterator_typed sol_dof_block_begin = dof_group.begin();
    typename dof_map_t::const_dof_iterator_typed sol_dof_block_end   = dof_group.begin();

    while ((sol_dof_block_end != dof_group.end()) && (cell_idx_in_metric < m_nb_blocks))
    {
      const mesh::MeshEntity sol_cell = sol_dof_block_end->mesh_entity();

      if ((m_first_cell_idx <= sol_cell.idx()) && (sol_cell.idx() <= m_last_cell_idx))
      {
        const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dof_block_end->tcell();
        const mesh::PointSetTag geo_pt_set_tag              = tcell_view.pt_set_id();
        const mesh::PointSetTagExt active_geo_cell_tag_ext(geo_pt_set_tag, P0,
                                                           mesh::CellTransform::NO_TRANS, 0u);
        const mesh::PointSetTagExt active_sol_cell_tag_ext(sol_cell.pt_set_id(), P0,
                                                           mesh::CellTransform::NO_TRANS, 0u);

        const mesh::sf::SFTag geo_sf_tag =
            sf_generator(geo_pt_set_tag.elem_shape(), geo_pt_set_tag.poly_order());
        const mesh::PointSetTag geo_quad_tag =
            quad_generator(geo_pt_set_tag.elem_shape(), geo_pt_set_tag.poly_order());
        const mesh::PointSetTagExt geo_quad_tag_ext(geo_quad_tag, P0, mesh::CellTransform::NO_TRANS,
                                                    0);
        const mesh::DiscreteElemKey geo_key(active_geo_cell_tag_ext, geo_sf_tag, geo_quad_tag_ext);

        m_geo_key_cache.push_back(geo_key);
      }
      cell_idx_in_metric++;
    } // If solution cell index is within range

    ++sol_dof_block_end;
  }
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
