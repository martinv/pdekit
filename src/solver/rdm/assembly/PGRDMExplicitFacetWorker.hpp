#ifndef PDEKIT_PG_RDM_Explicit_Facet_Worker_hpp
#define PDEKIT_PG_RDM_Explicit_Facet_Worker_hpp

#include "interpolation/FluxSpaceMetric.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryCache.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/SolutionCache.hpp"
#include "interpolation/SolutionSpaceMetric.hpp"
#include "mesh/point_set/QuadraturePermutation.hpp"

#include "solver/rdm/RDMethodConstData.hpp"
#include "solver/rdm/RDMethodScratchData.hpp"
#include "solver/rdm/RDTimeUpdate.hpp"
#include "solver/rdm/facetsplitters/FacetSchemeSelector.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

namespace detail
{

template <typename MeshConfig, typename Physics, typename FacetScheme>
class PGRDMExplicitFacetWorker
{
  private:
  using facet_scheme_type = typename internal::FacetSchemeSelector<Physics, FacetScheme>::type;

  using facet_method_data_type = typename facet_scheme_type::method_data;

  public:
  using f_space_cells  = interpolation::FunctionSpace<MeshConfig>;
  using f_space_facets = interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1>;

  /// Default constructor
  PGRDMExplicitFacetWorker();

  /// Default destructor
  ~PGRDMExplicitFacetWorker();

  void configure_cell_contour_spaces(const typename f_space_cells::ptr &geo_cell_contour_space,
                                     const typename f_space_cells::ptr &sol_cell_contour_space,
                                     const SFunc sf_type, const PointSetID quad_type,
                                     const Uint quadrature_order);

  void configure_facet_spaces(const mesh::Tria<MeshConfig> &cell_topology,
                              typename result_of::dof_map_t<MeshConfig> const &sol_cells,
                              const SFunc sf_type, const PointSetID quad_type,
                              const Uint quadrature_order);

  void iterate_over_facets(const mesh::Tria<MeshConfig> &cell_topology,
                           const typename result_of::dof_map_t<MeshConfig> &sol_dofs,
                           const interpolation::VectorMeshFunction<Real> &solution,
                           interpolation::VectorMeshFunction<Real> &nodal_residuals,
                           RDTimeUpdate &time_update, const Uint first_cell_idx,
                           const Uint last_cell_idx);

  private:
  enum
  {
    FACET_DIM           = MeshConfig::TDIM - _1D,
    EDGE_DIM            = _1D,
    FACET_DATA_TOPO_DIM = MeshConfig::TDIM - _1D
  };

  enum
  {
    NEQ = Physics::NEQ
  };

  /// Geometry cache and metric for facets
  using facet_geo_cache_type = interpolation::GeometryCache<MeshConfig::GDIM>;
  using facet_geo_metric_type =
      interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, FACET_DIM>;

  /// Solution cache and metric for facets
  using facet_sol_cache_type  = interpolation::SolutionCache;
  using facet_sol_metric_type = interpolation::SolutionSpaceMetric<MeshConfig, FACET_DIM>;

  /// Metric for fluxes on facets
  using facet_flux_metric_type = interpolation::FluxSpaceMetric<MeshConfig, Physics, FACET_DIM>;

  using scal_function_ptr = interpolation::ScalarMeshFunction<Real>::ptr;
  using vect_function_ptr = interpolation::VectorMeshFunction<Real>::ptr;

  /// DATA

  /// Pointer to geometry space on cell facets
  typename f_space_facets::ptr m_geo_facet_space;

  /// Pointer to solution space on cell facets
  typename f_space_facets::ptr m_sol_facet_space;

  /// METRIC CONTAINERS FOR DATA ON MESH SKELETON
  /// HERE WE OPERATE ON ENTITIES WITH codim = 1

  /// Geometry cache for facets, left side of cell interface
  facet_geo_cache_type m_geo_cache_facets_L;

  /// Geometry metric for facets, left side of cell interface
  facet_geo_metric_type m_geo_metric_facets_L;

  /// Geometry cache for facets, right side of cell interface
  facet_geo_cache_type m_geo_cache_facets_R;

  /// Geometry metric for facets, right side of cell interface
  facet_geo_metric_type m_geo_metric_facets_R;

  /// Solution cache for facets
  facet_sol_cache_type m_sol_cache_facets_L;
  facet_sol_cache_type m_sol_cache_facets_R;

  /// Interpolated values and derivatives of the solution u_h
  /// on the mesh skeleton (facets)
  facet_sol_metric_type m_sol_metric_facets_L;
  facet_sol_metric_type m_sol_metric_facets_R;

  /// Interpolated values of flux on the mesh facets
  facet_flux_metric_type m_flux_metric_facets_L;
  facet_flux_metric_type m_flux_metric_facets_R;

  /// Number of blocks to fill in the geometry cache
  Uint m_geo_facet_blk_size;

  /// Number of blocks to fill in the solution cache
  Uint m_sol_facet_blk_size;

  /// Number of cells to fill the cell solution space cache
  Uint m_cell_blk_size;

  /// Type of shape functions
  SFunc m_sf_type;

  /// Quadrature type
  PointSetID m_quad_type;

  /// Quadrature order
  Uint m_quad_order;

  /// A map which associates to each element type in solution space a
  /// corresponding reference element for fluxes

  common::DataMap<mesh::PointSetTagExt, RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>>
      m_method_data_facet_L;
  common::DataMap<mesh::PointSetTagExt, RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>>
      m_method_data_facet_R;

  /// Concrete RDS scheme
  typename internal::FacetSchemeSelector<Physics, FacetScheme>::type m_facet_scheme;
};

// ----------------------------------------------------------------------------
// DISCONTINUOUS RD METHOD - IMPLEMENTATION FOR
// VARIABLE BETA SCHEMES
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename FacetScheme>
PGRDMExplicitFacetWorker<MeshConfig, Physics, FacetScheme>::PGRDMExplicitFacetWorker()
    : m_sol_facet_space(nullptr), m_geo_facet_blk_size(0u), m_sol_facet_blk_size(0u),
      m_cell_blk_size(0u), m_sf_type(SFunc::Undefined), m_quad_type(PointSetID::Undefined),
      m_quad_order(0)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename FacetScheme>
PGRDMExplicitFacetWorker<MeshConfig, Physics, FacetScheme>::~PGRDMExplicitFacetWorker()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename FacetScheme>
void PGRDMExplicitFacetWorker<MeshConfig, Physics, FacetScheme>::configure_cell_contour_spaces(
    const typename f_space_cells::ptr &geo_cell_contour_space,
    const typename f_space_cells::ptr &sol_cell_contour_space, const SFunc sf_type,
    const PointSetID quad_type, const Uint quadrature_order)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename FacetScheme>
void PGRDMExplicitFacetWorker<MeshConfig, Physics, FacetScheme>::configure_facet_spaces(
    const mesh::Tria<MeshConfig> &cell_topology,
    typename result_of::dof_map_t<MeshConfig> const &sol_cells, const SFunc sf_type,
    const PointSetID quad_type, const Uint quadrature_order)
{
  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, sf_type, order, ModalBasis::Modal);
  };

  auto cell_quad_generator = [quadrature_order, quad_type](const ElemShape shape,
                                                           const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order, quad_type);
  };

  // ---------------------
  // GEOMETRY SPACE CONFIG
  // ---------------------
  m_geo_facet_space = std::make_shared<f_space_facets>();
  m_geo_facet_space->set_reference_fe_values_on_skeleton(cell_topology, sf_generator,
                                                         cell_quad_generator);

  m_geo_facet_blk_size = 150; // sol_cells.nb_skeleton_facets();

  m_sf_type    = sf_type;
  m_quad_type  = quad_type;
  m_quad_order = quadrature_order;

  // Allocation of cache and metric for facets
  m_geo_cache_facets_L.allocate(m_geo_facet_space->discrete_elements().cbegin(),
                                m_geo_facet_space->discrete_elements().cend(),
                                m_geo_facet_blk_size);
  m_geo_cache_facets_R.allocate(m_geo_facet_space->discrete_elements().cbegin(),
                                m_geo_facet_space->discrete_elements().cend(),
                                m_geo_facet_blk_size);
  m_geo_metric_facets_L.allocate_buffer(m_geo_facet_space->discrete_elements().cbegin(),
                                        m_geo_facet_space->discrete_elements().cend(),
                                        m_geo_facet_blk_size);
  m_geo_metric_facets_R.allocate_buffer(m_geo_facet_space->discrete_elements().cbegin(),
                                        m_geo_facet_space->discrete_elements().cend(),
                                        m_geo_facet_blk_size);

  // ---------------------
  // SOLUTION SPACE CONFIG
  // ---------------------
  m_sol_facet_space = std::make_shared<f_space_facets>();
  m_sol_facet_space->set_reference_fe_values_on_skeleton(cell_topology, sol_cells, sf_generator,
                                                         cell_quad_generator);

  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> const &facet_elem_type_map =
      m_sol_facet_space->reference_elements();

  m_sol_facet_blk_size = m_geo_facet_blk_size;

  // Cache and metric for solution on facets (mesh skeleton)
  // FIXME: theoretically, this could be remove as the same data might be
  // provided by m_sol_cache_cell_contour_{left,right} and
  // m_sol_metric_cell_contour_{left,right}
  m_sol_cache_facets_L.allocate(m_sol_facet_space->reference_elements().cbegin(),
                                m_sol_facet_space->reference_elements().cend(),
                                m_sol_facet_blk_size, NEQ);
  m_sol_cache_facets_R.allocate(m_sol_facet_space->reference_elements().cbegin(),
                                m_sol_facet_space->reference_elements().cend(),
                                m_sol_facet_blk_size, NEQ);

  const common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> &geo_facet_type_map =
      m_geo_facet_space->reference_elements();

  const typename common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>::const_iterator
      first_ref_elem = geo_facet_type_map.cbegin();

  const mesh::PointSetTag first_geo_facet_tag = first_ref_elem.key_value().std_region_tag();
  const Uint first_geo_facet_order            = first_geo_facet_tag.poly_order();

  m_sol_metric_facets_L.allocate_buffer(m_sol_facet_space->reference_elements().cbegin(),
                                        m_sol_facet_space->reference_elements().cend(),
                                        m_sol_facet_blk_size, NEQ);
  m_sol_metric_facets_R.allocate_buffer(m_sol_facet_space->reference_elements().cbegin(),
                                        m_sol_facet_space->reference_elements().cend(),
                                        m_sol_facet_blk_size, NEQ);

  m_flux_metric_facets_L.allocate_buffer(
      SFunc::Lagrange, first_geo_facet_order, m_sol_facet_space->reference_elements().cbegin(),
      m_sol_facet_space->reference_elements().cend(), m_sol_facet_blk_size);

  m_flux_metric_facets_R.allocate_buffer(
      SFunc::Lagrange, first_geo_facet_order, m_sol_facet_space->reference_elements().cbegin(),
      m_sol_facet_space->reference_elements().cend(), m_sol_facet_blk_size);

  m_method_data_facet_L.clear();
  m_method_data_facet_R.clear();

  // ------------------------------
  // Prepare method data for facets
  // ------------------------------

  for (common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>::const_iterator it =
           facet_elem_type_map.cbegin();
       it != facet_elem_type_map.cend(); ++it)
  {
    common::PtrHandle<RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>> method_data_L =
        m_method_data_facet_L.create(it.key_value());
    const interpolation::FEValues &fev = (*it.data_ptr());
    (*method_data_L).resize_variables(fev);

    common::PtrHandle<RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>> method_data_R =
        m_method_data_facet_R.create(it.key_value());
    (*method_data_R).resize_variables(fev);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename FacetScheme>
void PGRDMExplicitFacetWorker<MeshConfig, Physics, FacetScheme>::iterate_over_facets(
    const mesh::Tria<MeshConfig> &cell_topology,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &solution,
    interpolation::VectorMeshFunction<Real> &nodal_residuals, RDTimeUpdate &time_update,
    const Uint first_cell_idx, const Uint last_cell_idx)
{
  auto sf_generator = [this](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, m_sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [this](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, m_quad_order, m_quad_type);
  };

  using node_value_type = typename interpolation::VectorMeshFunction<Real>::entry_type;
  // using const_node_value_type = typename
  // interpolation::VectorMeshFunction<Real>::const_entry_type;

  // using cell_dofs_type = typename result_of::dof_map_t<MeshConfig>;

  // std::cout << "Total number of facets = " << cell_topology.skeleton_size()
  // << std::endl; std::cout << "Number of active facets = " <<
  // cell_topology.active_skeleton_size() << std::endl;

  const Uint nb_facets        = cell_topology.active_skeleton_size();
  const Uint nb_buffer_blocks = (nb_facets % m_geo_facet_blk_size == 0)
                                    ? nb_facets / m_geo_facet_blk_size
                                    : nb_facets / m_geo_facet_blk_size + 1;

  mesh::EntityDofRealign sol_facet_permutation_L, sol_facet_permutation_R;

  /*
  const typename f_space_facets::quad_generator_fcn &geo_trace_quad_gen =
      m_geo_facet_space->quad_generator();
  */

  const typename f_space_facets::quad_generator_fcn &sol_trace_quad_gen =
      m_sol_facet_space->quad_generator();

  // **********************************************************
  // PHASE 1:
  // **********************************************************
  // Loop over all facets, integrate normal fluxes on facets
  // and add the resulting integrals to element residuals
  // which are in the end obtained as
  // phi = \int_{\partial T} F\cdot n dS
  // by summing the facet contributions
  // At the same time, subtract the facet contributions across
  // the discontinuous interfaces and use them directly to
  // obtain and distribute facet residuals

  const interpolation::ScalarMeshFunction<Real> &nodal_dual_volume =
      time_update.nodal_dual_volume();

  for (Uint bb = 0; bb < nb_buffer_blocks; ++bb)
  {
    const Uint first_facet_in_block = bb * m_geo_facet_blk_size;
    const Uint last_facet_in_block  = std::min((bb + 1) * m_geo_facet_blk_size - 1, nb_facets - 1);

    m_geo_cache_facets_L.flush();
    m_geo_cache_facets_R.flush();

    m_sol_cache_facets_L.flush();
    m_sol_cache_facets_R.flush();

    for (Uint f = first_facet_in_block; f <= last_facet_in_block; ++f)
    {
      // std::cout << "f = " << f << std::endl;

      const mesh::TraceIncidences facet_blk =
          cell_topology.active_skeleton_entry(FACET_DIM, mesh::ActiveIdx(f));

      // Get the permutation sign of the first entity in the block
      // (position 0)
      const mesh::EntityRealignCode pcode_L = facet_blk.permutation(0).get().code();

      // If the entity on position 0 in incidence block needs to be
      // flipped, it is the entity on the rhs of this facet
      const Uint idx_in_block_L = (pcode_L.nb_flips() == 0) ? 0 : 1;

      // If the facet block has size 2, we have 2 entities forming the
      // face (i.e K+ and K-). In that case, idx_R is the 'other' index
      // (choosing from 0 and 1) than idx_L, and can be determined as
      // idx_R = (idx_L + 1) % 2 If the facet block has size 1, then idx_R
      // has to be the same as index left: equal to 0
      const Uint idx_in_block_R = (facet_blk.size() == 2) ? (idx_in_block_L + 1) % 2 : 0;

      // Get the permutation sign of the second entity in the block
      const mesh::EntityRealignCode pcode_R = facet_blk.permutation(idx_in_block_R).get().code();

      // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS
      // OF CELLS, NOT THEIR ACTIVE INDICES !!! IT IS FOR THIS REASON THAT
      // WE FIRST NEED TO RETRIEVE TO TOPOLOGICAL CELLS AND GET THEIR
      // ACTIVE INDICES BEFORE GETTING THE CORRECT CELLS FROM THE DOF
      // HANDLER
      const mesh::CellTopologyView<MeshConfig> tcell_L =
          cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_in_block_L)));
      const mesh::CellTopologyView<MeshConfig> tcell_R =
          cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_in_block_R)));

      // Build discrete element key for left facet (geometry)
      const mesh::PointSetTag geo_facet_reg_L =
          tcell_L.pt_set_id(FACET_DIM, facet_blk.local_id(idx_in_block_L));
      const ElemShape facet_shape_L = geo_facet_reg_L.elem_shape();
      const Uint geo_facet_order_L  = geo_facet_reg_L.poly_order();

      const mesh::PointSetTagExt geo_facet_reg_ext_L(geo_facet_reg_L, P0, pcode_L.adapt_op_id(),
                                                     pcode_L.local_pos_in_parent());
      const mesh::sf::SFTag geo_facet_sf_L     = sf_generator(facet_shape_L, geo_facet_order_L);
      const mesh::PointSetTag geo_facet_quad_L = quad_generator(facet_shape_L, geo_facet_order_L);
      const mesh::PointSetTagExt geo_facet_quad_ext_L(geo_facet_quad_L, P0, pcode_L.adapt_op_id(),
                                                      pcode_L.local_pos_in_parent());

      const mesh::DiscreteElemKey geo_key_L(geo_facet_reg_ext_L, geo_facet_sf_L,
                                            geo_facet_quad_ext_L);

      // Build discrete element key for right facet (geometry)
      const mesh::PointSetTag geo_facet_reg_R =
          tcell_R.pt_set_id(FACET_DIM, facet_blk.local_id(idx_in_block_R));
      const ElemShape facet_shape_R = geo_facet_reg_R.elem_shape();
      const Uint geo_facet_order_R  = geo_facet_reg_R.poly_order();

      const mesh::PointSetTagExt geo_facet_reg_ext_R(geo_facet_reg_R, P0, pcode_R.adapt_op_id(),
                                                     pcode_R.local_pos_in_parent());
      const mesh::sf::SFTag geo_facet_sf_R     = sf_generator(facet_shape_R, geo_facet_order_R);
      const mesh::PointSetTag geo_facet_quad_R = quad_generator(facet_shape_R, geo_facet_order_R);
      const mesh::PointSetTagExt geo_facet_quad_ext_R(geo_facet_quad_R, P0, pcode_R.adapt_op_id(),
                                                      pcode_R.local_pos_in_parent());

      const mesh::DiscreteElemKey geo_key_R(geo_facet_reg_ext_R, geo_facet_sf_R,
                                            geo_facet_quad_ext_R);

      // Get the cell and facet on the left- and right-hand side of
      // interface in SOLUTION space
      const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
      const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

      const mesh::MeshEntity sol_cell_L = sol_dofs.active_cell(mesh::ActiveIdx(active_cell_id_L));
      const mesh::MeshEntity sol_cell_R = sol_dofs.active_cell(mesh::ActiveIdx(active_cell_id_R));

      const mesh::MeshEntity sol_facet_L =
          sol_cell_L.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_L));
      const mesh::MeshEntity sol_facet_R =
          sol_cell_R.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_R));

      const mesh::PointSetTagExt sol_key_L(sol_facet_L.pt_set_id(), P0, pcode_L.adapt_op_id(),
                                           pcode_L.local_pos_in_parent());

      const mesh::PointSetTagExt sol_key_R(sol_facet_R.pt_set_id(), P0, pcode_R.adapt_op_id(),
                                           pcode_R.local_pos_in_parent());

      /*
      std::cout << "===================================================="
      << std::endl; std::cout << "Geo cell left = " << geo_cell_L <<
      std::endl; std::cout << "Geo cell right = " << geo_cell_R <<
      std::endl; std::cout << "Geo facet = " << geo_facet << std::endl;

      std::cout << "Solution cell left = " << sol_cell_L << std::endl;
      std::cout << "Solution facet left = " << sol_facet_L << std::endl;
      std::cout << "Solution cell right = " << sol_cell_R << std::endl;
      std::cout << "Solution facet right = " << sol_facet_R << std::endl;

      std::cout << "Facet has " << facet_blk.size() << " entries" <<
      std::endl; std::cout << "Permutation left = " <<
      facet_blk.permutation(idx_L).get().code()
                << "   " <<
      facet_blk.permutation(idx_L).get().type_id().as_string() <<
      std::endl; std::cout << "Permutation right = " <<
      facet_blk.permutation(idx_R).get().code()
                << "   " <<
      facet_blk.permutation(idx_R).get().type_id().as_string() <<
      std::endl; std::cout << std::endl;
      */

      const mesh::CellGeometry<MeshConfig::GDIM> geo_facet_coords_L =
          tcell_L.coordinates(FACET_DIM, facet_blk.local_id(idx_in_block_L));

      const mesh::CellGeometry<MeshConfig::GDIM> geo_facet_coords_R =
          tcell_R.coordinates(FACET_DIM, facet_blk.local_id(idx_in_block_R));

      m_geo_cache_facets_L.push_back_to_buffer(geo_facet_coords_L, geo_key_L);
      m_geo_cache_facets_R.push_back_to_buffer(geo_facet_coords_R, geo_key_R);

      sol_facet_permutation_L.change_type(sol_facet_L.pt_set_id(),
                                          facet_blk.permutation(idx_in_block_L).get().code());
      sol_facet_permutation_R.change_type(sol_facet_R.pt_set_id(),
                                          facet_blk.permutation(idx_in_block_R).get().code());

      m_sol_cache_facets_L.push_back_to_buffer(sol_facet_L, solution, sol_key_L);
      m_sol_cache_facets_R.push_back_to_buffer(sol_facet_R, solution, sol_key_R);

      /*
      m_sol_cache_facets_R.push_back_to_buffer(
          sol_facet_R, sol_facet_R_permutation, solution,
          mesh::DefaultStdRegMapKey(sol_facet_R.std_region_id(), 0u, P0));
      */

    } // Loop over facets in block

    // Metric data computation for values on mesh skeleton (facets)
    m_geo_metric_facets_L.empty_buffer();
    m_geo_metric_facets_R.empty_buffer();

    m_geo_metric_facets_L.evaluate(m_geo_cache_facets_L, interpolation::RebuildMetricIndex{true});
    m_geo_metric_facets_R.evaluate(m_geo_cache_facets_R, interpolation::RebuildMetricIndex{true});

    m_sol_metric_facets_L.empty_buffer();
    m_sol_metric_facets_R.empty_buffer();

    m_sol_metric_facets_L.evaluate(m_geo_metric_facets_L, m_sol_cache_facets_L,
                                   interpolation::ComputeMetricDerivs{true},
                                   interpolation::RebuildMetricIndex{true});
    m_sol_metric_facets_R.evaluate(m_geo_metric_facets_L, m_sol_cache_facets_R,
                                   interpolation::ComputeMetricDerivs{true},
                                   interpolation::RebuildMetricIndex{true});

    m_flux_metric_facets_L.empty_buffer();
    m_flux_metric_facets_R.empty_buffer();

    m_flux_metric_facets_L.evaluate(m_geo_cache_facets_L, m_geo_metric_facets_L,
                                    m_sol_cache_facets_L, m_sol_metric_facets_L,
                                    interpolation::RebuildMetricIndex{true});
    m_flux_metric_facets_R.evaluate(m_geo_cache_facets_R, m_geo_metric_facets_R,
                                    m_sol_cache_facets_R, m_sol_metric_facets_R,
                                    interpolation::RebuildMetricIndex{true});

    Uint cell_idx_in_metric = 0;

    mesh::QuadraturePermutation quad_p_L, quad_p_R;

    // Now loop over the same block of facets again and
    // - compute the facet residuals
    // - distribute the nodal contributions of facet residuals
    // - compute cell residuals by contour integration (accumulation
    //   of contributions from each cell)

    for (Uint f = first_facet_in_block; f <= last_facet_in_block; ++f)
    {
      const mesh::TraceIncidences facet_blk =
          cell_topology.active_skeleton_entry(FACET_DIM, mesh::ActiveIdx(f));

      // Get the permutation sign of the first entity in the block
      // (position 0)
      const mesh::EntityRealignCode permutation_code = facet_blk.permutation(0).get().code();

      // If the entity on position 0 in incidence block needs to be
      // flipped, it is the entity on the rhs of this facet
      const Uint idx_in_block_L = (permutation_code.nb_flips() == 0) ? 0 : 1;

      // If the facet block has size 2, we have 2 entities forming the
      // face (i.e K+ and K-). In that case, idx_R is the 'other' index
      // (choosing from 0 and 1) than idx_L, and can be determined as
      // idx_R = (idx_L + 1) % 2 If the facet block has size 1, then idx_R
      // has to be the same as index left: equal to 0
      const Uint idx_in_block_R = (facet_blk.size() == 2) ? (idx_in_block_L + 1) % 2 : 0;

      // const mesh::EntityPermutation p_L = facet_blk.permutation(idx_L);
      // const mesh::EntityPermutation p_R = facet_blk.permutation(idx_R);

      // const mesh::MeshEntity geo_cell =
      // geo_cells.active_cell(facet_blk.cell_id(idx_L)); const
      // mesh::MeshEntity geo_facet = geo_cell.sub_entity(FACET_DIM,
      // facet_blk.local_id(idx_L));

      // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS
      // OF CELLS, NOT THEIR ACTIVE INDICES !!!
      const mesh::CellTopologyView<MeshConfig> tcell_L =
          cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_in_block_L)));
      const mesh::CellTopologyView<MeshConfig> tcell_R =
          cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_in_block_R)));

      const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
      const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

      // Get the cell and facet on the left- and right-hand side of
      // interface in SOLUTION space
      const mesh::MeshEntity sol_cell_L = sol_dofs.active_cell(active_cell_id_L);
      const mesh::MeshEntity sol_cell_R = sol_dofs.active_cell(active_cell_id_R);

      const mesh::MeshEntity sol_facet_L =
          sol_cell_L.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_L));
      const mesh::MeshEntity sol_facet_R =
          sol_cell_R.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_R));

      /*
      std::cout << "L [" << active_cell_id_L << "] = ";
      for(Uint v = 0; v < sol_cell_L.nb_vert(); ++v)
      {
        std::cout << sol_cell_L.vertex(v) << " ";
      }
      std::cout << "R = [" << active_cell_id_R << "] = ";
      for(Uint v = 0; v < sol_cell_R.nb_vert(); ++v)
      {
        std::cout << sol_cell_R.vertex(v) << " ";
      }
      std::cout << std::endl;

      current_face.clear();
      for(Uint v = 0; v < sol_facet_R.nb_vert(); ++v)
      {
        current_face.insert(sol_facet_R.vertex(v));
      }

      if ( face_to_find == current_face )
      {
        std::cout << "FOUND PROBLEMATIC FACE" << std::endl;
      }
      */

      /*
      const mesh::DefaultStdRegMapKey
      sol_key_L(sol_facet_L.std_region_id(),
                                                local_sub_tag_L.local_pos_in_parent(),
      P0, local_sub_tag_L.adapt_op_id());

      const mesh::DefaultStdRegMapKey
      sol_key_R(sol_facet_R.std_region_id(),
                                                local_sub_tag_R.local_pos_in_parent(),
      P0, local_sub_tag_R.adapt_op_id());
      */

      const mesh::PointSetTag sol_facet_tag_L = sol_facet_L.pt_set_id();
      const mesh::PointSetTag sol_facet_tag_R = sol_facet_R.pt_set_id();

      sol_facet_permutation_L.change_type(sol_facet_tag_L,
                                          facet_blk.permutation(idx_in_block_L).get().code());
      sol_facet_permutation_R.change_type(sol_facet_tag_R,
                                          facet_blk.permutation(idx_in_block_R).get().code());

      const mesh::PointSetTag contour_quad_tag_L =
          sol_trace_quad_gen(sol_facet_tag_L.elem_shape(), sol_facet_tag_L.poly_order());
      const mesh::PointSetTag contour_quad_tag_R =
          sol_trace_quad_gen(sol_facet_tag_R.elem_shape(), sol_facet_tag_R.poly_order());

      quad_p_L.change_type(contour_quad_tag_L, facet_blk.local_id(idx_in_block_L),
                           facet_blk.permutation(idx_in_block_L).get().code());

      quad_p_R.change_type(contour_quad_tag_R, facet_blk.local_id(idx_in_block_R),
                           facet_blk.permutation(idx_in_block_R).get().code());

      // const math::MatrixBlock<Real> facet_node_coords =
      //    m_geo_cache_facets.cell_values(cell_idx_in_metric);

      const typename facet_geo_metric_type::cellwise_metric geo_facet_met_L =
          m_geo_metric_facets_L.cellwise_values(cell_idx_in_metric);

      const typename facet_geo_metric_type::cellwise_metric geo_facet_met_R =
          m_geo_metric_facets_R.cellwise_values(cell_idx_in_metric);

      const math::DenseConstMatView<Real> sol_facet_nodal_values_L =
          m_sol_cache_facets_L.cell_values(cell_idx_in_metric);

      const math::DenseConstMatView<Real> sol_facet_nodal_values_R =
          m_sol_cache_facets_R.cell_values(cell_idx_in_metric);

      const mesh::PointSetTagExt key_data_L =
          mesh::PointSetTagExt(sol_facet_L.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u);

      const mesh::PointSetTagExt key_data_R =
          mesh::PointSetTagExt(sol_facet_R.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u);

      common::PtrHandle<RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>> facet_method_data_L =
          m_method_data_facet_L.std_region_data(key_data_L);

      common::PtrHandle<RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>> facet_method_data_R =
          m_method_data_facet_R.std_region_data(key_data_R);

      const typename facet_sol_metric_type::cellwise_metric sol_facet_met_L =
          m_sol_metric_facets_L.cellwise_values(cell_idx_in_metric);

      const typename facet_sol_metric_type::cellwise_metric sol_facet_met_R =
          m_sol_metric_facets_R.cellwise_values(cell_idx_in_metric);

      const typename facet_flux_metric_type::cellwise_metric flux_facet_met_L =
          m_flux_metric_facets_L.cellwise_values(cell_idx_in_metric);

      const typename facet_flux_metric_type::cellwise_metric flux_facet_met_R =
          m_flux_metric_facets_R.cellwise_values(cell_idx_in_metric);

      m_facet_scheme.facet_nodal_contributions(
          geo_facet_met_L, geo_facet_met_R, sol_facet_nodal_values_L, sol_facet_nodal_values_R,
          sol_facet_met_L, quad_p_L, sol_facet_met_R, quad_p_R, flux_facet_met_L, flux_facet_met_R,
          *facet_method_data_L, *facet_method_data_R);

      if (idx_in_block_L != idx_in_block_R)
      {
        // ---------------------------------------------------------
        // Accumulate NODAL residuals contributed by the left facet
        // ---------------------------------------------------------
        math::DenseVecView<Real> elem_node_res_L = (*facet_method_data_L).elem_node_res();
        for (Uint n = 0; n < sol_facet_L.nb_vert(); ++n)
        {
          const Uint sol_vertex_id       = sol_facet_L.vertex(n);
          node_value_type nodal_residual = nodal_residuals.value(sol_vertex_id);
          const Real inv_node_volume     = 1. / nodal_dual_volume[sol_vertex_id];

          for (Uint eq = 0; eq < NEQ; ++eq)
          {
            nodal_residual[eq] += inv_node_volume * elem_node_res_L[n * NEQ + eq];
          }
        }

        // ---------------------------------------------------------
        // Accumulate NODAL residuals contributed by the right facet
        // ---------------------------------------------------------
        math::DenseVecView<Real> elem_node_res_R = (*facet_method_data_R).elem_node_res();
        for (Uint n = 0; n < sol_facet_R.nb_vert(); ++n)
        {
          /*
          const Uint sol_vertex_id =
              sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(n));
          */

          const Uint sol_vertex_id       = sol_facet_R.vertex(n);
          node_value_type nodal_residual = nodal_residuals.value(sol_vertex_id);
          const Real inv_node_volume     = 1. / nodal_dual_volume[sol_vertex_id];

          for (Uint eq = 0; eq < NEQ; ++eq)
          {
            nodal_residual[eq] += inv_node_volume * elem_node_res_R[n * NEQ + eq];
          }
        }
      } // If idx_L and idx_R are different (i.e. we are in interior face)

      // --------------------------------------------------------------
      // Accumulate to update coefficient contributed by the left facet
      // --------------------------------------------------------------

      for (Uint n = 0; n < sol_facet_L.nb_vert(); ++n)
      {
        const Real volume_L = nodal_dual_volume[sol_facet_L.vertex(n)];

        time_update.accumulate_nodal_wave_speed(
            sol_facet_L.vertex(n), volume_L * (*facet_method_data_L).m_elem_wave_speed[n]);
      }

      if (idx_in_block_L != idx_in_block_R)
      {

        // ---------------------------------------------------------------
        // Accumulate to update coefficient contributed by the right
        // facet This is done only in case the right facet exists (i.e.
        // is different from the left facet)
        // ---------------------------------------------------------------

        for (Uint n = 0; n < sol_facet_R.nb_vert(); ++n)
        {
          /*
          const Real volume_R = nodal_dual_volume[sol_facet_R.vertex(
              sol_facet_R_permutation.get().vertex(n))];

          time_update.accumulate_nodal_wave_speed(
              sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(n)),
              volume_R * (*facet_method_data_R).m_elem_wave_speed[n]);
          */

          const Real volume_R = nodal_dual_volume[sol_facet_R.vertex(n)];

          time_update.accumulate_nodal_wave_speed(
              sol_facet_R.vertex(n), volume_R * (*facet_method_data_R).m_elem_wave_speed[n]);
        }
      } // If idx_L and idx_R are different (i.e. we are in interior face)

      cell_idx_in_metric++;

    } // Loop over facets in block

  } // Loop over all facet blocks
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
