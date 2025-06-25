#ifndef PDEKIT_PG_RDM_Full_Jac_Facet_Assembly_Worker_hpp
#define PDEKIT_PG_RDM_Full_Jac_Facet_Assembly_Worker_hpp

#include "interpolation/FluxSpaceMetric.hpp"
#include "mesh/point_set/QuadraturePermutation.hpp"

#include "linear_system/TpetraFwd.hpp"
#include "mesh/KeyCache.hpp"
#include "solver/rdm/RDMethodScratchData.hpp"
#include "solver/rdm/RDTimeUpdate.hpp"
#include "solver/rdm/blending_coeff/RDBlendingCoeff.hpp"
#include "solver/rdm/cellsplitters/CellSchemeSelector.hpp"
#include "solver/rdm/facetsplitters/FacetSchemeSelector.hpp"

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
class PGRDMFullJacFacetAssemblyWorker
{
  private:
  using cell_scheme_type  = typename internal::CellSchemeSelector<Physics, CellScheme>::type;
  using facet_scheme_type = typename internal::FacetSchemeSelector<Physics, FacetScheme>::type;

  using cell_method_data_type  = typename cell_scheme_type::method_data;
  using facet_method_data_type = typename facet_scheme_type::method_data;

  public:
  /// TYPEDEFS

  using tria_t         = typename result_of::tria_t<MeshConfig>;
  using dof_map_t      = typename result_of::dof_map_t<MeshConfig>;
  using f_space_cells  = interpolation::FunctionSpace<MeshConfig>;
  using f_space_facets = interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1>;

  PGRDMFullJacFacetAssemblyWorker(const Uint worker_id);

  ~PGRDMFullJacFacetAssemblyWorker();

  void configure_cell_spaces(dof_map_t const &sol_cells, const Uint first_cell_idx,
                             const Uint last_cell_idx, const Uint nb_blocks, const SFunc sf_type,
                             const PointSetID quad_type, const Uint quadrature_order);

  void configure_facet_spaces(const tria_t &cell_topology, const dof_map_t &sol_cells,
                              const Uint first_facet_idx, const Uint last_facet_idx,
                              const SFunc sf_type, const PointSetID quad_type,
                              const Uint quadrature_order);

  void assemble_mat_and_rhs_facets_part(const tria_t &cell_topology, const dof_map_t &sol_cells,
                                        const interpolation::VectorMeshFunction<Real> &solution,
                                        RDTimeUpdate &time_update,
                                        const std::vector<bool> &is_Dirichlet_node,
                                        ls::TpetraCrsMatrix<Real> &mat,
                                        ls::TpetraMultiVector<Real> &rhs);

  void assemble_rhs_facets_part(const tria_t &cell_topology, const dof_map_t &sol_cells,
                                const interpolation::VectorMeshFunction<Real> &solution,
                                RDTimeUpdate &time_update,
                                const std::vector<bool> &is_Dirichlet_node,
                                ls::TpetraMultiVector<Real> &rhs);

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

  /// Geometry cache and metric for cells
  using cell_geo_cache_type  = interpolation::GeometryCache<MeshConfig::GDIM>;
  using cell_geo_metric_type = interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>;
  using cell_flux_metric_type =
      interpolation::FluxSpaceMetric<MeshConfig, Physics, MeshConfig::TDIM>;

  /// Solution cache and metric for cells
  using cell_sol_cache_type  = interpolation::SolutionCache;
  using cell_sol_metric_type = interpolation::SolutionSpaceMetric<MeshConfig, MeshConfig::TDIM>;

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

  /// Pre-fill the cache with discrete element keys for geometry
  void fill_key_cache(const tria_t &cell_topology, const dof_map_t &sol_cells);

  /// DATA
  /// Identifier of this worker
  const Uint m_worker_id;

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

  /// Cache for computed nodal residuals (facets)
  facet_sol_cache_type m_res_cache_facets_L;
  facet_sol_cache_type m_res_cache_facets_R;

  /// Interpolated values and derivatives of the solution u_h
  /// on the mesh skeleton (facets)
  facet_sol_metric_type m_sol_metric_facets_L;
  facet_sol_metric_type m_sol_metric_facets_R;

  /// Interpolated values of flux on the mesh facets
  facet_flux_metric_type m_flux_metric_facets_L;
  facet_flux_metric_type m_flux_metric_facets_R;

  /// Index of first cell on which we should iterate
  Uint m_first_cell_idx;

  /// Index of last cell on which we should iterate
  Uint m_last_cell_idx;

  /// Index of first facet on which we should iterate
  Uint m_first_facet_idx;

  /// Index of last facet on which we should iterate
  Uint m_last_facet_idx;

  /// Number of cells to fill the cell solution space cache
  Uint m_nb_blocks;

  /// Type of shape functions
  SFunc m_sf_type;

  /// Quadrature type
  PointSetID m_quad_type;

  /// Quadrature order
  Uint m_quad_order;

  /// Storage for discrete element keys
  mesh::KeyCache<mesh::DiscreteElemKey> m_geo_key_cache_L;
  mesh::KeyCache<mesh::DiscreteElemKey> m_geo_key_cache_R;

  /// Buffer for values that should be accumulated into system matrix
  std::vector<std::tuple<Uint, Uint, Real>> m_mat_buffer;

  /// Buffer for values that should be accumulated into RHS vector
  std::vector<std::tuple<Uint, Real>> m_rhs_buffer;

  /// Buffer for values that should be accumulated into 'time update' object
  std::vector<std::tuple<Uint, Real>> m_time_update_buffer;

  /// A map which associates to each element type in solution space a
  /// corresponding reference element for fluxes
  common::DataMap<mesh::PointSetTagExt, cell_method_data_type> m_method_data_cell;
  common::DataMap<mesh::PointSetTagExt, facet_method_data_type> m_method_data_facet_L;
  common::DataMap<mesh::PointSetTagExt, facet_method_data_type> m_method_data_facet_R;

  /// Concrete RDS scheme
  typename internal::CellSchemeSelector<Physics, CellScheme>::type m_cell_scheme;
  typename internal::FacetSchemeSelector<Physics, FacetScheme>::type m_facet_scheme;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
PGRDMFullJacFacetAssemblyWorker<MeshConfig, Physics, CellScheme,
                                FacetScheme>::PGRDMFullJacFacetAssemblyWorker(const Uint worker_id)
    : m_worker_id(worker_id), m_first_cell_idx(1), m_last_cell_idx(0), m_first_facet_idx(1),
      m_last_facet_idx(0), m_nb_blocks(1), m_sf_type(SFunc::Undefined),
      m_quad_type(PointSetID::Undefined), m_quad_order(0)
{
  m_mat_buffer.resize(0);
  m_rhs_buffer.resize(0);
  m_time_update_buffer.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
PGRDMFullJacFacetAssemblyWorker<MeshConfig, Physics, CellScheme,
                                FacetScheme>::~PGRDMFullJacFacetAssemblyWorker()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void PGRDMFullJacFacetAssemblyWorker<MeshConfig, Physics, CellScheme, FacetScheme>::
    configure_cell_spaces(dof_map_t const &sol_cells, const Uint first_cell_idx,
                          const Uint last_cell_idx, const Uint nb_blocks, const SFunc sf_type,
                          const PointSetID quad_type, const Uint quadrature_order)
{
  m_first_cell_idx = first_cell_idx;
  m_last_cell_idx  = last_cell_idx;

  // Use number of blocks given by m_cell_blk_size for cell cache
  m_nb_blocks = nb_blocks;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void PGRDMFullJacFacetAssemblyWorker<MeshConfig, Physics, CellScheme, FacetScheme>::
    configure_facet_spaces(const tria_t &cell_topology, const dof_map_t &sol_cells,
                           const Uint first_facet_idx, const Uint last_facet_idx,
                           const SFunc sf_type, const PointSetID quad_type,
                           const Uint quadrature_order)
{
  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
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

  // std::cout << "DEBUG {GEO}:" << std::endl;
  // m_geo_facet_space->print_discrete_element_types();

  m_first_facet_idx = first_facet_idx;
  m_last_facet_idx  = last_facet_idx;

  m_sf_type    = sf_type;
  m_quad_type  = quad_type;
  m_quad_order = quadrature_order;

  const Uint facet_block_size = 1;

  // Allocation of cache and metric for facets
  m_geo_cache_facets_L.allocate(m_geo_facet_space->discrete_elements().cbegin(),
                                m_geo_facet_space->discrete_elements().cend(), facet_block_size);
  m_geo_cache_facets_R.allocate(m_geo_facet_space->discrete_elements().cbegin(),
                                m_geo_facet_space->discrete_elements().cend(), facet_block_size);
  m_geo_metric_facets_L.allocate_buffer(m_geo_facet_space->discrete_elements().cbegin(),
                                        m_geo_facet_space->discrete_elements().cend(),
                                        facet_block_size);
  m_geo_metric_facets_R.allocate_buffer(m_geo_facet_space->discrete_elements().cbegin(),
                                        m_geo_facet_space->discrete_elements().cend(),
                                        facet_block_size);

  // ---------------------
  // SOLUTION SPACE CONFIG
  // ---------------------
  m_sol_facet_space = std::make_shared<f_space_facets>();
  m_sol_facet_space->set_reference_fe_values_on_skeleton(cell_topology, sol_cells, sf_generator,
                                                         cell_quad_generator);
  // std::cout << "DEBUG {SOL}:" << std::endl;
  // m_sol_facet_space->print_discrete_element_types();

  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> const &facet_elem_type_map =
      m_sol_facet_space->reference_elements();

  // Cache and metric for solution on facets (mesh skeleton)
  // FIXME: theoretically, this could be remove as the same data might be
  // provided by m_sol_cache_cell_contour_{left,right} and
  // m_sol_metric_cell_contour_{left,right}
  m_sol_cache_facets_L.allocate(m_sol_facet_space->reference_elements().cbegin(),
                                m_sol_facet_space->reference_elements().cend(), facet_block_size,
                                NEQ);
  m_sol_cache_facets_R.allocate(m_sol_facet_space->reference_elements().cbegin(),
                                m_sol_facet_space->reference_elements().cend(), facet_block_size,
                                NEQ);

  m_res_cache_facets_L.allocate(m_sol_facet_space->reference_elements().cbegin(),
                                m_sol_facet_space->reference_elements().cend(), facet_block_size,
                                NEQ);
  m_res_cache_facets_R.allocate(m_sol_facet_space->reference_elements().cbegin(),
                                m_sol_facet_space->reference_elements().cend(), facet_block_size,
                                NEQ);

  m_sol_metric_facets_L.allocate_buffer(m_sol_facet_space->reference_elements().cbegin(),
                                        m_sol_facet_space->reference_elements().cend(),
                                        facet_block_size, NEQ);
  m_sol_metric_facets_R.allocate_buffer(m_sol_facet_space->reference_elements().cbegin(),
                                        m_sol_facet_space->reference_elements().cend(),
                                        facet_block_size, NEQ);

  const common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> &geo_facet_type_map =
      m_geo_facet_space->reference_elements();

  const typename common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>::const_iterator
      first_ref_elem = geo_facet_type_map.cbegin();

  const mesh::PointSetTag first_geo_facet_tag = first_ref_elem.key_value().std_region_tag();
  const Uint first_geo_facet_order            = first_geo_facet_tag.poly_order();

  m_flux_metric_facets_L.allocate_buffer(
      SFunc::Lagrange, first_geo_facet_order, m_sol_facet_space->reference_elements().cbegin(),
      m_sol_facet_space->reference_elements().cend(), facet_block_size);

  m_flux_metric_facets_R.allocate_buffer(
      SFunc::Lagrange, first_geo_facet_order, m_sol_facet_space->reference_elements().cbegin(),
      m_sol_facet_space->reference_elements().cend(), facet_block_size);

  m_method_data_facet_L.clear();
  m_method_data_facet_R.clear();

  // ------------------------------
  // Prepare method data for facets
  // ------------------------------

  for (common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>::const_iterator it =
           facet_elem_type_map.cbegin();
       it != facet_elem_type_map.cend(); ++it)
  {
    common::PtrHandle<facet_method_data_type> method_data_L =
        m_method_data_facet_L.create(it.key_value());
    const interpolation::FEValues &fev = (*it.data_ptr());
    (*method_data_L).resize_variables(fev);

    common::PtrHandle<facet_method_data_type> method_data_R =
        m_method_data_facet_R.create(it.key_value());
    (*method_data_R).resize_variables(fev);
  }

  fill_key_cache(cell_topology, sol_cells);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void PGRDMFullJacFacetAssemblyWorker<MeshConfig, Physics, CellScheme, FacetScheme>::
    assemble_mat_and_rhs_facets_part(const tria_t &cell_topology, const dof_map_t &sol_cells,
                                     const interpolation::VectorMeshFunction<Real> &solution,
                                     RDTimeUpdate &time_update,
                                     const std::vector<bool> &is_Dirichlet_node,
                                     ls::TpetraCrsMatrix<Real> &mat,
                                     ls::TpetraMultiVector<Real> &rhs)
{
  mesh::EntityDofRealign sol_facet_L_permutation, sol_facet_R_permutation;
  mesh::QuadraturePermutation quad_p_L, quad_p_R;

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

  m_mat_buffer.resize(0);
  m_rhs_buffer.resize(0);
  m_time_update_buffer.resize(0);

  for (Uint f = m_first_facet_idx; f <= m_last_facet_idx; ++f)
  {
    // std::cout << "f = " << f << std::endl;

    const mesh::TraceIncidences facet_blk =
        cell_topology.active_skeleton_entry(FACET_DIM, mesh::ActiveIdx(f));

    // Get the permutation sign of the first entity in the block (position
    // 0)
    const mesh::EntityRealignCode pcode_L = facet_blk.permutation(0).get().code();

    // If the entity on position 0 in incidence block needs to be flipped,
    // it is the entity on the rhs of this facet
    const Uint idx_in_block_L = (pcode_L.nb_flips() == 0) ? 0 : 1;

    // If the facet block has size 2, we have 2 entities forming the face
    // (i.e K+ and K-). In that case, idx_R is the 'other' index (choosing
    // from 0 and 1) than idx_L, and can be determined as idx_R = (idx_L +
    // 1) % 2 If the facet block has size 1, then idx_R has to be the same
    // as index left: equal to 0
    const Uint idx_in_block_R = (facet_blk.size() == 2) ? (idx_in_block_L + 1) % 2 : 0;

    // Get the permutation sign of the first entity in the block (position
    // 0)
    const mesh::EntityRealignCode pcode_R = facet_blk.permutation(idx_in_block_R).get().code();

    // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS OF
    // CELLS, NOT THEIR ACTIVE INDICES !!! IT IS FOR THIS REASON THAT WE
    // FIRST NEED TO RETRIEVE TO TOPOLOGICAL CELLS AND GET THEIR ACTIVE
    // INDICES BEFORE GETTING THE CORRECT CELLS FROM THE DOF HANDLER
    const mesh::CellTopologyView<MeshConfig> tcell_L =
        cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_in_block_L)));
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_in_block_R)));

    const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
    const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

    const Uint geo_key_cache_pos_idx      = f - m_first_facet_idx;
    const mesh::DiscreteElemKey geo_key_L = m_geo_key_cache_L.key(geo_key_cache_pos_idx);
    const mesh::DiscreteElemKey geo_key_R = m_geo_key_cache_R.key(geo_key_cache_pos_idx);

    // Get the cell and facet on the left- and right-hand side of interface
    // in SOLUTION space
    const mesh::MeshEntity sol_cell_L = sol_cells.active_cell(mesh::ActiveIdx(active_cell_id_L));
    const mesh::MeshEntity sol_cell_R = sol_cells.active_cell(mesh::ActiveIdx(active_cell_id_R));

    const mesh::MeshEntity sol_facet_L =
        sol_cell_L.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_L));
    const mesh::MeshEntity sol_facet_R =
        sol_cell_R.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_R));

    const mesh::PointSetTagExt sol_key_L(sol_facet_L.pt_set_id(), P0, pcode_L.adapt_op_id(),
                                         pcode_L.local_pos_in_parent());

    const mesh::PointSetTagExt sol_key_R(sol_facet_R.pt_set_id(), P0, pcode_R.adapt_op_id(),
                                         pcode_R.local_pos_in_parent());

    /*
    std::cout << "====================================================" <<
    std::endl; std::cout << "Geo cell left = " << geo_cell_L << std::endl;
    std::cout << "Geo cell right = " << geo_cell_R << std::endl;
    std::cout << "Geo facet = " << geo_facet << std::endl;

    std::cout << "Solution cell left = " << sol_cell_L << std::endl;
    std::cout << "Solution facet left = " << sol_facet_L << std::endl;
    std::cout << "Solution cell right = " << sol_cell_R << std::endl;
    std::cout << "Solution facet right = " << sol_facet_R << std::endl;

    std::cout << "Facet has " << facet_blk.size() << " entries" <<
    std::endl; std::cout << "Permutation left = " <<
    facet_blk.permutation(idx_L).get().code()
              << "   " <<
    facet_blk.permutation(idx_L).get().type_id().as_string() << std::endl;
    std::cout << "Permutation right = " <<
    facet_blk.permutation(idx_R).get().code()
              << "   " <<
    facet_blk.permutation(idx_R).get().type_id().as_string() << std::endl;
    std::cout << std::endl;
    */

    m_geo_cache_facets_L.flush();
    m_geo_cache_facets_R.flush();

    m_sol_cache_facets_L.flush();
    m_sol_cache_facets_R.flush();

    const mesh::CellGeometry<MeshConfig::GDIM> geo_facet_coords_L =
        tcell_L.coordinates(FACET_DIM, facet_blk.local_id(idx_in_block_L));

    const mesh::CellGeometry<MeshConfig::GDIM> geo_facet_coords_R =
        tcell_R.coordinates(FACET_DIM, facet_blk.local_id(idx_in_block_R));

    m_geo_cache_facets_L.push_back_to_buffer(geo_facet_coords_L, geo_key_L);
    m_geo_cache_facets_R.push_back_to_buffer(geo_facet_coords_R, geo_key_R);

    sol_facet_L_permutation.change_type(sol_facet_L.pt_set_id(),
                                        facet_blk.permutation(idx_in_block_L).get().code());
    sol_facet_R_permutation.change_type(sol_facet_R.pt_set_id(),
                                        facet_blk.permutation(idx_in_block_R).get().code());

    m_sol_cache_facets_L.push_back_to_buffer(sol_facet_L, solution, sol_key_L);
    m_sol_cache_facets_R.push_back_to_buffer(sol_facet_R, solution, sol_key_R);

    /*
    m_sol_cache_facets_R.push_back_to_buffer(
        sol_facet_R, sol_facet_R_permutation, solution,
        mesh::DefaultStdRegMapKey(sol_facet_R.std_region_id(), 0u, P0));
    */

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
    m_sol_metric_facets_R.evaluate(m_geo_metric_facets_R, m_sol_cache_facets_R,
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

    const Uint cell_idx_in_metric = 0;
    const common::Range1D<Uint> facet_range(0, 0);

    const mesh::PointSetTag sol_tag_L = sol_cell_L.pt_set_id();
    const mesh::PointSetTag sol_tag_R = sol_cell_R.pt_set_id();

    const mesh::PointSetTag sol_facet_tag_L = sol_facet_L.pt_set_id();
    const mesh::PointSetTag sol_facet_tag_R = sol_facet_R.pt_set_id();

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

    common::PtrHandle<facet_method_data_type> facet_method_data_L =
        m_method_data_facet_L.std_region_data(key_data_L);

    common::PtrHandle<facet_method_data_type> facet_method_data_R =
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
      // Save the residuals to cache, this will be needed later
      // for finite differencing in order to obtain Jacobian
      // contributions by the facets
      // ---------------------------------------------------------

      math::DenseVecView<Real> facet_L_res = (*facet_method_data_L).elem_node_res();

      for (Uint nl = 0; nl < sol_facet_L.nb_vert(); ++nl)
      {
        const Real inv_volume_L = 1. / nodal_dual_volume[sol_facet_L.vertex(nl)];
        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          facet_L_res[nl * NEQ + eq] *= inv_volume_L;
        }
      }

      m_res_cache_facets_L.flush();
      m_res_cache_facets_L.push_vec_to_buffer(
          sol_facet_L, facet_L_res,
          mesh::PointSetTagExt(sol_facet_L.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

      math::DenseVecView<Real> facet_R_res = (*facet_method_data_R).elem_node_res();

      for (Uint nr = 0; nr < sol_facet_R.nb_vert(); ++nr)
      {
        /*
        const Real inv_volume_R =
            1. /
            nodal_dual_volume[sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(nr))];
        */

        const Real inv_volume_R = 1. / nodal_dual_volume[sol_facet_R.vertex(nr)];
        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          facet_R_res[nr * NEQ + eq] *= inv_volume_R;
        }
      }

      m_res_cache_facets_R.flush();
      m_res_cache_facets_R.push_vec_to_buffer(
          sol_facet_R, facet_R_res,
          mesh::PointSetTagExt(sol_facet_R.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

      // ---------------------------------------------------------
      // Accumulate NODAL residuals contributed by the left facet
      // ---------------------------------------------------------
      for (Uint n = 0; n < sol_facet_L.nb_vert(); ++n)
      {
        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          const Uint dof_idx = sol_facet_L.vertex(n) * NEQ + eq;
          if (!is_Dirichlet_node[dof_idx])
          {
            m_rhs_buffer.push_back(std::tuple<Uint, Real>(dof_idx, -facet_L_res[n * NEQ + eq]));
          }
        }
      }

      // ---------------------------------------------------------
      // Accumulate NODAL residuals contributed by the right facet
      // ---------------------------------------------------------
      for (Uint n = 0; n < sol_facet_R.nb_vert(); ++n)
      {
        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          /*
          const Uint dof_idx =
              sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(n))
          * NEQ
          + eq;
          */
          const Uint dof_idx = sol_facet_R.vertex(n) * NEQ + eq;
          if (!is_Dirichlet_node[dof_idx])
          {
            m_rhs_buffer.push_back(std::tuple<Uint, Real>(dof_idx, -facet_R_res[n * NEQ + eq]));
          }
        }
      }
    } // If idx_L and idx_R are different (i.e. we are in interior face)

    // --------------------------------------------------------------
    // Accumulate to update coefficient contributed by the left facet
    // --------------------------------------------------------------

    for (Uint n = 0; n < sol_facet_L.nb_vert(); ++n)
    {
      // update_coeff[sol_facet_L.vertex(n)] +=
      //    (*facet_method_data_L).m_elem_wave_speed[n];
      const Real volume_L = nodal_dual_volume[sol_facet_L.vertex(n)];

      /*
      time_update.accumulate_nodal_wave_speed(
          sol_facet_L.vertex(n), volume_L *
      (*facet_method_data_L).m_elem_wave_speed[n]);
      */

      m_time_update_buffer.push_back(std::tuple<Uint, Real>(
          sol_facet_L.vertex(n), volume_L * (*facet_method_data_L).m_elem_wave_speed[n]));
    }

    if (idx_in_block_L != idx_in_block_R)
    {
      // ---------------------------------------------------------------
      // Accumulate to update coefficient contributed by the right facet
      // This is done only in case the right facet exists (i.e. is
      // different from the left facet)
      // ---------------------------------------------------------------

      for (Uint n = 0; n < sol_facet_R.nb_vert(); ++n)
      {
        /*
        const Real volume_R =
            nodal_dual_volume[sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(n))];
        time_update.accumulate_nodal_wave_speed(
            sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(n)),
            volume_R * (*facet_method_data_R).m_elem_wave_speed[n]);
        */
        const Real volume_R = nodal_dual_volume[sol_facet_R.vertex(n)];
        /*
        time_update.accumulate_nodal_wave_speed(
            sol_facet_R.vertex(n),
            volume_R * (*facet_method_data_R).m_elem_wave_speed[n]);
        */

        m_time_update_buffer.push_back(std::tuple<Uint, Real>(
            sol_facet_R.vertex(n), volume_R * (*facet_method_data_R).m_elem_wave_speed[n]));
      }
    } // If idx_L and idx_R are different (i.e. we are in interior face)

    // **************************************************************************
    // 4) COMPUTE NUMERICAL JACOBIANS AND ACCUMULATE THEM TO SYSTEM MATRIX
    // **************************************************************************

    const Uint nb_dof_per_L_facet = sol_facet_L.nb_vert();
    const Uint nb_dof_per_R_facet = sol_facet_R.nb_vert();

    if (idx_in_block_L != idx_in_block_R)
    {
      const mesh::PointSetTagExt key_data_L =
          mesh::PointSetTagExt(sol_facet_L.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u);

      const mesh::PointSetTagExt key_data_R =
          mesh::PointSetTagExt(sol_facet_R.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u);

      common::PtrHandle<facet_method_data_type> facet_method_data_L =
          m_method_data_facet_L.std_region_data(key_data_L);

      common::PtrHandle<facet_method_data_type> facet_method_data_R =
          m_method_data_facet_R.std_region_data(key_data_R);

      const typename facet_geo_metric_type::cellwise_metric geo_facet_met_L =
          m_geo_metric_facets_L.cellwise_values(cell_idx_in_metric);

      const typename facet_geo_metric_type::cellwise_metric geo_facet_met_R =
          m_geo_metric_facets_R.cellwise_values(cell_idx_in_metric);

      // Cached residuals for left and right facet - these do not change
      // as we evaluate Jacobian contributions
      const math::DenseConstMatView<Real> res_nodal_values_L =
          m_res_cache_facets_L.cell_values(cell_idx_in_metric);

      const math::DenseConstMatView<Real> res_nodal_values_R =
          m_res_cache_facets_R.cell_values(cell_idx_in_metric);

      // **************************************************************************
      // 4a) PERTURB SOLUTION VALUES IN LEFT FACET AND COMPUTE DERIVATIVES
      //     OF RESIDUAL CONTRIBUTIONS TO LEFT FACET __AND__ RIGHT FACET
      // **************************************************************************

      // Metric terms coming from the right facet will not change while
      // we're perturbing the solution on the left facet, it is sufficient
      // to evaluate them here before looping over local DOFs of the left
      // facet
      m_sol_metric_facets_R.empty_buffer();
      m_sol_metric_facets_R.evaluate(m_geo_metric_facets_R, m_sol_cache_facets_R,
                                     interpolation::ComputeMetricDerivs{true},
                                     interpolation::RebuildMetricIndex{true});

      m_flux_metric_facets_R.empty_buffer();
      m_flux_metric_facets_R.evaluate(m_geo_cache_facets_R, m_geo_metric_facets_R,
                                      m_sol_cache_facets_R, m_sol_metric_facets_R,
                                      interpolation::RebuildMetricIndex{true});

      for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_L_facet; ++i_dof_in_elem)
      {
        for (Uint comp_u = 0; comp_u < NEQ; ++comp_u)
        {
          /*
          m_sol_cache_facets_L.perturb_values(mesh::PointSetTagExt(sol_facet_L.std_region_id(),
          P0, mesh::CellTransform::DO_NOTHING, 0u), i_dof_in_elem,
          comp_u);
          */

          m_sol_cache_facets_L.perturb_values(facet_range, i_dof_in_elem, comp_u);

          /*
          const math::DenseConstVecView<Real>
          unperturbed_node_in_L_facet =
              m_sol_cache_facets_L.unperturbed_values(mesh::PoinSetTagExt(
                  sol_facet_L.std_region_id(), P0,
          mesh::CellTransform::DO_NOTHING, 0u));
          */

          const math::DenseConstVecView<Real> unperturbed_node_in_L_facet =
              m_sol_cache_facets_L.unperturbed_values(facet_range);

          m_sol_metric_facets_L.empty_buffer();
          m_sol_metric_facets_L.evaluate(m_geo_metric_facets_L, m_sol_cache_facets_L,
                                         interpolation::ComputeMetricDerivs{true},
                                         interpolation::RebuildMetricIndex{true});

          m_flux_metric_facets_L.empty_buffer();
          m_flux_metric_facets_L.evaluate(m_geo_cache_facets_L, m_geo_metric_facets_L,
                                          m_sol_cache_facets_L, m_sol_metric_facets_L,
                                          interpolation::RebuildMetricIndex{true});

          /*
          m_sol_metric_facets_R.empty_buffer();
          m_sol_metric_facets_R.evaluate(m_geo_metric_facets_R,
          m_sol_cache_facets_R, true, true);
          */

          const math::DenseConstMatView<Real> sol_facet_nodal_values_L =
              m_sol_cache_facets_L.cell_values(cell_idx_in_metric);

          const math::DenseConstMatView<Real> sol_facet_nodal_values_R =
              m_sol_cache_facets_R.cell_values(cell_idx_in_metric);

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
              sol_facet_met_L, quad_p_L, sol_facet_met_R, quad_p_R, flux_facet_met_L,
              flux_facet_met_R, *facet_method_data_L, *facet_method_data_R);

          const Real inv_du = 1. / (sol_facet_nodal_values_L(i_dof_in_elem, comp_u) -
                                    unperturbed_node_in_L_facet[cell_idx_in_metric]);

          const Uint mat_col_idx = NEQ * sol_facet_L.vertex(i_dof_in_elem) + comp_u;

          // PERTURBATION OF LEFT FACET RESIDUAL WITH RESPECT TO
          // PERTURBATION OF LEFT FACET SOLUTION VALUE
          math::DenseVecView<Real> perturbed_res_nodal_values_L =
              (*facet_method_data_L).elem_node_res();

          for (Uint v = 0; v < sol_facet_L.nb_vert(); ++v)
          {
            const Uint left_vertex_id   = sol_facet_L.vertex(v);
            const Real inv_nodal_volume = 1. / nodal_dual_volume[left_vertex_id];

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              perturbed_res_nodal_values_L[v * NEQ + comp] *= inv_nodal_volume;
            }
          }

          for (Uint v = 0; v < sol_facet_L.nb_vert(); ++v)
          {
            const Uint left_vertex_id = sol_facet_L.vertex(v);

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              // For each component
              // When v == i_dof_in_elem, we are computing
              // derivatives of residuals corresponding to node v
              // with respect to components of u associated to v
              // (diagonal blocks) in the system matrix
              const Uint mat_row_idx = NEQ * left_vertex_id + comp;

              if (!is_Dirichlet_node[mat_row_idx])
              {
                const Real res_diff =
                    (perturbed_res_nodal_values_L[v * NEQ + comp] - res_nodal_values_L(v, comp)) *
                    inv_du;

                m_mat_buffer.push_back(
                    std::tuple<Uint, Uint, Real>(mat_row_idx, mat_col_idx, res_diff));

              } // If this is not Dirichlet node
            }   // Loop over equation components
          }     // Loop over vertices of left facet

          // PERTURBATION OF RIGHT FACET RESIDUAL WITH RESPECT TO
          // PERTURBATION OF LEFT FACET SOLUTION VALUE
          math::DenseVecView<Real> perturbed_res_nodal_values_R =
              (*facet_method_data_R).elem_node_res();

          for (Uint v = 0; v < sol_facet_R.nb_vert(); ++v)
          {
            /*
            const Uint right_vertex_id =
                sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(v));
            */
            const Uint right_vertex_id  = sol_facet_R.vertex(v);
            const Real inv_nodal_volume = 1. / nodal_dual_volume[right_vertex_id];

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              perturbed_res_nodal_values_R[v * NEQ + comp] *= inv_nodal_volume;
            }
          }

          for (Uint v = 0; v < sol_facet_R.nb_vert(); ++v)
          {
            /*
            const Uint right_vertex_id =
                sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(v));
            */
            const Uint right_vertex_id = sol_facet_R.vertex(v);

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              // For each component
              // When v == i_dof_in_elem, we are computing
              // derivatives of residuals corresponding to node v
              // with respect to components of u associated to v
              // (diagonal blocks) in the system matrix const Uint
              // mat_row_idx = NEQ * sol_facet_R.vertex(v) + comp;
              const Uint mat_row_idx = NEQ * right_vertex_id + comp;

              if (!is_Dirichlet_node[mat_row_idx])
              {
                const Real res_diff =
                    (perturbed_res_nodal_values_R[v * NEQ + comp] - res_nodal_values_R(v, comp)) *
                    inv_du;

                m_mat_buffer.push_back(
                    std::tuple<Uint, Uint, Real>(mat_row_idx, mat_col_idx, res_diff));

              } // If this is not Dirichlet node
            }   // Loop over equation components
          }     // Loop over vertices of left facet

          /*
          m_sol_cache_facets_L.remove_perturbation(mesh::PointSetTagExt(
              sol_facet_L.std_region_id(), P0,
          mesh::CellTransform::DO_NOTHING, 0u));
          */

          m_sol_cache_facets_L.remove_perturbation(facet_range);

        } // Loop over equation components

      } // Loop over local nodes of left facet

      // Metric terms coming from the left facet will not change while
      // we're perturbing the solution on the right facet, it is
      // sufficient to evaluate them here before looping over local DOFs
      // of the right facet
      m_sol_metric_facets_L.empty_buffer();
      m_sol_metric_facets_L.evaluate(m_geo_metric_facets_L, m_sol_cache_facets_L,
                                     interpolation::ComputeMetricDerivs{true},
                                     interpolation::RebuildMetricIndex{true});

      m_flux_metric_facets_L.empty_buffer();
      m_flux_metric_facets_L.evaluate(m_geo_cache_facets_L, m_geo_metric_facets_L,
                                      m_sol_cache_facets_L, m_sol_metric_facets_L,
                                      interpolation::RebuildMetricIndex{true});

      // **************************************************************************
      // 4b) PERTURB SOLUTION VALUES IN RIGHT FACET AND COMPUTE
      // DERIVATIVES
      //     OF RESIDUAL CONTRIBUTIONS TO LEFT FACET __AND__ RIGHT FACET
      // **************************************************************************

      for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_R_facet; ++i_dof_in_elem)
      {
        for (Uint comp_u = 0; comp_u < NEQ; ++comp_u)
        {
          /*
          m_sol_cache_facets_R.perturb_values(mesh::PointSetTagExt(sol_facet_R.std_region_id(),
          P0, mesh::CellTransform::DO_NOTHING, 0u), i_dof_in_elem,
          comp_u);
          */

          m_sol_cache_facets_R.perturb_values(facet_range, i_dof_in_elem, comp_u);

          /*
          const math::DenseConstVecView<Real>
          unperturbed_node_in_R_facet =
              m_sol_cache_facets_R.unperturbed_values(mesh::PointSetTagExt(
                  sol_facet_R.std_region_id(), P0,
          mesh::CellTransform::DO_NOTHING, 0u));
          */

          const math::DenseConstVecView<Real> unperturbed_node_in_R_facet =
              m_sol_cache_facets_R.unperturbed_values(facet_range);

          /*
          m_sol_metric_facets_L.empty_buffer();
          m_sol_metric_facets_L.evaluate(m_geo_metric_facets_L,
          m_sol_cache_facets_L, true, true);
          */

          m_sol_metric_facets_R.empty_buffer();
          m_sol_metric_facets_R.evaluate(m_geo_metric_facets_R, m_sol_cache_facets_R,
                                         interpolation::ComputeMetricDerivs{true},
                                         interpolation::RebuildMetricIndex{true});

          m_flux_metric_facets_R.empty_buffer();
          m_flux_metric_facets_R.evaluate(m_geo_cache_facets_R, m_geo_metric_facets_R,
                                          m_sol_cache_facets_R, m_sol_metric_facets_R,
                                          interpolation::RebuildMetricIndex{true});

          const math::DenseConstMatView<Real> sol_facet_nodal_values_L =
              m_sol_cache_facets_L.cell_values(cell_idx_in_metric);

          const math::DenseConstMatView<Real> sol_facet_nodal_values_R =
              m_sol_cache_facets_R.cell_values(cell_idx_in_metric);

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
              sol_facet_met_L, quad_p_L, sol_facet_met_R, quad_p_R, flux_facet_met_L,
              flux_facet_met_R, *facet_method_data_L, *facet_method_data_R);

          const Real inv_du = 1. / (sol_facet_nodal_values_R(i_dof_in_elem, comp_u) -
                                    unperturbed_node_in_R_facet[cell_idx_in_metric]);

          /*
          const Uint mat_col_idx =
              NEQ *
                  sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(i_dof_in_elem))
          + comp_u;
          */

          const Uint mat_col_idx = NEQ * sol_facet_R.vertex(i_dof_in_elem) + comp_u;

          // PERTURBATION OF LEFT FACET RESIDUAL WITH RESPECT TO
          // PERTURBATION OF RIGHT FACET SOLUTION VALUE
          math::DenseVecView<Real> perturbed_res_nodal_values_L =
              (*facet_method_data_L).elem_node_res();

          for (Uint v = 0; v < sol_facet_L.nb_vert(); ++v)
          {
            const Uint left_vertex_id   = sol_facet_L.vertex(v);
            const Real inv_nodal_volume = 1. / nodal_dual_volume[left_vertex_id];

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              perturbed_res_nodal_values_L[v * NEQ + comp] *= inv_nodal_volume;
            }
          }

          for (Uint v = 0; v < sol_facet_L.nb_vert(); ++v)
          {
            const Uint left_vertex_id = sol_facet_L.vertex(v);

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              // For each component
              // When v == i_dof_in_elem, we are computing
              // derivatives of residuals corresponding to node v
              // with respect to components of u associated to v
              // (diagonal blocks) in the system matrix
              const Uint mat_row_idx = NEQ * left_vertex_id + comp;

              if (!is_Dirichlet_node[mat_row_idx])
              {
                const Real res_diff =
                    (perturbed_res_nodal_values_L[v * NEQ + comp] - res_nodal_values_L(v, comp)) *
                    inv_du;

                m_mat_buffer.push_back(
                    std::tuple<Uint, Uint, Real>(mat_row_idx, mat_col_idx, res_diff));

              } // If this is not Dirichlet node
            }   // Loop over equation components
          }     // Loop over vertices of left facet

          // PERTURBATION OF RIGHT FACET RESIDUAL WITH RESPECT TO
          // PERTURBATION OF RIGHT FACET SOLUTION VALUE
          math::DenseVecView<Real> perturbed_res_nodal_values_R =
              (*facet_method_data_R).elem_node_res();

          for (Uint v = 0; v < sol_facet_R.nb_vert(); ++v)
          {
            /*
            const Uint right_vertex_id =
                sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(v));
            */

            const Uint right_vertex_id = sol_facet_R.vertex(v);

            const Real inv_nodal_volume = 1. / nodal_dual_volume[right_vertex_id];

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              perturbed_res_nodal_values_R[v * NEQ + comp] *= inv_nodal_volume;
            }
          }

          for (Uint v = 0; v < sol_facet_R.nb_vert(); ++v)
          {
            /*
            const Uint right_vertex_id =
                sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(v));
            */

            const Uint right_vertex_id = sol_facet_R.vertex(v);

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              // For each component
              // When v == i_dof_in_elem, we are computing
              // derivatives of residuals corresponding to node v
              // with respect to components of u associated to v
              // (diagonal blocks) in the system matrix const Uint
              // mat_row_idx = NEQ * sol_facet_R.vertex(v) + comp;
              const Uint mat_row_idx = NEQ * right_vertex_id + comp;

              if (!is_Dirichlet_node[mat_row_idx])
              {
                const Real res_diff =
                    (perturbed_res_nodal_values_R[v * NEQ + comp] - res_nodal_values_R(v, comp)) *
                    inv_du;

                m_mat_buffer.push_back(
                    std::tuple<Uint, Uint, Real>(mat_row_idx, mat_col_idx, res_diff));

              } // If this is not Dirichlet node
            }   // Loop over equation components
          }     // Loop over vertices of left facet

          /*
          m_sol_cache_facets_R.remove_perturbation(mesh::PointSetTagExt(
              sol_facet_R.std_region_id(), P0,
          mesh::CellTransform::DO_NOTHING, 0u));
          */

          m_sol_cache_facets_R.remove_perturbation(facet_range);

        } // Loop over equation components

      } // Loop over local nodes of right facet

    } // If left and right facet are different (hence the residual is not
      // zero)

    /*
    if (((f != m_first_facet_idx) && (f % m_nb_blocks == 0)) || (f ==
    m_last_facet_idx))
    {
      mat.add_values(m_mat_buffer);
      rhs.add_values(0, m_rhs_buffer);
      time_update.accumulate_wave_speeds(m_time_update_buffer);
      m_mat_buffer.resize(0);
      m_rhs_buffer.resize(0);
      m_time_update_buffer.resize(0);
    }
    */
  } // Loop over facets

  mat.add_values(m_mat_buffer);
  rhs.add_values(0, m_rhs_buffer);
  time_update.accumulate_wave_speeds(m_time_update_buffer);
  m_mat_buffer.resize(0);
  m_rhs_buffer.resize(0);
  m_time_update_buffer.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void PGRDMFullJacFacetAssemblyWorker<MeshConfig, Physics, CellScheme, FacetScheme>::
    assemble_rhs_facets_part(const tria_t &cell_topology, const dof_map_t &sol_cells,
                             const interpolation::VectorMeshFunction<Real> &solution,
                             RDTimeUpdate &time_update, const std::vector<bool> &is_Dirichlet_node,
                             ls::TpetraMultiVector<Real> &rhs)
{
  mesh::EntityDofRealign sol_facet_L_permutation, sol_facet_R_permutation;
  mesh::QuadraturePermutation quad_p_L, quad_p_R;

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

  m_rhs_buffer.resize(0);
  m_time_update_buffer.resize(0);

  for (Uint f = m_first_facet_idx; f <= m_last_facet_idx; ++f)
  {
    // std::cout << "f = " << f << std::endl;

    const mesh::TraceIncidences facet_blk =
        cell_topology.active_skeleton_entry(FACET_DIM, mesh::ActiveIdx(f));

    // Get the permutation sign of the first entity in the block (position
    // 0)
    const mesh::EntityRealignCode pcode_L = facet_blk.permutation(0).get().code();

    // If the entity on position 0 in incidence block needs to be flipped,
    // it is the entity on the rhs of this facet
    const Uint idx_in_block_L = (pcode_L.nb_flips() == 0) ? 0 : 1;

    // If the facet block has size 2, we have 2 entities forming the face
    // (i.e K+ and K-). In that case, idx_R is the 'other' index (choosing
    // from 0 and 1) than idx_L, and can be determined as idx_R = (idx_L +
    // 1) % 2 If the facet block has size 1, then idx_R has to be the same
    // as index left: equal to 0
    const Uint idx_in_block_R = (facet_blk.size() == 2) ? (idx_in_block_L + 1) % 2 : 0;

    // Get the permutation sign of the first entity in the block (position
    // 0)
    const mesh::EntityRealignCode pcode_R = facet_blk.permutation(idx_in_block_R).get().code();

    // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS OF
    // CELLS, NOT THEIR ACTIVE INDICES !!! IT IS FOR THIS REASON THAT WE
    // FIRST NEED TO RETRIEVE TO TOPOLOGICAL CELLS AND GET THEIR ACTIVE
    // INDICES BEFORE GETTING THE CORRECT CELLS FROM THE DOF HANDLER
    const mesh::CellTopologyView<MeshConfig> tcell_L =
        cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_in_block_L)));
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_in_block_R)));

    const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
    const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

    const Uint geo_key_cache_pos_idx      = f - m_first_facet_idx;
    const mesh::DiscreteElemKey geo_key_L = m_geo_key_cache_L.key(geo_key_cache_pos_idx);
    const mesh::DiscreteElemKey geo_key_R = m_geo_key_cache_R.key(geo_key_cache_pos_idx);

    // Get the cell and facet on the left- and right-hand side of interface
    // in SOLUTION space
    const mesh::MeshEntity sol_cell_L = sol_cells.active_cell(mesh::ActiveIdx(active_cell_id_L));
    const mesh::MeshEntity sol_cell_R = sol_cells.active_cell(mesh::ActiveIdx(active_cell_id_R));

    const mesh::MeshEntity sol_facet_L =
        sol_cell_L.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_L));
    const mesh::MeshEntity sol_facet_R =
        sol_cell_R.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_R));

    const mesh::PointSetTagExt sol_key_L(sol_facet_L.pt_set_id(), P0, pcode_L.adapt_op_id(),
                                         pcode_L.local_pos_in_parent());

    const mesh::PointSetTagExt sol_key_R(sol_facet_R.pt_set_id(), P0, pcode_R.adapt_op_id(),
                                         pcode_R.local_pos_in_parent());

    /*
    std::cout << "====================================================" <<
    std::endl; std::cout << "Geo cell left = " << geo_cell_L << std::endl;
    std::cout << "Geo cell right = " << geo_cell_R << std::endl;
    std::cout << "Geo facet = " << geo_facet << std::endl;

    std::cout << "Solution cell left = " << sol_cell_L << std::endl;
    std::cout << "Solution facet left = " << sol_facet_L << std::endl;
    std::cout << "Solution cell right = " << sol_cell_R << std::endl;
    std::cout << "Solution facet right = " << sol_facet_R << std::endl;

    std::cout << "Facet has " << facet_blk.size() << " entries" <<
    std::endl; std::cout << "Permutation left = " <<
    facet_blk.permutation(idx_L).get().code()
              << "   " <<
    facet_blk.permutation(idx_L).get().type_id().as_string() << std::endl;
    std::cout << "Permutation right = " <<
    facet_blk.permutation(idx_R).get().code()
              << "   " <<
    facet_blk.permutation(idx_R).get().type_id().as_string() << std::endl;
    std::cout << std::endl;
    */

    m_geo_cache_facets_L.flush();
    m_geo_cache_facets_R.flush();

    m_sol_cache_facets_L.flush();
    m_sol_cache_facets_R.flush();

    const mesh::CellGeometry<MeshConfig::GDIM> geo_facet_coords_L =
        tcell_L.coordinates(FACET_DIM, facet_blk.local_id(idx_in_block_L));

    const mesh::CellGeometry<MeshConfig::GDIM> geo_facet_coords_R =
        tcell_R.coordinates(FACET_DIM, facet_blk.local_id(idx_in_block_R));

    m_geo_cache_facets_L.push_back_to_buffer(geo_facet_coords_L, geo_key_L);
    m_geo_cache_facets_R.push_back_to_buffer(geo_facet_coords_R, geo_key_R);

    sol_facet_L_permutation.change_type(sol_facet_L.pt_set_id(),
                                        facet_blk.permutation(idx_in_block_L).get().code());
    sol_facet_R_permutation.change_type(sol_facet_R.pt_set_id(),
                                        facet_blk.permutation(idx_in_block_R).get().code());

    m_sol_cache_facets_L.push_back_to_buffer(sol_facet_L, solution, sol_key_L);
    m_sol_cache_facets_R.push_back_to_buffer(sol_facet_R, solution, sol_key_R);

    /*
    m_sol_cache_facets_R.push_back_to_buffer(
        sol_facet_R, sol_facet_R_permutation, solution,
        mesh::DefaultStdRegMapKey(sol_facet_R.std_region_id(), 0u, P0));
    */

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
    m_sol_metric_facets_R.evaluate(m_geo_metric_facets_R, m_sol_cache_facets_R,
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

    const Uint cell_idx_in_metric = 0;

    const mesh::PointSetTag sol_tag_L = sol_cell_L.pt_set_id();
    const mesh::PointSetTag sol_tag_R = sol_cell_R.pt_set_id();

    const mesh::PointSetTag sol_facet_tag_L = sol_facet_L.pt_set_id();
    const mesh::PointSetTag sol_facet_tag_R = sol_facet_R.pt_set_id();

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

    common::PtrHandle<facet_method_data_type> facet_method_data_L =
        m_method_data_facet_L.std_region_data(key_data_L);

    common::PtrHandle<facet_method_data_type> facet_method_data_R =
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
        const Real inv_volume_L = 1. / nodal_dual_volume[sol_facet_L.vertex(n)];

        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          elem_node_res_L[n * NEQ + eq] *= inv_volume_L;

          const Uint dof_idx = sol_facet_L.vertex(n) * NEQ + eq;

          if (!is_Dirichlet_node[dof_idx])
          {
            m_rhs_buffer.push_back(std::tuple<Uint, Real>(dof_idx, -elem_node_res_L[n * NEQ + eq]));
          }
        }
      }

      // ---------------------------------------------------------
      // Accumulate NODAL residuals contributed by the right facet
      // ---------------------------------------------------------
      math::DenseVecView<Real> elem_node_res_R = (*facet_method_data_R).elem_node_res();

      for (Uint n = 0; n < sol_facet_R.nb_vert(); ++n)
      {
        /*
        const Real inv_volume_R =
            1. /
            nodal_dual_volume[sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(n))];
        */

        const Real inv_volume_R = 1. / nodal_dual_volume[sol_facet_R.vertex(n)];

        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          elem_node_res_R[n * NEQ + eq] *= inv_volume_R;

          /*
          const Uint dof_idx =
              sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(n))
          * NEQ
          + eq;
          */
          const Uint dof_idx = sol_facet_R.vertex(n) * NEQ + eq;

          if (!is_Dirichlet_node[dof_idx])
          {
            m_rhs_buffer.push_back(std::tuple<Uint, Real>(dof_idx, -elem_node_res_R[n * NEQ + eq]));
          }
        }
      }
    } // If idx_L and idx_R are different (i.e. we are in interior face)

    // --------------------------------------------------------------
    // Accumulate to update coefficient contributed by the left facet
    // --------------------------------------------------------------

    for (Uint n = 0; n < sol_facet_L.nb_vert(); ++n)
    {
      // update_coeff[sol_facet_L.vertex(n)] +=
      //    (*facet_method_data_L).m_elem_wave_speed[n];
      const Real volume_L = nodal_dual_volume[sol_facet_L.vertex(n)];

      /*
      time_update.accumulate_nodal_wave_speed(
          sol_facet_L.vertex(n), volume_L *
      (*facet_method_data_L).m_elem_wave_speed[n]);
      */

      m_time_update_buffer.push_back(std::tuple<Uint, Real>(
          sol_facet_L.vertex(n), volume_L * (*facet_method_data_L).m_elem_wave_speed[n]));
    }

    if (idx_in_block_L != idx_in_block_R)
    {
      // ---------------------------------------------------------------
      // Accumulate to update coefficient contributed by the right facet
      // This is done only in case the right facet exists (i.e. is
      // different from the left facet)
      // ---------------------------------------------------------------

      for (Uint n = 0; n < sol_facet_R.nb_vert(); ++n)
      {
        /*
        const Real volume_R =
            nodal_dual_volume[sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(n))];
        time_update.accumulate_nodal_wave_speed(
            sol_facet_R.vertex(sol_facet_R_permutation.get().vertex(n)),
            volume_R * (*facet_method_data_R).m_elem_wave_speed[n]);
        */
        const Real volume_R = nodal_dual_volume[sol_facet_R.vertex(n)];

        /*
        time_update.accumulate_nodal_wave_speed(
            sol_facet_R.vertex(n),
            volume_R * (*facet_method_data_R).m_elem_wave_speed[n]);
        */

        m_time_update_buffer.push_back(std::tuple<Uint, Real>(
            sol_facet_R.vertex(n), volume_R * (*facet_method_data_R).m_elem_wave_speed[n]));
      }
    } // If idx_L and idx_R are different (i.e. we are in interior face)

    /*
    if (((f != m_first_facet_idx) && (f % m_nb_blocks == 0)) || (f ==
    m_last_facet_idx))
    {
      rhs.add_values(0, m_rhs_buffer);
      m_rhs_buffer.resize(0);

      time_update.accumulate_wave_speeds(m_time_update_buffer);
      m_time_update_buffer.resize(0);
    }
    */

    rhs.add_values(0, m_rhs_buffer);
    time_update.accumulate_wave_speeds(m_time_update_buffer);

    m_rhs_buffer.resize(0);
    m_time_update_buffer.resize(0);

  } // Loop over facets
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void PGRDMFullJacFacetAssemblyWorker<MeshConfig, Physics, CellScheme, FacetScheme>::fill_key_cache(
    const tria_t &cell_topology, const dof_map_t &sol_cells)
{
  m_geo_key_cache_L.clear();
  m_geo_key_cache_R.clear();

  auto sf_generator = [this](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, m_sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [this](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, m_quad_order, m_quad_type);
  };

  mesh::EntityDofRealign sol_facet_L_permutation, sol_facet_R_permutation;
  mesh::QuadraturePermutation quad_p_L, quad_p_R;

  for (Uint f = m_first_facet_idx; f <= m_last_facet_idx; ++f)
  {
    // std::cout << "f = " << f << std::endl;

    const mesh::TraceIncidences facet_blk =
        cell_topology.active_skeleton_entry(FACET_DIM, mesh::ActiveIdx(f));

    // Get the permutation sign of the first entity in the block (position
    // 0)
    const mesh::EntityRealignCode pcode_L = facet_blk.permutation(0).get().code();

    // If the entity on position 0 in incidence block needs to be flipped,
    // it is the entity on the rhs of this facet
    const Uint idx_in_block_L = (pcode_L.nb_flips() == 0) ? 0 : 1;

    // If the facet block has size 2, we have 2 entities forming the face
    // (i.e K+ and K-). In that case, idx_R is the 'other' index (choosing
    // from 0 and 1) than idx_L, and can be determined as idx_R = (idx_L +
    // 1) % 2 If the facet block has size 1, then idx_R has to be the same
    // as index left: equal to 0
    const Uint idx_in_block_R = (facet_blk.size() == 2) ? (idx_in_block_L + 1) % 2 : 0;

    // Get the permutation sign of the first entity in the block (position
    // 0)
    const mesh::EntityRealignCode pcode_R = facet_blk.permutation(idx_in_block_R).get().code();

    // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS OF
    // CELLS, NOT THEIR ACTIVE INDICES !!! IT IS FOR THIS REASON THAT WE
    // FIRST NEED TO RETRIEVE TO TOPOLOGICAL CELLS AND GET THEIR ACTIVE
    // INDICES BEFORE GETTING THE CORRECT CELLS FROM THE DOF HANDLER
    const mesh::CellTopologyView<MeshConfig> tcell_L =
        cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_in_block_L)));
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_in_block_R)));

    const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
    const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

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

    m_geo_key_cache_L.push_back(geo_key_L);
    m_geo_key_cache_R.push_back(geo_key_R);
  }
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
