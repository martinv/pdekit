#ifndef PDEKIT_D_RD_Method_Lax_Friedrichs_Implementation_hpp
#define PDEKIT_D_RD_Method_Lax_Friedrichs_Implementation_hpp

#include "interpolation/FluxSpaceMetric.hpp"
#include "mesh/point_set/QuadraturePermutation.hpp"
#include "solver/rdm/RDMethodScratchData.hpp"
#include "solver/rdm/RDTimeUpdate.hpp"
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

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
class DRDMethodLFImpl
{
  private:
  typedef typename internal::CellSchemeSelector<Physics, CellScheme>::type cell_scheme_type;
  typedef typename internal::FacetSchemeSelector<Physics, FacetScheme>::type facet_scheme_type;

  typedef typename cell_scheme_type::method_data cell_method_data_type;
  typedef typename facet_scheme_type::method_data facet_method_data_type;

  public:
  typedef interpolation::FunctionSpace<MeshConfig> f_space;

  /// Default constructor
  DRDMethodLFImpl();

  /// Default destructor
  ~DRDMethodLFImpl();

  /// Set the function space (reference elements)
  void configure_cell_spaces(const typename f_space::ptr &geo_cell_space,
                             const typename f_space::ptr &sol_cell_space);

  void configure_cell_contour_spaces(const typename f_space::ptr &geo_cell_contour_space,
                                     const typename f_space::ptr &sol_cell_contour_space);

  void configure_facet_spaces(const typename f_space::ptr &geo_facet_space,
                              const typename f_space::ptr &sol_facet_space);

  /// To complete one sweep over all entities in mesh, it is necessary to call
  /// 1) iterate_over_facets( ... )
  /// 2) iterate_over_cells( ... )
  /// It is important that the methods are called IN THAT ORDER (first
  /// 'iterate_over_facets' and then 'iterate_over_cells'
  void iterate_over_facets(const mesh::Tria<MeshConfig> &cell_topology,
                           const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
                           const typename result_of::dof_map_t<MeshConfig> &sol_dofs,
                           const typename f_space::vect_f &solution,
                           typename f_space::vect_f &nodal_residuals, RDTimeUpdate &time_update,
                           const Uint first_cell_idx, const Uint last_cell_idx);

  void iterate_over_cells(const mesh::Tria<MeshConfig> &cell_topology,
                          const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
                          const typename result_of::dof_map_t<MeshConfig> &sol_dofs,
                          const typename f_space::vect_f &solution,
                          typename f_space::vect_f &nodal_residuals, RDTimeUpdate &time_update,
                          const Uint first_cell_idx, const Uint last_cell_idx);

  void set_artificial_viscosity(
      const std::shared_ptr<interpolation::ScalarMeshFunction<Real>> &artificial_viscosity);

  private:
  enum
  {
    FACET_DIM           = MeshConfig::TDIM - _1D,
    EDGE_DIM            = _1D,
    FACET_DATA_TOPO_DIM = MeshConfig::TDIM - _1D
  };

  /// Geometry cache and metric for cells
  typedef interpolation::GeometryCache<MeshConfig::GDIM> cell_geo_cache_type;
  typedef interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> cell_geo_metric_type;
  typedef interpolation::FluxSpaceMetric<MeshConfig, Physics, MeshConfig::TDIM>
      cell_flux_metric_type;

  /// Solution cache and metric for cells
  typedef interpolation::SolutionCache cell_sol_cache_type;
  typedef interpolation::SolutionSpaceMetric<MeshConfig, MeshConfig::TDIM> cell_sol_metric_type;

  /// Geometry cache and metric for facets
  typedef interpolation::GeometryCache<MeshConfig::GDIM> facet_geo_cache_type;
  typedef interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, FACET_DIM>
      facet_geo_metric_type;

  /// Solution cache and metric for facets
  typedef interpolation::SolutionCache facet_sol_cache_type;
  typedef interpolation::SolutionSpaceMetric<MeshConfig, FACET_DIM> facet_sol_metric_type;

  typedef typename f_space::scal_f::ptr scal_function_ptr;
  typedef typename f_space::vect_f::ptr vect_function_ptr;

  /// DATA

  /// Pointer to the geometry space - cell interiors
  typename f_space::ptr m_geo_cell_space;

  /// Pointer to geometry space on cell facets
  typename f_space::ptr m_geo_facet_space;

  /// Pointer to the solution space - cell interiors
  typename f_space::ptr m_sol_cell_space;

  /// Pointer to solution space on cell facets
  typename f_space::ptr m_sol_facet_space;

  /// METRIC CONTAINERS FOR CELLS

  /// Geometry cache for cells
  cell_geo_cache_type m_geo_cache_cells;

  /// Geometry metric for cells
  cell_geo_metric_type m_geo_metric_cells;

  /// Solution cache for cells
  cell_sol_cache_type m_sol_cache_cells;

  /// Interpolated values and derivatives of the solution u_h
  cell_sol_metric_type m_sol_metric_cells;

  /// METRIC CONTAINERS FOR DATA ON MESH SKELETON
  /// HERE WE OPERATE ON ENTITIES WITH codim = 1

  /// Geometry cache for facets
  facet_geo_cache_type m_geo_cache_facets;

  /// Geometry metric for facets
  facet_geo_metric_type m_geo_metric_facets;

  /// Solution cache for facets
  facet_sol_cache_type m_sol_cache_facets_left;
  facet_sol_cache_type m_sol_cache_facets_right;

  /// Interpolated values and derivatives of the solution u_h
  /// on the mesh skeleton (facets)
  facet_sol_metric_type m_sol_metric_facets_left;
  facet_sol_metric_type m_sol_metric_facets_right;

  /// Number of blocks to fill in the geometry cache
  Uint m_geo_facet_blk_size;

  /// Number of blocks to fill in the solution cache
  Uint m_sol_facet_blk_size;

  /// Number of cells to fill the cell solution space cache
  Uint m_cell_blk_size;

  /// A map which associates to each element type in solution space a
  /// corresponding reference element for fluxes
  common::DataMap<mesh::PointSetTagExt, cell_method_data_type> m_method_data_cell;
  common::DataMap<mesh::PointSetTagExt, facet_method_data_type> m_method_data_left_facet;
  common::DataMap<mesh::PointSetTagExt, facet_method_data_type> m_method_data_right_facet;

  /// Pointers to the functions holding solution, residuals and update
  /// coefficients
  vect_function_ptr m_elem_residuals;

  /// Residual in one cell of the mesh
  typename Physics::FluxV m_active_cell_residual;

  /// Concrete RDS scheme
  cell_scheme_type m_cell_scheme;
  facet_scheme_type m_facet_scheme;

  /// Quadrature order in solution space
  Uint m_sol_quadrature_order;

  /// Type of quadrature in solution space for cell interior
  PointSetID m_sol_cell_quadrature_type;

  /// Type of quadrature in solution space for cell contours
  PointSetID m_sol_contour_quadrature_type;
};

// ----------------------------------------------------------------------------
// DISCONTINUOUS RD METHOD - IMPLEMENTATION FOR
// VARIABLE BETA SCHEMES
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
DRDMethodLFImpl<MeshConfig, Physics, CellScheme, FacetScheme>::DRDMethodLFImpl()
    : m_geo_cell_space(nullptr), m_geo_facet_space(nullptr), m_sol_cell_space(nullptr),
      m_sol_facet_space(nullptr), m_geo_facet_blk_size(0u), m_sol_facet_blk_size(0u),
      m_cell_blk_size(0u), m_sol_quadrature_order(P0),
      m_sol_cell_quadrature_type(PointSetID::Undefined),
      m_sol_contour_quadrature_type(PointSetID::Undefined)
{
  m_elem_residuals = vect_function_ptr(new typename f_space::vect_f("", "elem_residuals"));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
DRDMethodLFImpl<MeshConfig, Physics, CellScheme, FacetScheme>::~DRDMethodLFImpl()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodLFImpl<MeshConfig, Physics, CellScheme, FacetScheme>::configure_cell_spaces(
    const typename f_space::ptr &geo_cell_space, const typename f_space::ptr &sol_cell_space)
{
  // ---------------------
  // GEOMETRY SPACE CONFIG
  // ---------------------
  m_geo_cell_space = geo_cell_space;

  // Use number of blocks given by m_cell_blk_size for cell cache
  m_cell_blk_size = 150;

  m_geo_cache_cells.allocate(m_geo_cell_space->reference_elements(), m_cell_blk_size);
  m_geo_metric_cells.allocate_buffer(m_geo_cell_space->reference_elements(), m_cell_blk_size);

  // ---------------------
  // SOLUTION SPACE CONFIG
  // ---------------------
  m_sol_cell_space = sol_cell_space;

  m_sol_quadrature_order        = m_sol_cell_space->quadrature_order();
  m_sol_cell_quadrature_type    = m_sol_cell_space->quadrature_type();
  m_sol_contour_quadrature_type = m_sol_cell_quadrature_type;

  /*
  if (m_sol_cell_quadrature_type == PointSetID::Gauss)
  {
    m_sol_contour_quadrature_type = PointSetID::FaceGauss;
  }
  */

  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> const &sol_cell_type_map =
      m_sol_cell_space->reference_elements();

  // Cache for solution cells and data in cell volumes
  m_sol_cache_cells.allocate(sol_cell_type_map, m_cell_blk_size, Physics::NEQ);
  m_sol_metric_cells.allocate_buffer(sol_cell_type_map, m_cell_blk_size, Physics::NEQ);

  // ------------------------------
  // Prepare method data for cells
  // ------------------------------
  m_method_data_cell.clear();

  for (common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>::const_iterator it =
           sol_cell_type_map.cbegin();
       it != sol_cell_type_map.cend(); ++it)
  {
    const interpolation::FEValues &fev             = (*it.data_ptr());
    common::PtrHandle<cell_method_data_type> mdata = m_method_data_cell.create(it.key_value());
    (*mdata).resize_variables(fev);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodLFImpl<MeshConfig, Physics, CellScheme, FacetScheme>::configure_cell_contour_spaces(
    const typename f_space::ptr &geo_cell_contour_space,
    const typename f_space::ptr &sol_cell_contour_space)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodLFImpl<MeshConfig, Physics, CellScheme, FacetScheme>::configure_facet_spaces(
    const typename f_space::ptr &geo_facet_space, const typename f_space::ptr &sol_facet_space)
{
  // ---------------------
  // GEOMETRY SPACE CONFIG
  // ---------------------
  m_geo_facet_space = geo_facet_space;
  if (m_geo_facet_blk_size == 0)
    m_geo_facet_blk_size = 150; // geo_cells.nb_skeleton_facets();

  // Allocation of cache and metric for facets
  m_geo_cache_facets.allocate(m_geo_facet_space->reference_elements(), m_geo_facet_blk_size);
  m_geo_metric_facets.allocate_buffer(m_geo_facet_space->reference_elements(),
                                      m_geo_facet_blk_size);

  // ---------------------
  // SOLUTION SPACE CONFIG
  // ---------------------
  m_sol_facet_space = sol_facet_space;

  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> const &facet_elem_type_map =
      m_sol_facet_space->reference_elements();

  m_sol_facet_blk_size = m_geo_facet_blk_size;

  // Cache and metric for solution on facets (mesh skeleton)
  // FIXME: theoretically, this could be remove as the same data might be
  // provided by m_sol_cache_cell_contour_{left,right} and
  // m_sol_metric_cell_contour_{left,right}
  m_sol_cache_facets_left.allocate(m_sol_facet_space->reference_elements(), m_sol_facet_blk_size,
                                   Physics::NEQ);
  m_sol_cache_facets_right.allocate(m_sol_facet_space->reference_elements(), m_sol_facet_blk_size,
                                    Physics::NEQ);

  m_sol_metric_facets_left.allocate_buffer(m_sol_facet_space->reference_elements(),
                                           m_sol_facet_blk_size, Physics::NEQ);
  m_sol_metric_facets_right.allocate_buffer(m_sol_facet_space->reference_elements(),
                                            m_sol_facet_blk_size, Physics::NEQ);

  m_method_data_left_facet.clear();
  m_method_data_right_facet.clear();

  // ------------------------------
  // Prepare method data for facets
  // ------------------------------

  for (common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>::const_iterator it =
           facet_elem_type_map.cbegin();
       it != facet_elem_type_map.cend(); ++it)
  {
    common::PtrHandle<facet_method_data_type> method_data_left =
        m_method_data_left_facet.create(it.key_value());
    const interpolation::FEValues &fev = (*it.data_ptr());
    (*method_data_left).resize_variables(fev);

    common::PtrHandle<facet_method_data_type> method_data_right =
        m_method_data_right_facet.create(it.key_value());
    (*method_data_right).resize_variables(fev);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodLFImpl<MeshConfig, Physics, CellScheme, FacetScheme>::iterate_over_facets(
    const mesh::Tria<MeshConfig> &cell_topology,
    const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs,
    const typename f_space::vect_f &solution, typename f_space::vect_f &nodal_residuals,
    RDTimeUpdate &time_update, const Uint first_cell_idx, const Uint last_cell_idx)
{
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

  typedef typename interpolation::VectorMeshFunction<Real>::entry_type node_value_type;
  // typedef typename
  // interpolation::VectorMeshFunction<Real>::const_entry_type
  // const_node_value_type;

  // typedef typename result_of::dof_map_t<MeshConfig> cell_dofs_type;

  if (m_elem_residuals->nb_entries() != sol_dofs.nb_active_cells())
  {
    m_elem_residuals->resize(Physics::NEQ, sol_dofs.nb_active_cells());
  }
  m_elem_residuals->fill(0.0);

  const Uint nb_facets        = cell_topology.active_skeleton_size();
  const Uint nb_buffer_blocks = (nb_facets % m_geo_facet_blk_size == 0)
                                    ? nb_facets / m_geo_facet_blk_size
                                    : nb_facets / m_geo_facet_blk_size + 1;

  mesh::EntityDofRealign sol_facet_left_permutation, sol_facet_right_permutation;

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

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint bb = 0; bb < nb_buffer_blocks; ++bb)
  {
    const Uint first_facet_in_block = bb * m_geo_facet_blk_size;
    const Uint last_facet_in_block  = std::min((bb + 1) * m_geo_facet_blk_size - 1, nb_facets - 1);

    m_geo_cache_facets.flush();

    m_sol_cache_facets_left.flush();
    m_sol_cache_facets_right.flush();

    for (Uint f = first_facet_in_block; f <= last_facet_in_block; ++f)
    {
      // std::cout << "f = " << f << std::endl;

      const mesh::TraceIncidences facet_blk =
          cell_topology.active_skeleton_entry(MeshConfig::TDIM - 1, f);

      // Get the permutation sign of the first entity in the block
      // (position 0)
      const mesh::EntityRealignCode permutation_code = facet_blk.permutation(0).get().code();

      // If the entity on position 0 in incidence block needs to be
      // flipped, it is the entity on the rhs of this facet
      const Uint idx_left = (permutation_code.nb_flips() == 0) ? 0 : 1;

      // If the facet block has size 2, we have 2 entities forming the
      // face (i.e K+ and K-). In that case, idx_right is the 'other'
      // index (choosing from 0 and 1) than idx_left, and can be
      // determined as idx_right = (idx_left + 1) % 2 If the facet block
      // has size 1, then idx_right has to be the same as index left:
      // equal to 0
      const Uint idx_right = (facet_blk.size() == 2) ? (idx_left + 1) % 2 : 0;

      const mesh::MeshEntity geo_cell_left = geo_dofs.active_cell(facet_blk.cell_id(idx_left));
      // const mesh::MeshEntity geo_cell_right =
      // geo_dofs.active_cell(facet_blk.cell_id(idx_right));

      const mesh::MeshEntity geo_facet =
          geo_cell_left.sub_entity(FACET_DIM, facet_blk.local_id(idx_left));

      const mesh::CellTopologyView<MeshConfig> tcell_view_L =
          sol_dofs.tcell(mesh::ActiveIdx(facet_blk.cell_id(idx_left)));

      const mesh::MeshEntity sol_cell_left = sol_dofs.active_cell(facet_blk.cell_id(idx_left));
      const mesh::MeshEntity sol_facet_left =
          sol_cell_left.sub_entity(FACET_DIM, facet_blk.local_id(idx_left));

      const mesh::MeshEntity sol_cell_right = sol_dofs.active_cell(facet_blk.cell_id(idx_right));
      const mesh::MeshEntity sol_facet_right =
          sol_cell_right.sub_entity(FACET_DIM, facet_blk.local_id(idx_right));

      /*
      std::cout << "===================================================="
      << std::endl; std::cout << "Geo cell left = " << geo_cell_left <<
      std::endl; std::cout << "Geo cell right = " << geo_cell_right <<
      std::endl; std::cout << "Geo facet = " << geo_facet << std::endl;

      std::cout << "Solution cell left = " << sol_cell_left << std::endl;
      std::cout << "Solution facet left = " << sol_facet_left <<
      std::endl; std::cout << "Solution cell right = " << sol_cell_right
      << std::endl; std::cout << "Solution facet right = " <<
      sol_facet_right << std::endl;

      std::cout << "Facet has " << facet_blk.size() << " entries" <<
      std::endl; std::cout << "Permutation left = " <<
      facet_blk.permutation(idx_left).get().code()
                << "   " <<
      facet_blk.permutation(idx_left).get().type_id().as_string() <<
      std::endl; std::cout << "Permutation right = " <<
      facet_blk.permutation(idx_right).get().code()
                << "   " <<
      facet_blk.permutation(idx_right).get().type_id().as_string() <<
      std::endl; std::cout << std::endl;
      */

      const math::DenseConstMatView<Real> geo_facet_coords = loc_interpolator.transfer_coords(
          tcell_view_L.pt_set_id(FACET_DIM, facet_blk.local_id(idx_left)), geo_facet.pt_set_id(),
          tcell_view_L.coordinates(FACET_DIM, facet_blk.local_id(idx_left)));

      m_geo_cache_facets.push_back_to_buffer(
          geo_facet_coords,
          mesh::PointSetTagExt(geo_facet.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

      sol_facet_left_permutation.change_type(sol_facet_left.pt_set_id(),
                                             facet_blk.permutation(idx_left).get().code());
      sol_facet_right_permutation.change_type(sol_facet_right.pt_set_id(),
                                              facet_blk.permutation(idx_right).get().code());

      m_sol_cache_facets_left.push_back_to_buffer(
          sol_facet_left, solution,
          mesh::PointSetTagExt(sol_facet_left.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

      m_sol_cache_facets_right.push_back_to_buffer(
          sol_facet_right, solution,
          mesh::PointSetTagExt(sol_facet_right.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

    } // Loop over facets in block

    // Metric data computation for values on mesh skeleton (facets)
    m_geo_metric_facets.empty_buffer();
    m_geo_metric_facets.evaluate(m_geo_cache_facets, true);

    m_sol_metric_facets_left.empty_buffer();
    m_sol_metric_facets_right.empty_buffer();

    m_sol_metric_facets_left.evaluate(m_geo_metric_facets, m_sol_cache_facets_left, true);
    m_sol_metric_facets_right.evaluate(m_geo_metric_facets, m_sol_cache_facets_right, true);

    Uint cell_idx_in_metric = 0;

    mesh::QuadraturePermutation quad_p_left, quad_p_right;

    // Now loop over the same block of facets again and
    // - compute the facet residuals
    // - distribute the nodal contributions of facet residuals
    // - compute cell residuals by contour integration (accumulation
    //   of contributions from each cell)

    for (Uint f = first_facet_in_block; f <= last_facet_in_block; ++f)
    {
      const mesh::TraceIncidences facet_blk =
          cell_topology.active_skeleton_entry(MeshConfig::TDIM - 1, f);

      // Get the permutation sign of the first entity in the block
      // (position 0)
      const mesh::EntityRealignCode permutation_code = facet_blk.permutation(0).get().code();

      // If the entity on position 0 in incidence block needs to be
      // flipped, it is the entity on the rhs of this facet
      const Uint idx_left = (permutation_code.nb_flips() == 0) ? 0 : 1;

      // If the facet block has size 2, we have 2 entities forming the
      // face (i.e K+ and K-). In that case, idx_right is the 'other'
      // index (choosing from 0 and 1) than idx_left, and can be
      // determined as idx_right = (idx_left + 1) % 2 If the facet block
      // has size 1, then idx_right has to be the same as index left:
      // equal to 0
      const Uint idx_right = (facet_blk.size() == 2) ? (idx_left + 1) % 2 : 0;

      // const mesh::EntityPermutation p_left =
      // facet_blk.permutation(idx_left); const mesh::EntityPermutation
      // p_right = facet_blk.permutation(idx_right);

      // const mesh::MeshEntity geo_cell =
      // geo_cells.active_cell(facet_blk.cell_id(idx_left)); const
      // mesh::MeshEntity geo_facet = geo_cell.sub_entity(FACET_DIM,
      // facet_blk.local_id(idx_left));

      const mesh::MeshEntity sol_cell_left = sol_dofs.active_cell(facet_blk.cell_id(idx_left));
      const mesh::MeshEntity sol_facet_left =
          sol_cell_left.sub_entity(FACET_DIM, facet_blk.local_id(idx_left));

      const mesh::MeshEntity sol_cell_right = sol_dofs.active_cell(facet_blk.cell_id(idx_right));
      const mesh::MeshEntity sol_facet_right =
          sol_cell_right.sub_entity(FACET_DIM, facet_blk.local_id(idx_right));

      sol_facet_left_permutation.change_type(sol_facet_left.pt_set_id(),
                                             facet_blk.permutation(idx_left).get().code());
      sol_facet_right_permutation.change_type(sol_facet_right.pt_set_id(),
                                              facet_blk.permutation(idx_right).get().code());

      const mesh::PointSetTag sol_tag_left  = sol_cell_left.pt_set_id();
      const mesh::PointSetTag sol_tag_right = sol_cell_right.pt_set_id();

      /*
      const mesh::PointSetTag contour_quad_tag_left = mesh::PointSetTag(
          sol_tag_left.elem_shape(), m_sol_quadrature_order,
      m_sol_contour_quadrature_type);

      const mesh::PointSetTag contour_quad_tag_right = mesh::PointSetTag(
          sol_tag_right.elem_shape(), m_sol_quadrature_order,
      m_sol_contour_quadrature_type);
      */

      const mesh::PointSetTag contour_quad_tag_left =
          mesh::PointSetTag(sol_facet_left.pt_set_id().elem_shape(), m_sol_quadrature_order,
                            m_sol_contour_quadrature_type);

      const mesh::PointSetTag contour_quad_tag_right =
          mesh::PointSetTag(sol_facet_right.pt_set_id().elem_shape(), m_sol_quadrature_order,
                            m_sol_contour_quadrature_type);

      quad_p_left.change_type(contour_quad_tag_left, facet_blk.local_id(idx_left),
                              facet_blk.permutation(idx_left).get().code());

      quad_p_right.change_type(contour_quad_tag_right, facet_blk.local_id(idx_right),
                               facet_blk.permutation(idx_right).get().code());

      // const math::MatrixBlock<Real> facet_node_coords =
      //    m_geo_cache_facets.cell_values(cell_idx_in_metric);

      const typename facet_geo_metric_type::cellwise_metric geo_facet_met =
          m_geo_metric_facets.cellwise_values(cell_idx_in_metric);

      const math::DenseConstMatView<Real> sol_facet_nodal_values_left =
          m_sol_cache_facets_left.cell_values(cell_idx_in_metric);

      const math::DenseConstMatView<Real> sol_facet_nodal_values_right =
          m_sol_cache_facets_right.cell_values(cell_idx_in_metric);

      const mesh::PointSetTagExt key_data_left =
          mesh::PointSetTagExt(sol_facet_left.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u);

      const mesh::PointSetTagExt key_data_right =
          mesh::PointSetTagExt(sol_facet_right.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u);

      common::PtrHandle<facet_method_data_type> facet_method_data_left =
          m_method_data_left_facet.std_region_data(key_data_left);

      common::PtrHandle<facet_method_data_type> facet_method_data_right =
          m_method_data_right_facet.std_region_data(key_data_right);

      const typename facet_sol_metric_type::cellwise_metric sol_facet_met_left =
          m_sol_metric_facets_left.cellwise_values(cell_idx_in_metric);

      const typename facet_sol_metric_type::cellwise_metric sol_facet_met_right =
          m_sol_metric_facets_right.cellwise_values(cell_idx_in_metric);

      m_facet_scheme.facet_nodal_contributions(geo_facet_met, sol_facet_nodal_values_left,
                                               sol_facet_nodal_values_right, sol_facet_met_left,
                                               quad_p_left, sol_facet_met_right, quad_p_right,
                                               *facet_method_data_left, *facet_method_data_right);

      if (idx_left != idx_right)
      {
        // ---------------------------------------------------------
        // Accumulate NODAL residuals contributed by the left facet
        // ---------------------------------------------------------
        for (Uint n = 0; n < sol_facet_left.nb_vert(); ++n)
        {
          const Uint sol_vertex_id       = sol_facet_left.vertex(n);
          node_value_type nodal_residual = nodal_residuals.value(sol_vertex_id);
          const Real inv_node_volume     = 1. / nodal_dual_volume[sol_vertex_id];

          for (Uint eq = 0; eq < Physics::NEQ; ++eq)
          {
            nodal_residual[eq] +=
                inv_node_volume * (*facet_method_data_left).m_elem_node_res[n][eq];
          }
        }

        // ---------------------------------------------------------
        // Accumulate NODAL residuals contributed by the right facet
        // ---------------------------------------------------------
        for (Uint n = 0; n < sol_facet_right.nb_vert(); ++n)
        {
          const Uint sol_vertex_id       = sol_facet_right.vertex(n);
          node_value_type nodal_residual = nodal_residuals.value(sol_vertex_id);
          const Real inv_node_volume     = 1. / nodal_dual_volume[sol_vertex_id];

          for (Uint eq = 0; eq < Physics::NEQ; ++eq)
          {
            nodal_residual[eq] +=
                inv_node_volume * (*facet_method_data_right).m_elem_node_res[n][eq];
          }
        }
      } // If idx_left and idx_right are different (i.e. we are in
        // interior face)

      // --------------------------------------------------------------
      // Accumulate to CELL residuals contributed by the left facet
      // Accumulate to update coefficient contributed by the left facet
      // --------------------------------------------------------------
      node_value_type left_cell_residual = m_elem_residuals->value(sol_cell_left.idx());
      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        left_cell_residual[eq] += (*facet_method_data_left).m_flux_integral_in_elem[eq];
      }

      for (Uint n = 0; n < sol_facet_left.nb_vert(); ++n)
      {
        /*
        const Real volume_left =
        nodal_dual_volume[sol_facet_left.vertex(n)];

        time_update.accumulate_nodal_wave_speed(
            sol_facet_left.vertex(n), volume_left *
        (*facet_method_data_left).m_elem_wave_speed[n]);
        */

        time_update.accumulate_nodal_wave_speed(sol_facet_left.vertex(n),
                                                (*facet_method_data_left).m_elem_wave_speed[n]);
      }

      // ---------------------------------------------------------------
      // Accumulate to CELL residuals contributed by the right facet
      // Accumulate to update coefficient contributed by the right facet
      // This is done only in case the right facet exists (i.e. is
      // different from the left facet)
      // ---------------------------------------------------------------
      if (idx_left != idx_right)
      {
        node_value_type right_cell_residual = m_elem_residuals->value(sol_cell_right.idx());
        for (Uint eq = 0; eq < Physics::NEQ; ++eq)
        {
          right_cell_residual[eq] -= (*facet_method_data_right).m_flux_integral_in_elem[eq];
        }

        for (Uint n = 0; n < sol_facet_right.nb_vert(); ++n)
        {
          /*
          const Real volume_right =
          nodal_dual_volume[sol_facet_right.vertex(n)];

          time_update.accumulate_nodal_wave_speed(
              sol_facet_right.vertex(n),
              volume_right *
          (*facet_method_data_right).m_elem_wave_speed[n]);
          */

          time_update.accumulate_nodal_wave_speed(sol_facet_right.vertex(n),
                                                  (*facet_method_data_right).m_elem_wave_speed[n]);
        }
      } // If idx_left and idx_right are different (i.e. we are in
        // interior face)

      cell_idx_in_metric++;

    } // Loop over facets in block

  } // Loop over all facet blocks
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodLFImpl<MeshConfig, Physics, CellScheme, FacetScheme>::iterate_over_cells(
    const mesh::Tria<MeshConfig> &cell_topology,
    const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs,
    const typename f_space::vect_f &solution, typename f_space::vect_f &nodal_residuals,
    RDTimeUpdate &time_update, const Uint first_cell_idx, const Uint last_cell_idx)
{
  // **********************************************************
  // 'CELL' PHASE:
  // **********************************************************
  // Loop over all cells and distribute cell residuals that
  // were accumulated in phase 1

  typedef typename interpolation::VectorMeshFunction<Real>::entry_type node_value_type;

  const Uint nb_cells              = sol_dofs.nb_active_cells();
  const Uint nb_cell_buffer_blocks = (nb_cells % m_cell_blk_size == 0)
                                         ? nb_cells / m_cell_blk_size
                                         : nb_cells / m_cell_blk_size + 1;

  const interpolation::ScalarMeshFunction<Real> &nodal_dual_volume =
      time_update.nodal_dual_volume();

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint cell_blk = 0; cell_blk < nb_cell_buffer_blocks; ++cell_blk)
  {
    const Uint first_cell_in_blk = cell_blk * m_cell_blk_size;
    const Uint last_cell_in_blk  = std::min((cell_blk + 1) * m_cell_blk_size - 1, nb_cells - 1);

    m_geo_cache_cells.flush();
    m_sol_cache_cells.flush();

    // First loop over cells in block - push data in cache
    // TODO: maybe cells (MeshEntities) could be stored in a vector to be
    // reused in the second loop over cells in block?
    for (Uint c = first_cell_in_blk; c <= last_cell_in_blk; ++c)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));
      const mesh::MeshEntity geo_cell                     = geo_dofs.active_cell(c);
      const mesh::MeshEntity sol_cell                     = sol_dofs.active_cell(c);

      const math::DenseConstMatView<Real> geo_cell_coords = loc_interpolator.transfer_coords(
          tcell_view.pt_set_id(), geo_cell.pt_set_id(), tcell_view.coordinates());

      m_geo_cache_cells.push_back_to_buffer(
          geo_cell_coords,
          mesh::PointSetTagExt(geo_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

      m_sol_cache_cells.push_back_to_buffer(
          sol_cell, solution,
          mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

    } // Loop over cells in block

    m_geo_metric_cells.empty_buffer();
    m_geo_metric_cells.evaluate(m_geo_cache_cells, true);

    m_sol_metric_cells.empty_buffer();
    m_sol_metric_cells.evaluate(m_geo_metric_cells, m_sol_cache_cells, true);

    Uint cell_idx_in_metric = 0;

    // Second loop over cells in block - use the cellwise metric
    // to compute nodal contributions from element interiors
    for (Uint c = first_cell_in_blk; c <= last_cell_in_blk; ++c)
    {
      const mesh::MeshEntity sol_cell = sol_dofs.active_cell(c);

      const typename cell_geo_metric_type::cellwise_metric geo_cell_met =
          m_geo_metric_cells.cellwise_values(cell_idx_in_metric);

      const typename cell_sol_metric_type::cellwise_metric sol_cell_met =
          m_sol_metric_cells.cellwise_values(cell_idx_in_metric);

      const math::DenseConstMatView<Real> sol_nodal_values =
          m_sol_cache_cells.cell_values(cell_idx_in_metric);

      common::PtrHandle<cell_method_data_type> cell_method_data =
          m_method_data_cell.std_region_data(
              mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

      // Copy the cell residual obtained in phase 1 to the method data
      // variable called 'm_flux_integral_in_elem'
      node_value_type elem_residual = m_elem_residuals->value(sol_cell.idx());
      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        m_active_cell_residual[eq] = elem_residual[eq];
      }

      m_cell_scheme.compute_adv_residuals_no_flux_metric(
          sol_nodal_values, geo_cell_met, sol_cell_met, m_active_cell_residual, *cell_method_data);

      // --------------------------------------------------
      // Accumulate nodal residuals contributed by the cell
      // --------------------------------------------------
      for (Uint n = 0; n < sol_cell.nb_vert(); ++n)
      {
        const Uint sol_vertex_id      = sol_cell.vertex(n);
        node_value_type node_residual = nodal_residuals.value(sol_vertex_id);
        const Real inv_node_volume    = 1. / nodal_dual_volume[sol_vertex_id];

        for (Uint eq = 0; eq < Physics::NEQ; ++eq)
        {
          node_residual[eq] += inv_node_volume * (*cell_method_data).m_elem_node_res[n][eq];
        }

        time_update.accumulate_nodal_wave_speed(sol_cell.vertex(n),
                                                (*cell_method_data).m_elem_wave_speed[n]);
      }

      cell_idx_in_metric++;

    } // Loop over cells in block

  } // Loop over all cell blocks
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodLFImpl<MeshConfig, Physics, CellScheme, FacetScheme>::set_artificial_viscosity(
    const std::shared_ptr<interpolation::ScalarMeshFunction<Real>> &artificial_viscosity)
{
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
