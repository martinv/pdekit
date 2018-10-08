#ifndef PDEKIT_PG_RDM_Explicit_Facet_Worker_Tr_Diff_hpp
#define PDEKIT_PG_RDM_Explicit_Facet_Worker_Tr_Diff_hpp

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
class PGRDMExplicitFacetWorkerTrDiff
{
  private:
  typedef typename internal::FacetSchemeSelector<Physics, FacetScheme>::type facet_scheme_type;

  typedef typename facet_scheme_type::method_data facet_method_data_type;

  public:
  typedef interpolation::FunctionSpace<MeshConfig> f_space_cells;
  typedef interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1> f_space_facets;

  /// Default constructor
  PGRDMExplicitFacetWorkerTrDiff();

  /// Default destructor
  ~PGRDMExplicitFacetWorkerTrDiff();

  void configure_cell_contour_spaces(const typename f_space_cells::ptr &geo_cell_contour_space,
                                     const typename f_space_cells::ptr &sol_cell_contour_space);

  void configure_facet_spaces(const typename f_space_facets::ptr &geo_facet_space,
                              const typename f_space_facets::ptr &sol_facet_space);

  void iterate_over_facets(const mesh::Tria<MeshConfig> &cell_topology,
                           const typename result_of::dof_handler<MeshConfig>::type &geo_dofs,
                           const typename result_of::dof_handler<MeshConfig>::type &sol_dofs,
                           const interpolation::VectorMeshFunction<Real> &solution,
                           interpolation::VectorMeshFunction<Real> &nodal_residuals,
                           RDTimeUpdate &time_update, const Uint first_cell_idx,
                           const Uint last_cell_idx);

  private:
  enum
  {
    CELL_DIM            = MeshConfig::TDIM,
    FACET_DIM           = MeshConfig::TDIM - _1D,
    EDGE_DIM            = _1D,
    FACET_DATA_TOPO_DIM = MeshConfig::TDIM - _1D
  };

  enum
  {
    NEQ = Physics::NEQ
  };

  /// Geometry cache and metric for facets
  typedef interpolation::GeometryCache<MeshConfig> facet_geo_cache_type;
  typedef interpolation::GeometryMetric<MeshConfig, FACET_DIM> facet_geo_metric_type;

  /// Solution cache and metric for facets
  typedef interpolation::SolutionCache facet_sol_cache_type;
  typedef interpolation::SolutionSpaceMetric<MeshConfig, FACET_DIM> facet_sol_metric_type;

  /// Metric for fluxes on facets
  typedef interpolation::FluxSpaceMetric<MeshConfig, Physics, FACET_DIM> facet_flux_metric_type;

  typedef interpolation::ScalarMeshFunction<Real>::ptr scal_function_ptr;
  typedef interpolation::VectorMeshFunction<Real>::ptr vect_function_ptr;

  /// DATA

  /// Pointer to geometry space on cell facets
  typename f_space_facets::ptr m_geo_facet_space;

  /// Pointer to solution space on cell facets
  typename f_space_facets::ptr m_sol_facet_space;

  /// METRIC CONTAINERS FOR DATA ON MESH SKELETON
  /// HERE WE OPERATE ON ENTITIES WITH codim = 1

  /// Geometry cache for facets, left side of cell interface
  facet_geo_cache_type m_geo_cache_facets_left;

  /// Geometry metric for facets, left side of cell interface
  facet_geo_metric_type m_geo_metric_facets_left;

  /// Geometry cache for facets, right side of cell interface
  facet_geo_cache_type m_geo_cache_facets_right;

  /// Geometry metric for facets, right side of cell interface
  facet_geo_metric_type m_geo_metric_facets_right;

  /// Solution cache for facets
  facet_sol_cache_type m_sol_cache_facets_left;
  facet_sol_cache_type m_sol_cache_facets_right;

  /// Interpolated values and derivatives of the solution u_h
  /// on the mesh skeleton (facets)
  facet_sol_metric_type m_sol_metric_facets_left;
  facet_sol_metric_type m_sol_metric_facets_right;

  /// Interpolated values of flux on the mesh facets
  facet_flux_metric_type m_flux_metric_facets_left;
  facet_flux_metric_type m_flux_metric_facets_right;

  /// Number of blocks to fill in the geometry cache
  Uint m_geo_facet_blk_size;

  /// Number of blocks to fill in the solution cache
  Uint m_sol_facet_blk_size;

  /// Number of cells to fill the cell solution space cache
  Uint m_cell_blk_size;

  /// A map which associates to each element type in solution space a
  /// corresponding reference element for fluxes

  mesh::DataMap<RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>> m_method_data_left_facet;
  mesh::DataMap<RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>> m_method_data_right_facet;

  /// Concrete RDS scheme
  typename internal::FacetSchemeSelector<Physics, FacetScheme>::type m_facet_scheme;

  /// Quadrature order in solution space
  Uint m_sol_quadrature_order;

  /// Type of quadrature in solution space for cell interior
  PointSetID m_sol_cell_quadrature_type;

  /// Type of quadrature in solution space for cell contours
  PointSetID m_sol_contour_quadrature_type;

  /// Blending coefficient for nonlinear schemes
  scal_function_ptr m_blending_coeff;

  /// Pointer to the function that determines discontinuity sensor
  scal_function_ptr m_art_viscosity;
};

// ----------------------------------------------------------------------------
// DISCONTINUOUS RD METHOD - IMPLEMENTATION FOR
// VARIABLE BETA SCHEMES
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename FacetScheme>
PGRDMExplicitFacetWorkerTrDiff<MeshConfig, Physics,
                               FacetScheme>::PGRDMExplicitFacetWorkerTrDiff()
    : // m_geo_cell_space(nullptr), m_geo_facet_space(nullptr),
      // m_sol_cell_space(nullptr),
      m_sol_facet_space(nullptr), m_geo_facet_blk_size(0u), m_sol_facet_blk_size(0u),
      m_cell_blk_size(0u), m_sol_quadrature_order(P0),
      m_sol_cell_quadrature_type(PointSetID::Undefined),
      m_sol_contour_quadrature_type(PointSetID::Undefined), m_art_viscosity(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename FacetScheme>
PGRDMExplicitFacetWorkerTrDiff<MeshConfig, Physics, FacetScheme>::~PGRDMExplicitFacetWorkerTrDiff()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename FacetScheme>
void PGRDMExplicitFacetWorkerTrDiff<MeshConfig, Physics, FacetScheme>::
    configure_cell_contour_spaces(const typename f_space_cells::ptr &geo_cell_contour_space,
                                  const typename f_space_cells::ptr &sol_cell_contour_space)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename FacetScheme>
void PGRDMExplicitFacetWorkerTrDiff<MeshConfig, Physics, FacetScheme>::configure_facet_spaces(
    const typename f_space_facets::ptr &geo_facet_space,
    const typename f_space_facets::ptr &sol_facet_space)
{
  // ---------------------
  // GEOMETRY SPACE CONFIG
  // ---------------------
  m_geo_facet_space = geo_facet_space;

  m_geo_facet_blk_size = 150; // geo_cells.nb_skeleton_facets();

  // Allocation of cache and metric for facets
  m_geo_cache_facets_left.allocate(m_geo_facet_space->reference_elements(), m_geo_facet_blk_size);
  m_geo_cache_facets_right.allocate(m_geo_facet_space->reference_elements(), m_geo_facet_blk_size);
  m_geo_metric_facets_left.allocate_buffer(m_geo_facet_space->reference_elements(),
                                           m_geo_facet_blk_size);
  m_geo_metric_facets_right.allocate_buffer(m_geo_facet_space->reference_elements(),
                                            m_geo_facet_blk_size);

  // ---------------------
  // SOLUTION SPACE CONFIG
  // ---------------------
  m_sol_facet_space = sol_facet_space;

  m_sol_quadrature_order        = m_sol_facet_space->quadrature_order();
  m_sol_contour_quadrature_type = m_sol_facet_space->quadrature_type();

  mesh::DataMap<interpolation::FEValues> const &facet_elem_type_map =
      m_sol_facet_space->reference_elements();

  m_sol_facet_blk_size = m_geo_facet_blk_size;

  // Cache and metric for solution on facets (mesh skeleton)
  // FIXME: theoretically, this could be remove as the same data might be
  // provided by m_sol_cache_cell_contour_{left,right} and
  // m_sol_metric_cell_contour_{left,right}
  m_sol_cache_facets_left.allocate(m_sol_facet_space->reference_elements(), m_sol_facet_blk_size,
                                   NEQ);
  m_sol_cache_facets_right.allocate(m_sol_facet_space->reference_elements(), m_sol_facet_blk_size,
                                    NEQ);

  const mesh::DataMap<interpolation::FEValues> &geo_facet_type_map =
      m_geo_facet_space->reference_elements();

  const typename mesh::DataMap<interpolation::FEValues>::const_iterator first_ref_elem =
      geo_facet_type_map.cbegin();

  const mesh::PointSetTag first_geo_facet_tag = first_ref_elem.key_value().std_region_tag();
  const Uint first_geo_facet_order            = first_geo_facet_tag.poly_order();

  m_sol_metric_facets_left.allocate_buffer(m_sol_facet_space->reference_elements(),
                                           m_sol_facet_blk_size, NEQ);
  m_sol_metric_facets_right.allocate_buffer(m_sol_facet_space->reference_elements(),
                                            m_sol_facet_blk_size, NEQ);

  m_flux_metric_facets_left.allocate_buffer(SFunc::Lagrange, first_geo_facet_order,
                                            m_sol_facet_space->reference_elements(),
                                            m_sol_facet_blk_size);

  m_flux_metric_facets_right.allocate_buffer(SFunc::Lagrange, first_geo_facet_order,
                                             m_sol_facet_space->reference_elements(),
                                             m_sol_facet_blk_size);

  m_method_data_left_facet.clear();
  m_method_data_right_facet.clear();

  // ------------------------------
  // Prepare method data for facets
  // ------------------------------

  for (mesh::DataMap<interpolation::FEValues>::const_iterator it = facet_elem_type_map.cbegin();
       it != facet_elem_type_map.cend(); ++it)
  {
    common::PtrHandle<RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>> method_data_left =
        m_method_data_left_facet.create(it.key_value());
    const interpolation::FEValues &fev = (*it.data_ptr());
    (*method_data_left).resize_variables(fev);

    common::PtrHandle<RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>> method_data_right =
        m_method_data_right_facet.create(it.key_value());
    (*method_data_right).resize_variables(fev);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename FacetScheme>
void PGRDMExplicitFacetWorkerTrDiff<MeshConfig, Physics, FacetScheme>::iterate_over_facets(
    const mesh::Tria<MeshConfig> &cell_topology,
    const typename result_of::dof_handler<MeshConfig>::type &geo_dofs,
    const typename result_of::dof_handler<MeshConfig>::type &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &solution,
    interpolation::VectorMeshFunction<Real> &nodal_residuals, RDTimeUpdate &time_update,
    const Uint first_cell_idx, const Uint last_cell_idx)
{
  typedef typename interpolation::VectorMeshFunction<Real>::entry_type node_value_type;
  typedef typename interpolation::VectorMeshFunction<Real>::const_entry_type const_node_value_type;

  typedef typename result_of::dof_handler<MeshConfig>::type cell_dofs_type;

  // std::cout << "Total number of facets = " << cell_topology.skeleton_size()
  // << std::endl; std::cout << "Number of active facets = " <<
  // cell_topology.active_skeleton_size() << std::endl;

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

  for (Uint bb = 0; bb < nb_buffer_blocks; ++bb)
  {
    const Uint first_facet_in_block = bb * m_geo_facet_blk_size;
    const Uint last_facet_in_block  = std::min((bb + 1) * m_geo_facet_blk_size - 1, nb_facets - 1);

    m_geo_cache_facets_left.flush();
    m_geo_cache_facets_right.flush();

    m_sol_cache_facets_left.flush();
    m_sol_cache_facets_right.flush();

    for (Uint f = first_facet_in_block; f <= last_facet_in_block; ++f)
    {
      // std::cout << "f = " << f << std::endl;

      const mesh::TraceIncidences facet_blk =
          cell_topology.active_skeleton_entry(MeshConfig::TDIM - 1, f);

      // Get the permutation sign of the first entity in the block
      // (position 0)
      const mesh::EntityRealignCode pcode_L = facet_blk.permutation(0).get().code();

      // If the entity on position 0 in incidence block needs to be
      // flipped, it is the entity on the rhs of this facet
      const Uint idx_in_block_L = (pcode_L.nb_flips() == 0) ? 0 : 1;

      // If the facet block has size 2, we have 2 entities forming the
      // face (i.e K+ and K-). In that case, idx_right is the 'other'
      // index (choosing from 0 and 1) than idx_left, and can be
      // determined as idx_right = (idx_left + 1) % 2 If the facet block
      // has size 1, then idx_right has to be the same as index left:
      // equal to 0
      const Uint idx_in_block_R = (facet_blk.size() == 2) ? (idx_in_block_L + 1) % 2 : 0;

      // Get the permutation sign of the second entity in the block
      const mesh::EntityRealignCode pcode_R = facet_blk.permutation(idx_in_block_R).get().code();

      // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS
      // OF CELLS, NOT THEIR ACTIVE INDICES !!! IT IS FOR THIS REASON THAT
      // WE FIRST NEED TO RETRIEVE TO TOPOLOGICAL CELLS AND GET THEIR
      // ACTIVE INDICES BEFORE GETTING THE CORRECT CELLS FROM THE DOF
      // HANDLER
      const mesh::CellTopologyView<MeshConfig> tcell_L =
          cell_topology.cell(facet_blk.cell_id(idx_in_block_L));
      const mesh::CellTopologyView<MeshConfig> tcell_R =
          cell_topology.cell(facet_blk.cell_id(idx_in_block_R));

      const Uint active_cell_id_L = tcell_L.active_idx();
      const Uint active_cell_id_R = tcell_R.active_idx();

      // Get the cell and facet on the left- and right-hand side of
      // interface in GEOMETRY space
      const mesh::MeshEntity geo_cell_L = geo_dofs.active_cell(active_cell_id_L);
      const mesh::MeshEntity geo_cell_R = geo_dofs.active_cell(active_cell_id_R);

      const mesh::MeshEntity geo_facet_L =
          geo_cell_L.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_L));
      const mesh::MeshEntity geo_facet_R =
          geo_cell_R.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_R));

      const mesh::PointSetTagExt geo_key_L(geo_facet_L.pt_set_id(), P0, pcode_L.adapt_op_id(),
                                           pcode_L.local_pos_in_parent());

      const mesh::PointSetTagExt geo_key_R(geo_facet_R.pt_set_id(), P0, pcode_R.adapt_op_id(),
                                           pcode_R.local_pos_in_parent());

      // Get the cell and facet on the left- and right-hand side of
      // interface in SOLUTION space
      const mesh::MeshEntity sol_cell_L = sol_dofs.active_cell(active_cell_id_L);
      const mesh::MeshEntity sol_cell_R = sol_dofs.active_cell(active_cell_id_R);

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

      const mesh::DofCoordinates<MeshConfig> geo_facet_coords_L =
          geo_dofs.active_cell_coords(geo_facet_L);

      const mesh::DofCoordinates<MeshConfig> geo_facet_coords_R =
          geo_dofs.active_cell_coords(geo_facet_R);

      m_geo_cache_facets_left.push_back_to_buffer(geo_facet_coords_L, geo_key_L);
      m_geo_cache_facets_right.push_back_to_buffer(geo_facet_coords_R, geo_key_R);

      sol_facet_left_permutation.change_type(sol_facet_L.pt_set_id(),
                                             facet_blk.permutation(idx_in_block_L).get().code());
      sol_facet_right_permutation.change_type(sol_facet_R.pt_set_id(),
                                              facet_blk.permutation(idx_in_block_R).get().code());

      m_sol_cache_facets_left.push_back_to_buffer(sol_facet_L, solution, sol_key_L);
      m_sol_cache_facets_right.push_back_to_buffer(sol_facet_R, solution, sol_key_R);

      /*
      m_sol_cache_facets_right.push_back_to_buffer(
          sol_facet_right, sol_facet_right_permutation, solution,
          mesh::DefaultStdRegMapKey(sol_facet_right.std_region_id(), 0u,
      P0));
      */

    } // Loop over facets in block

    // Metric data computation for values on mesh skeleton (facets)
    m_geo_metric_facets_left.empty_buffer();
    m_geo_metric_facets_right.empty_buffer();

    m_geo_metric_facets_left.evaluate(m_geo_cache_facets_left, true);
    m_geo_metric_facets_right.evaluate(m_geo_cache_facets_right, true);

    m_sol_metric_facets_left.empty_buffer();
    m_sol_metric_facets_right.empty_buffer();

    m_sol_metric_facets_left.evaluate(m_geo_metric_facets_left, m_sol_cache_facets_left, true);
    m_sol_metric_facets_right.evaluate(m_geo_metric_facets_left, m_sol_cache_facets_right, true);

    m_flux_metric_facets_left.empty_buffer();
    m_flux_metric_facets_right.empty_buffer();

    m_flux_metric_facets_left.evaluate(m_geo_cache_facets_left, m_geo_metric_facets_left,
                                       m_sol_cache_facets_left, m_sol_metric_facets_left, true);
    m_flux_metric_facets_right.evaluate(m_geo_cache_facets_right, m_geo_metric_facets_right,
                                        m_sol_cache_facets_right, m_sol_metric_facets_right, true);

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
      const Uint idx_in_block_L = (permutation_code.nb_flips() == 0) ? 0 : 1;

      // If the facet block has size 2, we have 2 entities forming the
      // face (i.e K+ and K-). In that case, idx_right is the 'other'
      // index (choosing from 0 and 1) than idx_left, and can be
      // determined as idx_right = (idx_left + 1) % 2 If the facet block
      // has size 1, then idx_right has to be the same as index left:
      // equal to 0
      const Uint idx_in_block_R = (facet_blk.size() == 2) ? (idx_in_block_L + 1) % 2 : 0;

      // const mesh::EntityPermutation p_left =
      // facet_blk.permutation(idx_left); const mesh::EntityPermutation
      // p_right = facet_blk.permutation(idx_right);

      // const mesh::MeshEntity geo_cell =
      // geo_cells.active_cell(facet_blk.cell_id(idx_left)); const
      // mesh::MeshEntity geo_facet = geo_cell.sub_entity(FACET_DIM,
      // facet_blk.local_id(idx_left));

      // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS
      // OF CELLS, NOT THEIR ACTIVE INDICES !!!
      const mesh::CellTopologyView<MeshConfig> tcell_L =
          cell_topology.cell(facet_blk.cell_id(idx_in_block_L));
      const mesh::CellTopologyView<MeshConfig> tcell_R =
          cell_topology.cell(facet_blk.cell_id(idx_in_block_R));

      const Uint active_cell_id_L = tcell_L.active_idx();
      const Uint active_cell_id_R = tcell_R.active_idx();

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

      sol_facet_left_permutation.change_type(sol_facet_L.pt_set_id(),
                                             facet_blk.permutation(idx_in_block_L).get().code());
      sol_facet_right_permutation.change_type(sol_facet_R.pt_set_id(),
                                              facet_blk.permutation(idx_in_block_R).get().code());
      /*
      const mesh::PointSetTag sol_tag_left = sol_cell_L.std_region_id();
      const mesh::PointSetTag sol_tag_right = sol_cell_R.std_region_id();
      */

      /*
      const mesh::PointSetTag contour_quad_tag_left = mesh::PointSetTag(
          sol_tag_left.elem_shape(), m_sol_quadrature_order,
      m_sol_contour_quadrature_type);

      const mesh::PointSetTag contour_quad_tag_right = mesh::PointSetTag(
          sol_tag_right.elem_shape(), m_sol_quadrature_order,
      m_sol_contour_quadrature_type);
      */

      const mesh::PointSetTag contour_quad_tag_left =
          mesh::PointSetTag(sol_facet_L.pt_set_id().elem_shape(), m_sol_quadrature_order,
                            m_sol_contour_quadrature_type);

      const mesh::PointSetTag contour_quad_tag_right =
          mesh::PointSetTag(sol_facet_R.pt_set_id().elem_shape(), m_sol_quadrature_order,
                            m_sol_contour_quadrature_type);

      quad_p_left.change_type(contour_quad_tag_left, facet_blk.local_id(idx_in_block_L),
                              facet_blk.permutation(idx_in_block_L).get().code());

      quad_p_right.change_type(contour_quad_tag_right, facet_blk.local_id(idx_in_block_R),
                               facet_blk.permutation(idx_in_block_R).get().code());

      // const math::MatrixBlock<Real> facet_node_coords =
      //    m_geo_cache_facets.cell_values(cell_idx_in_metric);

      const typename facet_geo_metric_type::cellwise_metric geo_facet_met_L =
          m_geo_metric_facets_left.cellwise_values(cell_idx_in_metric);

      const typename facet_geo_metric_type::cellwise_metric geo_facet_met_R =
          m_geo_metric_facets_right.cellwise_values(cell_idx_in_metric);

      const math::DenseConstMatView<Real> sol_facet_nodal_values_left =
          m_sol_cache_facets_left.cell_values(cell_idx_in_metric);

      const math::DenseConstMatView<Real> sol_facet_nodal_values_right =
          m_sol_cache_facets_right.cell_values(cell_idx_in_metric);

      const mesh::PointSetTagExt key_data_left =
          mesh::PointSetTagExt(sol_facet_L.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u);

      const mesh::PointSetTagExt key_data_right =
          mesh::PointSetTagExt(sol_facet_R.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u);

      common::PtrHandle<RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>> facet_method_data_left =
          m_method_data_left_facet.std_region_data(key_data_left);

      common::PtrHandle<RDMethodScratchData<Physics, FACET_DATA_TOPO_DIM>> facet_method_data_right =
          m_method_data_right_facet.std_region_data(key_data_right);

      const typename facet_sol_metric_type::cellwise_metric sol_facet_met_left =
          m_sol_metric_facets_left.cellwise_values(cell_idx_in_metric);

      const typename facet_sol_metric_type::cellwise_metric sol_facet_met_right =
          m_sol_metric_facets_right.cellwise_values(cell_idx_in_metric);

      const typename facet_flux_metric_type::cellwise_metric flux_facet_met_left =
          m_flux_metric_facets_left.cellwise_values(cell_idx_in_metric);

      const typename facet_flux_metric_type::cellwise_metric flux_facet_met_right =
          m_flux_metric_facets_right.cellwise_values(cell_idx_in_metric);

      m_facet_scheme.facet_nodal_contributions(
          geo_facet_met_L, geo_facet_met_R, sol_facet_nodal_values_left,
          sol_facet_nodal_values_right, sol_facet_met_left, quad_p_left, sol_facet_met_right,
          quad_p_right, flux_facet_met_left, flux_facet_met_right, *facet_method_data_left,
          *facet_method_data_right);

      if (idx_in_block_L != idx_in_block_R)
      {
        // ---------------------------------------------------------
        // Accumulate NODAL residuals contributed by the left facet
        // ---------------------------------------------------------
        for (Uint n = 0; n < sol_facet_L.nb_vert(); ++n)
        {
          const Uint sol_vertex_id       = sol_facet_L.vertex(n);
          node_value_type nodal_residual = nodal_residuals.value(sol_vertex_id);
          const Real inv_node_volume     = 1. / nodal_dual_volume[sol_vertex_id];

          for (Uint eq = 0; eq < NEQ; ++eq)
          {
            nodal_residual[eq] +=
                inv_node_volume * (*facet_method_data_left).m_elem_node_res[n * NEQ + eq];
          }
        }

        // ---------------------------------------------------------
        // Accumulate NODAL residuals contributed by the right facet
        // ---------------------------------------------------------
        for (Uint n = 0; n < sol_facet_R.nb_vert(); ++n)
        {
          /*
          const Uint sol_vertex_id =
              sol_facet_right.vertex(sol_facet_right_permutation.get().vertex(n));
          */

          const Uint sol_vertex_id       = sol_facet_R.vertex(n);
          node_value_type nodal_residual = nodal_residuals.value(sol_vertex_id);
          const Real inv_node_volume     = 1. / nodal_dual_volume[sol_vertex_id];

          for (Uint eq = 0; eq < NEQ; ++eq)
          {
            nodal_residual[eq] +=
                inv_node_volume * (*facet_method_data_right).m_elem_node_res[n * NEQ + eq];
          }
        }
      } // If idx_left and idx_right are different (i.e. we are in
        // interior face)

      // --------------------------------------------------------------
      // Accumulate to update coefficient contributed by the left facet
      // --------------------------------------------------------------

      for (Uint n = 0; n < sol_facet_L.nb_vert(); ++n)
      {
        const Real volume_left = nodal_dual_volume[sol_facet_L.vertex(n)];

        time_update.accumulate_nodal_wave_speed(
            sol_facet_L.vertex(n), volume_left * (*facet_method_data_left).m_elem_wave_speed[n]);
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
          const Real volume_right =
          nodal_dual_volume[sol_facet_right.vertex(
              sol_facet_right_permutation.get().vertex(n))];

          time_update.accumulate_nodal_wave_speed(
              sol_facet_right.vertex(sol_facet_right_permutation.get().vertex(n)),
              volume_right *
          (*facet_method_data_right).m_elem_wave_speed[n]);
          */

          const Real volume_right = nodal_dual_volume[sol_facet_R.vertex(n)];

          time_update.accumulate_nodal_wave_speed(
              sol_facet_R.vertex(n),
              volume_right * (*facet_method_data_right).m_elem_wave_speed[n]);
        }
      } // If idx_left and idx_right are different (i.e. we are in
        // interior face)

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
