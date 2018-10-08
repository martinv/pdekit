#ifndef RDM_Weak_Boundary_Condition_hpp
#define RDM_Weak_Boundary_Condition_hpp

#include "solver/rdm/bc/RDMBCBase.hpp"
#include "solver/rdm/bc/RDMBCMetricData.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

// Base class to implement boundary conditions
// Physics ... physical model of given boundary condition
// DIM     ... space dimension of entity to which the boundary condition is
//             applied
// Example: BoundaryCondition<Dim::_3D,Dim::_2D> would be a boundary condition
// applied to a surface (2D entity) of a 3D mesh

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
class WeakBC : public RDMBCBase<MeshConfig, Physics, BcDim>
{
  public:
  using rdm_bc_base = RDMBCBase<MeshConfig, Physics, BcDim>;

  using mesh_config   = typename rdm_bc_base::mesh_config;
  using physics_type  = typename rdm_bc_base::physics_type;
  using mesh_type     = typename rdm_bc_base::mesh_type;
  using mesh_boundary = typename rdm_bc_base::mesh_boundary;
  using cell_dofs     = typename rdm_bc_base::cell_dofs;

  using f_space = typename rdm_bc_base::f_space;

  protected:
  /// TYPEDEFS

  using geo_cell_metric_type  = typename rdm_bc_base::geo_cell_metric_type;
  using sol_cell_metric_type  = typename rdm_bc_base::sol_cell_metric_type;
  using flux_cell_metric_type = typename rdm_bc_base::flux_cell_metric_type;

  public:
  /// Default constructor
  WeakBC();

  /// Constructor
  WeakBC(const std::string &name);

  /// Destructor
  ~WeakBC() override;

  /// Return the type of the boundary condition
  RDMBCType bc_type() const override;

  /// Set the function space (reference elements)
  void configure_mesh_data(const std::shared_ptr<const mesh_type> &mesh,
                           const common::PtrHandle<const cell_dofs> geo_dofs,
                           const common::PtrHandle<const cell_dofs> sol_dofs,
                           const std::string &domain_name) override;

  /// Configure the solution function
  void configure_spaces(const SFunc sf_type, const PointSetID quad_type,
                        const Uint quadrature_order) override;

  /// Apply the boundary condition
  void apply(RDTimeUpdate &time_update) override;

  /// Compute entries for global Jacobian (in implicit methods) using finite
  /// differencing
  void global_lsys_fin_diff_jacobian_entries(
      RDTimeUpdate &time_update, interpolation::SolutionCache &res_cache,
      interpolation::SolutionCache &res_cache_perturbed,
      std::vector<std::tuple<Uint, Uint, Real>> &mat_buffer,
      std::vector<std::tuple<Uint, Real>> &rhs_buffer) override;

  void global_lsys_rhs_entries(RDTimeUpdate &time_update, interpolation::SolutionCache &res_cache,
                               std::vector<std::tuple<Uint, Real>> &rhs_buffer) override;

  protected:
  /// TYPES AND TYPEDEFS
  using const_dof_iterator = typename mesh::BoundaryFacets<MeshConfig, BcDim>::const_dof_iterator;
  using const_dof_iterator_typed =
      typename mesh::BoundaryFacets<MeshConfig, BcDim>::const_dof_iterator_typed;
  using sol_elem_range_typed = common::IteratorRange<const_dof_iterator_typed>;
  const std::vector<sol_elem_range_typed> &sol_ranges() const;
  using bc_metric_data_t = RDMBCMetricData<MeshConfig, Physics, BcDim>;

  /// METHODS
  bc_metric_data_t const &metric_data() const;

  /// Compute the nodal residuals  on one cell
  void compute_on_element_weak(const Uint idx_in_metric,
                               geo_cell_metric_type const &cell_geo_metric,
                               sol_cell_metric_type const &cell_sol_metric,
                               flux_cell_metric_type const &cell_flux_metric);

  /// DATA

  /// Corrective residual for each node
  /// Dimensions: [NEQ x nb_nodes]
  common::DataMap<mesh::PointSetTagExt, math::DenseDVec<Real>> m_bface_residual;

  /// Update coefficient for each node in boundary element
  common::DataMap<mesh::PointSetTagExt, math::DenseDVec<Real>> m_bface_elem_update_coeff;

  /// The weak correction ('residual') integrated over one facet
  typename Physics::Properties::FluxV m_total_facet_residual;

  /// Sum of partial ('nodal', dof-related) residuals. This
  /// should be equal to m_total_facet_residual. If not,
  /// we need to scale the nodal residuals to satisfy this condition
  /// Otherwise the BC would not be conservative
  typename Physics::Properties::FluxV m_sum_dof_residuals;

  private:
  /// Apply the boundary condition given cached data and metric
  /// This supposes that the the geometry/solution cache and metric in the
  /// data have been processed a-priori!
  void compute_residuals(RDMBCMetricData<MeshConfig, Physics, BcDim> const &data,
                         const_dof_iterator_typed const &sol_iter_begin,
                         const_dof_iterator_typed const &sol_iter_end, RDTimeUpdate &time_update,
                         interpolation::SolutionCache &res_cache);

  /// Ranges of elements (for each element type separately)
  /// in the solution dof container
  std::vector<sol_elem_range_typed> m_sol_ranges;

  std::unique_ptr<bc_metric_data_t> m_metric_data;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::WeakBC()
    : RDMBCBase<MeshConfig, Physics, BcDim>(), m_metric_data(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::WeakBC(const std::string &name)
    : RDMBCBase<MeshConfig, Physics, BcDim>(name), m_metric_data(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::~WeakBC()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
RDMBCType WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::bc_type() const
{
  return BC_TYPE_WEAK;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
void WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::configure_mesh_data(
    const std::shared_ptr<const mesh_type> &mesh, const common::PtrHandle<const cell_dofs> geo_dofs,
    const common::PtrHandle<const cell_dofs> sol_dofs, const std::string &domain_name)
{
  // First call the method of the parent
  rdm_bc_base::configure_mesh_data(mesh, geo_dofs, sol_dofs, domain_name);

  mesh::BoundaryFacets<MeshConfig, BcDim> const &boundary =
      *(this->mesh()->all_boundaries().domain(this->domain_name()));

  using cell_dofs_type                = typename result_of::dof_map_t<MeshConfig>;
  cell_dofs_type const &sol_cell_dofs = *this->sol_dofs();

  boundary.all_bdry_dof_ranges(sol_cell_dofs, m_sol_ranges);

  for (const sol_elem_range_typed &sol_elem_range : m_sol_ranges)
  {
    const mesh::MeshEntity sol_cell = (*sol_elem_range.begin()).mesh_entity();

    common::PtrHandle<math::DenseDVec<Real>> residual_vec = m_bface_residual.create(
        mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));
    (*residual_vec).resize(sol_cell.nb_vert() * Physics::NEQ);

    common::PtrHandle<math::DenseDVec<Real>> elem_update_coeff = m_bface_elem_update_coeff.create(
        mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));
    (*elem_update_coeff).resize(sol_cell.nb_vert());
    (*elem_update_coeff).fill(0.0);
  }

  return;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
void WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::configure_spaces(const SFunc sf_type,
                                                                      const PointSetID quad_type,
                                                                      const Uint quadrature_order)
{
  // First call the method of the parent
  rdm_bc_base::configure_spaces(sf_type, quad_type, quadrature_order);

  const mesh::MeshBoundarySet<MeshConfig> &mesh_boundaries = this->mesh()->all_boundaries();
  const mesh::BoundaryFacets<MeshConfig, BcDim> &boundary =
      *(mesh_boundaries.domain(this->domain_name()));

  if (!m_metric_data)
  {
    const mesh::CellTopologyView<MeshConfig> first_cell_geo_cell = *(*this->mesh()).cbegin_cells();
    const Uint geo_cell_order = first_cell_geo_cell.pt_set_id().poly_order();

    m_metric_data = std::unique_ptr<bc_metric_data_t>(new bc_metric_data_t());
    m_metric_data->allocate_cache(geo_cell_order, *rdm_bc_base::geo_space(),
                                  *rdm_bc_base::sol_space(), boundary.nb_active_cells());
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
void WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::apply(RDTimeUpdate &time_update)
{
  interpolation::VectorMeshFunction<Real> &solution  = *(this->solution());
  interpolation::VectorMeshFunction<Real> &residuals = *(this->residuals());

  m_metric_data->m_geo_cache.flush();
  m_metric_data->m_geo_metric.empty_buffer();

  m_metric_data->fill_geo_cache((*this->mesh()), m_sol_ranges, (*rdm_bc_base::geo_space()));
  m_metric_data->m_geo_metric.evaluate(m_metric_data->m_geo_cache,
                                       interpolation::RebuildMetricIndex{true});

  // typedef typename result_of::dof_map_t<MeshConfig> cell_dofs_type;
  using node_value_type = typename interpolation::VectorMeshFunction<Real>::entry_type;
  // typedef typename
  // interpolation::VectorMeshFunction<Real>::const_entry_type
  // const_node_value_type;

  m_metric_data->m_sol_metric.empty_buffer();
  m_metric_data->m_sol_cache.flush();
  m_metric_data->m_flux_metric.empty_buffer();

  const_dof_iterator_typed sol_elem_it;

  // Fill solution cache
  for (sol_elem_range_typed const &sol_elem_range : m_sol_ranges)
  {
    for (sol_elem_it = sol_elem_range.begin(); sol_elem_it != sol_elem_range.end(); ++sol_elem_it)
    {
      const mesh::MeshEntity sol_cell = sol_elem_it->mesh_entity();
      m_metric_data->m_sol_cache.push_back_to_buffer(
          sol_cell, solution,
          mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));
    }
  }

  // Evaluate solution metric
  m_metric_data->m_sol_metric.evaluate(m_metric_data->m_geo_metric, m_metric_data->m_sol_cache,
                                       interpolation::ComputeMetricDerivs{true},
                                       interpolation::RebuildMetricIndex{true});
  m_metric_data->m_flux_metric.evaluate(m_metric_data->m_geo_cache, m_metric_data->m_geo_metric,
                                        m_metric_data->m_sol_cache, m_metric_data->m_sol_metric,
                                        interpolation::RebuildMetricIndex{true});

  // Loop over the computed metric data and evaluate boundary residuals
  Uint idx_in_metric = 0;

  const interpolation::ScalarMeshFunction<Real> &dual_nodal_volume =
      time_update.nodal_dual_volume();

  // Fill solution cache

  for (sol_elem_range_typed const &sol_elem_range : m_sol_ranges)
  {
    // Don't use a reference here. The following does compile but segfaults
    // with clang: const mesh::MeshEntity &sol_cell =
    // *(sol_elem_range.begin()); Use this instead:
    const mesh::MeshEntity sol_cell = (*sol_elem_range.begin()).mesh_entity();

    math::DenseDVec<Real> &bface_res = *(m_bface_residual.std_region_data(
        mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u)));

    math::DenseDVec<Real> &bface_update_coeff_vector = *(m_bface_elem_update_coeff.std_region_data(
        mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u)));

    for (sol_elem_it = sol_elem_range.begin(); sol_elem_it != sol_elem_range.end(); ++sol_elem_it)
    {
      const mesh::MeshEntity sol_cell = sol_elem_it->mesh_entity();

      bface_res.fill(0.0);

      bface_update_coeff_vector.fill(0.0);

      geo_cell_metric_type const cell_geo_metric =
          m_metric_data->m_geo_metric.cellwise_values(idx_in_metric);
      sol_cell_metric_type const cell_sol_metric =
          m_metric_data->m_sol_metric.cellwise_values(idx_in_metric);
      flux_cell_metric_type const cell_flux_metric =
          m_metric_data->m_flux_metric.cellwise_values(idx_in_metric);

      // Compute corrective residual on one boundary element
      compute_on_element_weak(idx_in_metric, cell_geo_metric, cell_sol_metric, cell_flux_metric);

      /*
      const mesh::IncidencePair idx_in_parent =
      sol_elem_it.idx_in_parent(); const Real parent_cell_volume =
      time_update.cell_volume(idx_in_parent.cell_idx);

      for (Uint n = 0; n < bface_update_coeff_vector.size(); ++n)
      {
        bface_update_coeff_vector[n] *= parent_cell_volume;
      }
      */

      for (Uint n = 0; n < sol_cell.nb_vert(); ++n)
      {
        node_value_type node_residual = residuals.value(sol_cell.vertex(n));

        // (*rdm_bc_base::m_rd_time_update)[sol_cell.vertex(n)] +=
        // bface_update_coeff_vector[n];
        // time_update.accumulate_nodal_wave_speed(sol_cell.vertex(n),
        // bface_update_coeff_vector[n]);
        const Real nodal_volume = dual_nodal_volume[sol_cell.vertex(n)];

        // IMPORTANT: THE DIMENSION OF UPDATE COEFFICIENT IS k+ = |S| *
        // lambda, WHERE lambda IS THE MAXIMUM EIGENVALUE
        // bface_update_coeff_vector[n] ONLY CONTAINS THE VALUES OF
        // lambda, HENCE THE WAVE SPEED, WHICH SHOULD HOLD [(k+)*|S|]
        // MUST BE COMPUTED AS (k+) * |S|^2

        time_update.accumulate_nodal_wave_speed(sol_cell.vertex(n),
                                                nodal_volume * bface_update_coeff_vector[n]);

        // ON THE OTHER HAND, IN THE MAIN SOLVER, THE GOVERNING EQUATION
        // ISGIVEN BY du/dt = 1./|S| * RDM_res, WITH RDM_res BEING THE
        // RDS ACCUMULATED RESIDUALS. HENCE THE BOUNDARY FACE RESIDUALS
        // (bface_res[n,eq] MUST BE DIVIDED BY NODAL VOLUME

        const Real inv_nodal_volume = 1. / nodal_volume;

        for (Uint eq = 0; eq < Physics::NEQ; ++eq)
        {
          node_residual[eq] += inv_nodal_volume * bface_res[n * Physics::NEQ + eq];
        }
      }

      idx_in_metric++;

    } // Loop over one type of boundary elements (solution space)

  } // Loop over boundary element ranges
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
void WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::compute_residuals(
    RDMBCMetricData<MeshConfig, Physics, BcDim> const &data,
    typename mesh::BoundaryFacets<MeshConfig, BcDim>::const_dof_iterator_typed const
        &sol_iter_begin,
    typename mesh::BoundaryFacets<MeshConfig, BcDim>::const_dof_iterator_typed const &sol_iter_end,
    RDTimeUpdate &time_update, interpolation::SolutionCache &res_cache)
{
  // const_dof_iterator bcell_geo_iter =
  const_dof_iterator_typed bcell_sol_iter = sol_iter_begin;

  // Since we are iterating over boundary cells of the same type, the
  // 'bface_res' data will always be the same (resized to the same number of
  // local cell nodes)
  const mesh::MeshEntity sol_cell  = bcell_sol_iter->mesh_entity();
  math::DenseDVec<Real> &bface_res = *(m_bface_residual.std_region_data(
      mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u)));

  // Since we are iterating over boundary cells of the same type, the
  // 'bface_update_coeff_vector' will always be the same (resized to the same
  // number of local cell nodes)
  math::DenseDVec<Real> &bface_update_coeff_vector = *(m_bface_elem_update_coeff.std_region_data(
      mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u)));

  const interpolation::ScalarMeshFunction<Real> &dual_nodal_volume =
      time_update.nodal_dual_volume();

  // This is a loop over all boundaries in the SOLUTION space mesh
  Uint idx_in_metric = 0;
  for (; bcell_sol_iter != sol_iter_end; ++bcell_sol_iter)
  {
    const mesh::MeshEntity sol_cell = bcell_sol_iter->mesh_entity();

    bface_res.fill(0.0);

    bface_update_coeff_vector.fill(0.0);

    geo_cell_metric_type const cell_geo_metric = data.m_geo_metric.cellwise_values(idx_in_metric);
    sol_cell_metric_type const cell_sol_metric = data.m_sol_metric.cellwise_values(idx_in_metric);
    flux_cell_metric_type const cell_flux_metric =
        data.m_flux_metric.cellwise_values(idx_in_metric);

    // Compute corrective residual on one boundary element
    // This call fills the data in 'bface_res' and
    // 'bface_update_coeff_vector'
    compute_on_element_weak(idx_in_metric, cell_geo_metric, cell_sol_metric, cell_flux_metric);

    for (Uint n = 0; n < sol_cell.nb_vert(); ++n)
    {
      // (*rdm_bc_base::m_rd_time_update)[sol_cell.vertex(n)] +=
      // bface_update_coeff_vector[n];
      // time_update.accumulate_nodal_wave_speed(sol_cell.vertex(n),
      // bface_update_coeff_vector[n]);
      const Real nodal_volume = dual_nodal_volume[sol_cell.vertex(n)];

      // IMPORTANT: THE DIMENSION OF UPDATE COEFFICIENT IS k+ = |S| *
      // lambda, WHERE lambda IS THE MAXIMUM EIGENVALUE
      // bface_update_coeff_vector[n] ONLY CONTAINS THE VALUES OF lambda,
      // HENCE THE WAVE SPEED, WHICH SHOULD HOLD [(k+)*|S|] MUST BE
      // COMPUTED AS (k+) * |S|^2
      //

      time_update.accumulate_nodal_wave_speed(sol_cell.vertex(n),
                                              nodal_volume * bface_update_coeff_vector[n]);
      // ON THE OTHER HAND, IN THE MAIN SOLVER, THE GOVERNING EQUATION IS
      // GIVEN BY du/dt = 1./|S| * RDM_res, WITH RDM_res BEING THE RDS
      // ACCUMULATED RESIDUALS. HENCE THE BOUNDARY FACE RESIDUALS
      // (bface_res[n,eq] MUST BE DIVIDED BY NODAL VOLUME

      const Real inv_nodal_volume = 1. / nodal_volume;

      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        bface_res[n * Physics::NEQ + eq] *= inv_nodal_volume;
      }
    }

    res_cache.push_vec_to_buffer(
        sol_cell, bface_res,
        mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

    idx_in_metric++;
  } // End of loop over boundary cells
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
void WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::global_lsys_fin_diff_jacobian_entries(
    RDTimeUpdate &time_update, interpolation::SolutionCache &residual_cache,
    interpolation::SolutionCache &residual_cache_perturbed,
    std::vector<std::tuple<Uint, Uint, Real>> &mat_buffer,
    std::vector<std::tuple<Uint, Real>> &rhs_buffer)
{
  interpolation::VectorMeshFunction<Real> &solution = *(this->solution());

  const_dof_iterator_typed sol_bdry_iter;

  m_metric_data->m_geo_cache.flush();
  m_metric_data->m_sol_cache.flush();

  residual_cache.flush();

  const auto geo_sf_generator   = rdm_bc_base::geo_space()->sf_generator();
  const auto geo_quad_generator = rdm_bc_base::geo_space()->quad_generator();

  for (sol_elem_range_typed const &bdry_elem_range : m_sol_ranges)
  {
    // ------------------------------------------
    // Weak BC processing
    // Ia) Fill the geometry cache
    // Ib) Fill solution cache and optionally source cache
    // ------------------------------------------

    for (sol_bdry_iter = bdry_elem_range.begin(); sol_bdry_iter != bdry_elem_range.end();
         ++sol_bdry_iter)
    {
      const mesh::PointSetTag geo_cell_type = sol_bdry_iter->geo_pt_set_id();
      const mesh::CellGeometry<MeshConfig::GDIM> geo_bdry_elem_coords =
          sol_bdry_iter->cell_geometry();
      const mesh::MeshEntity sol_bdry_elem = sol_bdry_iter->mesh_entity();

      const ElemShape cell_shape = geo_cell_type.elem_shape();
      const Uint geo_order       = geo_cell_type.poly_order();

      const mesh::PointSetTagExt geo_pt_set_ext(geo_cell_type, P0, mesh::CellTransform::NO_TRANS,
                                                0u);
      const mesh::sf::SFTag geo_sf         = geo_sf_generator(cell_shape, geo_order);
      const mesh::PointSetTag geo_quad_set = geo_quad_generator(cell_shape, geo_order);
      const mesh::PointSetTagExt geo_quad_set_ext(geo_quad_set, P0, mesh::CellTransform::NO_TRANS,
                                                  0u);
      const mesh::DiscreteElemKey geo_key(geo_pt_set_ext, geo_sf, geo_quad_set_ext);

      m_metric_data->m_geo_cache.push_back_to_buffer(geo_bdry_elem_coords, geo_key);

      m_metric_data->m_sol_cache.push_back_to_buffer(
          sol_bdry_elem, solution,
          mesh::PointSetTagExt(sol_bdry_elem.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));
    }

    m_metric_data->m_geo_metric.empty_buffer();
    m_metric_data->m_sol_metric.empty_buffer();

    m_metric_data->m_geo_metric.evaluate(m_metric_data->m_geo_cache,
                                         interpolation::RebuildMetricIndex{true});
    m_metric_data->m_sol_metric.evaluate(m_metric_data->m_geo_metric, m_metric_data->m_sol_cache,
                                         interpolation::ComputeMetricDerivs{true},
                                         interpolation::RebuildMetricIndex{true});

    m_metric_data->m_flux_metric.empty_buffer();
    m_metric_data->m_flux_metric.evaluate(m_metric_data->m_geo_cache, m_metric_data->m_geo_metric,
                                          m_metric_data->m_sol_cache, m_metric_data->m_sol_metric,
                                          interpolation::RebuildMetricIndex{true});

    // ------------------------------------------
    // Weak BC processing
    // II) Evaluate metric terms
    // ------------------------------------------

    compute_residuals(*m_metric_data, bdry_elem_range.begin(), bdry_elem_range.end(), time_update,
                      residual_cache);

    // --------------------------------------------------
    // Weak BC processing
    // III) Compute the residuals and store them in cache
    // --------------------------------------------------

    Uint cell_idx_in_metric = 0u;
    for (sol_bdry_iter = bdry_elem_range.begin(); sol_bdry_iter != bdry_elem_range.end();
         ++sol_bdry_iter)
    {
      const mesh::MeshEntity sol_bdry_elem = sol_bdry_iter->mesh_entity();
      const math::DenseConstMatView<Real> elem_nodal_residuals =
          residual_cache.cell_values(cell_idx_in_metric);

      for (Uint n = 0; n < sol_bdry_elem.nb_vert(); ++n)
      {
        for (Uint eq = 0; eq < Physics::NEQ; ++eq)
        {
          const Uint dof_idx = sol_bdry_elem.vertex(n) * Physics::NEQ + eq;
          rhs_buffer.push_back(std::tuple<Uint, Real>(dof_idx, -elem_nodal_residuals(n, eq)));
          //(*m_rhs).add_value(dof_idx, -elem_nodal_residuals(n, eq),
          // 0);
        }
      }

      cell_idx_in_metric++;
    } // Loop over boundary elements of one type

    // ---------------------------------------------------
    // Weak BC processing
    // IV) Compute numerical Jacobians and accumulate them
    //     to system matrix
    // ---------------------------------------------------

    const common::Range1D<Uint> elem_range(0, cell_idx_in_metric - 1);

    const mesh::MeshEntity first_bdry_elem = (*bdry_elem_range.begin()).mesh_entity();
    const Uint nb_dof_per_elem             = first_bdry_elem.nb_vert();
    const mesh::PointSetTag cell_type_id   = first_bdry_elem.pt_set_id();

    for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_elem; ++i_dof_in_elem)
    {
      for (Uint comp_u = 0; comp_u < Physics::NEQ; ++comp_u)
      {
        /*
        m_data->m_sol_cache.perturb_values(
            mesh::PointSetTagExt(cell_type_id, P0,
        mesh::CellTransform::DO_NOTHING, 0u), i_dof_in_elem, comp_u);
        */

        m_metric_data->m_sol_cache.perturb_values(elem_range, i_dof_in_elem, comp_u);

        /*
        const math::DenseConstVecView<Real> unperturbed_node_in_cell =
            m_data->m_sol_cache.unperturbed_values(
                mesh::PointSetTagExt(cell_type_id, P0,
        mesh::CellTransform::DO_NOTHING, 0u));
        */

        const math::DenseConstVecView<Real> unperturbed_node_in_cell =
            m_metric_data->m_sol_cache.unperturbed_values(elem_range);

        m_metric_data->m_sol_metric.empty_buffer();
        m_metric_data->m_sol_metric.evaluate(
            m_metric_data->m_geo_metric, m_metric_data->m_sol_cache,
            interpolation::ComputeMetricDerivs{true}, interpolation::RebuildMetricIndex{true});

        m_metric_data->m_flux_metric.empty_buffer();
        m_metric_data->m_flux_metric.evaluate(
            m_metric_data->m_geo_cache, m_metric_data->m_geo_metric, m_metric_data->m_sol_cache,
            m_metric_data->m_sol_metric, interpolation::RebuildMetricIndex{true});

        residual_cache_perturbed.flush();
        compute_residuals(*m_metric_data, bdry_elem_range.begin(), bdry_elem_range.end(),
                          time_update, residual_cache_perturbed);

        cell_idx_in_metric = 0u;
        for (sol_bdry_iter = bdry_elem_range.begin(); sol_bdry_iter != bdry_elem_range.end();
             ++sol_bdry_iter)
        {
          const mesh::MeshEntity sol_bdry_elem = sol_bdry_iter->mesh_entity();

          const math::DenseConstMatView<Real> sol_nodal_values =
              m_metric_data->m_sol_cache.cell_values(cell_idx_in_metric);

          const math::DenseConstMatView<Real> elem_nodal_residuals =
              residual_cache.cell_values(cell_idx_in_metric);

          const math::DenseConstMatView<Real> elem_nodal_residuals_perturbed =
              residual_cache_perturbed.cell_values(cell_idx_in_metric);

          const Real inv_du = 1. / (sol_nodal_values(i_dof_in_elem, comp_u) -
                                    unperturbed_node_in_cell[cell_idx_in_metric]);

          // Compute finite difference derivatives of residuals with
          // respect to du (which is the perturbation of component
          // comp_u of solution in node i_dof_in_elem)
          for (Uint v = 0; v < sol_bdry_elem.nb_vert(); ++v)
          {
            for (Uint comp = 0; comp < Physics::NEQ; ++comp)
            {
              // For each component
              // When v == i_dof_in_elem, we are computing
              // derivatives of residuals corresponding to node v
              // with respect to components of u associated to v
              // (diagonal blocks) in the system matrix
              const Uint global_row_idx = Physics::NEQ * sol_bdry_elem.vertex(v) + comp;

              // if (!is_Dirichlet_node[global_row_idx])
              // {
              const Uint global_col_idx =
                  Physics::NEQ * sol_bdry_elem.vertex(i_dof_in_elem) + comp_u;
              const Real residual_finite_diff = inv_du * (elem_nodal_residuals_perturbed(v, comp) -
                                                          elem_nodal_residuals(v, comp));
              mat_buffer.push_back(std::tuple<Uint, Uint, Real>(global_row_idx, global_col_idx,
                                                                residual_finite_diff));
            }
          }

          cell_idx_in_metric++;

        } // Loop over boundary elements of one type

        /*
        m_data->m_sol_cache.remove_perturbation(
            mesh::PointSetTagExt(cell_type_id, P0,
        mesh::CellTransform::DO_NOTHING, 0u));
        */

        m_metric_data->m_sol_cache.remove_perturbation(elem_range);

      } // Loop over equations components

    } // Loop over local dofs in one element

  } // Loop over all element types present in one boundary
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
void WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::global_lsys_rhs_entries(
    RDTimeUpdate &time_update, interpolation::SolutionCache &res_cache,
    std::vector<std::tuple<Uint, Real>> &rhs_buffer)
{
  interpolation::VectorMeshFunction<Real> &solution = *(this->solution());

  const_dof_iterator_typed sol_bdry_iter;

  m_metric_data->m_geo_cache.flush();
  m_metric_data->m_sol_cache.flush();
  res_cache.flush();

  const auto geo_sf_generator   = rdm_bc_base::geo_space()->sf_generator();
  const auto geo_quad_generator = rdm_bc_base::geo_space()->quad_generator();

  for (sol_elem_range_typed const &bdry_elem_range : m_sol_ranges)
  {
    // ------------------------------------------
    // Weak BC processing
    // Ia) Fill the geometry cache
    // Ib) Fill solution cache and optionally source cache
    // ------------------------------------------

    for (sol_bdry_iter = bdry_elem_range.begin(); sol_bdry_iter != bdry_elem_range.end();
         ++sol_bdry_iter)
    {
      const mesh::PointSetTag geo_cell_type = sol_bdry_iter->geo_pt_set_id();
      const mesh::CellGeometry<MeshConfig::GDIM> geo_bdry_elem_coords =
          sol_bdry_iter->cell_geometry();
      const mesh::MeshEntity sol_bdry_elem = sol_bdry_iter->mesh_entity();

      const ElemShape cell_shape = geo_cell_type.elem_shape();
      const Uint geo_order       = geo_cell_type.poly_order();

      const mesh::PointSetTagExt geo_pt_set_ext(geo_cell_type, P0, mesh::CellTransform::NO_TRANS,
                                                0u);
      const mesh::sf::SFTag geo_sf         = geo_sf_generator(cell_shape, geo_order);
      const mesh::PointSetTag geo_quad_set = geo_quad_generator(cell_shape, geo_order);
      const mesh::PointSetTagExt geo_quad_set_ext(geo_quad_set, P0, mesh::CellTransform::NO_TRANS,
                                                  0u);
      const mesh::DiscreteElemKey geo_key(geo_pt_set_ext, geo_sf, geo_quad_set_ext);

      m_metric_data->m_geo_cache.push_back_to_buffer(geo_bdry_elem_coords, geo_key);

      m_metric_data->m_sol_cache.push_back_to_buffer(
          sol_bdry_elem, solution,
          mesh::PointSetTagExt(sol_bdry_elem.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));
    }

    m_metric_data->m_geo_metric.empty_buffer();
    m_metric_data->m_sol_metric.empty_buffer();

    m_metric_data->m_geo_metric.evaluate(m_metric_data->m_geo_cache,
                                         interpolation::RebuildMetricIndex{true});
    m_metric_data->m_sol_metric.evaluate(m_metric_data->m_geo_metric, m_metric_data->m_sol_cache,
                                         interpolation::ComputeMetricDerivs{true},
                                         interpolation::RebuildMetricIndex{true});

    m_metric_data->m_flux_metric.empty_buffer();
    m_metric_data->m_flux_metric.evaluate(m_metric_data->m_geo_cache, m_metric_data->m_geo_metric,
                                          m_metric_data->m_sol_cache, m_metric_data->m_sol_metric,
                                          interpolation::RebuildMetricIndex{true});

    // ------------------------------------------
    // Weak BC processing
    // II) Evaluate metric terms
    // ------------------------------------------

    compute_residuals(*m_metric_data, bdry_elem_range.begin(), bdry_elem_range.end(), time_update,
                      res_cache);

    // --------------------------------------------------
    // Weak BC processing
    // III) Compute the residuals and store them in cache
    // --------------------------------------------------

    Uint cell_idx_in_metric = 0u;
    for (sol_bdry_iter = bdry_elem_range.begin(); sol_bdry_iter != bdry_elem_range.end();
         ++sol_bdry_iter)
    {
      const mesh::MeshEntity sol_bdry_elem = sol_bdry_iter->mesh_entity();
      const math::DenseConstMatView<Real> elem_nodal_residuals =
          res_cache.cell_values(cell_idx_in_metric);

      for (Uint n = 0; n < sol_bdry_elem.nb_vert(); ++n)
      {
        for (Uint eq = 0; eq < Physics::NEQ; ++eq)
        {
          const Uint dof_idx = sol_bdry_elem.vertex(n) * Physics::NEQ + eq;
          // (*m_rhs).add_value(dof_idx, -elem_nodal_residuals(n, eq),
          // 0);
          rhs_buffer.push_back(std::tuple<Uint, Real>(dof_idx, -elem_nodal_residuals(n, eq)));
        }
      }

      cell_idx_in_metric++;
    } // Loop over boundary elements of one type

  } // Loop over all element types present in one boundary
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
const std::vector<typename WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::sol_elem_range_typed>
    &WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::sol_ranges() const
{
  return m_sol_ranges;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
typename WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::bc_metric_data_t const &WeakBC<
    MeshConfig, Physics, BcDim, ConcreteBC>::metric_data() const
{
  return *m_metric_data;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim, typename ConcreteBC>
void WeakBC<MeshConfig, Physics, BcDim, ConcreteBC>::compute_on_element_weak(
    const Uint idx_in_metric, geo_cell_metric_type const &cell_geo_metric,
    sol_cell_metric_type const &cell_sol_metric, flux_cell_metric_type const &cell_flux_metric)
{
  static_cast<ConcreteBC *>(this)->compute_on_element_weak(idx_in_metric, cell_geo_metric,
                                                           cell_sol_metric, cell_flux_metric);

  math::DenseDVec<Real> &bface_res =
      *(m_bface_residual.std_region_data(cell_sol_metric.std_region_type()));

  // Rescale the residuals assigned to each dof/node so that the their sum
  // is equal to the total integrated residual over the boundary face.
  // This is necessary for bases that do not form partition of unity
  m_sum_dof_residuals.fill(0.0);

  const Uint nb_local_dof = cell_sol_metric.nb_dof_in_cell();
  ;
  for (Uint n = 0; n < nb_local_dof; ++n)
  {
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      m_sum_dof_residuals[eq] += bface_res[n * Physics::NEQ + eq];
    }
  }

  for (Uint eq = 0; eq < Physics::NEQ; ++eq)
  {
    if (std::abs(m_sum_dof_residuals[eq]) > 1.e-14)
    {
      const Real scale_factor = m_total_facet_residual[eq] / m_sum_dof_residuals[eq];
      for (Uint n = 0; n < nb_local_dof; ++n)
      {
        bface_res[n * Physics::NEQ + eq] *= scale_factor;
      }
    }
  }
}

// -----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
