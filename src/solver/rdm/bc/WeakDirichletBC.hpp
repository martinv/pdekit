#ifndef RDM_Weak_Dirichlet_BC_hpp
#define RDM_Weak_Dirichlet_BC_hpp

#include "solver/rdm/bc/RDMBCProjectionData.hpp"
#include "solver/rdm/bc/WeakBC.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, Uint BcDim = Physics::DIM - 1>
class WeakDirichletBC
    : public WeakBC<MeshConfig, Physics, BcDim, WeakDirichletBC<MeshConfig, Physics, BcDim>>
{
  public:
  using base = WeakBC<MeshConfig, Physics, BcDim, WeakDirichletBC<MeshConfig, Physics, BcDim>>;
  using mesh_config           = typename base::mesh_config;
  using physics_type          = typename base::physics_type;
  using mesh_type             = typename base::mesh_type;
  using cell_dofs             = typename base::cell_dofs;
  using geo_cell_metric_type  = typename base::geo_cell_metric_type;
  using sol_cell_metric_type  = typename base::sol_cell_metric_type;
  using flux_cell_metric_type = typename base::flux_cell_metric_type;

  using f_space = typename base::f_space;

  /// Default constructor
  WeakDirichletBC();

  /// Constructor
  WeakDirichletBC(const std::string &name);

  /// Destructor
  ~WeakDirichletBC() override;

  /// Set the function space (reference elements)
  void configure_mesh_data(const std::shared_ptr<const mesh_type> &mesh,
                           const common::PtrHandle<const cell_dofs> geo_dofs,
                           const common::PtrHandle<const cell_dofs> sol_dofs,
                           const std::string &domain_name) override;

  /// Configure the solution function
  void configure_spaces(const SFunc sf_type, const PointSetID quad_type,
                        const Uint quadrature_order) override;

  /// Set an analytical expression if needed for the bc
  void set_expression(Real (*expr_ptr)(
      const math::DenseConstVecView<Real> &,
      const interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint)) override;

  /// Compute the corrective residuals on one element
  void compute_on_element_weak(const Uint idx_in_metric,
                               geo_cell_metric_type const &cell_geo_metric,
                               sol_cell_metric_type const &cell_sol_metric,
                               flux_cell_metric_type const &cell_flux_metric);

  /// Apply the boundary condition
  // void apply();

  protected:
  using bc_projection_data_t     = RDMBCProjectionData<MeshConfig, BcDim>;
  using const_dof_iterator       = typename base::const_dof_iterator;
  using const_dof_iterator_typed = typename base::const_dof_iterator_typed;

  /// Pointer to function which computes the value of Dirichlet BC
  /// @param math::VectorBlock<Real> ... vector of point coordinates
  /// @param inteprolation::FunctionSpace<MeshConfig>::Function_t::column_type
  /// ... solution values
  ///        given point in case the Dirichlet BC depends on solution as well
  /// @param Uint ... component (when we're solving a system of equations)

  Real (*m_expression_ptr)(const math::DenseConstVecView<Real> &,
                           const interpolation::VectorMeshFunction<Real>::const_entry_type &,
                           const Uint);

  /// Corrective residual for each node
  /// Dimensions: [NEQ x nb_nodes]
  // mesh::StdRegionDataMap<math::DynamicVector<Real>> m_bface_residual;

  /// Update coefficient for each node in boundary element
  common::DataMap<mesh::PointSetTagExt, math::DenseDVec<Real>> m_bface_elem_update_coeff;
  std::unique_ptr<bc_projection_data_t> m_projection_data;

  typename Physics::Properties m_phys_properties;
  typename Physics::Properties::SolGradM m_solution_gradient;
  typename Physics::Properties::FluxV m_ghost_state;
  typename Physics::Properties::FluxV m_state_at_infty;
  typename Physics::Properties::FluxV m_qd_pt_flux_vector;
  typename Physics::Properties::FluxV m_ghost_flux_vector;

  /// Vector of eigenvalues of projected Jacobian
  typename Physics::Properties::FluxV m_phys_eigvals;

  /// Re-composed Jacobian matrix
  typename Physics::JM GhostJacobian;

  private:
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakDirichletBC<MeshConfig, Physics, BcDim>::WeakDirichletBC() : base(), m_projection_data(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakDirichletBC<MeshConfig, Physics, BcDim>::WeakDirichletBC(const std::string &name) : base(name)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakDirichletBC<MeshConfig, Physics, BcDim>::~WeakDirichletBC()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakDirichletBC<MeshConfig, Physics, BcDim>::configure_mesh_data(
    const std::shared_ptr<const mesh_type> &mesh, const common::PtrHandle<const cell_dofs> geo_dofs,
    const common::PtrHandle<const cell_dofs> sol_dofs, const std::string &domain_name)
{
  base::configure_mesh_data(mesh, geo_dofs, sol_dofs, domain_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakDirichletBC<MeshConfig, Physics, BcDim>::configure_spaces(const SFunc sf_type,
                                                                   const PointSetID quad_type,
                                                                   const Uint quadrature_order)
{
  base::configure_spaces(sf_type, quad_type, quadrature_order);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void WeakDirichletBC<MeshConfig, Physics, BcDIM>::set_expression(
    Real (*expr_ptr)(const math::DenseConstVecView<Real> &,
                     const interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint))
{
  m_expression_ptr = expr_ptr;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakDirichletBC<MeshConfig, Physics, BcDim>::compute_on_element_weak(
    const Uint idx_in_metric, geo_cell_metric_type const &cell_geo_metric,
    sol_cell_metric_type const &cell_sol_metric, flux_cell_metric_type const &cell_flux_metric)
{
#if 0
  if (!m_projection_data)
  {

    const mesh::MeshBoundarySet<MeshConfig> &mesh_boundaries = this->mesh()->all_boundaries();
    const mesh::BoundaryFacets<MeshConfig, BcDim> &boundary =
        *(mesh_boundaries.domain(this->domain_name()));

    const_dof_iterator geo_elem_it = boundary.cbegin(*this->geo_dofs());

    m_projection_data = std::unique_ptr<bc_projection_data_t>(new bc_projection_data_t());
    m_projection_data->setup_projection_lhs(geo_elem_it, this->sol_ranges());
    m_projection_data->setup_projection_rhs(geo_elem_it, this->sol_ranges(),
                                            this->metric_data().m_geo_metric,
                                            this->metric_data().m_sol_metric, m_expression_ptr);

    m_projection_data->project();
  }

#endif

  const Uint nb_qd_pts = cell_geo_metric.nb_qd_pts();
  const Uint nb_nodes  = cell_sol_metric.nb_dof_in_cell();

  const math::DenseConstVecView<Real> jacobians         = cell_geo_metric.jdet();
  const math::DenseDVec<Real> &weights                  = cell_geo_metric.pt_weights();
  const math::DenseConstMatView<Real> qd_pts_phys_coord = cell_geo_metric.interpolated_coords();
  const math::DenseConstMatView<Real> normals_at_qd_pts = cell_geo_metric.normals();

  const math::DenseDMat<Real> &V = cell_sol_metric.reference_sf_values();

  const math::DenseConstMatView<Real> sol_at_qd_pts = cell_sol_metric.field_values();

  // base::m_bface_residual.fill(0.0);

  // Integrate the corrective residual along the boundary face
  // Compute the update coefficients

  math::DenseDVec<Real> &bface_res =
      *(base::m_bface_residual.std_region_data(cell_sol_metric.std_region_type()));

  math::DenseDVec<Real> &bface_update_coeff_vector =
      *(base::m_bface_elem_update_coeff.std_region_data(cell_sol_metric.std_region_type()));
  bface_update_coeff_vector.fill(0.0);

  Real eig_max = 0.0;

  base::m_total_facet_residual.fill(0.0);

  for (Uint q = 0; q < nb_qd_pts; ++q)
  {
    const Real wq = jacobians[q] * weights[q];

    const math::DenseConstVecView<Real> state_at_qd_pt  = sol_at_qd_pts.row_transpose(q);
    const math::DenseConstVecView<Real> coord_at_qd_pt  = qd_pts_phys_coord.row_transpose(q);
    const math::DenseConstVecView<Real> normal_at_qd_pt = normals_at_qd_pts.row_transpose(q);

    Physics::compute_properties(qd_pts_phys_coord.row_transpose(q), state_at_qd_pt,
                                m_solution_gradient, m_phys_properties);

    Physics::build_K_mat(m_phys_properties, normal_at_qd_pt, m_phys_eigvals, GhostJacobian,
                         [](const Real eigvalue) { return (std::min(0.0, eigvalue)); });

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      m_state_at_infty[eq] = (*m_expression_ptr)(coord_at_qd_pt, state_at_qd_pt, eq);
    }

    eig_max = 1.e-9;

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      eig_max           = std::max(eig_max, std::abs(m_phys_eigvals[eq]));
      m_ghost_state[eq] = m_state_at_infty[eq] - state_at_qd_pt[eq];
    }

    /*
    for (Uint n = 0; n < bface_res.size(); ++n)
    {
      bface_update_coeff_vector[n] += wq * eig_max;
    }
    */

    for (Uint n = 0; n < nb_nodes; ++n)
    {
      bface_update_coeff_vector[n] = std::max(bface_update_coeff_vector[n], eig_max);
    }

    m_ghost_flux_vector = GhostJacobian * m_ghost_state;
    base::m_total_facet_residual += wq * m_ghost_flux_vector;

    for (Uint n = 0; n < nb_nodes; ++n)
    {
      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        bface_res[n * Physics::NEQ + eq] += V(q, n) * wq * m_ghost_flux_vector[eq];
      }
    }
  } // Loop over quadrature points
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
