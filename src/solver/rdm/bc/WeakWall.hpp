#ifndef PDEKIT_RDM_Weak_Wall_hpp
#define PDEKIT_RDM_Weak_Wall_hpp

#include "common/DataMap.hpp"
#include "solver/NumFlux.hpp"
#include "solver/rdm/bc/WeakBC.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, Uint BcDim = Physics::DIM - 1>
class WeakWall : public WeakBC<MeshConfig, Physics, BcDim, WeakWall<MeshConfig, Physics, BcDim>>
{
  private:
  using base = WeakBC<MeshConfig, Physics, BcDim, WeakWall<MeshConfig, Physics, BcDim>>;
  using geo_cell_metric_type  = typename base::geo_cell_metric_type;
  using sol_cell_metric_type  = typename base::sol_cell_metric_type;
  using flux_cell_metric_type = typename base::flux_cell_metric_type;

  public:
  /// Default constructor
  WeakWall();

  /// Constructor
  WeakWall(const std::string &name);

  /// Destructor
  ~WeakWall() override;

  /// Compute the corrective residuals on one element
  void compute_on_element_weak(const Uint idx_in_metric,
                               geo_cell_metric_type const &cell_geo_metric,
                               sol_cell_metric_type const &cell_sol_metric,
                               flux_cell_metric_type const &cell_flux_metric);

  /// Set a parameter for the boundary conditions
  void set_parameter(const std::string &param_name, const Real value) override;

  private:
  typename Physics::Properties m_phys_properties;
  typename Physics::Properties::SolGradM m_solution_gradient;
  typename Physics::Properties::FluxV m_flux_vector;
  typename Physics::Properties::FluxV m_ghost_flux_vector;
  // typename Physics::Properties::FluxV m_bface_residual;

  /// Values of flux Jacobian, it's inverse and the diagonal
  /// matrix of Jacobian eigenvalues
  typename Physics::JM Rv;
  typename Physics::JM Lv;

  /// Matrix of negative eigenvalues
  typename Physics::JM Dv;

  // Numerical flux for one of the implementations
  NumFluxLaxFriedrichs<Physics> m_num_flux;

  // Value of relaxation coefficient
  Real m_relax_coeff;
  bool m_use_relax_coeff;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakWall<MeshConfig, Physics, BcDim>::WeakWall()
    : base(), m_relax_coeff(1.0), m_use_relax_coeff(false)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakWall<MeshConfig, Physics, BcDim>::WeakWall(const std::string &name) : base(name)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakWall<MeshConfig, Physics, BcDim>::~WeakWall()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakWall<MeshConfig, Physics, BcDim>::compute_on_element_weak(
    const Uint idx_in_metric, geo_cell_metric_type const &cell_geo_metric,
    sol_cell_metric_type const &cell_sol_metric, flux_cell_metric_type const &cell_flux_metric)
{
  // base::m_bface_residual.fill(0.0); // This is done in the parent class

  const Uint nb_qd_pts = cell_geo_metric.nb_qd_pts();
  const Uint nb_nodes  = cell_sol_metric.nb_dof_in_cell();

  const math::DenseConstVecView<Real> jacobians         = cell_geo_metric.jdet();
  const math::DenseDVec<Real> &weights                  = cell_geo_metric.pt_weights();
  const math::DenseConstMatView<Real> qd_pts_phys_coord = cell_geo_metric.interpolated_coords();
  const math::DenseConstMatView<Real> normals_at_qd_pts = cell_geo_metric.normals();

  const math::DenseDMat<Real> &V = cell_sol_metric.reference_sf_values();

  const math::DenseConstMatView<Real> sol_at_qd_pts = cell_sol_metric.field_values();

  if (!m_use_relax_coeff)
  {
    const Uint poly_order = cell_sol_metric.std_region_type().std_region_tag().poly_order();
    // m_relax_coeff = 1. / std::max(1u, poly_order);
    m_relax_coeff = (poly_order <= 1) ? 1.0 : 1. / (3. * poly_order);
  }

  /*
  std::cout << "==========================================================="
  << std::endl; std::cout << "jacobians = " << std::endl << jacobians <<
  std::endl; std::cout << "weights =  " << std::endl << weights   <<
  std::endl; std::cout << "coords at quadrature points = " << std::endl <<
  qd_pts_phys_coord << std::endl;
  std::cout << "normals at quadrature points = " << std::endl <<
  normals_at_qd_pts << std::endl;
  std::cout << "Vandermonde matrix = " << std::endl << V << std::endl;
  std::cout << "Solution at quadrature points = " << std::endl <<
  sol_at_qd_pts
  << std::endl << std::endl;
  */

  // math::DynamicVector<Real> &bface_res =
  //    *(base::m_bface_residual.std_region_data(cell_sol_metric.std_region_type()));

  math::DenseDVec<Real> &bface_res = *(base::m_bface_residual.std_region_data(mesh::PointSetTagExt(
      cell_sol_metric.std_region_type().std_region_tag(), P0, mesh::CellTransform::NO_TRANS, 0u)));

  // math::DynamicVector<Real> &bface_update_coeff_vector =
  //     *(base::m_bface_elem_update_coeff.std_region_data(cell_sol_metric.std_region_type()));

  math::DenseDVec<Real> &bface_update_coeff_vector =
      *(base::m_bface_elem_update_coeff.std_region_data(
          mesh::PointSetTagExt(cell_sol_metric.std_region_type().std_region_tag(), P0,
                               mesh::CellTransform::NO_TRANS, 0u)));

  // bface_update_coeff_vector.fill(0.0); // This is done in the parent class

  // Integrate the corrective residual along the boundary face
  Real eig_max = 0.0;

  typename Physics::SolV u_ghost;
  base::m_total_facet_residual.fill(0.0);

  for (Uint q = 0; q < nb_qd_pts; ++q)
  {
    const Real wq = jacobians[q] * weights[q];

    Physics::compute_properties(qd_pts_phys_coord.row_transpose(q), sol_at_qd_pts.row_transpose(q),
                                m_solution_gradient, m_phys_properties);

    m_flux_vector.fill(0.0);

    // std::cout << "Coordinate = " << G.m_xq.row_transpose(q) << std::endl;
    // std::cout << "Normal = " << base::m_normals[q] << std::endl;

    const math::DenseConstVecView<Real> normal = normals_at_qd_pts.row_transpose(q);

    Physics::flux_jacobian_eigen_structure(m_phys_properties, normal, Rv, Lv, Dv);

    eig_max = 1.e-9;
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      eig_max = std::max(std::abs(Dv(eq, eq)), eig_max);
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

    // -------------------------------------------------------
    // First implementation of BC: F_wall = (0, p*nx, p*ny, 0)
    // -------------------------------------------------------

#if 0

    /*
    m_ghost_flux_vector.fill(0.0);

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      m_ghost_flux_vector[d + 1] = normal[d] * m_phys_properties.P;
    }

    Physics::flux(m_phys_properties, normal, m_flux_vector);

    // const Real inv_nb_nodes = 1./base::m_bface_residual.rows();

    for (Uint n = 0; n < nb_nodes; ++n)
    {
      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        bface_res[n * Physics::NEQ + eq] +=
            V(q, n) * wq * (m_ghost_flux_vector[eq] - m_flux_vector[eq]);
        // base::m_bface_residual(n,eq) += inv_nb_nodes * (
        // m_ghost_flux_vector[eq] - m_flux_vector[eq] );
      }
    }
    */

    m_ghost_flux_vector.fill(0.0);
    m_flux_vector.fill(0.0);

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      const math::DenseConstMatView<Real> flux_data_qd_pt = cell_flux_metric.flux_values(d);
      m_flux_vector += normal[d] * flux_data_qd_pt.row_transpose(q);

      m_ghost_flux_vector[d + 1] = normal[d] * m_phys_properties.P;
    }

    // Physics::flux(m_phys_properties, normal, m_flux_vector);

    // const Real inv_nb_nodes = 1./base::m_bface_residual.rows();

    base::m_total_facet_residual += wq * (m_ghost_flux_vector - m_flux_vector);

    for (Uint n = 0; n < nb_nodes; ++n)
    {
      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        bface_res[n * Physics::NEQ + eq] +=
            m_relax_coeff * V(q, n) * wq * (m_ghost_flux_vector[eq] - m_flux_vector[eq]);
        // base::m_bface_residual(n,eq) += inv_nb_nodes * (
        // m_ghost_flux_vector[eq] - m_flux_vector[eq] );
      }
    }

#endif

    // -------------------------------------------------------
    // Second implementation of BC - Stefano's version
    // -------------------------------------------------------

#if 1
    Physics::flux(m_phys_properties, normal, m_flux_vector);

    // u_ghost.fill(0.0);

    // Modify the physical properties so that
    // velocity becomes v <- v - (v.n)n

    u_ghost[0] = m_phys_properties.rho;

    Real v_dot_n = 0.0;
    Real v2_new  = 0.0;

    const Real v2_old = m_phys_properties.v2;

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      v_dot_n += m_phys_properties.V[d] * normal[d];
    }

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      // Normal component of velocity
      u_ghost[d + 1] = m_phys_properties.V[d] - v_dot_n * normal[d];
      v2_new += u_ghost[d + 1] * u_ghost[d + 1];
      // The conservative variable has to be multiplied by density
      u_ghost[d + 1] *= m_phys_properties.rho;
    }

    u_ghost[Physics::DIM + 1] =
        m_phys_properties.rhoE + 0.5 * m_phys_properties.rho * (v2_new - v2_old);

    Physics::compute_properties(qd_pts_phys_coord.row_transpose(q), u_ghost, m_solution_gradient,
                                m_phys_properties);
    Physics::flux(m_phys_properties, normal, m_ghost_flux_vector);

    base::m_total_facet_residual += wq * (m_ghost_flux_vector - m_flux_vector);

    for (Uint n = 0; n < nb_nodes; ++n)
    {
      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        bface_res[n * Physics::NEQ + eq] +=
            m_relax_coeff * V(q, n) * wq * (m_ghost_flux_vector[eq] - m_flux_vector[eq]);
      }
    }
#endif

    // -------------------------------------------------------
    // Third implementation of BC - numerical Flux
    // -------------------------------------------------------

#if 0
    Physics::flux(m_phys_properties, normal, m_flux_vector);

    // u_ghost.fill(0.0);

    // Modify the physical properties so that
    // velocity becomes v <- v - (v.n)n

    u_ghost[0] = m_phys_properties.rho;

    Real v_dot_n = 0.0;
    Real v2_new = 0.0;

    // const Real v2_old = m_phys_properties.v2;

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      v_dot_n += m_phys_properties.V[d] * normal[d];
    }

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      // Normal component of velocity
      u_ghost[d + 1] = m_phys_properties.V[d] - v_dot_n * normal[d];
      v2_new += u_ghost[d + 1] * u_ghost[d + 1];
      // The conservative variable has to be multiplied by density
      u_ghost[d + 1] *= m_phys_properties.rho;
    }

    u_ghost[Physics::DIM + 1] = m_phys_properties.rhoE;

    m_num_flux.compute(qd_pts_phys_coord.row_transpose(q), qd_pts_phys_coord.row_transpose(q),
                       normal, sol_at_qd_pts.row_transpose(q), u_ghost, m_ghost_flux_vector);

    base::m_total_facet_residual += wq * (m_ghost_flux_vector - m_flux_vector);

    for (Uint n = 0; n < nb_nodes; ++n)
    {
      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        bface_res[n * Physics::NEQ + eq] +=
            V(q, n) * wq * (m_ghost_flux_vector[eq] - m_flux_vector[eq]);
      }
    }
#endif

  } // Loop over integration points
}

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakWall<MeshConfig, Physics, BcDim>::set_parameter(const std::string &param_name,
                                                         const Real value)
{
  if (param_name == "relax_coeff")
  {
    m_relax_coeff     = value;
    m_use_relax_coeff = true;
  }
  else if (param_name == "automatic_relax_coeff")
  {
    m_relax_coeff     = 1.0;
    m_use_relax_coeff = false;
  }
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
