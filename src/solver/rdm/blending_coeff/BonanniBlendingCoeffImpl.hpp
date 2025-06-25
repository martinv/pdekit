#ifndef PDEKIT_Solver_RDM_Bonanni_Blending_Coeff_hpp
#define PDEKIT_Solver_RDM_Bonanni_Blending_Coeff_hpp

#include "common/DataMap.hpp"
#include "interpolation/SolutionSpaceMetric.hpp"
#include "math/DenseDVec.hpp"
#include "math/MathConstants.hpp"
#include "math/unary_ops/VectorNorm.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"
#include "solver/rdm/blending_coeff/BlendingCoeffImplBase.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

namespace detail
{

template <typename MeshConfig, typename Physics>
class BonanniBlendingCoeffImpl : public BlendingCoeffImplBase<MeshConfig>
{
  public:
  /// Default constructor
  BonanniBlendingCoeffImpl();

  /// Default destructor
  ~BonanniBlendingCoeffImpl() override;

  /// Resize internal data
  void setup(typename result_of::dof_map_t<MeshConfig> const &sol_dofs) override;

  /// Evaluate the blending coefficient
  void calculate_impl(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                      const interpolation::VectorMeshFunction<Real> &u,
                      interpolation::ScalarMeshFunction<Real> &function) override;

  /// Clear all internal data
  void clear() override;

  /// Set a parameter
  void set_param(const std::string &param_name, const Real param_value) override;

  private:
  /// DATA
  bool m_use_ho_theta_scaling;

  /// Map of Finite Element values describing element geometry
  std::vector<mesh::DiscreteElemKey> m_geo_fe_values_map;

  /// Map of Finite Element values for solution space
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> m_sol_fe_values_map;

  /// Geometry cache
  interpolation::GeometryCache<MeshConfig::GDIM> m_geo_cache;

  /// Geometry metric
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> m_geo_metric;

  /// Cache for pressure field
  interpolation::SolutionCache m_press_cache;

  /// Metric for pressure field
  interpolation::SolutionSpaceMetric<MeshConfig> m_press_metric;

  /// Cache for velocity values
  interpolation::SolutionCache m_vel_cache;

  /// Metric for velocity values
  interpolation::SolutionSpaceMetric<MeshConfig> m_vel_metric;

  /// Pressure values in one element
  math::DenseDVec<Real> m_press_in_elem;

  /// Velocity values in one element
  math::DenseDVec<Real> m_vel_in_elem;

  /// Characteristic element size
  std::vector<Real> m_lc;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
BonanniBlendingCoeffImpl<MeshConfig, Physics>::BonanniBlendingCoeffImpl()
    : m_use_ho_theta_scaling(false)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
BonanniBlendingCoeffImpl<MeshConfig, Physics>::~BonanniBlendingCoeffImpl()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void BonanniBlendingCoeffImpl<MeshConfig, Physics>::setup(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs)
{
  mesh::StdPointSet quad;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell                     = sol_dofs.active_cell(mesh::ActiveIdx(c));

    // Type of element we are processing
    const mesh::PointSetTag geo_elem_type_tag = tcell_view.pt_set_id();
    const mesh::PointSetTag sol_elem_type_tag = sol_cell.pt_set_id();

    // Corresponding quadrature - we need to integrate pressure gradient in
    // each element.
    const Uint quad_order =
        std::max(sol_elem_type_tag.poly_order(), geo_elem_type_tag.poly_order());

    // Finite Element values for geometry element
    const mesh::PointSetTagExt geo_elem_type_tag_ext(geo_elem_type_tag, quad_order,
                                                     mesh::CellTransform::NO_TRANS, 0u);
    const mesh::sf::SFTag geo_sf(geo_elem_type_tag.elem_shape(), SFunc::Lagrange,
                                 geo_elem_type_tag.poly_order(), ModalBasis::Modal);
    const mesh::PointSetTag geo_eval_pt_set(geo_elem_type_tag.elem_shape(), quad_order,
                                            PointSetID::Gauss);
    const mesh::PointSetTagExt geo_eval_pt_set_ext(geo_eval_pt_set, quad_order,
                                                   mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey geo_key(geo_elem_type_tag_ext, geo_sf, geo_eval_pt_set_ext);

    mesh::add_unique_discr_elem_key(m_geo_fe_values_map, geo_key);

    const mesh::PointSetTagExt sol_elem_key =
        mesh::PointSetTagExt(sol_elem_type_tag, quad_order, mesh::CellTransform::NO_TRANS, 0u);

    quad.change_type(geo_eval_pt_set);

    common::PtrHandle<interpolation::FEValues> sol_fe_values =
        m_sol_fe_values_map.std_region_data(sol_elem_key);

    if (sol_fe_values.is_null())
    {
      // Solution shape function for this element
      const mesh::sf::SFTag sol_sf_tag =
          mesh::sf::SFTag(sol_elem_type_tag.elem_shape(), SFunc::Lagrange,
                          sol_elem_type_tag.poly_order(), ModalBasis::Modal);

      sol_fe_values = m_sol_fe_values_map.create(sol_elem_key);
      (*sol_fe_values).configure(sol_elem_type_tag, sol_sf_tag);
      (*sol_fe_values).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }
  } // Loop over all active cells

  m_geo_cache.allocate(m_geo_fe_values_map.cbegin(), m_geo_fe_values_map.cend(), 1u);
  m_geo_metric.allocate_buffer(m_geo_fe_values_map.cbegin(), m_geo_fe_values_map.cend(), 1u);

  m_press_cache.allocate(m_sol_fe_values_map.cbegin(), m_sol_fe_values_map.cend(), 1u, 1u);
  m_press_metric.allocate_buffer(m_sol_fe_values_map.cbegin(), m_sol_fe_values_map.cend(), 1u, 1u);

  // This cache can hold only one block, but with Physics::DIM fields (2 or 3
  // velocity components)
  m_vel_cache.allocate(m_sol_fe_values_map.cbegin(), m_sol_fe_values_map.cend(), 1u, Physics::DIM);
  m_vel_metric.allocate_buffer(m_sol_fe_values_map.cbegin(), m_sol_fe_values_map.cend(), 1u,
                               Physics::DIM);

  m_lc.resize(sol_dofs.nb_active_cells());
  m_lc.assign(sol_dofs.nb_active_cells(), 0.0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void BonanniBlendingCoeffImpl<MeshConfig, Physics>::calculate_impl(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &u,
    interpolation::ScalarMeshFunction<Real> &function)
{
  typename Physics::SolGradM u_grad_in_node;
  typename Physics::Properties phys_props;

  Real elem_volume;
  typename Physics::CoordV v_mean_in_elem;
  typename Physics::CoordV p_grad_in_elem;

  typename Physics::CoordV v_mean_in_domain;
  v_mean_in_domain.fill(0.0);
  Real domain_volume = 0.0;

  Real p_min = std::numeric_limits<Real>::max();
  Real p_max = -p_min;

  m_lc.assign(sol_dofs.nb_active_cells(), 0.0);

  // Vector of pressure gradients in one element
  // Each entry corresponds to derivatives of p with respect to one coordinate
  // variable The derivatives are stored in ConstMatrixBlock<Real> holding
  // derivative values in each quadrature point Example: dp[X1][3] holds
  // derivatives of p with respect to X1 = Y in fourt quadrature point
  std::vector<math::DenseConstMatView<Real>> dp(Physics::DIM);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    m_geo_cache.flush();
    m_geo_metric.empty_buffer();

    m_press_cache.flush();
    m_press_metric.empty_buffer();

    m_vel_cache.flush();
    m_vel_metric.empty_buffer();

    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell                     = sol_dofs.active_cell(mesh::ActiveIdx(c));

    // Type of element we are processing
    const mesh::PointSetTag geo_elem_type_tag = tcell_view.pt_set_id();
    const mesh::PointSetTag sol_elem_type_tag = sol_cell.pt_set_id();

    // Corresponding quadrature - we need to integrate pressure gradient in
    // each element.
    const Uint quad_order =
        std::max(sol_elem_type_tag.poly_order(), geo_elem_type_tag.poly_order());

    const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell_view.coordinates();

    const math::DenseConstMatView<Real> sol_cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), sol_cell.pt_set_id(), tcell_view.coordinates());

    const mesh::PointSetTagExt geo_elem_type_tag_ext(geo_elem_type_tag, quad_order,
                                                     mesh::CellTransform::NO_TRANS, 0u);
    const mesh::sf::SFTag geo_sf(geo_elem_type_tag.elem_shape(), SFunc::Lagrange,
                                 geo_elem_type_tag.poly_order(), ModalBasis::Modal);
    const mesh::PointSetTag geo_eval_pt_set(geo_elem_type_tag.elem_shape(), quad_order,
                                            PointSetID::Gauss);
    const mesh::PointSetTagExt geo_eval_pt_set_ext(geo_eval_pt_set, quad_order,
                                                   mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey geo_key(geo_elem_type_tag_ext, geo_sf, geo_eval_pt_set_ext);

    m_geo_cache.push_back_to_buffer(geo_cell_coords, geo_key);

    const mesh::PointSetTagExt sol_elem_key =
        mesh::PointSetTagExt(sol_elem_type_tag, quad_order, mesh::CellTransform::NO_TRANS, 0u);

    m_press_in_elem.resize(sol_cell.nb_vert());
    m_vel_in_elem.resize(sol_cell.nb_vert() * Physics::DIM);

    for (Uint v = 0; v < sol_cell.nb_vert(); ++v)
    {
      const math::DenseConstVecView<Real> node_coords = sol_cell_coords.row_transpose(v);
      const interpolation::VectorMeshFunction<Real>::const_entry_type u_in_node =
          u.const_value(sol_cell.vertex(v));
      Physics::compute_properties(node_coords, u_in_node, u_grad_in_node, phys_props);

      m_press_in_elem[v] = phys_props.P;

      p_min = std::min(p_min, phys_props.P);
      p_max = std::max(p_max, phys_props.P);

      for (Uint d = 0; d < Physics::DIM; ++d)
      {
        m_vel_in_elem[v * Physics::DIM + d] = phys_props.V[d];
      }
    }

    m_press_cache.push_vec_to_buffer(sol_cell, m_press_in_elem, sol_elem_key);
    m_vel_cache.push_vec_to_buffer(sol_cell, m_vel_in_elem, sol_elem_key);

    m_geo_metric.evaluate(m_geo_cache, interpolation::RebuildMetricIndex{true});
    m_press_metric.evaluate(m_geo_metric, m_press_cache, interpolation::ComputeMetricDerivs{true},
                            interpolation::RebuildMetricIndex{true});
    m_vel_metric.evaluate(m_geo_metric, m_vel_cache, interpolation::ComputeMetricDerivs{false},
                          interpolation::RebuildMetricIndex{true});

    const typename interpolation::GeometryMetric<
        MeshConfig::GDIM, MeshConfig::TDIM>::cellwise_metric cell_geo_metric =
        m_geo_metric.cellwise_values(0);
    const typename interpolation::SolutionSpaceMetric<MeshConfig>::cellwise_metric
        cell_pressure_metric = m_press_metric.cellwise_values(0);
    const typename interpolation::SolutionSpaceMetric<MeshConfig>::cellwise_metric
        cell_velocity_metric = m_vel_metric.cellwise_values(0);

    // Get jacobian values in quadrature points of current element
    const math::DenseConstVecView<Real> jq = cell_geo_metric.jdet();
    // Get integration weights for current element
    const math::DenseDVec<Real> &wq = cell_geo_metric.pt_weights();

    const math::DenseConstMatView<Real> vel_in_qd_pt = cell_velocity_metric.field_values();

    /*
    std::cout << "---------------------------------------------" <<
    std::endl; std::cout << "Coordinates = " << std::endl
              << sol_cell_coords << std::endl;
    m_press_metric.print_fill_status();
    std::cout << "Pressure values = " << std::endl;
    for (Uint i = 0; i < m_press_in_elem.size(); ++i)
    {
      std::cout << m_press_in_elem[i][0] << " ";
    }
    std::cout << std::endl;

    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      std::cout << "Pressure gradients: dp/dx" << d << ":" << std::endl;
      std::cout << cell_pressure_metric.field_derivatives(d) << std::endl;
    }
    */

    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      dp[d] = cell_pressure_metric.field_derivatives(d);
    }

    elem_volume = 0.0;
    v_mean_in_elem.fill(0.0);
    p_grad_in_elem.fill(0.0);

    // Compute the mean velocity in one element and gradient of pressure in
    // one element
    for (Uint q = 0; q < cell_geo_metric.nb_qd_pts(); ++q)
    {
      const Real jw_q = jq[q] * wq[q];
      elem_volume += jw_q;

      for (Uint d = 0; d < Physics::DIM; ++d)
      {
        v_mean_in_elem[d] += jw_q * vel_in_qd_pt(q, d);
        p_grad_in_elem[d] += jw_q * dp[d](q, 0);
      }
    }

    // Accumulate to domain volume
    domain_volume += elem_volume;

    // Initialize value of theta for given cell
    function[c] = 0.0;

    // function[c] = math::norm_e2(p_grad_in_elem) / elem_volume;

    // For the moment, store in the resulting vector of thetas the value of
    // int(grad p) DOT v_mean
    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      v_mean_in_domain[d] += v_mean_in_elem[d];
      v_mean_in_elem[d] /= elem_volume;
      p_grad_in_elem[d] /= elem_volume;
      function[c] += v_mean_in_elem[d] * p_grad_in_elem[d];
    }

    const Real vel_norm = math::norm_e2(v_mean_in_elem);

    // Multiply theta in each cell by ||v||
    function[c] *= vel_norm;

    m_lc[c] = std::pow(elem_volume, 1. / Physics::DIM);

  } // Loop over all active cells

  // Divide mean velocity in domain by total domain volume
  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    v_mean_in_domain[d] /= domain_volume;
  }

  const Real delta_pv2 =
      (p_max - p_min) * math::norm_e2(v_mean_in_domain) * math::norm_e2(v_mean_in_domain);

  if (m_use_ho_theta_scaling)
  {

    for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
    {
      // Take the '+' operator of the function
      const Real sc               = std::max(function[c] / delta_pv2, 0.0);
      const mesh::MeshEntity cell = sol_dofs.active_cell(mesh::ActiveIdx(c));
      const Real p_order          = cell.pt_set_id().poly_order();

      function[c] = std::min(1., sc * sc * std::pow(m_lc[c], p_order));
      // function[c] = std::min(1., std::pow(sc * m_lc[c], p_order));
    }
  }
  else
  {
    for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
    {
      // Take the '+' operator of the function
      const Real sc = std::max(function[c] / delta_pv2, 0.0);
      function[c]   = std::min(1., sc * sc * m_lc[c]);
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void BonanniBlendingCoeffImpl<MeshConfig, Physics>::clear()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void BonanniBlendingCoeffImpl<MeshConfig, Physics>::set_param(const std::string &param_name,
                                                              const Real param_value)
{
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
