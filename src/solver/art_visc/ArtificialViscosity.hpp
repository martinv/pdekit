#ifndef PDEKIT_Solver_Artificial_Viscosity_hpp
#define PDEKIT_Solver_Artificial_Viscosity_hpp

#include <array>

#include "common/DataMap.hpp"
#include "interpolation/GeometryCache.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/mesh_function/MeshFunctionTools.hpp"
#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/MathConstants.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/local_topology/TraceEntityTuple.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"

namespace pdekit
{

namespace solver
{

template <typename MeshConfig>
class ArtificialViscosity
{
  public:
  /// Default constructor
  ArtificialViscosity();

  /// Delete copy constructor
  ArtificialViscosity(const ArtificialViscosity &other) = delete;

  /// Default destructor
  ~ArtificialViscosity();

  /// Delete assignment operator
  ArtificialViscosity &operator=(const ArtificialViscosity &rhs) = delete;

  /// Clear all internal data
  void clear();

  /// Prepare the internal data to perform per-element projection
  /// @param reference_elements ... set of element types in the __SOLUTION__
  /// space
  void prepare_indicator_operators(typename result_of::dof_map_t<MeshConfig> const &sol_dofs);

  /// Compute the indicator values in mesh nodes
  void compute_indicator_normalized(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                                    const interpolation::VectorMeshFunction<Real> &u,
                                    interpolation::ScalarMeshFunction<Real> &indicator,
                                    const Real kappa = 5.0);

  /// Compute the indicator values in mesh nodes
  void compute_indicator(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                         const interpolation::VectorMeshFunction<Real> &u,
                         interpolation::ScalarMeshFunction<Real> &indicator,
                         const Real kappa = 5.0);

  /// Compute artificial viscosity in mesh vertices
  void compute_artificial_viscosity_nodal(mesh::Tria<MeshConfig> &input_mesh,
                                          typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                                          const interpolation::VectorMeshFunction<Real> &u,
                                          interpolation::ScalarMeshFunction<Real> &art_visc,
                                          const Real kappa = 5.0, const Real factor = 1.0);

  /// Compute cellwise artificial viscosity
  void compute_artificial_viscosity(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                                    const interpolation::VectorMeshFunction<Real> &u,
                                    interpolation::ScalarMeshFunction<Real> &art_visc,
                                    const Real kappa = 5.0, const Real factor = 1.0);

  private:
  struct IndicatorScratchData
  {
    /// Mass matrix for local L2 projection
    math::DenseDMat<Real> M;

    /// Inverse of the mass matrix
    math::DenseDMat<Real> Minv;

    /// Vandermonde matrix for modal basis where
    /// the columns corresponding to highest-order modes
    /// are zeroed out
    math::DenseDMat<Real> V_modal_p_m_1;

    /// Vandermonde matrix to convert from modal to nodal
    /// expansion coefficients in element dofs
    math::DenseDMat<Real> V_modal_to_nodal_dof;

    /// Vandermonde matrix to convert from nodal to modal
    /// expansion coefficients in element dofs
    math::DenseDMat<Real> V_nodal_to_modal_dof;

    /// Interpolation matrix to transfer modal solution from element
    /// nodes to quadrature points
    math::DenseDMat<Real> V_modal_dof_to_quad;

    /// Interpolation matrix to interpolate solution from element
    /// nodes to quadrature points
    math::DenseDMat<Real> V_nodal_dof_to_quad;

    /// Local RHS during L2 projection
    math::DenseDVec<Real> rhs;

    /// Nodal solution values in one element (input)
    math::DenseDVec<Real> u_nodal_dof;

    /// Solution values (nodal or modal) interpolated to quadrature
    /// points of one element
    math::DenseDVec<Real> u_nodal_qd_pt;

    /// Solution values (modal) interpolated to quadrature
    /// points of one element with the lower-order modal
    /// Vandermonde matrix (highest-order modes switched off)
    math::DenseDVec<Real> u_quad_p_m_1;

    /// Modal expansion coefficients in one element (output)
    math::DenseDVec<Real> u_modal_dof;
  };

  /// Map of Finite Element values describing element geometry
  std::vector<mesh::DiscreteElemKey> m_geo_fe_values_map;

  /// Map of Finite Element values with modal basis of polynomial order p
  /// Modal basis is evaluated at quadrature points in element
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> m_modal_fe_values_qd_pt_map;

  /// Map of Finite Element values with modal basis of polynomial order p
  /// Modal basis is evaluated at degrees of freedom in element
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> m_modal_fe_values_dof_map;

  /// Map of scratch data for computation of local L2 projection
  common::DataMap<mesh::PointSetTagExt, IndicatorScratchData> m_scratch_data_map;

  /// Geometry cache
  interpolation::GeometryCache<MeshConfig::GDIM> m_geo_cache;

  /// Geometry metric
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> m_geo_metric;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
ArtificialViscosity<MeshConfig>::ArtificialViscosity()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
ArtificialViscosity<MeshConfig>::~ArtificialViscosity()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void ArtificialViscosity<MeshConfig>::clear()
{
  m_geo_fe_values_map.clear();
  m_modal_fe_values_qd_pt_map.clear();
  m_modal_fe_values_dof_map.clear();
  m_scratch_data_map.clear();
  m_geo_cache.clear();
  m_geo_metric.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void ArtificialViscosity<MeshConfig>::prepare_indicator_operators(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs)
{
  mesh::StdRegion std_region;
  mesh::sf::ShapeFunction shape_function;
  mesh::StdPointSet quad;
  math::DenseDVec<Real> weights;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell                     = sol_dofs.active_cell(mesh::ActiveIdx(c));

    // Type of element we are processing
    const mesh::PointSetTag geo_elem_type_tag = tcell_view.pt_set_id();
    const mesh::PointSetTag sol_elem_type_tag = sol_cell.pt_set_id();

    // Modal shape function set that will be used on this element type
    const mesh::sf::SFTag geo_sf_tag =
        mesh::sf::SFTag(geo_elem_type_tag.elem_shape(), SFunc::Lagrange,
                        geo_elem_type_tag.poly_order(), ModalBasis::Modal);

    const mesh::sf::SFTag modal_sf_tag =
        mesh::sf::SFTag(sol_elem_type_tag.elem_shape(), SFunc::Modal,
                        sol_elem_type_tag.poly_order(), ModalBasis::Modal);
    shape_function.change_type(sol_elem_type_tag, modal_sf_tag);

    // Corresponding quadrature. We want to compute the mass matrix for
    // local L2 projection, hence the quadrature order should be 2 x
    // polynomial order of the shape functions on the element: we will
    // integrate the terms of the type phi_i x phi_j
    const Uint quad_order =
        std::max(2 * sol_elem_type_tag.poly_order(), geo_elem_type_tag.poly_order());

    // Finite Element values for geometry element
    const mesh::PointSetTagExt geo_pt_set_ext(geo_elem_type_tag, quad_order,
                                              mesh::CellTransform::NO_TRANS, 0u);
    const mesh::PointSetTag geo_eval_pt_set(geo_elem_type_tag.elem_shape(), quad_order,
                                            PointSetID::Gauss);
    const mesh::PointSetTagExt geo_eval_pt_set_ext(geo_eval_pt_set, quad_order,
                                                   mesh::CellTransform::NO_TRANS, 0u);

    const mesh::DiscreteElemKey geo_key(geo_pt_set_ext, geo_sf_tag, geo_eval_pt_set_ext);

    mesh::add_unique_discr_elem_key(m_geo_fe_values_map, geo_key);

    quad.change_type(geo_eval_pt_set);

    // Finite Element values for given modal expansion given by
    // 'modal_sf_tag'
    common::PtrHandle<interpolation::FEValues> modal_fe_values_qd_pt =
        m_modal_fe_values_qd_pt_map.std_region_data(
            mesh::PointSetTagExt(sol_elem_type_tag, quad_order, mesh::CellTransform::NO_TRANS, 0u));

    if (modal_fe_values_qd_pt.is_null())
    {
      modal_fe_values_qd_pt = m_modal_fe_values_qd_pt_map.create(
          mesh::PointSetTagExt(sol_elem_type_tag, quad_order, mesh::CellTransform::NO_TRANS, 0u));

      (*modal_fe_values_qd_pt).configure(sol_elem_type_tag, modal_sf_tag);
      (*modal_fe_values_qd_pt).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
      // std::cout << "V modal = " << std::endl <<
      // (*fe_values).Vandermonde() << std::endl;

      // Values of modal functions in dofs (Vandermonde matrix of primal
      // basis)
      common::PtrHandle<interpolation::FEValues> modal_fe_values_dof =
          m_modal_fe_values_dof_map.create(mesh::PointSetTagExt(sol_elem_type_tag, quad_order,
                                                                mesh::CellTransform::NO_TRANS, 0u));

      std_region.change_type(sol_elem_type_tag);
      (*modal_fe_values_dof).configure(sol_elem_type_tag, modal_sf_tag);
      weights.resize(std_region.get().nb_nodes());
      weights.fill(1.0 / std_region.get().nb_nodes());
      (*modal_fe_values_dof).fill_Vandermonde(std_region.get().coordinates(), weights);

      // std::cout << "Modal Vandermonde matrix in DOFs:" << std::endl;
      // std::cout << (*modal_fe_values_dof).Vandermonde() << std::endl;
    }

    common::PtrHandle<IndicatorScratchData> scratch_data = m_scratch_data_map.std_region_data(
        mesh::PointSetTagExt(sol_elem_type_tag, quad_order, mesh::CellTransform::NO_TRANS, 0u));
    if (scratch_data.is_null())
    {
      scratch_data = m_scratch_data_map.create(
          mesh::PointSetTagExt(sol_elem_type_tag, quad_order, mesh::CellTransform::NO_TRANS, 0u));

      const Uint nb_elem_nodes = (*modal_fe_values_qd_pt).Vandermonde().cols();

      (*scratch_data).M.resize(nb_elem_nodes, nb_elem_nodes);
      (*scratch_data).Minv.resize(nb_elem_nodes, nb_elem_nodes);

      const math::DenseDMat<Real> &V_modal_qd_pt = (*modal_fe_values_qd_pt).Vandermonde();
      math::DenseDMat<Real> &V_modal_p_m_1       = (*scratch_data).V_modal_p_m_1;

      V_modal_p_m_1.resize(V_modal_qd_pt.rows(), V_modal_qd_pt.cols());
      V_modal_p_m_1 = V_modal_qd_pt;

      const math::DenseConstVecView<bool> is_ho_mode =
          shape_function.get().is_leading_expansion_term();

      for (Uint m = 0; m < is_ho_mode.size(); ++m)
      {
        if (is_ho_mode[m])
        {
          for (Uint i = 0; i < V_modal_p_m_1.rows(); ++i)
          {
            V_modal_p_m_1(i, m) = 0.0;
          }
        }
      }

      // std::cout << "V_modal = " << std::endl << V_modal << std::endl;
      // std::cout << "V_modal(p-1) = " << std::endl << V_modal_p_m_1 <<
      // std::endl;

      // Prepare matrix for transformation from modal to nodal space
      (*scratch_data).V_modal_to_nodal_dof.resize(nb_elem_nodes, nb_elem_nodes);
      shape_function.get().compute_ref_values(std_region.get().coordinates(),
                                              (*scratch_data).V_modal_to_nodal_dof);

      // Prepare matrix for transformation from nodal to modal space
      (*scratch_data).V_nodal_to_modal_dof.resize(nb_elem_nodes, nb_elem_nodes);
      (*scratch_data).V_modal_to_nodal_dof.inv((*scratch_data).V_nodal_to_modal_dof);

      // Prepare matrix for transformation from nodes to quadrature points
      // (IN MODAL SPACE)
      (*scratch_data).V_modal_dof_to_quad.resize(quad.get().size(), nb_elem_nodes);
      shape_function.get().compute_ref_values(quad.get().coordinates(),
                                              (*scratch_data).V_modal_dof_to_quad);

      // Prepare matrix for transformation from nodes to quadrature points
      // (IN NODAL SPACE)
      (*scratch_data).V_nodal_dof_to_quad.resize(quad.get().size(), nb_elem_nodes);

      const mesh::sf::SFTag nodal_sf_tag =
          mesh::sf::SFTag(sol_elem_type_tag.elem_shape(), SFunc::Lagrange,
                          sol_elem_type_tag.poly_order(), ModalBasis::Modal);
      shape_function.change_type(sol_elem_type_tag, nodal_sf_tag);
      shape_function.get().compute_ref_values(quad.get().coordinates(),
                                              (*scratch_data).V_nodal_dof_to_quad);

      // std::cout << "V node to quad = " << std::endl <<
      // (*scratch_data).V_node_to_quad << std::endl;

      (*scratch_data).rhs.resize(nb_elem_nodes);
      (*scratch_data).u_nodal_dof.resize(nb_elem_nodes);
      (*scratch_data).u_nodal_qd_pt.resize(quad.get().size());
      (*scratch_data).u_quad_p_m_1.resize(quad.get().size());
      (*scratch_data).u_modal_dof.resize(nb_elem_nodes);
    }
  }

  m_geo_cache.allocate(m_geo_fe_values_map.cbegin(), m_geo_fe_values_map.cend(), 1u);
  m_geo_metric.allocate_buffer(m_geo_fe_values_map.cbegin(), m_geo_fe_values_map.cend(), 1u);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void ArtificialViscosity<MeshConfig>::compute_indicator_normalized(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &u,
    interpolation::ScalarMeshFunction<Real> &indicator, const Real kappa)
{
  compute_artificial_viscosity(sol_dofs, u, indicator, kappa);

  Real indicator_max = 0.0;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    indicator_max = std::max(indicator_max, indicator[c]);
  }

  if (indicator_max < 1.e-14)
  {
    return;
  }

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    indicator[c] /= indicator_max;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void ArtificialViscosity<MeshConfig>::compute_indicator(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &u,
    interpolation::ScalarMeshFunction<Real> &indicator, const Real kappa)
{
  compute_artificial_viscosity(sol_dofs, u, indicator, kappa);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void ArtificialViscosity<MeshConfig>::compute_artificial_viscosity(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &u,
    interpolation::ScalarMeshFunction<Real> &art_visc, const Real kappa, const Real factor)
{
  Real cell_volume;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell                     = sol_dofs.active_cell(mesh::ActiveIdx(c));

    const Uint quad_order =
        std::max(2 * sol_cell.pt_set_id().poly_order(), tcell_view.pt_set_id().poly_order());

    m_geo_cache.flush();

    const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell_view.coordinates();

    const mesh::PointSetTag geo_elem_type_tag = tcell_view.pt_set_id();
    const mesh::PointSetTagExt geo_pt_set_ext(geo_elem_type_tag, quad_order,
                                              mesh::CellTransform::NO_TRANS, 0u);
    const mesh::sf::SFTag geo_sf_tag =
        mesh::sf::SFTag(geo_elem_type_tag.elem_shape(), SFunc::Lagrange,
                        geo_elem_type_tag.poly_order(), ModalBasis::Modal);
    const mesh::PointSetTag geo_eval_pt_set(geo_elem_type_tag.elem_shape(), quad_order,
                                            PointSetID::Gauss);
    const mesh::PointSetTagExt geo_eval_pt_set_ext(geo_eval_pt_set, quad_order,
                                                   mesh::CellTransform::NO_TRANS, 0u);

    const mesh::DiscreteElemKey geo_key(geo_pt_set_ext, geo_sf_tag, geo_eval_pt_set_ext);

    m_geo_cache.push_back_to_buffer(geo_cell_coords, geo_key);

    m_geo_metric.empty_buffer();
    m_geo_metric.evaluate(m_geo_cache, interpolation::RebuildMetricIndex{true});

    IndicatorScratchData &sd = *(m_scratch_data_map.std_region_data(
        mesh::PointSetTagExt(sol_cell.pt_set_id(), quad_order, mesh::CellTransform::NO_TRANS, 0u)));

    typedef typename interpolation::VectorMeshFunction<Real>::const_entry_type nodal_value_type;

    for (Uint n = 0; n < sol_cell.nb_vert(); ++n)
    {
      const nodal_value_type u_node = u.const_value(sol_cell.vertex(n));
      sd.u_nodal_dof[n]             = u_node[0];
    }

    const typename interpolation::GeometryMetric<
        MeshConfig::GDIM, MeshConfig::TDIM>::cellwise_metric cell_geo_metric =
        m_geo_metric.cellwise_values(0);

    const math::DenseConstVecView<Real> jq = cell_geo_metric.jdet();
    const math::DenseDVec<Real> &wq        = cell_geo_metric.pt_weights();

    const Uint nb_cell_nodes = sol_cell.nb_vert();

    const interpolation::FEValues &modal_fe_values = (*m_modal_fe_values_qd_pt_map.std_region_data(
        mesh::PointSetTagExt(sol_cell.pt_set_id(), quad_order, mesh::CellTransform::NO_TRANS, 0u)));
    const math::DenseDMat<Real> &V_modal           = modal_fe_values.Vandermonde();

    cell_volume = 0.0;

    // Interpolate values from nodes to quadrature points
    sd.u_nodal_qd_pt = sd.V_nodal_dof_to_quad * sd.u_nodal_dof;

    sd.M.fill(0.0);
    sd.rhs.fill(0.0);

    for (Uint q = 0; q < cell_geo_metric.nb_qd_pts(); ++q)
    {
      const Real wj = jq[q] * wq[q];
      cell_volume += wj;

      for (Uint i = 0; i < nb_cell_nodes; ++i)
      {
        for (Uint j = 0; j < nb_cell_nodes; ++j)
        {
          sd.M(i, j) += wj * V_modal(q, i) * V_modal(q, j);
        }
        sd.rhs[i] += wj * sd.u_nodal_qd_pt[q] * V_modal(q, i);
      }
    } // Loop over quadrature points

    // Invert the mass matrix for local L2 projection
    sd.M.inv(sd.Minv);

    // Compute the modal expansion coefficients in dofs
    sd.u_modal_dof = sd.Minv * sd.rhs;

    // std::cout << "Nodal values in dofs:" << std::endl << sd.u_quad <<
    // std::endl; std::cout << "Modal expansion coefficients: " << std::endl
    // << sd.u_modal_coeffs << std::endl;

    // These values are values of u in quadrature points computed the modal
    // expansion This should yield the same result as
    //
    // sd.u_quad = sd.V_node_to_quad * sd.u_nodal;
    //
    // which was computed already above. The line below is therefore
    // reduntant and could be commented out
    sd.u_nodal_qd_pt = V_modal * sd.u_modal_dof; // Was alread computed and could be commented out
    // std::cout << "Modal values in quadrature points:" << std::endl <<
    // sd.u_quad << std::endl;

    sd.u_quad_p_m_1 = sd.V_modal_p_m_1 * sd.u_modal_dof;

    // std::cout << "uq(p) = " << std::endl << sd.u_quad << std::endl;
    // std::cout << "uq(p-1) = " << std::endl << sd.u_quad_p_m_1 <<
    // std::endl << std::endl;

    // Compute the actual indicator
    Real nom   = 0.0;
    Real denom = 0.0;
    // std::cout << "j = ";
    // std::cout.precision(10);
    for (Uint q = 0; q < cell_geo_metric.nb_qd_pts(); ++q)
    {
      const Real wj = jq[q] * wq[q];

      nom += wj * (sd.u_nodal_qd_pt[q] - sd.u_quad_p_m_1[q]) *
             (sd.u_nodal_qd_pt[q] - sd.u_quad_p_m_1[q]);
      denom += wj * sd.u_nodal_qd_pt[q] * sd.u_nodal_qd_pt[q];
      // std::cout << " " << jq[q];
    }

    // std::cout << std::endl;

    // std::cout << "[" << nom << "," << denom << "]" << std::endl;

    // Element-wise value
    const Real Se = (denom <= 1.e-10) ? 0.0 : nom / denom;
    const Real se = std::log(Se) / std::log(10.);

    const Uint p_sol = sol_cell.pt_set_id().poly_order();
    const Real h     = std::pow(cell_volume, 1. / MeshConfig::GDIM);

    const Real eps0 = h / p_sol;
    const Real s0   = 1. / (p_sol * p_sol * p_sol * p_sol);
    // const Real kappa = 0.3 * s0;
    // const Real kappa = 5.0;

    /*
    std::cout << "=========================" << std::endl;
    std::cout << "cell " << c << std::endl;
    std::cout << "Se = " << Se << std::endl;
    std::cout << "se = " << se << std::endl;
    std::cout << "p = " << p_sol << std::endl;
    std::cout << "h = " << h << std::endl;
    std::cout << "eps0 = " << eps0 << std::endl;
    std::cout << "s0 = " << s0 << std::endl;
    std::cout << "kappa = " << kappa << std::endl;
    */

    if (se < (s0 - kappa))
    {
      art_visc[c] = 0.0;
    }
    else if (se > (s0 + kappa))
    {
      art_visc[c] = factor * eps0;
    }
    else
    {
      art_visc[c] = factor * 0.5 * eps0 * (1. + std::sin(math::pi * (se - s0) / (2. * kappa)));
    }

  } // Loop over cells in the dof handler
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void ArtificialViscosity<MeshConfig>::compute_artificial_viscosity_nodal(
    mesh::Tria<MeshConfig> &input_mesh, typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &u,
    interpolation::ScalarMeshFunction<Real> &art_visc, const Real kappa, const Real factor)
{
  // -----------------------------------------------
  // Compute the artifical viscosity values in cells
  // -----------------------------------------------

  interpolation::ScalarMeshFunction<Real> art_visc_cellwise("", "art_visc_cellwise");
  art_visc_cellwise.resize(sol_dofs.nb_active_cells());

  compute_artificial_viscosity(sol_dofs, u, art_visc_cellwise, kappa, factor);

  // -------------------------------------------
  // Move the values from cells to cell vertices
  // -------------------------------------------
  interpolation::MeshFunctionTools::interp_from_cells_to_nodes_p1_max(input_mesh, sol_dofs,
                                                                      art_visc_cellwise, art_visc);
}

// ----------------------------------------------------------------------------

} // namespace solver

} // namespace pdekit

#endif
