#ifndef PDEKIT_Interpolation_Operator_Residual_Cellwise_hpp
#define PDEKIT_Interpolation_Operator_Residual_Cellwise_hpp

#include "interpolation/FluxSpaceMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "mesh/MeshConfig.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
class OperatorResidualCellwise
{
  public:
  /// Default constructor
  OperatorResidualCellwise();

  /// Default destructor
  ~OperatorResidualCellwise();

  /// Prepare internal variables by inspecting the dof storage for geometry
  /// and solution
  void setup(typename result_of::dof_map_t<MeshConfig> const &sol_dofs);

  /// Compute the operator residual in every dof of the mesh
  /// @note For the moment, the 'operator residual' is just divergence of
  /// fluxes
  void evaluate(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                interpolation::VectorMeshFunction<Real> const &solution,
                interpolation::VectorMeshFunction<Real> &operator_res);

  private:
  /// TYPEDEFS
  typedef interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric_type;
  typedef interpolation::FluxSpaceMetric<MeshConfig, Physics> flux_metric_type;

  /// DATA
  /// Cache for geometry elements
  interpolation::GeometryCache<MeshConfig::GDIM> m_geo_cache;

  /// Metric of geometry elements
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> m_geo_metric;

  /// Cache for solution elements
  interpolation::SolutionCache m_sol_cache;

  /// Cache for solution metric
  interpolation::SolutionSpaceMetric<MeshConfig> m_sol_metric;

  /// Metric for interpolation of fluxes
  flux_metric_type m_flux_metric;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
OperatorResidualCellwise<MeshConfig, Physics>::OperatorResidualCellwise()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
OperatorResidualCellwise<MeshConfig, Physics>::~OperatorResidualCellwise()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void OperatorResidualCellwise<MeshConfig, Physics>::setup(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs)
{
  std::vector<mesh::DiscreteElemKey> geo_cell_map;
  common::DataMap<mesh::PointSetTagExt, FEValues> sol_cell_map;

  mesh::StdPointSet quad;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell                     = sol_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::PointSetTag geo_tag = tcell_view.pt_set_id();
    const mesh::PointSetTag sol_tag = sol_cell.pt_set_id();

    // The GEOMETRY element key also knows what is the polynomial order of
    // the corresponding SOLUTION element - this is useful in adapted mesh,
    // when multiple solution elements (with different polynomial order)
    // might have the same geometrical support element
    //
    // REMARK: this also takes care of situation when the geometry and
    // solution elements have different nodal distributions (equidistant,
    // warpblend), because the geometry FEValues are generated using shape
    // functions associated with the DOFs of the geometry element, but
    // EVALUATED AT THE DOFS OF THE SOLUTION ELEMENT, which can have
    // different nodal distribution!

    const Uint p_order = std::max(1u, std::max(2 * geo_tag.poly_order(), 2 * sol_tag.poly_order()));

    const mesh::PointSetTagExt geo_tag_ext =
        mesh::PointSetTagExt(geo_tag, p_order, mesh::CellTransform::NO_TRANS, 0);

    const mesh::PointSetTagExt sol_tag_ext =
        mesh::PointSetTagExt(sol_tag, p_order, mesh::CellTransform::NO_TRANS, 0);

    const mesh::sf::SFTag geo_basis_tag = mesh::sf::SFTag(geo_tag.elem_shape(), SFunc::Lagrange,
                                                          geo_tag.poly_order(), ModalBasis::Modal);

    // The Vandermonde matrix is computed in the DOFs of the SOLUTION
    // element
    const mesh::PointSetTag eval_pt_set_tag_geo(geo_tag.elem_shape(), p_order, PointSetID::Gauss);
    const mesh::PointSetTagExt eval_pt_set_tag_geo_ext(eval_pt_set_tag_geo, p_order,
                                                       mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey geo_key(geo_tag_ext, geo_basis_tag, eval_pt_set_tag_geo_ext);

    bool geo_key_found = false;
    for (const auto &cell : geo_cell_map)
    {
      if (cell == geo_key)
      {
        geo_key_found = true;
        break;
      }
    }
    if (!geo_key_found)
    {
      geo_cell_map.push_back(geo_key);
    }

    common::PtrHandle<interpolation::FEValues> fe_values_sol_ptr =
        sol_cell_map.std_region_data(sol_tag_ext);

    if (fe_values_sol_ptr.is_null())
    {
      fe_values_sol_ptr = sol_cell_map.create(sol_tag_ext);
      (*fe_values_sol_ptr)
          .configure(sol_cell.pt_set_id(),
                     mesh::sf::SFTag(sol_tag.elem_shape(), SFunc::Lagrange, sol_tag.poly_order(),
                                     ModalBasis::Modal));

      const mesh::PointSetTag quad_tag(sol_tag.elem_shape(), p_order, PointSetID::Gauss);
      quad.change_type(quad_tag);

      // Fill Vandermonde matrix
      (*fe_values_sol_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());

      // (*fe_values_sol_ptr).print();
    }
  } // Loop over cells

  // Clear old data from all cache and metric objects
  m_geo_cache.clear();
  m_geo_metric.clear();
  m_sol_cache.clear();
  m_sol_metric.clear();
  m_flux_metric.clear();

  // Setup cache and metric objects for new data
  const Uint nb_blocks = 1;
  m_geo_cache.allocate(geo_cell_map.cbegin(), geo_cell_map.cend(), nb_blocks);
  m_geo_metric.allocate_buffer(geo_cell_map.cbegin(), geo_cell_map.cend(), nb_blocks);

  m_sol_cache.allocate(sol_cell_map.cbegin(), sol_cell_map.cend(), nb_blocks, Physics::NEQ);
  m_sol_metric.allocate_buffer(sol_cell_map.cbegin(), sol_cell_map.cend(), nb_blocks, Physics::NEQ);

  const mesh::CellTopologyView<MeshConfig> first_tcell_view = sol_dofs.tcell(mesh::ActiveIdx(0));

  m_flux_metric.allocate_buffer(SFunc::Lagrange, first_tcell_view.pt_set_id().poly_order(),
                                sol_cell_map.cbegin(), sol_cell_map.cend(), nb_blocks);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void OperatorResidualCellwise<MeshConfig, Physics>::evaluate(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
    interpolation::VectorMeshFunction<Real> const &solution,
    interpolation::VectorMeshFunction<Real> &operator_res)
{
  operator_res.resize(Physics::NEQ, sol_dofs.nb_active_cells());

  const Uint cell_idx_in_metric = 0;

  math::DenseDVec<Real> operator_res_in_cell(Physics::NEQ);
  math::DenseDVec<Real> operator_res_at_qd_pt(Physics::NEQ);

  std::array<math::DenseConstMatView<Real>, Physics::DIM> grad_F;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell                     = sol_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::PointSetTag geo_tag = tcell_view.pt_set_id();
    const mesh::PointSetTag sol_tag = sol_cell.pt_set_id();

    const Uint p_order = std::max(1u, std::max(2 * geo_tag.poly_order(), 2 * sol_tag.poly_order()));

    const mesh::PointSetTagExt geo_tag_ext =
        mesh::PointSetTagExt(geo_tag, p_order, mesh::CellTransform::NO_TRANS, 0);

    const mesh::sf::SFTag geo_basis_tag = mesh::sf::SFTag(geo_tag.elem_shape(), SFunc::Lagrange,
                                                          geo_tag.poly_order(), ModalBasis::Modal);

    // The Vandermonde matrix is computed in the DOFs of the SOLUTION
    // element
    const mesh::PointSetTag eval_pt_set_tag_geo(geo_tag.elem_shape(), p_order, PointSetID::Gauss);
    const mesh::PointSetTagExt eval_pt_set_tag_geo_ext(eval_pt_set_tag_geo, p_order,
                                                       mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey geo_key(geo_tag_ext, geo_basis_tag, eval_pt_set_tag_geo_ext);

    const mesh::PointSetTagExt sol_key =
        mesh::PointSetTagExt(sol_tag, p_order, mesh::CellTransform::NO_TRANS, 0);

    const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell_view.coordinates();

    m_geo_cache.flush();
    m_geo_cache.push_back_to_buffer(geo_cell_coords, geo_key);

    m_sol_cache.flush();
    m_sol_cache.push_back_to_buffer(sol_cell, solution, sol_key);

    m_geo_metric.empty_buffer();
    m_sol_metric.empty_buffer();
    m_flux_metric.empty_buffer();

    m_geo_metric.evaluate(m_geo_cache, RebuildMetricIndex{true});
    m_sol_metric.evaluate(m_geo_metric, m_sol_cache, ComputeMetricDerivs{true},
                          RebuildMetricIndex{true});
    m_flux_metric.evaluate(m_geo_cache, m_geo_metric, m_sol_cache, m_sol_metric,
                           RebuildMetricIndex{true});

    const typename geo_metric_type::cellwise_metric geo_cell_met =
        m_geo_metric.cellwise_values(cell_idx_in_metric);

    const typename flux_metric_type::cellwise_metric flux_cell_met =
        m_flux_metric.cellwise_values(cell_idx_in_metric);

    math::DenseDVec<Real> const &weight_at_qd_pt     = geo_cell_met.pt_weights();
    const math::DenseConstVecView<Real> det_at_qd_pt = geo_cell_met.jdet();

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      grad_F[d] = flux_cell_met.flux_derivatives(d);
    }

    operator_res_in_cell.fill(0.0);

    // Real cell_volume = 0.0;

    for (Uint q = 0; q < geo_cell_met.nb_qd_pts(); ++q)
    {
      operator_res_at_qd_pt.fill(0.0);
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        operator_res_at_qd_pt += grad_F[dim].row_transpose(q);
      }

      const Real wq = weight_at_qd_pt[q] * det_at_qd_pt[q];
      // cell_volume += wq;

      operator_res_in_cell += wq * operator_res_at_qd_pt;
    }

    interpolation::VectorMeshFunction<Real>::entry_type op_res_value = operator_res.value(c);

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      // op_res_value[eq] = 1. / cell_volume * operator_res_in_cell[eq];
      op_res_value[eq] = operator_res_in_cell[eq];
    }
  } // Loop over cells in dof handler
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
