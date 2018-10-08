#ifndef PDEKIT_Solver_Postprocessing_hpp
#define PDEKIT_Solver_Postprocessing_hpp

#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/SolutionSpaceMetric.hpp"
#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/adaptation/LocalInterpolator.hpp"

namespace pdekit
{

namespace solver
{

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class PostprocessingUtils
{
  public:
  PostprocessingUtils() = default;

  PostprocessingUtils(const PostprocessingUtils &other) = default;

  ~PostprocessingUtils() = default;

  PostprocessingUtils &operator=(const PostprocessingUtils &rhs) = default;

  static Real compute_function_norm(
      const typename pdekit::result_of::dof_map_t<MeshConfig> &geo_cell_dofs,
      const typename pdekit::result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
      const interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM> &geo_space,
      const interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM> &sol_space,
      const interpolation::ScalarMeshFunction<Real> &f);

  static Real compute_function_norm(
      const typename pdekit::result_of::dof_map_t<MeshConfig> &geo_cell_dofs,
      const typename pdekit::result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
      const interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM> &geo_space,
      const interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM> &sol_space,
      const interpolation::VectorMeshFunction<Real> &f);

  template <typename PhysModel>
  static void compute_Euler_Ma(
      const typename pdekit::result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
      const interpolation::VectorMeshFunction<Real> &solution,
      interpolation::ScalarMeshFunction<Real> &f);

  template <typename PhysModel>
  static void compute_Euler_entropy(
      const typename pdekit::result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
      const interpolation::VectorMeshFunction<Real> &solution, const Real rho_in, const Real p_in,
      interpolation::ScalarMeshFunction<Real> &f);

  template <typename PhysModel, typename Evaluator>
  static void compute_Euler_quantity(
      const typename pdekit::result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
      const interpolation::VectorMeshFunction<Real> &solution, const Evaluator &eval,
      interpolation::ScalarMeshFunction<Real> &f);
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Real PostprocessingUtils<MeshConfig>::compute_function_norm(
    const typename pdekit::result_of::dof_map_t<MeshConfig> &geo_cell_dofs,
    const typename pdekit::result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
    const interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM> &geo_space,
    const interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM> &sol_space,
    const interpolation::ScalarMeshFunction<Real> &f)
{
  interpolation::GeometryCache<MeshConfig::GDIM> geo_cache;
  geo_cache.allocate(geo_space.discrete_elements().cbegin(), geo_space.discrete_elements().cend(),
                     1u);

  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric;
  geo_metric.allocate_buffer(geo_space.discrete_elements().cbegin(),
                             geo_space.discrete_elements().cend(), 1u);

  interpolation::SolutionCache sol_cache;
  // This should be rather
  // sol_cache.allocate(sol_space.reference_elements(), 1u, f.nb_fields());
  sol_cache.allocate(sol_space.reference_elements().cbegin(), sol_space.reference_elements().cend(),
                     1u, 1u);

  interpolation::SolutionSpaceMetric<MeshConfig> sol_metric;
  // This should be rather
  // sol_metric.allocate_buffer(sol_space.reference_elements(), 1u,
  // f.nb_fields());
  sol_metric.allocate_buffer(sol_space.reference_elements().cbegin(),
                             sol_space.reference_elements().cend(), 1u, 1u);

  Real sum_cell_contributions = 0.0;

  const auto geo_sf_generator   = geo_space.sf_generator();
  const auto geo_quad_generator = geo_space.quad_generator();

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < geo_cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity geo_cell = geo_cell_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell = sol_cell_dofs.active_cell(mesh::ActiveIdx(c));

    geo_cache.flush();

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), geo_cell.pt_set_id(), tcell_view.coordinates());

    const mesh::PointSetTag geo_pt_set = geo_cell.pt_set_id();
    const ElemShape elem_shape         = geo_pt_set.elem_shape();
    const Uint geo_order               = geo_pt_set.poly_order();

    const mesh::PointSetTagExt geo_pt_set_ext(geo_pt_set, P0, mesh::CellTransform::NO_TRANS, 0u);
    const mesh::sf::SFTag geo_sf = geo_sf_generator(elem_shape, geo_order);

    const mesh::PointSetTag geo_quad = geo_quad_generator(elem_shape, geo_order);
    const mesh::PointSetTagExt geo_quad_ext(geo_quad, P0, mesh::CellTransform::NO_TRANS, 0u);

    const mesh::DiscreteElemKey geo_key(geo_pt_set_ext, geo_sf, geo_quad_ext);

    geo_cache.push_back_to_buffer(cell_coords, geo_key);

    sol_cache.flush();
    sol_cache.push_back_to_buffer(
        sol_cell, f,
        mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

    geo_metric.empty_buffer();
    geo_metric.evaluate(geo_cache, interpolation::RebuildMetricIndex{true});

    sol_metric.empty_buffer();
    sol_metric.evaluate(geo_metric, sol_cache, interpolation::ComputeMetricDerivs{true},
                        interpolation::RebuildMetricIndex{true});

    typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>::cellwise_metric
        cell_geo_met = geo_metric.cellwise_values(0);
    typename interpolation::SolutionSpaceMetric<MeshConfig>::cellwise_metric cell_sol_met =
        sol_metric.cellwise_values(0);

    // Coordinates of quadrature points in physical space
    // const math::ConstMatrixBlock<Real> qd_pts =
    // cell_geo_met.interpolated_coords();

    // Vector of integration weights in quadrature points
    const math::DenseDVec<Real> &wq = cell_geo_met.pt_weights();

    // Vector of transformation jacobians in quadrature points
    const math::DenseConstVecView<Real> jq = cell_geo_met.jdet();

    // Values of solution in quadrature points
    const math::DenseConstMatView<Real> sol_q = cell_sol_met.field_values();

    for (Uint q = 0; q < cell_geo_met.nb_qd_pts(); ++q)
    {
      const Real value2 = sol_q(q, 0) * sol_q(q, 0);
      sum_cell_contributions += wq[q] * jq[q] * value2;
    }
  } // Loop over cells
  return std::sqrt(sum_cell_contributions);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Real PostprocessingUtils<MeshConfig>::compute_function_norm(
    const typename pdekit::result_of::dof_map_t<MeshConfig> &geo_cell_dofs,
    const typename pdekit::result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
    const interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM> &geo_space,
    const interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM> &sol_space,
    const interpolation::VectorMeshFunction<Real> &f)
{
  interpolation::GeometryCache<MeshConfig::GDIM> geo_cache;
  geo_cache.allocate(geo_space.discrete_elements().cbegin(), geo_space.discrete_elements().cend(),
                     1u);

  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric;
  geo_metric.allocate_buffer(geo_space.discrete_elements().cbegin(),
                             geo_space.discrete_elements().cend(), 1u);

  interpolation::SolutionCache sol_cache;
  // This should be rather
  // sol_cache.allocate(sol_space.reference_elements(), 1u, f.nb_fields());
  sol_cache.allocate(sol_space.reference_elements().cbegin(), sol_space.reference_elements().cend(),
                     1u, 1u);

  interpolation::SolutionSpaceMetric<MeshConfig> sol_metric;
  // This should be rather
  // sol_metric.allocate_buffer(sol_space.reference_elements(), 1u,
  // f.nb_fields());
  sol_metric.allocate_buffer(sol_space.reference_elements().cbegin(),
                             sol_space.reference_elements().cend(), 1u, 1u);

  Real sum_cell_contributions = 0.0;

  const auto geo_sf_generator   = geo_space.sf_generator();
  const auto geo_quad_generator = geo_space.quad_generator();

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < geo_cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity geo_cell = geo_cell_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell = sol_cell_dofs.active_cell(mesh::ActiveIdx(c));

    geo_cache.flush();

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), geo_cell.pt_set_id(), tcell_view.coordinates());

    const mesh::PointSetTag geo_pt_set = geo_cell.pt_set_id();
    const ElemShape elem_shape         = geo_pt_set.elem_shape();
    const Uint geo_order               = geo_pt_set.poly_order();

    const mesh::PointSetTagExt geo_pt_set_ext(geo_pt_set, P0, mesh::CellTransform::NO_TRANS, 0u);
    const mesh::sf::SFTag geo_sf = geo_sf_generator(elem_shape, geo_order);

    const mesh::PointSetTag geo_quad = geo_quad_generator(elem_shape, geo_order);
    const mesh::PointSetTagExt geo_quad_ext(geo_quad, P0, mesh::CellTransform::NO_TRANS, 0u);

    const mesh::DiscreteElemKey geo_key(geo_pt_set_ext, geo_sf, geo_quad_ext);

    geo_cache.push_back_to_buffer(cell_coords, geo_key);

    sol_cache.flush();
    sol_cache.push_back_to_buffer(
        sol_cell, f,
        mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));

    geo_metric.empty_buffer();
    geo_metric.evaluate(geo_cache, interpolation::RebuildMetricIndex{true});

    sol_metric.empty_buffer();
    sol_metric.evaluate(geo_metric, sol_cache, interpolation::ComputeMetricDerivs{true},
                        interpolation::RebuildMetricIndex{true});

    typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>::cellwise_metric
        cell_geo_met = geo_metric.cellwise_values(0);
    typename interpolation::SolutionSpaceMetric<MeshConfig>::cellwise_metric cell_sol_met =
        sol_metric.cellwise_values(0);

    // Coordinates of quadrature points in physical space
    // const math::ConstMatrixBlock<Real> qd_pts =
    // cell_geo_met.interpolated_coords();

    // Vector of integration weights in quadrature points
    const math::DenseDVec<Real> &wq = cell_geo_met.pt_weights();

    // Vector of transformation jacobians in quadrature points
    const math::DenseConstVecView<Real> jq = cell_geo_met.jdet();

    // Values of solution in quadrature points
    const math::DenseConstMatView<Real> sol_q = cell_sol_met.field_values();

    for (Uint q = 0; q < cell_geo_met.nb_qd_pts(); ++q)
    {
      const Real value2 = sol_q(q, 0) * sol_q(q, 0);
      sum_cell_contributions += wq[q] * jq[q] * value2;
    }
  } // Loop over cells
  return std::sqrt(sum_cell_contributions);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename PhysModel>
void PostprocessingUtils<MeshConfig>::compute_Euler_Ma(
    const typename pdekit::result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
    const interpolation::VectorMeshFunction<Real> &solution,
    interpolation::ScalarMeshFunction<Real> &f)
{
  PhysModel physical_model;
  typename PhysModel::Properties properties;
  typename PhysModel::Properties::SolGradM gradient_matrix;

  mesh::adapt::LocalInterpolator loc_interpolator;

  typename PhysModel::Properties::CoordV one_node_coord;

  for (Uint c = 0; c < sol_cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity active_cell = sol_cell_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> sol_dof_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), active_cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < active_cell.nb_vert(); ++n)
    {
      for (Uint d = 0; d < PhysModel::Properties::DIM; ++d)
      {
        one_node_coord[d] = sol_dof_coords(n, d);
      }

      physical_model.compute_properties(one_node_coord, solution.const_value(active_cell.vertex(n)),
                                        gradient_matrix, properties);
      f[active_cell.vertex(n)] = properties.Ma;
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename PhysModel>
void PostprocessingUtils<MeshConfig>::compute_Euler_entropy(
    const typename pdekit::result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
    const interpolation::VectorMeshFunction<Real> &solution, const Real rho_in, const Real p_in,
    interpolation::ScalarMeshFunction<Real> &f)
{
  PhysModel physical_model;
  typename PhysModel::Properties properties;
  typename PhysModel::Properties::SolGradM gradient_matrix;

  mesh::adapt::LocalInterpolator loc_interpolator;

  typename PhysModel::Properties::CoordV one_node_coord;

  for (Uint c = 0; c < sol_cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity active_cell = sol_cell_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> sol_dof_coords = loc_interpolator.transfer_coords(
        tcell_view.cell_type(), active_cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < active_cell.nb_vert(); ++n)
    {
      for (Uint d = 0; d < PhysModel::Properties::DIM; ++d)
      {
        one_node_coord[d] = sol_dof_coords(n, d);
      }

      physical_model.compute_properties(one_node_coord, solution.const_value(active_cell.vertex(n)),
                                        gradient_matrix, properties);
      f[active_cell.vertex(n)] =
          (properties.P / p_in) / (std::pow(properties.rho / rho_in, 1.4)) - 1.;
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename PhysModel, typename Evaluator>
void PostprocessingUtils<MeshConfig>::compute_Euler_quantity(
    const typename pdekit::result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
    const interpolation::VectorMeshFunction<Real> &solution, const Evaluator &eval,
    interpolation::ScalarMeshFunction<Real> &f)
{
  PhysModel physical_model;
  typename PhysModel::Properties properties;
  typename PhysModel::Properties::SolGradM gradient_matrix;

  mesh::adapt::LocalInterpolator loc_interpolator;

  typename PhysModel::Properties::CoordV one_node_coord;

  for (Uint c = 0; c < sol_cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity active_cell = sol_cell_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> sol_dof_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), active_cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < active_cell.nb_vert(); ++n)
    {
      for (Uint d = 0; d < PhysModel::Properties::DIM; ++d)
      {
        one_node_coord[d] = sol_dof_coords(n, d);
      }

      physical_model.compute_properties(one_node_coord, solution.const_value(active_cell.vertex(n)),
                                        gradient_matrix, properties);
      f[active_cell.vertex(n)] = eval(properties);
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void extract_boundary_data(const mesh::Tria<MeshConfig> &mesh,
                           const typename result_of::dof_map_t<MeshConfig> &dofs_in,
                           const interpolation::VectorMeshFunction<Real> &data_in,
                           const std::string &boundary_name,
                           interpolation::VectorMeshFunction<Real> &coords,
                           interpolation::VectorMeshFunction<Real> &data_out)

{
  const typename result_of::mesh_boundary_set_t<MeshConfig> &boundary_set = mesh.all_boundaries();
  const std::shared_ptr<mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdry_domain =
      boundary_set.domain(boundary_name);

  // Count number of entities on boundary
  Uint nb_bdry_dofs = 0;

  for (Uint bc = 0; bc < bdry_domain->nb_active_cells(); ++bc)
  {
    const mesh::MeshEntity bcell = bdry_domain->active_cell(dofs_in, mesh::ActiveIdx(bc));
    nb_bdry_dofs += bcell.nb_vert();
  }

  coords.resize(MeshConfig::GDIM, nb_bdry_dofs);
  data_out.resize(data_in.nb_fields(), nb_bdry_dofs);

  nb_bdry_dofs = 0;

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint bc = 0; bc < bdry_domain->nb_active_cells(); ++bc)
  {
    const mesh::CellTopologyView<MeshConfig> adj_tcell_view =
        bdry_domain->adjacent_tcell(mesh::ActiveIdx(bc));
    const mesh::MeshEntity bcell = bdry_domain->active_cell(dofs_in, mesh::ActiveIdx(bc));

    const Uint loc_id = bdry_domain->local_id(mesh::ActiveIdx(bc));

    const math::DenseConstMatView<Real> bcell_coords = loc_interpolator.transfer_coords(
        adj_tcell_view.pt_set_id(MeshConfig::TDIM - 1, loc_id), bcell.pt_set_id(),
        adj_tcell_view.coordinates(MeshConfig::TDIM - 1, loc_id));

    for (Uint v = 0; v < bcell.nb_vert(); ++v)
    {
      const math::DenseConstVecView<Real> node_coords = bcell_coords.row_transpose(v);
      coords.insert_value(nb_bdry_dofs, node_coords);
      data_out.insert_value(nb_bdry_dofs, data_in.const_value(bcell.vertex(v)));
      nb_bdry_dofs++;
    }
  }
}

// ----------------------------------------------------------------------------

} // namespace solver

} // namespace pdekit

#endif
