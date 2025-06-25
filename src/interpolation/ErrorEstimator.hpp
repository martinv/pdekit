#ifndef PDEKIT_Interpolation_Error_Estimator_hpp
#define PDEKIT_Interpolation_Error_Estimator_hpp

#include "interpolation/GeometryMetric.hpp"
#include "interpolation/SolutionSpaceMetric.hpp"

namespace pdekit
{

namespace interpolation
{

template <typename MeshConfig>
class ErrorEstimator
{
  public:
  /// Constructor
  ErrorEstimator(typename mesh::Tria<MeshConfig>::const_shared_ptr mesh);

  /// Compute error in nodes
  template <typename ExactSolution, typename T>
  void compute_pointwise_L2_error(const ExactSolution &exact_sol,
                                  const typename result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
                                  const VectorMeshFunction<T> &f) const;

  /// Compute L1 norm of error
  template <typename ExactSolution, typename T>
  void compute_L1_error(const ExactSolution &exact_sol, const mesh::Tria<MeshConfig> &mesh_in,
                        const typename result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
                        const VectorMeshFunction<T> &f) const;

  /// Compute L2 norm of error
  template <typename ExactSolution, typename T>
  void compute_L2_error(const ExactSolution &exact_sol, const mesh::Tria<MeshConfig> &mesh_in,
                        const typename result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
                        const VectorMeshFunction<T> &f) const;

  /// Compute infinity norm of error
  template <typename ExactSolution, typename T>
  void compute_infty_error(const ExactSolution &exact_sol, const mesh::Tria<MeshConfig> &mesh_in,
                           const typename result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
                           const VectorMeshFunction<T> &f) const;

  private:
  // Type of norm in which error is computer
  enum ErrorNorm
  {
    POINTWISE_L2 = 1,
    L1           = 2,
    L2           = 3,
    INFTY        = 4
  };

  /// The dummy template parameter is necessary because C++ does not allow
  /// full specialization of class nested in a template class
  /// Only partial specializations are allowed
  template <Uint ErrNorm, Uint Dummy>
  struct ErrorPointContribution;

  // --------------------------------------------------------------------------

  template <Uint Dummy>
  struct ErrorPointContribution<L1, Dummy>
  {
    inline static void contrib_at_one_point(const Real wj, const Real exact, const Real numerical,
                                            Real &current_norm)
    {
      current_norm += wj * std::abs(exact - numerical);
    }
  };

  // --------------------------------------------------------------------------

  template <Uint Dummy>
  struct ErrorPointContribution<L2, Dummy>
  {
    inline static void contrib_at_one_point(const Real wj, const Real exact, const Real numerical,
                                            Real &current_norm)
    {
      current_norm += wj * (exact - numerical) * (exact - numerical);
    }
  };

  // --------------------------------------------------------------------------

  template <Uint Dummy>
  struct ErrorPointContribution<INFTY, Dummy>
  {
    inline static void contrib_at_one_point(const Real wj, const Real exact, const Real numerical,
                                            Real &current_norm)
    {
      current_norm = std::max(current_norm, std::abs(exact - numerical));
    }
  };

  // --------------------------------------------------------------------------

  /// The actual implementation of norm computation
  template <typename ExactSolution, typename T, Uint NormType>
  void compute_error_impl(const ExactSolution &exact_sol, const mesh::Tria<MeshConfig> &mesh_in,
                          const typename result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
                          const VectorMeshFunction<T> &f, Real &sum_cell_contributions) const;

  template <typename DofIterator>
  void detect_discrete_elems(
      const common::IteratorRange<DofIterator> &dofs,
      const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
      const std::function<mesh::PointSetTag(const ElemShape shape, const Uint tcell_order,
                                            const Uint dofs_order)> &quad_generator,
      const bool apply_to_geometry, const bool use_filter,
      std::vector<mesh::DiscreteElemKey> &discrete_elems,
      common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> &fe_values) const;

  /// Mesh on which this error is computed
  typename mesh::Tria<MeshConfig> const &m_mesh;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
ErrorEstimator<MeshConfig>::ErrorEstimator(typename mesh::Tria<MeshConfig>::const_shared_ptr mesh)
    : m_mesh(*mesh)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename ExactSolution, typename T>
void ErrorEstimator<MeshConfig>::compute_pointwise_L2_error(
    const ExactSolution &exact_sol, const typename result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
    const VectorMeshFunction<T> &f) const
{

  //  mesh::CellConnectivity const & connectivity = m_mesh.topology()(dim);
  const Uint nb_nodes = sol_cell_dofs.nb_nodes();

  Real sum_squares = 0.0;

  std::vector<bool> node_processed(sol_cell_dofs.nb_nodes());
  node_processed.assign(sol_cell_dofs.nb_nodes(), false);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < sol_cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = sol_cell_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_cell_dofs.tcell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      if (!node_processed[cell.vertex(n)])
      {
        const math::DenseConstVecView<Real> point = cell_coords.row_transpose(n);

        typename VectorMeshFunction<T>::const_entry_type const f_node =
            f.const_value(cell.vertex(n));

        const Real delta = f_node[0] - exact_sol.value(point);
        sum_squares += delta * delta;

        node_processed[cell.vertex(n)] = true;
      }
    }
  }

  const Real h     = 1.0 / std::pow(nb_nodes, 1. / MeshConfig::TDIM);
  const Real error = std::pow(sum_squares, 1. / MeshConfig::TDIM) * h;
  std::cout.precision(14);

  std::cout << "#=============================================================="
               "======================="
            << std::endl;
  std::cout << "# n_nodes   h=1/NDOF^(1/dim)      abs(log h)      pointwise "
               "err.      log(pointwise)"
            << std::endl;

  std::cout << " " << nb_nodes << "       " << h << "   " << std::abs(std::log(h)) << "   " << error
            << "   ";
  std::cout << std::log(error) << std::endl;
  std::cout << "#=============================================================="
               "======================="
            << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename ExactSolution, typename T>
void ErrorEstimator<MeshConfig>::compute_L1_error(
    const ExactSolution &exact_sol, const mesh::Tria<MeshConfig> &mesh_in,
    const typename result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
    const VectorMeshFunction<T> &f) const
{
  Real l1_error = 0.0;
  compute_error_impl<ExactSolution, T, L1>(exact_sol, mesh_in, sol_cell_dofs, f, l1_error);

  const Uint nb_nodes = sol_cell_dofs.nb_nodes();

  const Real h = 1.0 / std::pow(nb_nodes, 1. / MeshConfig::TDIM);

  std::cout.precision(14);

  std::cout << "#=============================================================="
               "===================="
            << std::endl;
  std::cout << "# n_nodes   h=1/NDOF^(1/dim)      abs(log h)         L1 "
               "error         log(L1 error)"
            << std::endl;

  std::cout << " " << nb_nodes << "       " << h << "   " << std::abs(std::log(h)) << "   "
            << l1_error << "   ";
  std::cout << std::log(l1_error) << std::endl;
  std::cout << "#=============================================================="
               "===================="
            << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename ExactSolution, typename T>
void ErrorEstimator<MeshConfig>::compute_L2_error(
    const ExactSolution &exact_sol, const mesh::Tria<MeshConfig> &mesh_in,
    const typename result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
    const VectorMeshFunction<T> &f) const
{
  Real sum_error_squares = 0;
  compute_error_impl<ExactSolution, T, L2>(exact_sol, mesh_in, sol_cell_dofs, f, sum_error_squares);

  const Uint nb_nodes = sol_cell_dofs.nb_nodes();

  const Real h        = 1.0 / std::pow(nb_nodes, 1. / MeshConfig::TDIM);
  const Real l2_error = std::pow(sum_error_squares, 1. / 2.);

  std::cout.precision(14);

  std::cout << "#=============================================================="
               "===================="
            << std::endl;
  std::cout << "# n_nodes   h=1/NDOF^(1/dim)      abs(log h)         L2 "
               "error         log(L2 error)"
            << std::endl;

  std::cout << " " << nb_nodes << "       " << h << "   " << std::abs(std::log(h)) << "   "
            << l2_error << "   ";
  std::cout << std::log(l2_error) << std::endl;
  std::cout << "#=============================================================="
               "===================="
            << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename ExactSolution, typename T>
void ErrorEstimator<MeshConfig>::compute_infty_error(
    const ExactSolution &exact_sol, const mesh::Tria<MeshConfig> &mesh_in,
    const typename result_of::dof_map_t<MeshConfig> &sol_cell_dofs,
    const VectorMeshFunction<T> &f) const
{
  Real infty_norm = 0;
  compute_error_impl<ExactSolution, T, INFTY>(exact_sol, mesh_in, sol_cell_dofs, f, infty_norm);

  const Uint nb_nodes = sol_cell_dofs.nb_nodes();

  const Real h = 1.0 / std::pow(nb_nodes, 1. / MeshConfig::TDIM);

  std::cout.precision(14);

  std::cout << "#=============================================================="
               "===================="
            << std::endl;
  std::cout << "# n_nodes   h=1/NDOF^(1/dim)      abs(log h)       infty "
               "error         log(infty error)"
            << std::endl;

  std::cout << " " << nb_nodes << "       " << h << "   " << std::abs(std::log(h)) << "   "
            << infty_norm << "   ";
  std::cout << std::log(infty_norm) << std::endl;
  std::cout << "#=============================================================="
               "===================="
            << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename ExactSolution, typename T, Uint NormType>
void ErrorEstimator<MeshConfig>::compute_error_impl(
    const ExactSolution &exact_sol, const mesh::Tria<MeshConfig> &mesh_in,
    const typename result_of::dof_map_t<MeshConfig> &sol_cell_dofs, const VectorMeshFunction<T> &f,
    Real &sum_cell_contributions) const
{
  using dof_type_t = result_of::dof_map_t<MeshConfig>;
  const common::IteratorRange<typename dof_type_t::const_dof_iterator> dofs_range(
      sol_cell_dofs.cbegin(), sol_cell_dofs.cend());

  auto sf_generator = [](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };
  auto quad_generator = [](const ElemShape shape, const Uint tcell_order, const Uint dofs_order) {
    return mesh::PointSetTag(shape, std::max(tcell_order, dofs_order), PointSetID::Gauss);
  };

  std::vector<mesh::DiscreteElemKey> discrete_elements_mesh;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> fe_values_mesh;

  std::vector<mesh::DiscreteElemKey> discrete_elements_dofs;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> fe_values_dofs;

  detect_discrete_elems(dofs_range, sf_generator, quad_generator, true, true,
                        discrete_elements_mesh, fe_values_mesh);
  detect_discrete_elems(dofs_range, sf_generator, quad_generator, false, true,
                        discrete_elements_dofs, fe_values_dofs);

  interpolation::GeometryCache<MeshConfig::GDIM> geo_cache;
  geo_cache.allocate(discrete_elements_mesh.cbegin(), discrete_elements_mesh.cend(), 1u);

  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric;
  geo_metric.allocate_buffer(discrete_elements_mesh.cbegin(), discrete_elements_mesh.cend(), 1u);

  interpolation::SolutionCache sol_cache;
  // This should be rather
  sol_cache.allocate(fe_values_dofs.cbegin(), fe_values_dofs.cend(), 1u, 1u);

  interpolation::SolutionSpaceMetric<MeshConfig> sol_metric;
  // This should be rather
  sol_metric.allocate_buffer(fe_values_dofs.cbegin(), fe_values_dofs.cend(), 1u, 1u);

  sum_cell_contributions = 0.0;

  for (Uint c = 0; c < mesh_in.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = mesh_in.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell                = sol_cell_dofs.active_cell(mesh::ActiveIdx(c));

    const ElemShape cell_shape = tcell.cell_shape();
    const Uint cell_geo_order  = tcell.pt_set_id().poly_order();

    mesh::PointSetTagExt cell_type_tag_ext(tcell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS,
                                           0u);
    mesh::sf::SFTag const sf_type_tag = sf_generator(cell_shape, cell_geo_order);

    const mesh::PointSetTag quad_tag =
        quad_generator(cell_shape, cell_geo_order, sol_cell.pt_set_id().poly_order());
    const mesh::PointSetTagExt quad_tag_ext(quad_tag, P0, mesh::CellTransform::NO_TRANS, 0u);

    const mesh::DiscreteElemKey cell_key(cell_type_tag_ext, sf_type_tag, quad_tag_ext);

    const mesh::CellGeometry<MeshConfig::GDIM> cell_coords = tcell.coordinates();

    geo_cache.flush();
    geo_cache.push_back_to_buffer(cell_coords, cell_key);

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
    const math::DenseConstMatView<Real> qd_pts = cell_geo_met.interpolated_coords();

    // Vector of integration weights in quadrature points
    const math::DenseDVec<Real> &wq = cell_geo_met.pt_weights();

    // Vector of transformation jacobians in quadrature points
    const math::DenseConstVecView<Real> jq = cell_geo_met.jdet();

    // Values of solution in quadrature points
    const math::DenseConstMatView<Real> sol_q = cell_sol_met.field_values();

    for (Uint q = 0; q < cell_geo_met.nb_qd_pts(); ++q)
    {
      ErrorPointContribution<NormType, 0u>::contrib_at_one_point(
          wq[q] * jq[q], exact_sol.value(qd_pts.row_transpose(q)), sol_q(q, 0),
          sum_cell_contributions);
    }
  } // Loop over cells
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename DofIterator>
void ErrorEstimator<MeshConfig>::detect_discrete_elems(
    const common::IteratorRange<DofIterator> &dofs,
    const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
    const std::function<mesh::PointSetTag(const ElemShape shape, const Uint tcell_order,
                                          const Uint dofs_order)> &quad_generator,
    const bool apply_to_geometry, const bool use_filter,
    std::vector<mesh::DiscreteElemKey> &discrete_elems,
    common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> &fe_values) const
{
  discrete_elems.clear();
  fe_values.clear();

  mesh::StdPointSet quadrature;

  for (DofIterator dof_iter = dofs.begin(); dof_iter != dofs.end(); ++dof_iter)
  {
    const mesh::PointSetTag mesh_std_reg_tag = dof_iter->geo_pt_set_id();
    const mesh::PointSetTag dof_std_reg_tag  = dof_iter->pt_set_id();

    const ElemShape shape = dof_std_reg_tag.elem_shape();
    const Uint generator_poly_deg =
        apply_to_geometry ? mesh_std_reg_tag.poly_order() : dof_std_reg_tag.poly_order();

    mesh::sf::SFTag const sf_type_tag = sf_generator(shape, generator_poly_deg);

    const mesh::PointSetTag quad_tag =
        quad_generator(shape, mesh_std_reg_tag.poly_order(), dof_std_reg_tag.poly_order());

    quadrature.change_type(quad_tag);

    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      const mesh::PointSetTag generator_tag =
          apply_to_geometry ? mesh_std_reg_tag : dof_std_reg_tag;
      const mesh::PointSetTagExt ext_tag(generator_tag, P0, mesh::CellTransform::NO_TRANS, ie);
      common::PtrHandle<FEValues> elem_ptr = fe_values.std_region_data(ext_tag);

      if (elem_ptr.is_null())
      {
        elem_ptr = fe_values.create(ext_tag);
        (*elem_ptr).configure(generator_tag, sf_type_tag);

        (*elem_ptr).fill_Vandermonde(quadrature.get().coordinates(ie), quadrature.get().weights(ie),
                                     use_filter);
      }
    }
    // Loop that generates new discrete elements
    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      const mesh::PointSetTag std_reg_generator_tag =
          apply_to_geometry ? mesh_std_reg_tag : dof_std_reg_tag;
      const mesh::PointSetTagExt std_reg_generator_tag_ext(std_reg_generator_tag, P0,
                                                           mesh::CellTransform::NO_TRANS, ie);
      const mesh::PointSetTagExt quad_tag_ext(quad_tag, P0, mesh::CellTransform::NO_TRANS, 0);

      const mesh::DiscreteElemKey candidate(std_reg_generator_tag_ext, sf_type_tag, quad_tag_ext);

      mesh::add_unique_discr_elem_key(discrete_elems, candidate);
    }

  } // Iterator loop over dofs
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
