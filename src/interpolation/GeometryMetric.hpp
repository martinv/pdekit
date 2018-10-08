#ifndef PDEKIT_Interpolation_Geometry_Metric_hpp
#define PDEKIT_Interpolation_Geometry_Metric_hpp

#include <set>

#include "common/DataMap.hpp"
#include "common/Meta.hpp"
#include "interpolation/CellGeoMetric.hpp"
#include "interpolation/FEValues.hpp"
#include "interpolation/GeometryCache.hpp"
#include "interpolation/MetricFlags.hpp"
#include "mesh/point_set/QuadratureAdaptTransformAlgoFactory.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim = TopoDim>
class GeometryMetric
{
  private:
  enum
  {
    evaluate_boundary_normals = GeoDim > EvalDim ? 1 : 0
  };

  struct MetricData
  {
    enum
    {
      GDIM = GeoDim
    };
    enum
    {
      TDIM = EvalDim
    };

    /// Type of this element
    mesh::DiscreteElemKey key;

    /// Finite element values - SF evaluated at quadrature points,
    /// quadrature point coordinates and weights etc.
    interpolation::FEValues fe_values;

    /// How many elements have been filled
    Uint nb_elem_filled;

    /// Values after interpolation (OUTPUT)
    math::DenseDMat<Real> m_coord_values;

    /// Array of matrices of shape function derivatives in physical space
    std::array<math::DenseDMat<Real>, TopoDim> m_sf_derivatives;

    /// Array of matrices of derivatives dx/dxi
    std::array<math::DenseDMat<Real>, TopoDim> m_coord_derivatives;

    /// Array of Jacobi matrices in all integration points
    std::vector<math::DenseSMat<Real, TopoDim, EvalDim>> invJ;

    /// Array of jacobians (determinants of Jacobi matrices)
    std::vector<Real> det_j;

    /// Matrix of normals in quadrature points
    math::DenseDMat<Real> m_normals;

    /// One normal on a facet in reference space
    math::DenseSVec<Real, GeoDim> m_ref_normal_on_restriction;
  };

  public:
  using cellwise_metric =
      typename common::SelectType<evaluate_boundary_normals, CellGeoMetricWithNormals<MetricData>,
                                  CellGeoMetric<MetricData>>::type;

  using const_std_reg_iterator = std::set<mesh::PointSetTag>::const_iterator;

  /// Default constructor
  GeometryMetric();

  /// Default destructor
  ~GeometryMetric();

  /// Return the number of elements that fit in buffer
  Uint max_nb_blocks_in_buffer() const;

  /// Set the size of buffer
  template <typename DiscreteElemKeyIterator>
  void allocate_buffer(const DiscreteElemKeyIterator keys_begin,
                       const DiscreteElemKeyIterator keys_end, Uint const nb_blocks);

  /// Empty the data in the buffer
  void empty_buffer();

  /// Clear all internal data
  void clear();

  /// Return nb. values that have been pushed to buffer so far
  Uint nb_values_in_buffer() const;

  /// Collocation - interpolate the values
  void evaluate(GeometryCache<GeoDim> const &gc, const RebuildMetricIndex rebuild_metric_idx);

  /// Get one block of values
  cellwise_metric const cellwise_values(const Uint idx) const;

  /// Print all types this metric can process
  void print_types() const;

  private:
  /// TYPES

  using metric_data_map_t = common::DataMap<mesh::DiscreteElemKey, MetricData>;

  /// This class actually knows how to evaluate the derivatives, jacobians ...
  /// It is a helper class visible only to GeometryMetric. Client code
  /// of GeometryMetric will never need to use it directly.

  template <Uint GEODIM, Uint TDIM, Uint Dummy>
  class MetricComputer
  {
  };

  /// DATA

  /// Maximum number of cell entries the buffer can take
  Uint m_nb_blocks;

  /// A map of buffers: for each element type we have one buffer
  metric_data_map_t m_metric_data_map;

  /// Linear index of metric data, consists of a vector of pairs
  /// [PointSetTag,Uint]
  /// Value m_mdata_index[i] is a pair [FIRST,SECOND] and corresponds to
  /// metric data which can be found in
  /// m_metric_data_map.std_region_data(FIRST), cell block on position
  /// [SECOND] in this particular MetricData
  std::vector<std::pair<common::PtrHandle<MetricData>, Uint>> m_mdata_index;
};

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
GeometryMetric<GeoDim, TopoDim, EvalDim>::GeometryMetric()
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
GeometryMetric<GeoDim, TopoDim, EvalDim>::~GeometryMetric()
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
Uint GeometryMetric<GeoDim, TopoDim, EvalDim>::max_nb_blocks_in_buffer() const
{
  return m_nb_blocks;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
template <typename DiscreteElemKeyIterator>
void GeometryMetric<GeoDim, TopoDim, EvalDim>::allocate_buffer(
    const DiscreteElemKeyIterator keys_begin, const DiscreteElemKeyIterator keys_end,
    Uint const nb_blocks)
{
  /// NOTE: for discontinuous method, there is no need to allocate
  /// m_coord_buffer and copy data into it, but instead m_coord_buffer could
  /// be just a MatrixBlock<Real>. To make this class work with both
  /// continuous and discontinuous fields, each metric_data should have an
  /// array of 'raw' values and m_coord_buffer should become a
  /// MatrixBlock<Real>. Then the procedure would be as follows: I)  For
  /// continous method, 'raw' values would allocate sufficient memory
  ///     and let m_coord_buffer be a proxy to this block
  /// II) For discontinous method, no 'raw' values would be allocated, but
  /// only m_coord_buffer would serve as proxy to wherever the data is stored.
  /// Since the data is already discontinous, it does not need to be copied,
  /// but the work can be directly done on it by means of the proxy.

  m_metric_data_map.clear();

  m_nb_blocks = nb_blocks;

  DiscreteElemKeyIterator it = keys_begin;
  Uint n_elem_types          = 0;

  for (; it != keys_end; ++it)
  {
    n_elem_types++;
  }

  m_mdata_index.reserve(n_elem_types * nb_blocks);
  m_mdata_index.resize(0);

  const Uint nb_coord_per_node = GeoDim;

  mesh::StdRegion std_region;
  mesh::StdPointSet point_set;

  it = keys_begin;

  for (; it != keys_end; ++it)
  {
    common::PtrHandle<MetricData> met_data = m_metric_data_map.create(*it);

    (*met_data).key = *it;
    (*met_data).fe_values.configure(it->support().std_region_tag(), it->basis());

    point_set.change_type(it->eval_pts().std_region_tag());

    const mesh::CellTransform adapt_op_id = it->eval_pts().cell_transform_id();
    const Uint local_id                   = it->eval_pts().local_id();

    // if (adapt_op_id == mesh::CellTransform::NO_TRANS)
    if (!mesh::CellTransformTraits::is_refinement(adapt_op_id))
    {
      (*met_data).fe_values.fill_Vandermonde(point_set.get().coordinates(local_id),
                                             point_set.get().weights());
    }
    else
    {
      QuadratureAdaptTransformAlgoFactory::instance_type &quad_adapt_algo_fact =
          QuadratureAdaptTransformAlgoFactory::instance();
      const QuadratureAdaptTransformAlgoFactory::instance_type::const_product_base_ptr
          algo_transform = quad_adapt_algo_fact.create(adapt_op_id);

      math::DenseDMat<Real> transf_quad_coords;
      math::DenseDVec<Real> transf_quad_weights;

      algo_transform->compute_transformed_coords(point_set.get().coordinates(), local_id,
                                                 transf_quad_coords);
      algo_transform->compute_transformed_weights(point_set.get().weights(), local_id,
                                                  transf_quad_weights);

      (*met_data).fe_values.fill_Vandermonde(transf_quad_coords, transf_quad_weights);
    }

    (*met_data).nb_elem_filled = 0;

    // These variables are here just to make
    // the code below shorter to write and easier to read
    const Uint nb_dof_per_elem    = (*met_data).fe_values.nb_nodes();
    const Uint nb_qd_pts_per_elem = (*met_data).fe_values.nb_qd_pts();

    (*met_data).m_coord_values.resize(nb_qd_pts_per_elem, nb_blocks * nb_coord_per_node);
    (*met_data).m_coord_values.fill(0.0);

    for (Uint d = 0; d < TopoDim; ++d)
    {
      (*met_data).m_sf_derivatives[d].resize(nb_qd_pts_per_elem, nb_blocks * nb_dof_per_elem);
      (*met_data).m_sf_derivatives[d].fill(0.0);

      (*met_data).m_coord_derivatives[d].resize(nb_qd_pts_per_elem, nb_blocks * nb_coord_per_node);
      (*met_data).m_coord_derivatives[d].fill(0.0);
    }

    (*met_data).invJ.resize(nb_qd_pts_per_elem * m_nb_blocks);
    (*met_data).det_j.resize(nb_qd_pts_per_elem * m_nb_blocks);

    if (evaluate_boundary_normals || ((*met_data).key.eval_pts().cell_transform_id() ==
                                      mesh::CellTransform::RESTRICT_TO_CODIM_1))
    {
      (*met_data).m_normals.resize(nb_qd_pts_per_elem, nb_blocks * nb_coord_per_node);
    }

    if ((*met_data).key.eval_pts().cell_transform_id() == mesh::CellTransform::RESTRICT_TO_CODIM_1)
    {
      std_region.change_type((*met_data).key.support().std_region_tag());
      (*met_data).m_ref_normal_on_restriction = std_region.get().facet_normal(local_id);
    }
  }
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
void GeometryMetric<GeoDim, TopoDim, EvalDim>::empty_buffer()
{
  m_mdata_index.resize(0);

  for (typename metric_data_map_t::iterator it = m_metric_data_map.begin();
       it != m_metric_data_map.end(); ++it)
  {
    (*it.data_ptr()).nb_elem_filled = 0;
  }
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
void GeometryMetric<GeoDim, TopoDim, EvalDim>::clear()
{
  m_nb_blocks = 0;

  m_metric_data_map.clear();

  m_mdata_index.clear();
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
Uint GeometryMetric<GeoDim, TopoDim, EvalDim>::nb_values_in_buffer() const
{
  return m_mdata_index.size();
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
void GeometryMetric<GeoDim, TopoDim, EvalDim>::evaluate(GeometryCache<GeoDim> const &gc,
                                                        const RebuildMetricIndex rebuild_metric_idx)
{
  if (rebuild_metric_idx)
  {
    // Rebuild the data index for fast access to metric data
    m_mdata_index.resize(0);

    for (Uint i = 0; i < gc.nb_values_in_buffer(); ++i)
    {
      common::PtrHandle<MetricData> md = m_metric_data_map.std_region_data(gc.key(i));
      (*md).nb_elem_filled++;
      m_mdata_index.push_back(
          std::pair<common::PtrHandle<MetricData>, Uint>(md, gc.position_in_block(i)));
    }
  }

  for (typename metric_data_map_t::iterator it = m_metric_data_map.begin();
       it != m_metric_data_map.end(); ++it)
  {
    if ((*it.data_ptr()).nb_elem_filled > 0)
    {
      math::DenseConstMatView<Real> const nodal_coords = gc.buffer_data((*it.data_ptr()).key);
      MetricComputer<GeoDim, EvalDim, 0>::evaluate_derivatives_and_jacobians(nodal_coords,
                                                                             (*it.data_ptr()));
    }
  }
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
typename GeometryMetric<GeoDim, TopoDim, EvalDim>::cellwise_metric const GeometryMetric<
    GeoDim, TopoDim, EvalDim>::cellwise_values(const Uint idx) const
{
  cellwise_metric cm(m_mdata_index[idx].first, m_mdata_index[idx].second);
  return cm;
}

// ----------------------------------------------------------------------------
// Specializations of GeometryMetric::MetricComputer for some combinations
// of dimensions  - here 2D volume metric (GDIM = 2, TDIM = 2)
// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
template <Uint Dummy>
class GeometryMetric<GeoDim, TopoDim, EvalDim>::MetricComputer<_2D, _2D, Dummy>
{
  public:
  static void evaluate_derivatives_and_jacobians(
      const math::DenseConstMatView<Real> &coords,
      typename GeometryMetric<GeoDim, TopoDim, EvalDim>::MetricData &data)
  {

    math::DenseSMat<Real, _2D, _2D> J;
    // Derivatives of one shape function in physical and
    // reference space
    math::DenseSVec<Real, _2D> dPhys, dRef;

    // Step 1: compute values and derivatives in interpolation points,
    // values of derivatives only in reference space

    data.m_coord_values = data.fe_values.Vandermonde() * coords;

    math::DenseDMat<Real> &coord_deriv_XI0 = data.m_coord_derivatives[XI0];
    math::DenseDMat<Real> &coord_deriv_XI1 = data.m_coord_derivatives[XI1];

    // First use the 'derivatives' matrices to compute derivatives with
    // respect to reference space variables (xi,eta,zeta). This is because
    // dV_XI0 and dV_XI1 are derivative matrices in reference space
    math::DenseDMat<Real> const &dV_XI0 = data.fe_values.deriv_Vandermonde(XI0);
    math::DenseDMat<Real> const &dV_XI1 = data.fe_values.deriv_Vandermonde(XI1);

    coord_deriv_XI0 = dV_XI0 * coords;
    coord_deriv_XI1 = dV_XI1 * coords;

    // Step 2: compute derivatives in physical space

    math::DenseDMat<Real> &sf_deriv_XI0 = data.m_sf_derivatives[XI0];
    math::DenseDMat<Real> &sf_deriv_XI1 = data.m_sf_derivatives[XI1];

    const Uint nb_dof_per_elem    = data.fe_values.nb_nodes();
    const Uint nb_qd_pts_per_elem = data.fe_values.nb_qd_pts();

    for (Uint c = 0; c < data.nb_elem_filled; ++c)
    {
      const Uint coord_deriv_offset = GeoDim * c;

      const Uint sf_deriv_offset = nb_dof_per_elem * c;

      for (Uint q = 0; q < nb_qd_pts_per_elem; ++q)
      {
        J(XI0, X0) = coord_deriv_XI0(q, coord_deriv_offset + X0); // dX/dXi
        J(XI0, X1) = coord_deriv_XI0(q, coord_deriv_offset + X1); // dY/dXi

        J(XI1, X0) = coord_deriv_XI1(q, coord_deriv_offset + X0); // dX/dEta
        J(XI1, X1) = coord_deriv_XI1(q, coord_deriv_offset + X1); // dY/dEta

        data.det_j[nb_qd_pts_per_elem * c + q] = J(XI0, X0) * J(XI1, X1) - J(XI0, X1) * J(XI1, X0);
        const Real inv_det                     = 1.0 / data.det_j[nb_qd_pts_per_elem * c + q];

        math::DenseSMat<Real, _2D, _2D> &invJ = data.invJ[nb_qd_pts_per_elem * c + q];
        invJ(0, 0)                            = inv_det * J(1, 1);
        invJ(1, 1)                            = inv_det * J(0, 0);
        invJ(0, 1)                            = -inv_det * J(0, 1);
        invJ(1, 0)                            = -inv_det * J(1, 0);

        for (Uint n = 0; n < nb_dof_per_elem; ++n)
        {
          dRef[XI0] = dV_XI0(q, n);
          dRef[XI1] = dV_XI1(q, n);

          // dPhys = invJ * dRef;
          dPhys[X0] = invJ(0, 0) * dRef[XI0] + invJ(0, 1) * dRef[XI1];
          dPhys[X1] = invJ(1, 0) * dRef[XI0] + invJ(1, 1) * dRef[XI1];

          sf_deriv_XI0(q, sf_deriv_offset + n) = dPhys[X0];
          sf_deriv_XI1(q, sf_deriv_offset + n) = dPhys[X1];
        }
      } // Loop over quadrature points
    }   // Loop over cells

    if (data.key.eval_pts().cell_transform_id() == mesh::CellTransform::RESTRICT_TO_CODIM_1)
    {
      for (Uint c = 0; c < data.nb_elem_filled; ++c)
      {
        const Uint offset = c * GeoDim;

        for (Uint q = 0; q < nb_qd_pts_per_elem; ++q)
        {
          const math::DenseSMat<Real, _2D, _2D> &invJ = data.invJ[nb_qd_pts_per_elem * c + q];

          /*
          data.m_normals(q, offset + X0) = invJ(0, 0) *
          data.m_ref_normal_on_restriction[XI0] + invJ(1, 0) *
          data.m_ref_normal_on_restriction[XI1];

          data.m_normals(q, offset + X1) = invJ(0, 1) *
          data.m_ref_normal_on_restriction[XI0] + invJ(1, 1) *
          data.m_ref_normal_on_restriction[XI1];
          */
          data.m_normals(q, offset + X0) = invJ(0, 0) * data.m_ref_normal_on_restriction[XI0] +
                                           invJ(0, 1) * data.m_ref_normal_on_restriction[XI1];

          data.m_normals(q, offset + X1) = invJ(1, 0) * data.m_ref_normal_on_restriction[XI0] +
                                           invJ(1, 1) * data.m_ref_normal_on_restriction[XI1];

          const Real inv_norm =
              1. / std::sqrt(data.m_normals(q, offset + X0) * data.m_normals(q, offset + X0) +
                             data.m_normals(q, offset + X1) * data.m_normals(q, offset + X1));

          data.m_normals(q, offset + X0) *= inv_norm;
          data.m_normals(q, offset + X1) *= inv_norm;
        }
      }
    }

  } // Static method 'evaluate_derivatives and jacobians'
};

// ----------------------------------------------------------------------------
// Specializations of GeometryMetric::MetricComputer for some combinations
// of dimensions  - here 2D volume metric (GDIM = 3, TDIM = 3)
// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
template <Uint Dummy>
class GeometryMetric<GeoDim, TopoDim, EvalDim>::MetricComputer<_3D, _3D, Dummy>
{
  public:
  static void evaluate_derivatives_and_jacobians(
      const math::DenseConstMatView<Real> &coords,
      typename GeometryMetric<GeoDim, TopoDim, EvalDim>::MetricData &data)
  {

    math::DenseSMat<Real, _3D, _3D> J;
    // Derivatives of one shape function in physical and
    // reference space
    math::DenseSVec<Real, _3D> dPhys, dRef;

    // Step 1: compute values and derivatives in interpolation points,
    // values of derivatives only in reference space

    data.m_coord_values = data.fe_values.Vandermonde() * coords;

    math::DenseDMat<Real> &coord_deriv_XI0 = data.m_coord_derivatives[XI0];
    math::DenseDMat<Real> &coord_deriv_XI1 = data.m_coord_derivatives[XI1];
    math::DenseDMat<Real> &coord_deriv_XI2 = data.m_coord_derivatives[XI2];

    // First use the 'derivatives' matrices to compute derivatives with
    // respect to reference space variables (xi,eta,zeta). This is because
    // dV_XI0 and dV_XI1 are derivative matrices in reference space
    math::DenseDMat<Real> const &dV_XI0 = data.fe_values.deriv_Vandermonde(XI0);
    math::DenseDMat<Real> const &dV_XI1 = data.fe_values.deriv_Vandermonde(XI1);
    math::DenseDMat<Real> const &dV_XI2 = data.fe_values.deriv_Vandermonde(XI2);

    coord_deriv_XI0 = dV_XI0 * coords;
    coord_deriv_XI1 = dV_XI1 * coords;
    coord_deriv_XI2 = dV_XI2 * coords;

    // Step 2: compute derivatives in physical space
    math::DenseDMat<Real> &sf_deriv_XI0 = data.m_sf_derivatives[XI0];
    math::DenseDMat<Real> &sf_deriv_XI1 = data.m_sf_derivatives[XI1];
    math::DenseDMat<Real> &sf_deriv_XI2 = data.m_sf_derivatives[XI2];

    const Uint nb_dof_per_elem    = data.fe_values.nb_nodes();
    const Uint nb_qd_pts_per_elem = data.fe_values.nb_qd_pts();

    for (Uint c = 0; c < data.nb_elem_filled; ++c)
    {
      const Uint coord_deriv_offset = GeoDim * c;

      const Uint sf_deriv_offset = nb_dof_per_elem * c;

      for (Uint q = 0; q < nb_qd_pts_per_elem; ++q)
      {
        J(XI0, X0) = coord_deriv_XI0(q, coord_deriv_offset + X0); // dX/dXi
        J(XI0, X1) = coord_deriv_XI0(q, coord_deriv_offset + X1); // dY/dXi
        J(XI0, X2) = coord_deriv_XI0(q, coord_deriv_offset + X2); // dZ/dXi

        J(XI1, X0) = coord_deriv_XI1(q, coord_deriv_offset + X0); // dX/dEta
        J(XI1, X1) = coord_deriv_XI1(q, coord_deriv_offset + X1); // dY/dEta
        J(XI1, X2) = coord_deriv_XI1(q, coord_deriv_offset + X2); // dZ/dEta

        J(XI2, X0) = coord_deriv_XI2(q, coord_deriv_offset + X0); // dX/dZeta
        J(XI2, X1) = coord_deriv_XI2(q, coord_deriv_offset + X1); // dY/dZeta
        J(XI2, X2) = coord_deriv_XI2(q, coord_deriv_offset + X2); // dZ/dZeta

        data.det_j[nb_qd_pts_per_elem * c + q] = J.det();
        // const Real inv_det = 1.0/data.det_j[data.nb_qd_pts_per_elem *
        // c + q];

        math::DenseSMat<Real, _3D, _3D> &invJ = data.invJ[nb_qd_pts_per_elem * c + q];
        J.inv(invJ);

        for (Uint n = 0; n < nb_dof_per_elem; ++n)
        {
          dRef[XI0] = dV_XI0(q, n);
          dRef[XI1] = dV_XI1(q, n);
          dRef[XI2] = dV_XI2(q, n);

          // dPhys = invJ * dRef;
          dPhys[X0] = invJ(0, 0) * dRef[XI0] + invJ(0, 1) * dRef[XI1] + invJ(0, 2) * dRef[XI2];
          dPhys[X1] = invJ(1, 0) * dRef[XI0] + invJ(1, 1) * dRef[XI1] + invJ(1, 2) * dRef[XI2];
          dPhys[X2] = invJ(2, 0) * dRef[XI0] + invJ(2, 1) * dRef[XI1] + invJ(2, 2) * dRef[XI2];

          sf_deriv_XI0(q, sf_deriv_offset + n) = dPhys[X0];
          sf_deriv_XI1(q, sf_deriv_offset + n) = dPhys[X1];
          sf_deriv_XI2(q, sf_deriv_offset + n) = dPhys[X2];
        }
      } // Loop over quadrature points
    }   // Loop over cells

    if (data.key.eval_pts().cell_transform_id() == mesh::CellTransform::RESTRICT_TO_CODIM_1)
    {
      for (Uint c = 0; c < data.nb_elem_filled; ++c)
      {
        const Uint offset = c * GeoDim;

        for (Uint q = 0; q < nb_qd_pts_per_elem; ++q)
        {
          const math::DenseSMat<Real, _3D, _3D> &invJ = data.invJ[nb_qd_pts_per_elem * c + q];

          data.m_normals(q, offset + X0) = invJ(0, 0) * data.m_ref_normal_on_restriction[XI0] +
                                           invJ(0, 1) * data.m_ref_normal_on_restriction[XI1] +
                                           invJ(0, 2) * data.m_ref_normal_on_restriction[XI2];

          data.m_normals(q, offset + X1) = invJ(1, 0) * data.m_ref_normal_on_restriction[XI0] +
                                           invJ(1, 1) * data.m_ref_normal_on_restriction[XI1] +
                                           invJ(1, 2) * data.m_ref_normal_on_restriction[XI2];

          data.m_normals(q, offset + X2) = invJ(2, 0) * data.m_ref_normal_on_restriction[XI0] +
                                           invJ(2, 1) * data.m_ref_normal_on_restriction[XI1] +
                                           invJ(2, 2) * data.m_ref_normal_on_restriction[XI2];

          const Real inv_norm =
              1. / std::sqrt(data.m_normals(q, offset + X0) * data.m_normals(q, offset + X0) +
                             data.m_normals(q, offset + X1) * data.m_normals(q, offset + X1) +
                             data.m_normals(q, offset + X2) * data.m_normals(q, offset + X2));

          data.m_normals(q, offset + X0) *= inv_norm;
          data.m_normals(q, offset + X1) *= inv_norm;
          data.m_normals(q, offset + X2) *= inv_norm;
        }
      }
    }

  } // Static method 'evaluate_derivatives and jacobians'
};

// ----------------------------------------------------------------------------
// Specializations of GeometryMetric::MetricComputer for some combinations
// of dimensions  - here 1D line metric in 2D space (GDIM = 2, DIM = 1)
// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
template <Uint Dummy>
class GeometryMetric<GeoDim, TopoDim, EvalDim>::MetricComputer<_2D, _1D, Dummy>
{
  public:
  static void evaluate_derivatives_and_jacobians(
      const math::DenseConstMatView<Real> &coords,
      typename GeometryMetric<GeoDim, TopoDim, EvalDim>::MetricData &data)
  {
    math::DenseSVec<Real, _2D> normal;

    // Step 1: compute values and derivatives in interpolation points,
    // values of derivatives only
    //         in reference space
    // If no data have been filled, there's nothing to compute

    data.m_coord_values = data.fe_values.Vandermonde() * coords;

    math::DenseDMat<Real> &coord_deriv_XI0 = data.m_coord_derivatives[XI0];

    // First use the 'derivatives' matrices to compute derivatives with
    // respect to reference space variables (xi,eta,zeta).
    coord_deriv_XI0 = data.fe_values.deriv_Vandermonde(XI0) * coords;

    // Step 2: compute derivatives in physical space
    const Uint nb_qd_pts_per_elem = data.fe_values.nb_qd_pts();

    for (Uint c = 0; c < data.nb_elem_filled; ++c)
    {
      const Uint offset = c * GeoDim;

      for (Uint q = 0; q < nb_qd_pts_per_elem; ++q)
      {
        normal[X0] = coord_deriv_XI0(q, offset + X1);  // dY/dXi
        normal[X1] = -coord_deriv_XI0(q, offset + X0); // dX/dXi

        const Real norm     = std::sqrt(normal[X0] * normal[X0] + normal[X1] * normal[X1]);
        const Real inv_norm = 1.0 / norm;

        data.m_normals(q, offset + X0) = inv_norm * normal[X0];
        data.m_normals(q, offset + X1) = inv_norm * normal[X1];

        data.det_j[nb_qd_pts_per_elem * c + q] = norm;

      } // Loop over quadrature points
    }   // Loop over cells

  } // Static method 'evaluate_derivatives and jacobians'
};

// ----------------------------------------------------------------------------
// Specializations of GeometryMetric::MetricComputer for some combinations
// of dimensions  - here 2D surface metric in 3D space (GDIM = 3, DIM = 2)
// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
template <Uint Dummy>
class GeometryMetric<GeoDim, TopoDim, EvalDim>::MetricComputer<_3D, _2D, Dummy>
{
  public:
  static void evaluate_derivatives_and_jacobians(
      const math::DenseConstMatView<Real> &coords,
      typename GeometryMetric<GeoDim, TopoDim, EvalDim>::MetricData &data)
  {

    math::DenseSVec<Real, _3D> normal;

    // Step 1: compute values and derivatives in interpolation points,
    // values of derivatives only in reference space If no data have been
    // filled, there's nothing to compute

    data.m_coord_values = data.fe_values.Vandermonde() * coords;

    math::DenseDMat<Real> &coord_deriv_XI0 = data.m_coord_derivatives[XI0];
    math::DenseDMat<Real> &coord_deriv_XI1 = data.m_coord_derivatives[XI1];

    // First use the 'derivatives' matrices to compute derivatives with
    // respect to reference space variables (xi,eta,zeta).
    coord_deriv_XI0 = data.fe_values.deriv_Vandermonde(XI0) * coords;
    coord_deriv_XI1 = data.fe_values.deriv_Vandermonde(XI1) * coords;

    // Step 2: compute derivatives in physical space
    const Uint nb_qd_pts_per_elem = data.fe_values.nb_qd_pts();

    for (Uint c = 0; c < data.nb_elem_filled; ++c)
    {
      const Uint offset = c * GeoDim;

      for (Uint q = 0; q < nb_qd_pts_per_elem; ++q)
      {
        normal[X0] = coord_deriv_XI0(q, offset + X1) * coord_deriv_XI1(q, offset + X2) -
                     coord_deriv_XI0(q, offset + X2) * coord_deriv_XI1(q, offset + X1);

        normal[X1] = coord_deriv_XI0(q, offset + X2) * coord_deriv_XI1(q, offset + X0) -
                     coord_deriv_XI0(q, offset + X0) * coord_deriv_XI1(q, offset + X2);

        normal[X2] = coord_deriv_XI0(q, offset + X0) * coord_deriv_XI1(q, offset + X1) -
                     coord_deriv_XI0(q, offset + X1) * coord_deriv_XI1(q, offset + X0);

        const Real norm =
            std::sqrt(normal[X] * normal[X] + normal[Y] * normal[Y] + normal[Z] * normal[Z]);
        const Real inv_norm = 1.0 / norm;

        data.det_j[nb_qd_pts_per_elem * c + q] = norm;

        data.m_normals(q, offset + X0) = inv_norm * normal[X0];
        data.m_normals(q, offset + X1) = inv_norm * normal[X1];
        data.m_normals(q, offset + X2) = inv_norm * normal[X2];

      } // Loop over quadrature points
    }   // Loop over cells

  } // Static method 'evaluate_derivatives and jacobians'
};

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TopoDim, Uint EvalDim>
void GeometryMetric<GeoDim, TopoDim, EvalDim>::print_types() const
{
  for (typename metric_data_map_t::const_iterator it = m_metric_data_map.cbegin();
       it != m_metric_data_map.cend(); ++it)
  {
    const mesh::DiscreteElemKey key = it.key_value();
    std::cout << "[" << key << "], fill status: " << (*it.data_ptr()).nb_elem_filled << "/"
              << m_nb_blocks << " items" << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
