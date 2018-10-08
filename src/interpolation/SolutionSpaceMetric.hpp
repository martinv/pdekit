#ifndef PDEKIT_Interpolation_Solution_Space_Metric_hpp
#define PDEKIT_Interpolation_Solution_Space_Metric_hpp

#include "common/DataMap.hpp"
#include "common/Meta.hpp"
#include "interpolation/CellSolutionMetric.hpp"
#include "interpolation/FEValues.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/MetricFlags.hpp"
#include "interpolation/SolutionCache.hpp"
#include "math/DenseMatView.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM = MeshConfig::TDIM>
class SolutionSpaceMetric
{
  private:
  // enum { evaluate_boundary_normals = MeshConfig::GDIM > DIM ? 1 : 0  };

  /// TYPES
  struct MetricData
  {
    enum
    {
      GDIM = MeshConfig::GDIM
    };
    enum
    {
      TDIM = DIM
    };

    /// Type of this element
    mesh::PointSetTagExt std_region_type;

    /// Finite element values - SF evaluated at quadrature points,
    /// quadrature point coordinates and weights etc.
    interpolation::FEValues fe_values;

    /// How many elements have been filled
    Uint nb_elem_filled;

    /// Values after interpolation (OUTPUT)
    math::DenseDMat<Real> m_values;

    /// Array of matrices of derivatives in PHYSICAL space
    std::array<math::DenseDMat<Real>, DIM> m_derivatives;
  };

  public:
  /*
  typedef typename common::SelectType<evaluate_boundary_normals,
                                      CellwiseMetricWithNormals<MeshConfig,DIM>,
                                      CellSolutionMetric<MeshConfig,DIM>
  >::type cellwise_metric;
  */

  using cellwise_metric = CellSolutionMetric<MetricData>;

  /// Default constructor
  SolutionSpaceMetric();

  /// Default destructor
  ~SolutionSpaceMetric();

  /// Return the number of elements that fit in buffer
  Uint max_nb_blocks_in_buffer() const;

  /// Return the number of fields per dof
  Uint nb_fields() const;

  /// Set the size of buffer
  template <typename FeValsIterator>
  void allocate_buffer(const FeValsIterator fe_vals_begin, const FeValsIterator fe_vals_end,
                       Uint const nb_blocks, Uint const nb_fields);

  /// Empty the data in the buffer
  void empty_buffer();

  /// Clear all internal data
  void clear();

  /// Return nb. values that have been pushed to buffer so far
  Uint nb_values_in_buffer() const;

  /// Collocation - interpolate the values and their derivatives
  void evaluate(GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &gmetric,
                SolutionCache const &sc, const ComputeMetricDerivs compute_derivatives,
                const RebuildMetricIndex rebuild_index);

  /// Get one block of values
  cellwise_metric const cellwise_values(const Uint idx) const;

  /// Print some info about data contained in the solution space metric
  void print_fill_status() const;

  private:
  /// TYPES

  /// This class actually knows how to evaluate the derivatives, jacobians ...
  /// It is a helper class visible only to FunctionSpaceMetric. Client code
  /// of FunctionSpaceMetric will never need to use it directly.

  template <Uint GEODIM, Uint TDIM, Uint Dummy>
  class MetricComputer
  {
  };

  /// DATA

  /// Maximum number of cell entries the buffer can take
  Uint m_nb_blocks;

  /// How many fields per dof do we consider.
  /// Example: 2D Euler equations will require 4 fields per dof
  /// (density, momentum x, momentum y, energy)
  /// Another example: space coordinates (in 3D, for example), can be seen
  /// as 3 fields per dof: each node has x,y and z coordinates
  Uint m_nb_fields;

  /// A map of buffers: for each element type we have one buffer
  common::DataMap<mesh::PointSetTagExt, MetricData> m_metric_data_map;

  /// Linear index of metric data, consists of a vector of pairs
  /// [PointSetTag,Uint]
  /// Value m_mdata_index[i] is a pair [FIRST,SECOND] and corresponds to
  /// metric data which can be found in
  /// m_metric_data_map.std_region_data(FIRST), cell block on position
  /// [SECOND] in this particular MetricData
  std::vector<std::pair<common::PtrHandle<MetricData>, Uint>> m_mdata_index;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
SolutionSpaceMetric<MeshConfig, DIM>::SolutionSpaceMetric()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
SolutionSpaceMetric<MeshConfig, DIM>::~SolutionSpaceMetric()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
Uint SolutionSpaceMetric<MeshConfig, DIM>::max_nb_blocks_in_buffer() const
{
  return m_nb_blocks;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
Uint SolutionSpaceMetric<MeshConfig, DIM>::nb_fields() const
{
  return m_nb_fields;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
template <typename FeValsIterator>
void SolutionSpaceMetric<MeshConfig, DIM>::allocate_buffer(const FeValsIterator fe_vals_begin,
                                                           const FeValsIterator fe_vals_end,
                                                           Uint const nb_blocks,
                                                           Uint const nb_fields)
{
  /// NOTE: for discontinuous method, there is no need to allocate m_buffer
  /// and copy data into it, but instead m_buffer could be just a
  /// MatrixBlock<Real>. To make this class work with both continuous and
  /// discontinuous fields, each metric_data should have an array of 'raw'
  /// values and m_buffer should become a MatrixBlock<Real>. Then the
  /// procedure would be as follows: I) For continous method, 'raw' values
  /// would allocate sufficient memory
  ///     and let m_buffer be a proxy to this block
  /// II) For discontinous method, no 'raw' values would be allocated, but
  /// only
  ///     m_buffer would serve as proxy to wherever the data is stored. Since
  ///     the data is already discontinous, it does not need to be copied, but
  ///     the work can be directly done on it by means of the proxy.

  m_metric_data_map.clear();

  m_nb_blocks = nb_blocks;
  m_nb_fields = nb_fields;

  FeValsIterator it = fe_vals_begin;
  Uint n_elem_types = 0;

  for (; it != fe_vals_end; ++it)
  {
    n_elem_types++;
  }

  m_mdata_index.reserve(n_elem_types * nb_blocks);
  m_mdata_index.resize(0);

  it = fe_vals_begin;

  for (; it != fe_vals_end; ++it)
  {
    common::PtrHandle<MetricData> met_data = m_metric_data_map.create(it.key_value());

    (*met_data).std_region_type = it.key_value();

    FEValues::copy((*it.data_ptr()), (*met_data).fe_values);

    (*met_data).nb_elem_filled = 0;

    const Uint nb_qd_pts_per_elem = (*met_data).fe_values.nb_qd_pts();

    (*met_data).m_values.resize(nb_qd_pts_per_elem, nb_blocks * nb_fields);
    (*met_data).m_values.fill(0.0);

    for (Uint d = 0; d < DIM; ++d)
    {
      (*met_data).m_derivatives[d].resize(nb_qd_pts_per_elem, nb_blocks * nb_fields);
      (*met_data).m_derivatives[d].fill(0.0);
    }

    //(*met_data).invJ.resize((*met_data).nb_qd_pts_per_elem * m_nb_blocks);
    //(*met_data).det_j.resize((*met_data).nb_qd_pts_per_elem *
    // m_nb_blocks);

    /*
    if ( evaluate_boundary_normals )
    {
      (*met_data).m_normals.resize((*met_data).nb_qd_pts_per_elem,nb_blocks*m_nb_fields);
      // Note that here MeshConfig::GDIM is actually the number of fields!
    We should also be able
      // to write:
      //(*met_data).m_normals.resize((*met_data).nb_qd_pts_per_elem,nb_blocks*MeshConfig::GDIM);
    }
    */
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void SolutionSpaceMetric<MeshConfig, DIM>::empty_buffer()
{
  for (typename common::DataMap<mesh::PointSetTagExt, MetricData>::iterator it =
           m_metric_data_map.begin();
       it != m_metric_data_map.end(); ++it)
  {
    (*it.data_ptr()).nb_elem_filled = 0;
  }
  m_mdata_index.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void SolutionSpaceMetric<MeshConfig, DIM>::clear()
{
  m_nb_blocks = 0;

  m_nb_fields = 0;

  m_metric_data_map.clear();

  m_mdata_index.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
Uint SolutionSpaceMetric<MeshConfig, DIM>::nb_values_in_buffer() const
{
  return m_mdata_index.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void SolutionSpaceMetric<MeshConfig, DIM>::evaluate(
    GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &gmetric, SolutionCache const &sc,
    const ComputeMetricDerivs compute_derivatives, const RebuildMetricIndex rebuild_index)
{
  // 1) Rebuild index if necessary
  if (rebuild_index)
  {
    // Rebuild the data index for fast access to metric data
    m_mdata_index.resize(0);

    for (Uint i = 0; i < sc.nb_values_in_buffer(); ++i)
    {
      common::PtrHandle<MetricData> md = m_metric_data_map.std_region_data(sc.std_region_type(i));
      (*md).nb_elem_filled++;
      m_mdata_index.push_back(
          std::pair<common::PtrHandle<MetricData>, Uint>(md, sc.position_in_block(i)));
    }
  }

  // 2) Compute field values in quadrature points
  for (typename common::DataMap<mesh::PointSetTagExt, MetricData>::iterator it =
           m_metric_data_map.begin();
       it != m_metric_data_map.end(); ++it)
  {
    if ((*it.data_ptr()).nb_elem_filled != 0)
    {
      math::DenseConstMatView<Real> const buffer = sc.buffer_data((*it.data_ptr()).std_region_type);

      math::DenseDMat<Real> &values  = (*it.data_ptr()).m_values;
      math::DenseDMat<Real> const &V = (*it.data_ptr()).fe_values.Vandermonde();

      values = V * buffer;
    }
  }

  // 3) Compute field derivatives in quadrature points
  if (compute_derivatives)
  {
    MetricComputer<MeshConfig::GDIM, DIM, 0>::compute_derivatives(gmetric, sc, *this);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
typename SolutionSpaceMetric<MeshConfig, DIM>::cellwise_metric const SolutionSpaceMetric<
    MeshConfig, DIM>::cellwise_values(const Uint idx) const
{
  cellwise_metric cm(m_mdata_index[idx].first, m_mdata_index[idx].second, m_nb_fields);
  return cm;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void SolutionSpaceMetric<MeshConfig, DIM>::print_fill_status() const
{
  std::cout << "=== Solution space metric status ===" << std::endl;
  std::cout << "\tCan store " << m_metric_data_map.size() << " element types:" << std::endl;
  for (typename common::DataMap<mesh::PointSetTagExt, MetricData>::const_iterator it =
           m_metric_data_map.cbegin();
       it != m_metric_data_map.cend(); ++it)
  {
    std::cout << "\t\tElem type " << (*it.data_ptr()).std_region_type.std_region_tag().as_string()
              << std::endl;
    std::cout << "\t\tNb of elements filled: " << (*it.data_ptr()).nb_elem_filled << "/"
              << m_nb_blocks << std::endl;
    std::cout << "\tFE values:" << std::endl;
    (*it.data_ptr()).fe_values.print();
  }
  std::cout << "\tData index length: " << m_mdata_index.size() << std::endl;
}

// ----------------------------------------------------------------------------
// Specializations of FunctionSpaceMetric::MetricComputer for some combinations
// of dimensions  - here 2D volume metric (GDIM = 2, TDIM = 2)
// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
template <Uint Dummy>
class SolutionSpaceMetric<MeshConfig, DIM>::MetricComputer<_2D, _2D, Dummy>
{
  public:
  static void compute_derivatives(
      GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &support,
      SolutionCache const &sc, SolutionSpaceMetric<MeshConfig, DIM> &solution)
  {
    math::DenseSVec<Real, _2D> dPhys,
        dRef; // Derivatives of one shape function
              // in physical and reference space

    using geo_cell_metric =
        typename GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM>::cellwise_metric;
    using mdata_type = typename SolutionSpaceMetric<MeshConfig, DIM>::MetricData;

    // Step 1: compute values and derivatives in interpolation points,
    // values of derivatives only in reference space
    for (typename common::DataMap<mesh::PointSetTagExt, mdata_type>::iterator it =
             solution.m_metric_data_map.begin();
         it != solution.m_metric_data_map.end(); ++it)
    {
      mdata_type &data = (*it.data_ptr());

      // If no data have been filled, there's nothing to compute
      if (data.nb_elem_filled > 0)
      {
        // math::DenseDMat<Real> &derivatives_XI0 =
        // data.m_derivatives[XI0]; math::DenseDMat<Real>
        // &derivatives_XI1 = data.m_derivatives[XI1];

        math::DenseMatView<Real> derivatives_XI0 = data.m_derivatives[XI0].block(
            0, 0, data.fe_values.nb_qd_pts(), data.nb_elem_filled * solution.m_nb_fields);
        math::DenseMatView<Real> derivatives_XI1 = data.m_derivatives[XI1].block(
            0, 0, data.fe_values.nb_qd_pts(), data.nb_elem_filled * solution.m_nb_fields);

        // First use the 'derivatives' matrices to compute derivatives
        // with respect to reference space variables (xi,eta,zeta). This
        // is because dV_XI0 and dV_XI1 are derivative matrices in
        // reference space
        math::DenseConstMatView<Real> const buffer = sc.buffer_data(it.key_value());

        derivatives_XI0 = data.fe_values.deriv_Vandermonde(XI0) * buffer;
        derivatives_XI1 = data.fe_values.deriv_Vandermonde(XI1) * buffer;
      }
    }

    // Step 2: compute derivatives in physical space by reusing the
    // jacobians from the 'support' metric

    for (Uint c = 0; c < solution.m_mdata_index.size(); ++c)
    {
      // geo_cell_metric const geo_metric =
      // support.cellwise_values(solution.m_mdata_index[c].second);
      geo_cell_metric const geo_metric = support.cellwise_values(c);
      // std::cout << "Nb. nodes in geo = " << geo_metric.nb_dof_in_cell()
      // << std::endl;

      // const math::VectorBlock<Real> jacobians_at_qd_pt =
      // geo_metric.jdet();

      mdata_type &data_solution = (*solution.m_mdata_index[c].first);

      // math::DenseDMat<Real> &derivatives_XI0 =
      // data_solution.m_derivatives[XI0]; math::DenseDMat<Real>
      // &derivatives_XI1 = data_solution.m_derivatives[XI1];

      // const Uint offset_sol = solution.m_nb_fields * c;
      const Uint offset_sol = solution.m_nb_fields * solution.m_mdata_index[c].second;

      //      math::MatrixBlock<Real> derivatives_XI0 =
      //      data_solution.m_derivatives[XI0].block(
      //          0, 0, data_solution.fe_values.nb_qd_pts(),
      //          data_solution.nb_elem_filled *
      // solution.m_nb_fields);
      //      math::MatrixBlock<Real> derivatives_XI1 =
      //      data_solution.m_derivatives[XI1].block(
      //          0, 0, data_solution.fe_values.nb_qd_pts(),
      //          data_solution.nb_elem_filled *
      // solution.m_nb_fields);

      math::DenseMatView<Real> derivatives_XI0 = data_solution.m_derivatives[XI0].block(
          0, offset_sol, data_solution.fe_values.nb_qd_pts(), solution.m_nb_fields);
      math::DenseMatView<Real> derivatives_XI1 = data_solution.m_derivatives[XI1].block(
          0, offset_sol, data_solution.fe_values.nb_qd_pts(), solution.m_nb_fields);

      const Uint nb_qd_pts_per_elem = data_solution.fe_values.nb_qd_pts();
      for (Uint q = 0; q < nb_qd_pts_per_elem; ++q)
      {
        const math::DenseConstMatView<Real> invJ = geo_metric.inv_jacobi(q);

        // data_solution.det_j[data_solution.nb_qd_pts_per_elem * c + q]
        // = jacobians_at_qd_pt[q];

        for (Uint f = 0; f < solution.m_nb_fields; ++f)
        {
          /*
          dRef[XI0] = derivatives_XI0(q, offset_sol + f);
          dRef[XI1] = derivatives_XI1(q, offset_sol + f);
          */

          dRef[XI0] = derivatives_XI0(q, f);
          dRef[XI1] = derivatives_XI1(q, f);

          // dPhys = invJ * dRef;
          dPhys[X0] = invJ(0, 0) * dRef[XI0] + invJ(0, 1) * dRef[XI1];
          dPhys[X1] = invJ(1, 0) * dRef[XI0] + invJ(1, 1) * dRef[XI1];

          /*
          derivatives_XI0(q, offset_sol + f) = dPhys[X0];
          derivatives_XI1(q, offset_sol + f) = dPhys[X1];
          */

          derivatives_XI0(q, f) = dPhys[X0];
          derivatives_XI1(q, f) = dPhys[X1];
        }
      } // Loop over quadrature points

    } // Loop over all cells in buffer

  } // Static method 'evaluate_derivatives and jacobians'
};

// ----------------------------------------------------------------------------
// Specializations of FunctionSpaceMetric::MetricComputer for some combinations
// of dimensions  - here 3D volume metric (GDIM = 3, TDIM = 3)
// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
template <Uint Dummy>
class SolutionSpaceMetric<MeshConfig, DIM>::MetricComputer<_3D, _3D, Dummy>
{
  public:
  static void compute_derivatives(
      GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &support,
      SolutionCache const &sc, SolutionSpaceMetric<MeshConfig, DIM> &solution)
  {
    math::DenseSVec<Real, _3D> dPhys,
        dRef; // Derivatives of one shape function
              // in physical and reference space

    using geo_cell_metric =
        typename GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM>::cellwise_metric;
    using mdata_type = typename SolutionSpaceMetric<MeshConfig, DIM>::MetricData;

    // Step 1: compute values and derivatives in interpolation points,
    // values of derivatives only in reference space
    for (typename common::DataMap<mesh::PointSetTagExt, mdata_type>::iterator it =
             solution.m_metric_data_map.begin();
         it != solution.m_metric_data_map.end(); ++it)
    {
      mdata_type &data = (*it.data_ptr());

      // If no data have been filled, there's nothing to compute
      if (data.nb_elem_filled > 0)
      {
        math::DenseMatView<Real> derivatives_XI0 = data.m_derivatives[XI0].block(
            0, 0, data.fe_values.nb_qd_pts(), data.nb_elem_filled * solution.m_nb_fields);
        math::DenseMatView<Real> derivatives_XI1 = data.m_derivatives[XI1].block(
            0, 0, data.fe_values.nb_qd_pts(), data.nb_elem_filled * solution.m_nb_fields);
        math::DenseMatView<Real> derivatives_XI2 = data.m_derivatives[XI2].block(
            0, 0, data.fe_values.nb_qd_pts(), data.nb_elem_filled * solution.m_nb_fields);

        // First use the 'derivatives' matrices to compute derivatives
        // with respect to reference space variables (xi,eta,zeta). This
        // is because dV_XI0 and dV_XI1 are derivative matrices in
        // reference space
        math::DenseConstMatView<Real> const buffer = sc.buffer_data(it.key_value());

        derivatives_XI0 = data.fe_values.deriv_Vandermonde(XI0) * buffer;
        derivatives_XI1 = data.fe_values.deriv_Vandermonde(XI1) * buffer;
        derivatives_XI2 = data.fe_values.deriv_Vandermonde(XI2) * buffer;
      }
    }

    // Step 2: compute derivatives in physical space by reusing the
    // jacobians from the 'support' metric
    for (Uint c = 0; c < solution.m_mdata_index.size(); ++c)
    {
      // geo_cell_metric const geo_metric =
      // support.cellwise_values(solution.m_mdata_index[c].second);
      geo_cell_metric const geo_metric = support.cellwise_values(c);

      // const math::VectorBlock<Real> jacobians_at_qd_pt =
      // geo_metric.jdet();

      mdata_type &data_solution = (*solution.m_mdata_index[c].first);

      // math::DenseDMat<Real> &derivatives_XI0 =
      // data_solution.m_derivatives[XI0]; math::DenseDMat<Real>
      // &derivatives_XI1 = data_solution.m_derivatives[XI1];
      // math::DenseDMat<Real> &derivatives_XI2 =
      // data_solution.m_derivatives[XI2];

      // const Uint offset_sol = solution.m_nb_fields * c;
      const Uint offset_sol = solution.m_nb_fields * solution.m_mdata_index[c].second;

      math::DenseMatView<Real> derivatives_XI0 = data_solution.m_derivatives[XI0].block(
          0, offset_sol, data_solution.fe_values.nb_qd_pts(), solution.m_nb_fields);
      math::DenseMatView<Real> derivatives_XI1 = data_solution.m_derivatives[XI1].block(
          0, offset_sol, data_solution.fe_values.nb_qd_pts(), solution.m_nb_fields);
      math::DenseMatView<Real> derivatives_XI2 = data_solution.m_derivatives[XI2].block(
          0, offset_sol, data_solution.fe_values.nb_qd_pts(), solution.m_nb_fields);

      const Uint nb_qd_pts_per_elem = data_solution.fe_values.nb_qd_pts();
      for (Uint q = 0; q < nb_qd_pts_per_elem; ++q)
      {
        const math::DenseConstMatView<Real> invJ = geo_metric.inv_jacobi(q);

        // data_solution.det_j[data_solution.nb_qd_pts_per_elem * c + q]
        // = jacobians_at_qd_pt[q];

        for (Uint f = 0; f < solution.m_nb_fields; ++f)
        {
          dRef[XI0] = derivatives_XI0(q, f);
          dRef[XI1] = derivatives_XI1(q, f);
          dRef[XI2] = derivatives_XI2(q, f);

          // dPhys = invJ * dRef;
          dPhys[X0] = invJ(0, 0) * dRef[XI0] + invJ(0, 1) * dRef[XI1] + invJ(0, 2) * dRef[XI2];
          dPhys[X1] = invJ(1, 0) * dRef[XI0] + invJ(1, 1) * dRef[XI1] + invJ(1, 2) * dRef[XI2];
          dPhys[X2] = invJ(2, 0) * dRef[XI0] + invJ(2, 1) * dRef[XI1] + invJ(2, 2) * dRef[XI2];

          derivatives_XI0(q, f) = dPhys[X0];
          derivatives_XI1(q, f) = dPhys[X1];
          derivatives_XI2(q, f) = dPhys[X2];
        }
      } // Loop over quadrature points

    } // Loop over all cells in buffer

  } // Static method 'evaluate_derivatives and jacobians'
};

// ----------------------------------------------------------------------------
// Specializations of FunctionSpaceMetric::MetricComputer for some combinations
// of dimensions  - here 1D line metric in 2D space (GDIM = 2, DIM = 1)
// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
template <Uint Dummy>
class SolutionSpaceMetric<MeshConfig, DIM>::MetricComputer<_2D, _1D, Dummy>
{
  public:
  static void compute_derivatives(
      GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &support,
      SolutionCache const &sc, SolutionSpaceMetric<MeshConfig, DIM> &solution)
  {
    // math::StaticVector<Real,_2D> normal;

    using mdata_type = typename SolutionSpaceMetric<MeshConfig, DIM>::MetricData;

    // Step 1: compute values and derivatives in interpolation points,
    // values of derivatives only
    //         in reference space
    for (typename common::DataMap<mesh::PointSetTagExt, mdata_type>::iterator it =
             solution.m_metric_data_map.begin();
         it != solution.m_metric_data_map.end(); ++it)
    {
      mdata_type &data = (*it.data_ptr());

      // If no data have been filled, there's nothing to compute
      if (data.nb_elem_filled == 0)
        return;

      // First use the 'derivatives' matrices to compute derivatives with
      // respect to reference space variables (xi,eta,zeta). This is
      // because dV_XI0 and dV_XI1 are derivative matrices in reference
      // space
      math::DenseConstMatView<Real> const buffer = sc.buffer_data(it.key_value());

      data.m_derivatives[XI0] = data.fe_values.deriv_Vandermonde(XI0) * buffer;
    }

    /*
    // Step 2: compute derivatives in physical space by reusing the
    jacobians from the 'support' metric for(Uint c = 0; c <
    support.m_mdata_index.size(); ++c)
    {
      mdata_type const& data_support = (*support.m_mdata_index[c].first);
      mdata_type & data_solution = (*solution.m_mdata_index[c].first);

      const Uint nb_qd_pt = data_support.nb_qd_pts_per_elem;
      const Uint cell_idx_support =  support.m_mdata_index[c].second;
      const Uint cell_idx_solution = solution.m_mdata_index[c].second;

      for(Uint q = 0; q < data_support.nb_qd_pts_per_elem; ++q)
      {
        data_solution.det_j[cell_idx_solution*nb_qd_pt+q] =
    data_support.det_j[cell_idx_support*nb_qd_pt + q];
        data_solution.m_normals(q,cell_idx_solution*solution.m_nb_fields+X0)
    = data_support.m_normals(q,cell_idx_support*support.m_nb_fields+X0);
        data_solution.m_normals(q,cell_idx_solution*solution.m_nb_fields+X1)
    = data_support.m_normals(q,cell_idx_support*support.m_nb_fields+X1); }
    // Loop over quadrature points

    } // Loop over all cells in buffer
    */

  } // Static method 'evaluate_derivatives and jacobians'
};

// ----------------------------------------------------------------------------
// Specializations of FunctionSpaceMetric::MetricComputer for some combinations
// of dimensions  - here 2D surface metric in 3D space (GDIM = 3, DIM = 2)
// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
template <Uint Dummy>
class SolutionSpaceMetric<MeshConfig, DIM>::MetricComputer<_3D, _2D, Dummy>
{
  public:
  static void compute_derivatives(
      GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &support,
      SolutionCache const &sc, SolutionSpaceMetric<MeshConfig, DIM> &solution)
  {
    // math::StaticVector<Real, _3D> normal;

    using mdata_type = typename SolutionSpaceMetric<MeshConfig, DIM>::MetricData;

    // Step 1: compute values and derivatives in interpolation points,
    // values of derivatives only
    //         in reference space
    for (typename common::DataMap<mesh::PointSetTagExt, mdata_type>::iterator it =
             solution.m_metric_data_map.begin();
         it != solution.m_metric_data_map.end(); ++it)
    {
      mdata_type &data = (*it.data_ptr());

      // If no data have been filled, there's nothing to compute
      if (data.nb_elem_filled == 0)
        return;

      // First use the 'derivatives' matrices to compute derivatives with
      // respect
      // to reference space variables (xi,eta,zeta). This is because
      // dV_XI0 and dV_XI1 are derivative matrices in reference space
      math::DenseConstMatView<Real> const buffer = sc.buffer_data(it.key_value());

      data.m_derivatives[XI0] = data.fe_values.deriv_Vandermonde(XI0) * buffer;
      data.m_derivatives[XI1] = data.fe_values.deriv_Vandermonde(XI1) * buffer;
    }

    /*
    // Step 2: compute derivatives in physical space by reusing the
    jacobians from the 'support' metric for(Uint c = 0; c <
    support.m_mdata_index.size(); ++c)
    {
      mdata_type const& data_support = (*support.m_mdata_index[c].first);
      mdata_type & data_solution = (*solution.m_mdata_index[c].first);

      const Uint nb_qd_pt = data_support.nb_qd_pts_per_elem;
      const Uint cell_idx_support =  support.m_mdata_index[c].second;
      const Uint cell_idx_solution = solution.m_mdata_index[c].second;

      //const Uint offset_sol = solution.m_nb_fields * c;

      for(Uint q = 0; q < data_support.nb_qd_pts_per_elem; ++q)
      {
        data_solution.det_j[cell_idx_solution*nb_qd_pt+q] =
    data_support.det_j[cell_idx_support*nb_qd_pt + q];
        data_solution.m_normals(q,cell_idx_solution*solution.m_nb_fields+X0)
    = data_support.m_normals(q,cell_idx_support*support.m_nb_fields+X0);
        data_solution.m_normals(q,cell_idx_solution*solution.m_nb_fields+X1)
    = data_support.m_normals(q,cell_idx_support*support.m_nb_fields+X1);
        data_solution.m_normals(q,cell_idx_solution*solution.m_nb_fields+X2)
    = data_support.m_normals(q,cell_idx_support*support.m_nb_fields+X2);

      } // Loop over quadrature points

    } // Loop over all cells in buffer
    */

  } // Static method 'evaluate_derivatives and jacobians'
};

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
