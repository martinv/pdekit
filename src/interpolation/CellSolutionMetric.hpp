#ifndef PDEKIT_Interpolation_Cell_Solution_Metric_hpp
#define PDEKIT_Interpolation_Cell_Solution_Metric_hpp

#include "common/PtrHandle.hpp"
#include "interpolation/FEValues.hpp"
#include "math/DenseConstMatView.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"
#include "mesh/std_region/PointSetTagExt.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename MetricData>
class CellSolutionMetric
{

  public:
  /// Default constructor
  CellSolutionMetric();

  /// Constructor
  CellSolutionMetric(common::PtrHandle<MetricData> const data, const Uint idx,
                     const Uint nb_fields);

  /// Copy constructor
  CellSolutionMetric(const CellSolutionMetric &other_met);

  /// Assignment operator
  CellSolutionMetric &operator=(CellSolutionMetric const &met_rhs);

  /// Get the standard region type of the cell to which this cell metric
  /// corresponds
  const mesh::PointSetTagExt std_region_type() const;

  /// Get the number of dof
  Uint nb_dof_in_cell() const;

  /// Get the number of fields
  Uint nb_fields() const;

  /// Get the number of quadrature points
  Uint nb_qd_pts() const;

  /// Values of shape functions in reference space
  math::DenseDMat<Real> const &reference_sf_values() const;

  /// Get the matrix of interpolated values
  const math::DenseConstMatView<Real> field_values() const;

  /// Get the matrix of derivatives at interpolation points
  const math::DenseConstMatView<Real> field_derivatives(const Uint dim) const;

  protected:
  /// Reference to given metric data
  // MetricData const &m_mdata;
  common::PtrHandle<MetricData const> m_mdata;

  /// Cell index in that metric data
  Uint m_c_idx;

  /// Number of fields
  Uint m_nb_fields;
};

// ----------------------------------------------------------------------------

template <typename MetricData>
CellSolutionMetric<MetricData>::CellSolutionMetric()
    : m_mdata(nullptr), m_c_idx(0u), m_nb_fields(0u)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
CellSolutionMetric<MetricData>::CellSolutionMetric(common::PtrHandle<MetricData> const data,
                                                   const Uint idx, const Uint nb_fields)
    : m_mdata(data), m_c_idx(idx), m_nb_fields(nb_fields)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
CellSolutionMetric<MetricData>::CellSolutionMetric(const CellSolutionMetric &other_met)
    : m_mdata(other_met.m_mdata), m_c_idx(other_met.m_c_idx), m_nb_fields(other_met.m_nb_fields)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
CellSolutionMetric<MetricData> &CellSolutionMetric<MetricData>::operator=(
    CellSolutionMetric const &met_rhs)
{
  m_mdata     = met_rhs.m_mdata;
  m_c_idx     = met_rhs.m_c_idx;
  m_nb_fields = met_rhs.m_nb_fields;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline const mesh::PointSetTagExt CellSolutionMetric<MetricData>::std_region_type() const
{
  return (*m_mdata).std_region_type;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline Uint CellSolutionMetric<MetricData>::nb_dof_in_cell() const
{
  return (*m_mdata).fe_values.nb_nodes();
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline Uint CellSolutionMetric<MetricData>::nb_fields() const
{
  return m_nb_fields;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline Uint CellSolutionMetric<MetricData>::nb_qd_pts() const
{
  return (*m_mdata).fe_values.nb_qd_pts();
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline math::DenseDMat<Real> const &CellSolutionMetric<MetricData>::reference_sf_values() const
{
  return (*m_mdata).fe_values.Vandermonde();
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline const math::DenseConstMatView<Real> CellSolutionMetric<MetricData>::field_values() const
{
  const Uint nb_qd_pts = (*m_mdata).fe_values.nb_qd_pts();
  math::DenseConstMatView<Real> field_values(
      (*m_mdata).m_values.const_block(0u, m_c_idx * m_nb_fields, nb_qd_pts, m_nb_fields));
  return field_values;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline const math::DenseConstMatView<Real> CellSolutionMetric<MetricData>::field_derivatives(
    const Uint dim) const
{
  const Uint nb_qd_pts = (*m_mdata).fe_values.nb_qd_pts();
  math::DenseConstMatView<Real> field_derivatives(
      (*m_mdata).m_derivatives[dim].const_block(0u, m_c_idx * m_nb_fields, nb_qd_pts, m_nb_fields));

  return field_derivatives;
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
