#ifndef PDEKIT_Interpolation_CellGeoMetric_hpp
#define PDEKIT_Interpolation_CellGeoMetric_hpp

#include "common/PtrHandle.hpp"
#include "interpolation/FEValues.hpp"
#include "math/DenseConstMatView.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename MetricData>
class CellGeoMetric
{

  public:
  /// Default constructor
  CellGeoMetric();

  /// Constructor
  CellGeoMetric(common::PtrHandle<MetricData> const data, const Uint idx);

  /// Copy constructor
  CellGeoMetric(const CellGeoMetric &other_met);

  /// Assignment operator
  CellGeoMetric &operator=(const CellGeoMetric &met_rhs);

  /// Get the number of dof
  Uint nb_dof_in_cell() const;

  /// Get the number of fields
  Uint nb_fields() const;

  /// Get the number of quadrature points
  Uint nb_qd_pts() const;

  /// Point weights in reference space
  math::DenseDVec<Real> const &pt_weights() const;

  /// Values of shape functions in reference space
  math::DenseDMat<Real> const &reference_sf_values() const;

  /// Derivatives of shape functions in reference space
  math::DenseDMat<Real> const &reference_sf_derivatives(const Uint dim) const;

  /// Get the matrix of interpolated values
  const math::DenseConstMatView<Real> interpolated_coords() const;

  /// Get the matrix of derivatives at interpolation points
  const math::DenseConstMatView<Real> coord_transf_derivatives(const Uint dim) const;

  /// Get the matrix of shape function derivatives at interpolation points
  const math::DenseConstMatView<Real> sf_derivatives(const Uint dim) const;

  /// Get the inverse of the jacobi matrix
  const math::DenseConstMatView<Real> inv_jacobi(const Uint qd_pt) const;

  /// Return an array of jacobi determinants
  const math::DenseConstVecView<Real> jdet() const;

  /// Returns normals at quadrature points
  math::DenseConstMatView<Real> const normals() const;

  protected:
  /// Reference to given metric data
  // MetricData const &m_mdata;
  common::PtrHandle<MetricData const> m_mdata;

  /// Cell index in that metric data
  Uint m_c_idx;
};

// ----------------------------------------------------------------------------

template <typename MetricData>
CellGeoMetric<MetricData>::CellGeoMetric() : m_mdata(nullptr), m_c_idx(0u)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
CellGeoMetric<MetricData>::CellGeoMetric(common::PtrHandle<MetricData> const data, const Uint idx)
    : m_mdata(data), m_c_idx(idx)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
CellGeoMetric<MetricData>::CellGeoMetric(const CellGeoMetric &other_met)
    : m_mdata(other_met.m_mdata), m_c_idx(other_met.m_c_idx)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
CellGeoMetric<MetricData> &CellGeoMetric<MetricData>::operator=(const CellGeoMetric &met_rhs)
{
  m_mdata = met_rhs.m_mdata;
  m_c_idx = met_rhs.m_c_idx;
  return *this;
}

// ----------------------------------------------------------------------------

/// Get the number of dof
template <typename MetricData>
inline Uint CellGeoMetric<MetricData>::nb_dof_in_cell() const
{
  return (*m_mdata).fe_values.nb_nodes();
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline Uint CellGeoMetric<MetricData>::nb_fields() const
{
  return MetricData::GDIM;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline Uint CellGeoMetric<MetricData>::nb_qd_pts() const
{
  return (*m_mdata).fe_values.nb_qd_pts();
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline math::DenseDVec<Real> const &CellGeoMetric<MetricData>::pt_weights() const
{
  return (*m_mdata).fe_values.qw();
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline math::DenseDMat<Real> const &CellGeoMetric<MetricData>::reference_sf_values() const
{
  return (*m_mdata).fe_values.Vandermonde();
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline math::DenseDMat<Real> const &CellGeoMetric<MetricData>::reference_sf_derivatives(
    const Uint dim) const
{
  return (*m_mdata).fe_values.deriv_Vandermonde(dim);
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline const math::DenseConstMatView<Real> CellGeoMetric<MetricData>::interpolated_coords() const
{
  const Uint nb_qd_pts = (*m_mdata).fe_values.nb_qd_pts();
  math::DenseConstMatView<Real> field_values((*m_mdata).m_coord_values.const_block(
      0u, m_c_idx * MetricData::GDIM, nb_qd_pts, MetricData::GDIM));
  return field_values;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline const math::DenseConstMatView<Real> CellGeoMetric<MetricData>::coord_transf_derivatives(
    const Uint dim) const
{
  const Uint nb_qd_pts = (*m_mdata).fe_values.nb_qd_pts();
  math::DenseConstMatView<Real> field_derivatives((*m_mdata).m_coord_derivatives[dim].const_block(
      0u, m_c_idx * MetricData::GDIM, nb_qd_pts, MetricData::GDIM));
  return field_derivatives;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline const math::DenseConstMatView<Real> CellGeoMetric<MetricData>::sf_derivatives(
    const Uint dim) const
{
  const Uint nb_qd_pts       = (*m_mdata).fe_values.nb_qd_pts();
  const Uint nb_dof_per_elem = (*m_mdata).fe_values.nb_nodes();

  math::DenseConstMatView<Real> sf_derivatives((*m_mdata).m_sf_derivatives[dim].const_block(
      0u, m_c_idx * nb_dof_per_elem, nb_qd_pts, nb_dof_per_elem));
  return sf_derivatives;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline const math::DenseConstMatView<Real> CellGeoMetric<MetricData>::inv_jacobi(
    const Uint qd_pt) const
{
  const Uint nb_qd_pts = (*m_mdata).fe_values.nb_qd_pts();
  math::DenseConstMatView<Real> inv_jacobi((*m_mdata).invJ[m_c_idx * nb_qd_pts + qd_pt].const_block(
      0u, 0u, MetricData::GDIM, MetricData::GDIM));
  return inv_jacobi;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline const math::DenseConstVecView<Real> CellGeoMetric<MetricData>::jdet() const
{
  const Uint nb_qd_pts = (*m_mdata).fe_values.nb_qd_pts();
  math::DenseConstVecView<Real> jdet((*m_mdata).det_j.data() + m_c_idx * nb_qd_pts, nb_qd_pts);
  return jdet;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
math::DenseConstMatView<Real> const CellGeoMetric<MetricData>::normals() const
{
  const Uint nb_qd_pts = (*m_mdata).fe_values.nb_qd_pts();
  math::DenseConstMatView<Real> face_normals((*m_mdata).m_normals.const_block(
      0u, m_c_idx * MetricData::GDIM, nb_qd_pts, MetricData::GDIM));
  return face_normals;
}

// ----------------------------------------------------------------------------
// Cellwise metric to be used on boundaries
// Inherits from the CellGeoMetric above, but in addition provides boundary
// normals
// ----------------------------------------------------------------------------

template <typename MetricData>
class CellGeoMetricWithNormals : public CellGeoMetric<MetricData>
{
  private:
  typedef CellGeoMetric<MetricData> Base;

  public:
  /// Default constructor
  CellGeoMetricWithNormals();

  /// Constructor
  CellGeoMetricWithNormals(common::PtrHandle<MetricData> const data, const Uint idx);

  /// Copy constructor
  CellGeoMetricWithNormals(const CellGeoMetricWithNormals &other_met);

  /// Assignment operator
  CellGeoMetricWithNormals &operator=(const CellGeoMetricWithNormals &met_rhs);

  /// Returns normals at quadrature points
  math::DenseConstMatView<Real> const normals() const;

  private:
};

// ----------------------------------------------------------------------------

template <typename MetricData>
CellGeoMetricWithNormals<MetricData>::CellGeoMetricWithNormals() : CellGeoMetric<MetricData>()
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
CellGeoMetricWithNormals<MetricData>::CellGeoMetricWithNormals(
    common::PtrHandle<MetricData> const data, const Uint idx)
    : CellGeoMetric<MetricData>(data, idx)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
CellGeoMetricWithNormals<MetricData>::CellGeoMetricWithNormals(
    const CellGeoMetricWithNormals &other_met)
    : CellGeoMetric<MetricData>(other_met)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
CellGeoMetricWithNormals<MetricData> &CellGeoMetricWithNormals<MetricData>::operator=(
    const CellGeoMetricWithNormals &met_rhs)
{
  Base::m_mdata = Base::met_rhs.m_mdata;
  Base::m_c_idx = Base::met_rhs.m_c_idx;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline math::DenseConstMatView<Real> const CellGeoMetricWithNormals<MetricData>::normals() const
{
  const Uint nb_qd_pts = (*Base::m_mdata).fe_values.nb_qd_pts();
  math::DenseConstMatView<Real> face_normals(
      (*Base::m_mdata)
          .m_normals.const_block(0u, Base::m_c_idx * MetricData::GDIM, nb_qd_pts,
                                 MetricData::GDIM));
  return face_normals;
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
