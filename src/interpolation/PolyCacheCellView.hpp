#ifndef PDEKIT_Interpolation_Poly_Cache_Cell_View_hpp
#define PDEKIT_Interpolation_Poly_Cache_Cell_View_hpp

#include "common/PtrHandle.hpp"
#include "interpolation/FEValues.hpp"
#include "math/DenseDMat.hpp"
#include "mesh/DiscreteElemKey.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename MetricData>
class PolyCacheCellView
{

  public:
  /// Default constructor
  PolyCacheCellView();

  /// Constructor
  PolyCacheCellView(common::PtrHandle<MetricData const> const data);

  /// Copy constructor
  PolyCacheCellView(const PolyCacheCellView &other_met);

  /// Assignment operator
  PolyCacheCellView &operator=(PolyCacheCellView const &met_rhs);

  /// Get the standard region type of the cell to which this cell metric
  /// corresponds
  const mesh::DiscreteElemKey discrete_elem_key() const;

  /// Get the number of dof
  Uint nb_dof_in_cell() const;

  /// Get the number of quadrature points
  Uint nb_qd_pts() const;

  /// Values of shape functions in reference space
  const math::DenseDMat<Real> &reference_sf_values() const;

  /// Get the matrix of derivatives at interpolation points
  const math::DenseDMat<Real> &reference_sf_derivatives(const Uint dim) const;

  protected:
  /// Reference to given metric data
  common::PtrHandle<MetricData const> m_mdata;
};

// ----------------------------------------------------------------------------

template <typename MetricData>
PolyCacheCellView<MetricData>::PolyCacheCellView() : m_mdata(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
PolyCacheCellView<MetricData>::PolyCacheCellView(common::PtrHandle<MetricData const> const data)
    : m_mdata(data)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
PolyCacheCellView<MetricData>::PolyCacheCellView(const PolyCacheCellView &other_met)
    : m_mdata(other_met.m_mdata)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData>
PolyCacheCellView<MetricData> &PolyCacheCellView<MetricData>::operator=(
    PolyCacheCellView const &met_rhs)
{
  m_mdata = met_rhs.m_mdata;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline const mesh::DiscreteElemKey PolyCacheCellView<MetricData>::discrete_elem_key() const
{
  return (*m_mdata).pdata_type;
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline Uint PolyCacheCellView<MetricData>::nb_dof_in_cell() const
{
  return (*m_mdata).ref_fe_values.nb_nodes();
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline Uint PolyCacheCellView<MetricData>::nb_qd_pts() const
{
  return (*m_mdata).ref_fe_values.nb_qd_pts();
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline math::DenseDMat<Real> const &PolyCacheCellView<MetricData>::reference_sf_values() const
{
  return (*m_mdata).ref_fe_values.Vandermonde();
}

// ----------------------------------------------------------------------------

template <typename MetricData>
inline const math::DenseDMat<Real> &PolyCacheCellView<MetricData>::reference_sf_derivatives(
    const Uint dim) const
{
  return (*m_mdata).ref_fe_values.deriv_Vandermonde(dim);
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
