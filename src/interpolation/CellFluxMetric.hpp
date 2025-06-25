#ifndef PDEKIT_Interpolation_Cell_Flux_Metric_hpp
#define PDEKIT_Interpolation_Cell_Flux_Metric_hpp

#include "common/PtrHandle.hpp"
#include "interpolation/FEValues.hpp"
#include "math/DenseConstMatView.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
class CellFluxMetric
{

  public:
  /// Default constructor
  CellFluxMetric();

  /// Constructor
  CellFluxMetric(common::PtrHandle<MetricData> const data, const Uint idx);

  /// Copy constructor
  CellFluxMetric(const CellFluxMetric &other_met);

  /// Assignment operator
  CellFluxMetric &operator=(const CellFluxMetric &met_rhs);

  /// Get the number of dof
  Uint nb_dof_in_solution_cell() const;

  /// Get the number of dof in the flux element
  Uint nb_dof_in_flux_cell() const;

  /// Get the number of fields
  Uint nb_fields() const;

  /// Get the number of quadrature points
  Uint nb_qd_pts() const;

  /// Get the matrix of interpolated values
  const math::DenseConstMatView<Real> flux_values(const Uint dim) const;

  /// Get the matrix of derivatives at interpolation points
  const math::DenseConstMatView<Real> flux_derivatives(const Uint dim) const;

  protected:
  /// Reference to given metric data
  // MetricData const &m_mdata;
  common::PtrHandle<MetricData const> m_mdata;

  /// Cell index in that metric data
  Uint m_c_idx;
};

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
CellFluxMetric<MetricData, Physics>::CellFluxMetric() : m_mdata(nullptr), m_c_idx(0u)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
CellFluxMetric<MetricData, Physics>::CellFluxMetric(common::PtrHandle<MetricData> const data,
                                                    const Uint idx)
    : m_mdata(data), m_c_idx(idx)
{
}

// ----------------------------------------------------------------------------
template <typename MetricData, typename Physics>
CellFluxMetric<MetricData, Physics>::CellFluxMetric(const CellFluxMetric &other_met)
    : m_mdata(other_met.m_mdata), m_c_idx(other_met.m_c_idx)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
CellFluxMetric<MetricData, Physics> &CellFluxMetric<MetricData, Physics>::operator=(
    const CellFluxMetric &met_rhs)
{
  m_mdata = met_rhs.m_mdata;
  m_c_idx = met_rhs.m_c_idx;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
inline Uint CellFluxMetric<MetricData, Physics>::nb_dof_in_solution_cell() const
{
  return (*m_mdata).nb_dof_per_elem;
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
inline Uint CellFluxMetric<MetricData, Physics>::nb_dof_in_flux_cell() const
{
  return (*m_mdata).nb_flux_dof_per_elem;
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
inline Uint CellFluxMetric<MetricData, Physics>::nb_fields() const
{
  return Physics::NEQ;
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
inline Uint CellFluxMetric<MetricData, Physics>::nb_qd_pts() const
{
  return (*m_mdata).nb_qd_pts_per_elem;
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
inline const math::DenseConstMatView<Real> CellFluxMetric<MetricData, Physics>::flux_values(
    const Uint dim) const
{
  const Uint nb_qd_pts = (*m_mdata).nb_qd_pts_per_elem;
  math::DenseConstMatView<Real> flux_values((*m_mdata).m_flux_in_qd_pts[dim].const_block(
      0u, m_c_idx * Physics::NEQ, nb_qd_pts, Physics::NEQ));
  return flux_values;
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
inline const math::DenseConstMatView<Real> CellFluxMetric<MetricData, Physics>::flux_derivatives(
    const Uint dim) const
{
  const Uint nb_qd_pts = (*m_mdata).nb_qd_pts_per_elem;
  math::DenseConstMatView<Real> flux_derivatives((*m_mdata).m_flux_deriv_in_qd_pts[dim].const_block(
      0u, m_c_idx * Physics::NEQ, nb_qd_pts, Physics::NEQ));

  return flux_derivatives;
}

// ----------------------------------------------------------------------------
// Cellwise metric to be with vector valued fields
// Inherits from the CellwiseMetric above, but provides extra data
// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
class FacetFluxMetric : public CellFluxMetric<MetricData, Physics>
{

  public:
  /// Default constructor
  FacetFluxMetric();

  /// Constructor
  FacetFluxMetric(common::PtrHandle<MetricData> const data, const Uint idx);

  /// Copy constructor
  FacetFluxMetric(const FacetFluxMetric &other_met);

  /// Assignment operator
  FacetFluxMetric &operator=(const FacetFluxMetric &met_rhs);

  /// Get the matrix of interpolated derivatives - this method is overwritten
  /// to return an empty block, because flux derivatives are not available
  /// on facets
  inline const math::DenseConstMatView<Real> flux_derivatives(const Uint dim) const;

  private:
  typedef CellFluxMetric<MetricData, Physics> Base;
};

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
FacetFluxMetric<MetricData, Physics>::FacetFluxMetric() : CellFluxMetric<MetricData, Physics>()
{
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
FacetFluxMetric<MetricData, Physics>::FacetFluxMetric(common::PtrHandle<MetricData> const data,
                                                      const Uint idx)
    : CellFluxMetric<MetricData, Physics>(data, idx)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
FacetFluxMetric<MetricData, Physics>::FacetFluxMetric(const FacetFluxMetric &other_met)
    : CellFluxMetric<MetricData, Physics>(other_met)
{
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
FacetFluxMetric<MetricData, Physics> &FacetFluxMetric<MetricData, Physics>::operator=(
    const FacetFluxMetric &met_rhs)
{
  Base::m_mdata = Base::met_rhs.m_mdata;
  Base::m_c_idx = Base::met_rhs.m_c_idx;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename MetricData, typename Physics>
const math::DenseConstMatView<Real> FacetFluxMetric<MetricData, Physics>::flux_derivatives(
    const Uint dim) const
{
  return math::DenseConstMatView<Real>();
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
