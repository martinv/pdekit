#ifndef PDEKIT_Solver_FE_Metric_Data_hpp
#define PDEKIT_Solver_FE_Metric_Data_hpp

#include <array>

#include "interpolation/FEValues.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"

namespace pdekit
{

namespace solver
{

/// ============================================================================
///       CONTAINER FOR THE DATA IN REFERENCE SPACE
///       CONCERNING THE GEOMETRY OF THE ELEMENT
/// ============================================================================

/// GEODIM - geometrical dimension
/// TOPODIM - topological dimension
/// Example: triangle on 3D surface has geometrical dimension 3 (it's in 3D
/// space), but topological dimension
/// only 2 (it is a surface element)
/// A line in 2D has geometrical dimension 2 and topological dimension 1
/// A line in 3D has geometrical dimension 3 and topological dimension 1
/// Topological dimension for a given shape is always the same, no matter what
/// the geometrical dimension is

template <Uint GEODIM, Uint TDIM>
class FEMetricData
{
  public:
  /// METHODS

  /// Constructor
  FEMetricData();

  /// Destructor
  ~FEMetricData();

  /// This method resizes all member variables
  /// Depending on what's the finite element passed to it
  void resize_variables(interpolation::FEValues const &fe);

  /// finite element for geometry (support)
  interpolation::FEValues m_fe;

  /// Vandermonde matrices
  math::DenseDMat<Real> V;
  std::array<math::DenseDMat<Real>, TDIM> dV;

  /// Number of nodes in element
  Uint nb_nodes;

  /// Number of quadrature points for this element
  Uint nb_qd_pts;

  /// Nodes of one cell
  mesh::MeshEntity m_elem_nodes;

  /// Nodal coordinates of one cell
  math::DenseDMat<Real> m_node_coord;

  /// Quadrature coordinates of one cell
  math::DenseDMat<Real> m_quad_coord;

  /// Coordinates of quadrature points in physical space
  math::DenseDMat<Real> m_xq;

  /// Derivatives of the transformation from physical to
  /// reference space with respect to xi at all quadrature points
  /// Has 2(3) columns: each row is composed of (dx/dxi dy/dxi (dz/dxi)) for
  /// given
  /// quadrature point
  std::array<math::DenseDMat<Real>, TDIM> m_dxqd;

  /// Transformation Jacobian in each quadrature point
  math::DenseDVec<Real> m_j;
};

// ===========================================================================

template <Uint GEODIM, Uint TDIM>
FEMetricData<GEODIM, TDIM>::FEMetricData()
{
}

// ===========================================================================

template <Uint GEODIM, Uint TDIM>
FEMetricData<GEODIM, TDIM>::~FEMetricData()
{
}

// ============================================================================

template <Uint GEODIM, Uint TDIM>
void FEMetricData<GEODIM, TDIM>::resize_variables(interpolation::FEValues const &fe)
{
  interpolation::FEValues::copy(fe, m_fe);

  V.resize(fe.Vandermonde().rows(), fe.Vandermonde().cols());
  for (Uint r = 0; r < V.rows(); ++r)
  {
    for (Uint c = 0; c < V.cols(); ++c)
    {
      V(r, c) = fe.Vandermonde()(r, c);
    }
  }

  for (Uint dim = 0; dim < TDIM; ++dim)
  {
    dV[dim].resize(fe.deriv_Vandermonde(dim).rows(), fe.deriv_Vandermonde(dim).cols());
    for (Uint r = 0; r < dV[dim].rows(); ++r)
    {
      for (Uint c = 0; c < dV[dim].cols(); ++c)
      {
        dV[dim](r, c) = fe.deriv_Vandermonde(dim)(r, c);
      }
    }
  }

  nb_nodes  = m_fe.nb_nodes();
  nb_qd_pts = m_fe.qp().rows();

  m_node_coord.resize(nb_nodes, GEODIM);

  m_quad_coord.resize(m_fe.qp().rows(), m_fe.qp().cols());
  for (Uint i = 0; i < m_fe.qp().rows(); ++i)
  {
    for (Uint j = 0; j < m_fe.qp().cols(); ++j)
    {
      m_quad_coord(i, j) = m_fe.qp()(i, j);
    }
  }

  m_xq.resize(nb_qd_pts, GEODIM);

  for (Uint d = 0; d < TDIM; ++d)
  {
    m_dxqd[d].resize(nb_qd_pts, GEODIM);
  }
  m_j.resize(nb_qd_pts);
}

// ============================================================================

} // namespace solver

} // namespace pdekit

#endif
