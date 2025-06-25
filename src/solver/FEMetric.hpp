#ifndef PDEKIT_Solver_FE_Metric_hpp
#define PDEKIT_Solver_FE_Metric_hpp

#include <array>

#include "interpolation/FEValues.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"
#include "mesh/MeshEntity.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace solver
{

/// ============================================================================
///       CLASS FOR COMPUTING INTERPOLATED VALUES IN FINITE ELEMENT
///       COMPUTES SHAPE FUNCTION DERIVATIVES IN PHYSICAL SPACE AS WELL
/// ============================================================================

/// GEODIM - geometrical dimension
/// TDIM - topological dimension
/// Example: triangle on 3D surface has geometrical dimension 3 (it's in 3D
/// space), but topological dimension
/// only 2 (it is a surface element)
/// A line in 2D has geometrical dimension 2 and topological dimension 1
/// A line in 3D has geometrical dimension 3 and topological dimension 1
/// Topological dimension for a given shape is always the same, no matter what
/// the geometrical dimension is

template <Uint GEODIM, Uint TDIM>
class FEMetric
{
  public:
  /// METHODS

  /// Constructor
  FEMetric();

  /// Destructor
  ~FEMetric();

  /// Compute the [x,y,z] coordinates in quadrature points, i.e. compute
  /// [x(qd),y(qd),z(qd)]. In addition,
  /// compute the derivatives of the mapping (x,y) = (x(xi,eta),y(xi,eta),
  /// i.e. the mapping which determines the position of point (x,y) in
  /// physical space corresponding to another point (xi,eta) in reference
  /// space.
  /// @param V   ... matrix of size [(nb. of quadrature points) x (nb. of
  /// shape functions)]
  ///                that contains the values of all shape functions in all
  ///                quadrature points One row has the values of all shape
  ///                functions in one quadrature point
  /// @param dV  ... array of matrices where each matrix is of size
  ///                [(nb. of quadrature points) x (nb. of shape functions)]
  ///                and contains derivatives of all shape functions with
  ///                respect to one variable in reference space
  ///        dV[xi] has all values of d(phi)/d(xi), dV[eta] has all values
  ///        d(phi)/d(eta)
  /// @param node_coord ... matrix of nodal coordinates. The size is [(nb. of
  /// nodes) x dimension]
  ///                Each row has [x,y,z] coordinates of one node of given
  ///                element
  /// @param interpolated_cord ... resulting values (x,y,z) in all quadrature
  /// points.
  ///                              Each row has (x,y,z) in one quadrature
  ///                              point.
  /// @param deriv_interpolated_coord ... derivatives of the
  ///        mapping (xi,eta) -> (x(xi,eta),y(xi,eta)).
  ///        The result is an array of 3 matrices, each of them has as many
  ///        rows as we have quadrature points and 2 colums (in 2D) or 3
  ///        columns (in 3D). The first matrix has the following structure:
  ///        |dx/dxi_0, dy/dxi_0, dz/dxi_0|
  ///                                        |dx/dxi_1, dy/dxi_1, dz/dxi_1|
  ///                                                     . . .
  ///                                        |dx/dxi_n, dy/dxi_n, dz/dxi_n|
  ///        The second matrix is the same, but holds derivatives with respect
  ///        to eta The array 'deriv_interpolated_coord' serves as input for
  ///        the next method, called compute_derivatives_in_phys_space

  void interpolate(const math::DenseDMat<Real> &V,
                   const std::array<math::DenseDMat<Real>, TDIM> &dV,
                   const math::DenseDMat<Real> &node_coord,
                   math::DenseDMat<Real> &interpolated_coord,
                   std::array<math::DenseDMat<Real>, TDIM> &deriv_interpolated_coord) const;

  /// Given matrices of derivatives with respect to xi,eta,(zeta) in reference
  /// space,
  /// transform the derivatives to physical space and store in dphys
  /// @param jacobian ... derivatives of the transformation from physical to
  /// reference space.
  ///                   - jacobian[0] (or jacobian[xi]) is a matrix which has
  ///                   as
  ///                     many rows as there are quadrature points (or points
  ///                     in reference space where the derivatives should be
  ///                     computed) and 2  columns (in 2D) or 3 columns (in
  ///                     3D) Therefore each row of jacobian[xi] contains
  ///                     [dx/dxi,dy/dxi,dz/dxi] in one point.
  ///                   - Similarly, each row of jacobian[eta] contains
  ///                   [dx/deta,dy/deta,dz/deta]
  ///                     in one point.
  ///                     If GEODIM = 3, there is third matrix jacobian[zeta]
  ///
  /// @param dref     ... array of Vandermonde matrices of derivatives.
  ///                   - dref[xi] is a matrix of size
  ///                     [(nb. of quadrature points) x (nb. of shape
  ///                     functions)], where each row contains all derivatives
  ///                     of shape functions with respect to xi
  ///                   - dref[eta] contains all derivatives of shape
  ///                   functions
  ///                     with respect to eta
  ///                   - dref[zeta] contains all derivatives of shape
  ///                   functions
  ///                     with respect to zeta
  ///
  /// @param dphys    ... array of derivatives of shape functions in physical
  /// space.
  ///                   - dphys[X] is a matrix of size
  ///                     [nb. of quadrature points) x (nb. of shape
  ///                     functions)] that contains all derivatives of shape
  ///                     functions with respect to X. Each row has all
  ///                     derivatives of shape functions with respect to X in
  ///                     one point
  ///                   - dphys[Y] contains all derivatives of shape functions
  ///                   with respect to Y
  ///                   - analogously for dphys[Z] (if it is present, i.e. if
  ///                   GEODIM = 3)

  inline void compute_derivatives_in_phys_space(
      const std::array<math::DenseDMat<Real>, GEODIM> &jacobian,
      const std::array<math::DenseDMat<Real>, GEODIM> &dref,
      std::array<math::DenseDMat<Real>, GEODIM> &dphys);

  /// The same as above, but in addition, store the jacobians of the
  /// transformation in each quadrature point

  inline void compute_derivatives_in_phys_space(
      const std::array<math::DenseDMat<Real>, GEODIM> &jacobian,
      const std::array<math::DenseDMat<Real>, GEODIM> &dref,
      std::array<math::DenseDMat<Real>, GEODIM> &dphys, math::DenseDVec<Real> &jac_at_qd_pt);

  private:
  /// Derivatives of ONE shape function in reference and physical space
  /// computed at ONE point
  math::DenseSVec<Real, TDIM> m_dphi_ref;
  math::DenseSVec<Real, GEODIM> m_dphi_phys;

  /// Jacobian matrix of transformation physical->reference space
  /// and its inverse at one point
  math::DenseSMat<Real, TDIM, GEODIM> m_J;
  math::DenseSMat<Real, TDIM, GEODIM> m_J_inverse;
};

// ===========================================================================

template <Uint GEODIM, Uint TDIM>
FEMetric<GEODIM, TDIM>::FEMetric()
{
}

// ===========================================================================

template <Uint GEODIM, Uint TDIM>
FEMetric<GEODIM, TDIM>::~FEMetric()
{
}

// ============================================================================

template <Uint GEODIM, Uint TDIM>
void FEMetric<GEODIM, TDIM>::interpolate(
    const math::DenseDMat<Real> &V, const std::array<math::DenseDMat<Real>, TDIM> &dV,
    const math::DenseDMat<Real> &node_coord, math::DenseDMat<Real> &interpolated_coord,
    std::array<math::DenseDMat<Real>, TDIM> &deriv_interpolated_coord) const
{
  /// Coordinates of quadrature points in physical space
  /// [nb_qd_pts x GEODIM] = [nb_qd_pts x nb_nodes] * [nb_nodes * GEODIM]
  interpolated_coord = V * node_coord;

  /// Derivatives dx/dxi of the transformation physical -> reference space
  for (Uint dim = 0; dim < TDIM; ++dim)
  {
    deriv_interpolated_coord[dim] = dV[dim] * node_coord;
  }
}

// ============================================================================

template <Uint GEODIM, Uint TDIM>
void FEMetric<GEODIM, TDIM>::compute_derivatives_in_phys_space(
    const std::array<math::DenseDMat<Real>, GEODIM> &jacobian,
    const std::array<math::DenseDMat<Real>, GEODIM> &dref,
    std::array<math::DenseDMat<Real>, GEODIM> &dphys)
{
  for (Uint q = 0; q < jacobian[X].rows(); ++q)
  {
    for (Uint dim = 0; dim < GEODIM; ++dim)
    {
      m_J.insert_row(dim, jacobian[dim].row(q));
    }

    m_J.inv(m_J_inverse);

    for (Uint n = 0; n < dref[X].cols(); ++n)
    {
      for (Uint dim = 0; dim < GEODIM; ++dim)
      {
        m_dphi_ref[dim] = (dref[dim])(q, n);
      }

      m_dphi_phys = m_J_inverse * m_dphi_ref;

      for (Uint dim = 0; dim < GEODIM; ++dim)
      {
        dphys[dim](q, n) = m_dphi_phys[dim];
      }
    }
  }
}

// ============================================================================

template <Uint GEODIM, Uint TDIM>
void FEMetric<GEODIM, TDIM>::compute_derivatives_in_phys_space(
    const std::array<math::DenseDMat<Real>, GEODIM> &jacobian,
    const std::array<math::DenseDMat<Real>, GEODIM> &dref,
    std::array<math::DenseDMat<Real>, GEODIM> &dphys, math::DenseDVec<Real> &jac_at_qd_pt)
{
  for (Uint q = 0; q < jacobian[X].rows(); ++q)
  {
    for (Uint dim = 0; dim < GEODIM; ++dim)
    {
      m_J.insert_row(dim, jacobian[dim].row(q));
    }

    m_J.inv(m_J_inverse);

    jac_at_qd_pt[q] = m_J.det();

    for (Uint n = 0; n < dref[X].cols(); ++n)
    {
      for (Uint dim = 0; dim < GEODIM; ++dim)
      {
        m_dphi_ref[dim] = (dref[dim])(q, n);
      }

      m_dphi_phys = m_J_inverse * m_dphi_ref;

      for (Uint dim = 0; dim < GEODIM; ++dim)
      {
        dphys[dim](q, n) = m_dphi_phys[dim];
      }
    }
  }
}

// ============================================================================

// Specialization for 2D case: invert Jacobian matrix manually
template <>
void FEMetric<_2D, _2D>::compute_derivatives_in_phys_space(
    const std::array<math::DenseDMat<Real>, _2D> &jacobian,
    const std::array<math::DenseDMat<Real>, _2D> &dref,
    std::array<math::DenseDMat<Real>, _2D> &dphys, math::DenseDVec<Real> &jac_at_qd_pt)
{
  for (Uint q = 0; q < jacobian[X].rows(); ++q)
  {
    for (Uint dim = 0; dim < _2D; ++dim)
    {
      m_J.insert_row(dim, jacobian[dim].const_row(q));
    }

    jac_at_qd_pt[q]    = m_J(0, 0) * m_J(1, 1) - m_J(0, 1) * m_J(1, 0);
    const Real inv_jac = 1.0 / jac_at_qd_pt[q];

    for (Uint n = 0; n < dref[X].cols(); ++n)
    {
      m_dphi_ref[X0] = dref[X0](q, n);
      m_dphi_ref[X1] = dref[X1](q, n);

      m_dphi_phys[X] = inv_jac * (m_J(1, 1) * m_dphi_ref[XI0] - m_J(0, 1) * m_dphi_ref[XI1]);
      m_dphi_phys[Y] = inv_jac * (-m_J(1, 0) * m_dphi_ref[XI0] + m_J(0, 0) * m_dphi_ref[XI1]);

      dphys[X0](q, n) = m_dphi_phys[X0];
      dphys[X1](q, n) = m_dphi_phys[X1];
    }
  }
}

// ============================================================================

// Remark: specalization for 3D does not really bring any speedup ... (for RDM)
#if 0
// Specialization for 3D case: invert Jacobian matrix manually
template<>
void FEMetric<_3D,_3D>::compute_derivatives_in_phys_space(
                                  const std::array<math::DynamicMatrix<Real>,_3D>& jacobian,
                                  const std::array<math::DynamicMatrix<Real>,_3D>& dref,
                                        std::array<math::DynamicMatrix<Real>,_3D>& dphys,
                                                            math::DynamicVector<Real> & jac_at_qd_pt)
{
  for(Uint q = 0; q < jacobian[X].rows(); ++q)
  {
    for(Uint dim = 0; dim < _3D; ++dim)
    {
      m_J.insert_row(dim,jacobian[dim].row(q));
    }

    jac_at_qd_pt[q] = m_J(0,0)* ( m_J(1,1)*m_J(2,2)-m_J(2,1)*m_J(1,2) )
                    - m_J(0,1)* ( m_J(1,0)*m_J(2,2)-m_J(2,0)*m_J(1,2) )
                    + m_J(0,2)* ( m_J(1,0)*m_J(2,1)-m_J(2,0)*m_J(1,1) );

    const Real inv_jac = 1.0/jac_at_qd_pt[q];

    for(Uint n = 0; n < dref[X].cols(); ++n)
    {
      m_dphi_ref[X0] = dref[X0](q,n);
      m_dphi_ref[X1] = dref[X1](q,n);
      m_dphi_ref[X2] = dref[X2](q,n);

      m_dphi_phys[X0] = inv_jac * (   ( m_J(1,1)*m_J(2,2) - m_J(2,1)*m_J(1,2) ) * m_dphi_ref[XI0]
                                    - ( m_J(0,1)*m_J(2,2) - m_J(2,1)*m_J(0,2) ) * m_dphi_ref[XI1]
                                    + ( m_J(0,1)*m_J(1,2) - m_J(1,1)*m_J(0,2) ) * m_dphi_ref[XI2] );

      m_dphi_phys[X1] = inv_jac * ( - ( m_J(1,0)*m_J(2,2) - m_J(2,0)*m_J(1,2) ) * m_dphi_ref[XI0]
                                    + ( m_J(0,0)*m_J(2,2) - m_J(2,0)*m_J(0,2) ) * m_dphi_ref[XI1]
                                    - ( m_J(0,0)*m_J(1,2) - m_J(1,0)*m_J(0,2) ) * m_dphi_ref[XI2] );

      m_dphi_phys[X2] = inv_jac * (   ( m_J(1,0)*m_J(2,1) - m_J(2,0)*m_J(1,1) ) * m_dphi_ref[XI0]
                                    - ( m_J(0,0)*m_J(2,1) - m_J(2,0)*m_J(0,1) ) * m_dphi_ref[XI1]
                                    + ( m_J(0,0)*m_J(1,1) - m_J(1,0)*m_J(0,1) ) * m_dphi_ref[XI2] );

      dphys[X0](q,n) = m_dphi_phys[X0];
      dphys[X1](q,n) = m_dphi_phys[X1];
      dphys[X2](q,n) = m_dphi_phys[X2];
    }
  }
}
#endif

// ============================================================================

} // namespace solver

} // namespace pdekit

#endif
