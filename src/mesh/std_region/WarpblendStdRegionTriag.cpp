#include "mesh/std_region/WarpblendStdRegionTriag.hpp"

#include "math/polynomials/JacobiPolynomial.hpp"
#include "math/polynomials/PolyLib.hpp"
#include "mesh/std_region/EquidistStdRegionTriag.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

/// A triangle 3 edges, 1 face and  0 volumes (but we store one default
/// 'invalid' entity as its volume)

const std::array<Uint, 4> WarpblendStdRegionTriag::TopologyStorage = {1, 3, 1, 1};

// ----------------------------------------------------------------------------

const math::DenseSVec<Real, 15> WarpblendStdRegionTriag::alpopt =
    math::values_list(0.0000)(0.0000)(1.4152)(0.1001)(0.2751)(0.9800)(1.0999)(1.2832)(1.3648)(
        1.4773)(1.4959)(1.5743)(1.5770)(1.6223)(1.6258);

// ----------------------------------------------------------------------------

Uint WarpblendStdRegionTriag::nb_dof(const Uint poly_order)
{
  return (poly_order + 1) * (poly_order + 2) / 2;
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTriag::fill_topology_relations(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences)
{
  EquidistStdRegionTriag::fill_topology_relations(incidences);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTriag::fill_edge_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  EquidistStdRegionTriag::fill_edge_to_node_connectivity(poly_order, entity_vec);

  // Override the type of the entities
  for (Uint e = 0; e < 3; ++e)
  {
    entity_vec[e]->change_point_set_type(PointSetID::Warpblend);
  }
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTriag::fill_face_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  EquidistStdRegionTriag::fill_face_to_node_connectivity(poly_order, entity_vec);

  // Override the topology type of the face
  StdRegionEntity &face = *entity_vec[0];
  face.change_point_set_type(PointSetID::Warpblend);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTriag::fill_volume_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  StdRegionEntity &current_entity = *entity_vec[0];
  // This entity will not be resized - it has length 0 by default

  const common::ArrayView<const Uint, _1D, Uint> vert_ids_view;
  const common::ArrayView<const Uint, _1D, Uint> p1_vert_flags_view;

  current_entity.construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined),
                           vert_ids_view, p1_vert_flags_view, 0);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTriag::fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord)
{
  // The coordinates have to be written in such a manner that first come
  // vertex nodes, then edge nodes. For the remaining interior nodes, we
  // repeat the procedure recursively. Imagine we remove the edges of the
  // triangle. We will be left with a smaller interior triangle. We apply the
  // same numbering to it again: first vertices, then internal edge nodes.
  // After, we peel its edges off and so on until no nodes are left

  const Uint nb_nodes = (poly_order + 1) * (poly_order + 2) / 2;
  coord.resize(nb_nodes, _2D);

  // First vertex  (index 0)
  coord(0, KSI) = -1.0;
  coord(0, ETA) = -1.0;
  // Second vertex (index 1)
  coord(1, KSI) = 1.0;
  coord(1, ETA) = -1.0;
  // Third vertex  (index 2)
  coord(2, KSI) = -1.0;
  coord(2, ETA) = 1.0;

  if (poly_order == P2)
  {
    coord(3, KSI) = 0.0;
    coord(3, ETA) = -1.0;
    coord(4, KSI) = 0.0;
    coord(4, ETA) = 0.0;
    coord(5, KSI) = -1.0;
    coord(5, ETA) = 0.0;
  }

  if (poly_order > P2)
  {
    // We will use barycentric coordinates L0, L1, L2: x = X0*L0 + X1*L1 +
    // X2*L2
    //                                                 y = Y1*L0 + Y1*L1 +
    //                                                 Y2*L2
    // Given (x,y), barycentric coordinates can be computed from the
    // relations above:
    //              L0 = 1/((X0-X2)*(Y1-Y2)-(Y0-Y2)*(X1-X2)) * [
    //              (Y1-Y2)*(x-X2)
    // - (X1-X2)*(y-Y2)]
    //              L1 = 1/((X0-X2)*(Y1-Y2)-(Y0-Y2)*(X1-X2)) *
    //              [-(Y0-Y2)*(x-X2)
    // + (X0-X2)*(y-Y2)]
    //              L2 = 1 - L0 - L1
    // For the reference triangle with vertices [X0,Y0] = [-1,-1]
    //                                          [X1,Y1] = [ 1,-1]
    //                                          [X2,Y2] = [-1, 1],
    // the equations above can be written as
    //              L0 = 1/4 * [ -2*(x-X2) - 2*(y-Y2)] = 1/4 * [ -2*(x+1) -
    // 2*(y-1)]
    //              L1 = 1/4 * [  2*(x-X2) ]           = 1/4 * [  2*(x+1) ]
    //              L2 = 1 - L0 - L1
    //

    Real alpha;

    if (poly_order < 16)
    {
      alpha = alpopt[poly_order - 1];
    }
    else
    {
      alpha = 5. / 3.;
    }

    // Barycentric coordinates
    math::DenseDVec<Real> L0(nb_nodes);
    math::DenseDVec<Real> L1(nb_nodes);
    math::DenseDVec<Real> L2(nb_nodes);
    math::DenseDVec<Real> blend0(nb_nodes);
    math::DenseDVec<Real> blend1(nb_nodes);
    math::DenseDVec<Real> blend2(nb_nodes);
    math::DenseDVec<Real> warp_factor0(nb_nodes);
    math::DenseDVec<Real> warp_factor1(nb_nodes);
    math::DenseDVec<Real> warp_factor2(nb_nodes);

    // std::cout << "Coordinates of warpblend nodes before optimization:" <<
    // std::endl;

    // Barycentric coordinates of first vertex
    L0[0] = 1.0;
    L1[0] = 0.0;
    L2[0] = 0.0;
    // Barycentric coordinates of second vertex
    L0[1] = 0.0;
    L1[1] = 1.0;
    L2[1] = 0.0;
    // Barycentric coordinates of third vertex
    L0[2] = 0.0;
    L1[2] = 0.0;
    L2[2] = 1.0;

    const Real dx = 2. / poly_order;
    const Real dy = 2. / poly_order;

    // Node coordinates of the first edge
    Uint inode = 3;
    for (Uint i = 1; i < poly_order; ++i)
    {
      const Real x = -1. + i * dx; // x-coordinate
      const Real y = -1.0;         // y-coordinate

      L0[inode] = 0.25 * ((-2.) * (x + 1.) - 2. * (y - 1.));
      L1[inode] = 0.25 * (2. * (x + 1.));
      L2[inode] = 1. - L0[inode] - L1[inode];
      inode++;
    }

    // Node coordinates of the second edge
    for (Uint i = 1; i < poly_order; ++i)
    {
      const Real x = 1. - i * dx;
      const Real y = -1. + i * dy;

      L0[inode] = 0.25 * ((-2.) * (x + 1.) - 2. * (y - 1.));
      L1[inode] = 0.25 * (2. * (x + 1.));
      L2[inode] = 1. - L0[inode] - L1[inode];
      inode++;
    }

    // Node coordinates of the third edge
    for (Uint i = 1; i < poly_order; ++i)
    {
      const Real x = -1.;
      const Real y = 1. - i * dy;

      L0[inode] = 0.25 * ((-2.) * (x + 1.) - 2. * (y - 1.));
      L1[inode] = 0.25 * (2. * (x + 1.));
      L2[inode] = 1. - L0[inode] - L1[inode];
      inode++;
    }

    // Offset denotes how many layers of nodes from the boundary of the main
    // triangle are we (how many contours have we 'peeled off' so far)
    Uint offset = 1;
    // Polynomial order to which the remaining interior nodes correspond
    Int current_poly_order = poly_order - 3;

    Real x, y;

    while (current_poly_order >= 0)
    {
      // First corner node, this can alos be center node (e.g. in P3
      // triangle) Therefore we might have to stop the filling after this
      // node.
      x = -1. + offset * dx;
      y = -1. + offset * dy;

      L0[inode] = 0.25 * ((-2.) * (x + 1.) - 2. * (y - 1.));
      L1[inode] = 0.25 * (2. * (x + 1.));
      L2[inode] = 1. - L0[inode] - L1[inode];

      inode++;
      if (current_poly_order == _0D)
        break;

      // Second corner
      x = 1. - (2 * offset) * dx;
      y = -1. + offset * dy;

      L0[inode] = 0.25 * ((-2.) * (x + 1.) - 2. * (y - 1.));
      L1[inode] = 0.25 * (2. * (x + 1.));
      L2[inode] = 1. - L0[inode] - L1[inode];

      inode++;

      // Third corner
      x = -1. + offset * dx;
      y = 1. - (2 * offset) * dy;

      L0[inode] = 0.25 * ((-2.) * (x + 1.) - 2. * (y - 1.));
      L1[inode] = 0.25 * (2. * (x + 1.));
      L2[inode] = 1. - L0[inode] - L1[inode];

      inode++;

      if (current_poly_order == _1D)
        break;

      // Print high-order edge nodes:

      // Node coordinates of the first sub-triag edge
      for (Int i = 1; i < current_poly_order; ++i)
      {
        x = -1. + (offset + i) * dx; // x-coordinate
        y = -1. + offset * dy;       // y-coordinate

        L0[inode] = 0.25 * ((-2.) * (x + 1.) - 2. * (y - 1.));
        L1[inode] = 0.25 * (2. * (x + 1.));
        L2[inode] = 1. - L0[inode] - L1[inode];

        inode++;
      }

      // Node coordinates of the second edge
      for (Int i = 1; i < current_poly_order; ++i)
      {
        x = 1. - (2 * offset + i) * dx;
        y = -1. + (offset + i) * dy;

        L0[inode] = 0.25 * ((-2.) * (x + 1.) - 2. * (y - 1.));
        L1[inode] = 0.25 * (2. * (x + 1.));
        L2[inode] = 1. - L0[inode] - L1[inode];

        inode++;
      }

      // Node coordinates of the third edge
      for (Int i = 1; i < current_poly_order; ++i)
      {
        x = -1. + offset * dx;
        y = 1. - (2 * offset + i) * dy;

        L0[inode] = 0.25 * ((-2.) * (x + 1.) - 2. * (y - 1.));
        L1[inode] = 0.25 * (2. * (x + 1.));
        L2[inode] = 1. - L0[inode] - L1[inode];

        inode++;
      }

      current_poly_order -= 3;
      offset++;
    }

    // Cartesian coordinates in EQUILATERAL triangle corresponding
    // to the barycentric coordinates L0, L1, L2
    math::DenseDVec<Real> x0equi(nb_nodes);
    math::DenseDVec<Real> x1equi(nb_nodes);

    for (Uint n = 0; n < nb_nodes; ++n)
    {
      x0equi[n] = -L0[n] + L1[n];
      x1equi[n] = (-L0[n] - L1[n] + 2. * L2[n]) / std::sqrt(3.0);
    }

    // Compute blending function at each node for each edge
    for (Uint i = 0; i < nb_nodes; ++i)
    {
      blend0[i] = 4. * L0[i] * L1[i];
      blend1[i] = 4. * L2[i] * L1[i];
      blend2[i] = 4. * L2[i] * L0[i];
    }

    // Amount of warp for each node, for each edge
    math::DenseDVec<Real> diffL(nb_nodes);

    for (Uint i = 0; i < nb_nodes; ++i)
    {
      // diffL[i] = L2[i] - L1[i];
      diffL[i] = L1[i] - L0[i];
    }

    warpfactor(poly_order, diffL, warp_factor0);

    for (Uint i = 0; i < nb_nodes; ++i)
    {
      // diffL[i] = L0[i] - L2[i];
      diffL[i] = L2[i] - L1[i];
    }

    warpfactor(poly_order, diffL, warp_factor1);

    for (Uint i = 0; i < nb_nodes; ++i)
    {
      // diffL[i] = L1[i] - L0[i];
      diffL[i] = L0[i] - L2[i];
    }

    warpfactor(poly_order, diffL, warp_factor2);

    // Combine blend & warp
    for (Uint i = 0; i < nb_nodes; ++i)
    {
      warp_factor0[i] = blend0[i] * warp_factor0[i] * (1. + (alpha * L2[i]) * (alpha * L2[i]));
      warp_factor1[i] = blend1[i] * warp_factor1[i] * (1. + (alpha * L0[i]) * (alpha * L0[i]));
      warp_factor2[i] = blend2[i] * warp_factor2[i] * (1. + (alpha * L1[i]) * (alpha * L1[i]));
    }

    // Accumulate deformations associated with each edge
    for (Uint i = 0; i < nb_nodes; ++i)
    {
      x0equi[i] = x0equi[i] + 1. * warp_factor0[i] + std::cos(2. * PI / 3.) * warp_factor1[i] +
                  std::cos(4. * PI / 3.) * warp_factor2[i];
      x1equi[i] = x1equi[i] + 0. * warp_factor0[i] + std::sin(2. * PI / 3.) * warp_factor1[i] +
                  std::sin(4. * PI / 3.) * warp_factor2[i];
    }

    // Go from (x0equi,x1equi) in equilateral triangle to coordinates in
    // standard triangle

    for (Uint i = 0; i < nb_nodes; ++i)
    {
      L2[i] = (std::sqrt(3.0) * x1equi[i] + 1.0) / 3.0;
      L0[i] = (-3.0 * x0equi[i] - std::sqrt(3.0) * x1equi[i] + 2.0) / 6.0;
      L1[i] = (3.0 * x0equi[i] - std::sqrt(3.0) * x1equi[i] + 2.0) / 6.0;

      coord(i, XI0) = -L0[i] + L1[i] - L2[i];
      coord(i, XI1) = -L0[i] - L1[i] + L2[i];
    }

  } // If polynomial order > 2
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTriag::fill_facet_normals(math::DenseDMat<Real> &normals)
{
  EquidistStdRegionTriag::fill_facet_normals(normals);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTriag::warpfactor(const Uint poly_order,
                                         const math::DenseDVec<Real> &xi_coord,
                                         math::DenseDVec<Real> &warp_factor)
{
  // Purpose  : Compute scaled warp function at order N based on rout
  // interpolation nodes

  // Compute LGL and equidistant node distribution
  // Compute Gauss-Lobatto distribution of nodes on biunit interval
  const Uint nb_pts_on_interval = poly_order + 1;
  math::DenseDVec<Real> LGL(nb_pts_on_interval);
  math::DenseDVec<Real> x_equidist(nb_pts_on_interval);
  math::DenseDVec<Real> tmp(nb_pts_on_interval - 2);

  math::jac_zeros(nb_pts_on_interval - 2, tmp, 1, 1);

  LGL[0] = -1.0;
  LGL[1] = 1.0;

  for (Uint i = 0; i < (nb_pts_on_interval - 2); ++i)
  {
    LGL[i + 2] = tmp[i];
  }

  // std::cout << "LGL = " << LGL << std::endl;

  // Compute equidistant distribution of nodes on biunit interval
  const Real dx = 2. / (nb_pts_on_interval - 1);

  x_equidist[0] = -1.0;
  x_equidist[1] = 1.0;

  for (Uint i = 1; i < (nb_pts_on_interval - 1); ++i)
  {
    x_equidist[i + 1] = -1. + i * dx;
  }

  // Compute values of Lagrange shape functions in given point set
  // Starting from equidistant point set as reference
  // We will compute
  // 1) a matrix of Jacobi polynomials in equidistant nodes on biunit interval
  // called V_equidist
  // 2) a matrix of Jacobi polynomials up to order 'poly_order' in the nodes
  // 'xi_coord'. This matrix
  //    will be called V_jacobi
  // 3) finally, the matrix of Lagrange polynomial values in 'xi_coord',
  // called V_lagrange, will be given by
  //
  //    V_lagrange = V_jacobi * inverse(V_equidist)
  //

  math::DenseDMat<Real> V_equidist(poly_order + 1, poly_order + 1);
  math::DenseDMat<Real> V_equidist_inv(poly_order + 1, poly_order + 1);

  const Uint nb_xi_pts = xi_coord.size();
  // std::cout << "nb_xi_pts = " << nb_xi_pts << std::endl;

  math::DenseDMat<Real> V_jacobi(nb_xi_pts, poly_order + 1);
  math::DenseDMat<Real> V_lagrange(nb_xi_pts, poly_order + 1);

  // Generate the Vandermonde matrix on the equidistant point set
  math::JacobiPolynomial jp;

  // std::cout << "Filling V_equidist ... " << std::endl;

  for (Uint i = 0; i < (poly_order + 1); ++i)
  {
    for (Uint j = 0; j < (poly_order + 1); ++j)
    {
      V_equidist(i, j) = jp(j, 0, 0, x_equidist[i]);
    }
  }

  V_equidist.inv(V_equidist_inv);

  // Generate the (transpose of the) Vandermonde matrix on the xi points
  // given as input to this function
  for (Uint i = 0; i < nb_xi_pts; ++i)
  {
    for (Uint j = 0; j < (poly_order + 1); ++j)
    {
      V_jacobi(i, j) = jp(j, 0.0, 0.0, xi_coord[i]);
    }
  }

  V_lagrange = V_jacobi * V_equidist_inv;

  // Clip off very small values
  for (Uint r = 0; r < V_lagrange.rows(); ++r)
  {
    for (Uint c = 0; c < V_lagrange.cols(); ++c)
    {
      if (std::abs(V_lagrange(r, c)) < 1.e-14)
        V_lagrange(r, c) = 0.0;
    }
  }

  // Now we will re-use the vector LGL so that it stores the values of
  // LGL-x_equidist
  // This means that we have to subtract x_equidist[i] from existing value
  // LGL[i] for each i

  for (Uint pt = 0; pt < LGL.size(); ++pt)
  {
    LGL[pt] -= x_equidist[pt];
  }

  warp_factor.resize(nb_pts_on_interval);

  warp_factor = V_lagrange * LGL;

  for (Uint i = 0; i < warp_factor.size(); ++i)
  {
    if (std::abs(xi_coord[i]) < (1. - 1.e-10))
    {
      const Real sf  = 1.0 - xi_coord[i] * xi_coord[i];
      warp_factor[i] = warp_factor[i] / sf;
    }
  }
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTriag::fill_permutation(const Uint poly_order,
                                               const EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
  EquidistStdRegionTriag::fill_permutation(poly_order, permutation_code, permutation_vec);
} // End of method

// ----------------------------------------------------------------------------

void WarpblendStdRegionTriag::fill_rotation_permutation(const Uint poly_order,
                                                        std::vector<Uint> &permutation_vec)
{
  EquidistStdRegionTriag::fill_rotation_permutation(poly_order, permutation_vec);
} // End of method

// ----------------------------------------------------------------------------

void WarpblendStdRegionTriag::fill_flip_permutation(const Uint poly_order,
                                                    std::vector<Uint> &permutation_vec)
{
  EquidistStdRegionTriag::fill_flip_permutation(poly_order, permutation_vec);
} // End of method

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
