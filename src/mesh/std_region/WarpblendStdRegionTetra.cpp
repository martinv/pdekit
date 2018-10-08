#include "mesh/std_region/WarpblendStdRegionTetra.hpp"

#include "math/DenseConstVecView.hpp"
#include "math/DenseSMat.hpp"
#include "math/polynomials/PolyLib.hpp"
#include "math/unary_ops/VectorNorm.hpp"
#include "mesh/std_region/EquidistStdRegionTetra.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

/// A tetra has 1 volume, 4 triangular faces and 6 edges
const std::array<Uint, 4> WarpblendStdRegionTetra::TopologyStorage = {1, 6, 4, 1};

// ----------------------------------------------------------------------------

const math::DenseSVec<Real, 15> WarpblendStdRegionTetra::alpopt =
    math::values_list(0.)(0.)(0.)(0.1002)(1.1332)(1.5608)(1.3413)(1.2577)(1.1603)(1.10153)(0.6080)(
        0.4523)(0.8856)(0.8717)(0.9655);

// ----------------------------------------------------------------------------

/// A TETRAHEDRON OF POLYNOMIAL ORDER P HAS (P+1)*(P+2)*(P+3)/6 NODES
Uint WarpblendStdRegionTetra::nb_dof(const Uint poly_order)
{
  return ((poly_order + 1) * (poly_order + 2) * (poly_order + 3) / 6);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTetra::fill_topology_relations(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences)
{
  EquidistStdRegionTetra::fill_topology_relations(incidences);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTetra::fill_edge_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  EquidistStdRegionTetra::fill_edge_to_node_connectivity(poly_order, entity_vec);

  // Override the type of the entities
  for (Uint e = 0; e < 6; ++e)
  {
    entity_vec[e]->change_point_set_type(PointSetID::Warpblend);
  }
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTetra::fill_face_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  EquidistStdRegionTetra::fill_face_to_node_connectivity(poly_order, entity_vec);

  // Override the topology type of the face
  for (Uint e = 0; e < 4; ++e)
  {
    entity_vec[e]->change_point_set_type(PointSetID::Warpblend);
  }

} // end of method

// ----------------------------------------------------------------------------

void WarpblendStdRegionTetra::fill_volume_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  EquidistStdRegionTetra::fill_volume_to_node_connectivity(poly_order, entity_vec);

  // Set the topology type of the volume
  StdRegionEntity &volume = *entity_vec[0];
  volume.change_point_set_type(PointSetID::Warpblend);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTetra::fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord)
{
  EquidistStdRegionTetra::fill_coordinates(poly_order, coord);

  /*
  coord.resize(nb_nodes, _3D);

  Uint sk = 0;

  for (Uint n = 0; n < (poly_order + 1); ++n)
  {
    for (Uint m = 0; m < (poly_order + 1 - n); ++m)
    {
      for (Uint q = 0; q < (poly_order + 1 - n - m); ++q)
      {
        coord(sk, X0) = -1. + q * 2. / poly_order;
        coord(sk, X1) = -1. + m * 2. / poly_order;
        coord(sk, X2) = -1. + n * 2. / poly_order;
        sk++;
      }
    }
  }
  */

  // If the polynomial order is less than 3, the Warpblend
  // and equidistant nodes are identical, we do not need to continue
  if (poly_order < P3)
  {
    return;
  }

  const Real tol = 1.e-10;

  const Real alpha = poly_order < 15 ? alpopt[poly_order - 1] : 1.0;

  const math::DenseSVec<Real, _3D> v1 =
      math::values_list(-1.)(-1. / std::sqrt(3.))(-1. / std::sqrt(6.));
  const math::DenseSVec<Real, _3D> v2 =
      math::values_list(1.)(-1. / std::sqrt(3.))(-1. / std::sqrt(6.));
  const math::DenseSVec<Real, _3D> v3 =
      math::values_list(0.)(2. / std::sqrt(3))(-1. / std::sqrt(6.));
  const math::DenseSVec<Real, _3D> v4 = math::values_list(0.)(0.)(3. / std::sqrt(6.));

  // Orthogonal axis tangents on faces 1-4
  math::DenseSMat<Real, 4, _3D> t1, t2;
  for (Uint d = 0; d < _3D; ++d)
  {
    t1(0, d) = v2[d] - v1[d];
    t1(1, d) = v2[d] - v1[d];
    t1(2, d) = v3[d] - v2[d];
    t1(3, d) = v3[d] - v1[d];
    t2(0, d) = v3[d] - 0.5 * (v1[d] + v2[d]);
    t2(1, d) = v4[d] - 0.5 * (v1[d] + v2[d]);
    t2(2, d) = v4[d] - 0.5 * (v2[d] + v3[d]);
    t2(3, d) = v4[d] - 0.5 * (v1[d] + v3[d]);
  }

  // Normalize the tangent vectors
  for (Uint f = 0; f < 4; ++f)
  {
    const Real inv_norm_t1 = 1.0 / math::norm_e2(t1.const_row_transp(f));
    const Real inv_norm_t2 = 1.0 / math::norm_e2(t2.const_row_transp(f));

    for (Uint d = 0; d < _3D; ++d)
    {
      t1(f, d) *= inv_norm_t1;
      t2(f, d) *= inv_norm_t2;
    }
  }

  // Warp and blend for each face
  const Uint nb_nodes = ((poly_order + 1) * (poly_order + 2) * (poly_order + 3) / 6);

  // Form undeformed coordinates
  math::DenseDVec<Real> La(nb_nodes), Lb(nb_nodes), Lc(nb_nodes), Ld(nb_nodes);
  math::DenseDMat<Real> tmp_XYZ(nb_nodes, _3D);
  tmp_XYZ.fill(0.0);

  for (Uint n = 0; n < nb_nodes; ++n)
  {
    const Real L1 = 0.5 * (1. + coord(n, XI2));
    const Real L2 = 0.5 * (1. + coord(n, XI1));
    const Real L3 = -0.5 * (1. + coord(n, XI0) + coord(n, XI1) + coord(n, XI2));
    const Real L4 = 0.5 * (1. + coord(n, XI0));

    for (Uint d = 0; d < _3D; ++d)
    {
      tmp_XYZ(n, d) = L3 * v1[d] + L4 * v2[d] + L2 * v3[d] + L1 * v4[d];
    }
  }

  math::DenseDMat<Real> shift(nb_nodes, _3D);
  shift.fill(0.0);

  math::DenseDVec<Real> warp1(nb_nodes);
  math::DenseDVec<Real> warp2(nb_nodes);
  warp1.fill(0.0);
  warp2.fill(0.0);

  math::DenseDVec<Real> volume_blend(nb_nodes);

  // Main loop over faces
  for (Uint f = 0; f < 4; ++f)
  {
    switch (f)
    {
      case 0:
      {
        for (Uint n = 0; n < nb_nodes; ++n)
        {
          const Real L1 = 0.5 * (1. + coord(n, XI2));
          const Real L2 = 0.5 * (1. + coord(n, XI1));
          const Real L3 = -0.5 * (1. + coord(n, XI0) + coord(n, XI1) + coord(n, XI2));
          const Real L4 = 0.5 * (1. + coord(n, XI0));

          La[n] = L1;
          Lb[n] = L2;
          Lc[n] = L3;
          Ld[n] = L4;
        }
        break;
      }
      case 1:
      {
        for (Uint n = 0; n < nb_nodes; ++n)
        {
          const Real L1 = 0.5 * (1. + coord(n, XI2));
          const Real L2 = 0.5 * (1. + coord(n, XI1));
          const Real L3 = -0.5 * (1. + coord(n, XI0) + coord(n, XI1) + coord(n, XI2));
          const Real L4 = 0.5 * (1. + coord(n, XI0));

          La[n] = L2;
          Lb[n] = L1;
          Lc[n] = L3;
          Ld[n] = L4;
        }
        break;
      }
      case 2:
      {
        for (Uint n = 0; n < nb_nodes; ++n)
        {
          const Real L1 = 0.5 * (1. + coord(n, XI2));
          const Real L2 = 0.5 * (1. + coord(n, XI1));
          const Real L3 = -0.5 * (1. + coord(n, XI0) + coord(n, XI1) + coord(n, XI2));
          const Real L4 = 0.5 * (1. + coord(n, XI0));

          La[n] = L3;
          Lb[n] = L1;
          Lc[n] = L4;
          Ld[n] = L2;
        }
        break;
      }
      case 3:
      {
        for (Uint n = 0; n < nb_nodes; ++n)
        {
          const Real L1 = 0.5 * (1. + coord(n, XI2));
          const Real L2 = 0.5 * (1. + coord(n, XI1));
          const Real L3 = -0.5 * (1. + coord(n, XI0) + coord(n, XI1) + coord(n, XI2));
          const Real L4 = 0.5 * (1. + coord(n, XI0));

          La[n] = L4;
          Lb[n] = L1;
          Lc[n] = L3;
          Ld[n] = L2;
        }
        break;
      }
    };

    warp_shift_face_3D(poly_order, alpha, Lb, Lc, Ld, warp1, warp2);

    // compute volume blending
    for (Uint i = 0; i < nb_nodes; ++i)
    {
      volume_blend[i] = Lb[i] * Lc[i] * Ld[i];
    }

    // Modify linear blend
    for (Uint i = 0; i < nb_nodes; ++i)
    {
      const Real denom = (Lb[i] + 0.5 * La[i]) * (Lc[i] + 0.5 * La[i]) * (Ld[i] + 0.5 * La[i]);

      if (denom > tol)
      {
        volume_blend[i] = (1. + std::pow(alpha * La[i], 2.)) * volume_blend[i] / denom;
      }
    }

    // compute warp & blend
    // shift = shift+(blend.*warp1)*t1(face,:) + (blend.*warp2)*t2(face,:);
    for (Uint i = 0; i < nb_nodes; ++i)
    {
      const Real bw1 = volume_blend[i] * warp1[i];
      const Real bw2 = volume_blend[i] * warp2[i];
      for (Uint d = 0; d < _3D; ++d)
      {
        shift(i, d) += bw1 * t1(f, d) + bw2 * t2(f, d);
      }
    }

    // fix face warp
    for (Uint i = 0; i < nb_nodes; ++i)
    {
      const Uint tol_sum = (Lb[i] > tol) + (Lc[i] > tol) + (Ld[i] > tol);
      if ((La[i] < tol) && (tol_sum < 3))
      {
        for (Uint d = 0; d < _3D; ++d)
        {
          shift(i, d) = warp1[i] * t1(f, d) + warp2[i] * t2(f, d);
        }
      }
    }

  } // Loop over faces

  // Add the shift to the node coordinates
  for (Uint i = 0; i < nb_nodes; ++i)
  {
    for (Uint d = 0; d < _3D; ++d)
    {
      tmp_XYZ(i, d) += shift(i, d);
    }
  }

  // Convert back to the coordinate system of reference tetrahedron: compute
  // barycentric coordinates for each node and use them to determine
  // the position of the node in the reference tetrahedron
  math::DenseSMat<Real, _3D, _3D> M_barycentric, inv_M_barycentric;
  math::DenseSVec<Real, _3D> lambda, pt_coord_barycentric;

  const math::DenseSVec<Real, 3> v1_equidist = math::values_list(-1.)(-1.)(-1.);
  const math::DenseSVec<Real, 3> v2_equidist = math::values_list(1.)(-1.)(-1.);
  const math::DenseSVec<Real, 3> v3_equidist = math::values_list(-1.)(1.)(-1.);
  const math::DenseSVec<Real, 3> v4_equidist = math::values_list(-1.)(-1.)(1.);

  // Intentionally hardcode the coordinates of the first four nodes,
  // which HAVE to be set exactly to the coordinates of vertices of reference
  // tetrahedron
  coord(0, X0) = -1.;
  coord(0, X1) = -1.;
  coord(0, X2) = -1.;

  coord(1, X0) = 1.;
  coord(1, X1) = -1.;
  coord(1, X2) = -1.;

  coord(2, X0) = -1.;
  coord(2, X1) = 1.;
  coord(2, X2) = -1.;

  coord(3, X0) = -1.;
  coord(3, X1) = -1.;
  coord(3, X2) = 1.;

  for (Uint i = 4; i < nb_nodes; ++i)
  {
    M_barycentric(X0, 0) = v1[X0] - v4[X0];
    M_barycentric(X0, 1) = v2[X0] - v4[X0];
    M_barycentric(X0, 2) = v3[X0] - v4[X0];

    M_barycentric(X1, 0) = v1[X1] - v4[X1];
    M_barycentric(X1, 1) = v2[X1] - v4[X1];
    M_barycentric(X1, 2) = v3[X1] - v4[X1];

    M_barycentric(X2, 0) = v1[X2] - v4[X2];
    M_barycentric(X2, 1) = v2[X2] - v4[X2];
    M_barycentric(X2, 2) = v3[X2] - v4[X2];

    pt_coord_barycentric[X0] = tmp_XYZ(i, X0) - v4[X0];
    pt_coord_barycentric[X1] = tmp_XYZ(i, X1) - v4[X1];
    pt_coord_barycentric[X2] = tmp_XYZ(i, X2) - v4[X2];

    M_barycentric.inv(inv_M_barycentric);

    lambda = inv_M_barycentric * pt_coord_barycentric;

    for (Uint d = 0; d < _3D; ++d)
    {
      coord(i, d) = lambda[0] * v1_equidist[d] + lambda[1] * v2_equidist[d] +
                    lambda[2] * v3_equidist[d] +
                    (1. - lambda[0] - lambda[1] - lambda[2]) * v4_equidist[d];
    }
  }
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTetra::fill_facet_normals(math::DenseDMat<Real> &normals)
{
  EquidistStdRegionTetra::fill_facet_normals(normals);
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTetra::fill_permutation(const Uint poly_order,
                                               const EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTetra::warp_shift_face_3D(const Uint poly_order, const Real alpha,
                                                 const math::DenseDVec<Real> &L1,
                                                 const math::DenseDVec<Real> &L2,
                                                 const math::DenseDVec<Real> &L3,
                                                 math::DenseDVec<Real> &warp1,
                                                 math::DenseDVec<Real> &warp2)
{
  warp1.fill(0.0);
  warp2.fill(0.0);

  // 1) Compute GAuss-Lobatto-Legendre node distribution
  math::DenseDVec<Real> tmp(poly_order - 1);
  math::DenseDVec<Real> gauss_x(poly_order + 1);
  // math::zwgj(gauss_x, weights, poly_order-1, 1, 1);
  math::jac_zeros(poly_order - 1, tmp, 1, 1);
  gauss_x[0]          = 1.0;
  gauss_x[poly_order] = -1.0;
  // Reverse the order of the nodes - multiply by -1
  for (Uint i = 0; i < (poly_order - 1); ++i)
  {
    gauss_x[i + 1] = -tmp[i];
  }

  // 2) Compute blending function at each node for each edge
  const Uint nb_nodes = ((poly_order + 1) * (poly_order + 2) * (poly_order + 3) / 6);

  math::DenseDVec<Real> blend1(nb_nodes);
  math::DenseDVec<Real> blend2(nb_nodes);
  math::DenseDVec<Real> blend3(nb_nodes);

  for (Uint n = 0; n < nb_nodes; ++n)
  {
    blend1[n] = L2[n] * L3[n];
    blend2[n] = L1[n] * L3[n];
    blend3[n] = L1[n] * L2[n];
  }

  // 3) amount of warp for each node, for each edge
  math::DenseDVec<Real> warpfactor1(nb_nodes);
  math::DenseDVec<Real> warpfactor2(nb_nodes);
  math::DenseDVec<Real> warpfactor3(nb_nodes);

  tmp.resize(nb_nodes);

  tmp = L3 - L2;
  evalwarp(poly_order, gauss_x, tmp, warpfactor1);
  for (Uint i = 0; i < warpfactor1.size(); ++i)
  {
    warpfactor1[i] *= 4.0;
  }

  tmp = L1 - L3;
  evalwarp(poly_order, gauss_x, tmp, warpfactor2);
  for (Uint i = 0; i < warpfactor2.size(); ++i)
  {
    warpfactor2[i] *= 4.0;
  }

  tmp = L2 - L1;
  evalwarp(poly_order, gauss_x, tmp, warpfactor3);
  for (Uint i = 0; i < warpfactor3.size(); ++i)
  {
    warpfactor3[i] *= 4.0;
  }

  // 4) combine blend & warp
  math::DenseDVec<Real> face_warp1(nb_nodes);
  math::DenseDVec<Real> face_warp2(nb_nodes);
  math::DenseDVec<Real> face_warp3(nb_nodes);

  for (Uint i = 0; i < nb_nodes; ++i)
  {
    face_warp1[i] = blend1[i] * warpfactor1[i] * (1. + std::pow(alpha * L1[i], 2));
    face_warp2[i] = blend2[i] * warpfactor2[i] * (1. + std::pow(alpha * L2[i], 2));
    face_warp3[i] = blend3[i] * warpfactor3[i] * (1. + std::pow(alpha * L3[i], 2));
  }

  // 5) evaluate shift in equilateral triangle
  for (Uint i = 0; i < nb_nodes; ++i)
  {
    warp1[i] = 1. * face_warp1[i] + std::cos(2. * PI / 3.) * face_warp2[i] +
               std::cos(4. * PI / 3.) * face_warp3[i];
    warp2[i] = 0. * face_warp1[i] + std::sin(2. * PI / 3.) * face_warp2[i] +
               std::sin(4. * PI / 3.) * face_warp3[i];
  }
}

// ----------------------------------------------------------------------------

void WarpblendStdRegionTetra::evalwarp(const Uint poly_order, const math::DenseDVec<Real> &xnodes,
                                       const math::DenseDVec<Real> &xout,
                                       math::DenseDVec<Real> &warp)
{
  // Purpose: compute one-dimensional edge warping function

  warp.resize(xout.size());

  math::DenseDVec<Real> xeq(poly_order + 1);
  math::DenseDVec<Real> d(xout.size());

  for (Uint i = 0; i < (poly_order + 1); ++i)
  {
    xeq[i] = -1 + 2. * (poly_order - i) / poly_order;
  }

  for (Uint i = 0; i < (poly_order + 1); ++i)
  {
    d.fill(0.0);

    for (Uint k = 0; k < d.size(); ++k)
    {
      d[k] = xnodes[i] - xeq[i];
    }

    for (Uint j = 1; j < poly_order; ++j)
    {
      if (i != j)
      {
        for (Uint k = 0; k < d.size(); ++k)
        {
          d[k] *= (xout[k] - xeq[j]) / (xeq[i] - xeq[j]);
        }
      }
    }

    if (i != 0)
    {
      for (Uint k = 0; k < d.size(); ++k)
      {
        d[k] /= -(xeq[i] - xeq[0]);
      }
    }

    if (i != poly_order)
    {
      for (Uint k = 0; k < d.size(); ++k)
      {
        d[k] /= (xeq[i] - xeq[poly_order]);
      }
    }

    for (Uint k = 0; k < d.size(); ++k)
    {
      warp[k] += d[k];
    }
  } // Loop over i
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
