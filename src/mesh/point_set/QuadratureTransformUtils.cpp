#include "mesh/point_set/QuadratureTransformUtils.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseVecView.hpp"
#include "math/unary_ops/VectorNorm.hpp"

namespace pdekit
{

namespace mesh
{

namespace detail
{

// ----------------------------------------------------------------------------
// Some useful functions to generate quadrature points from another quadrature
// ----------------------------------------------------------------------------

void fill_line_quadrature_permutation(const Uint nb_qd_pts,
                                      const mesh::EntityRealignCode &permutation_code,
                                      std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  // First set permutation_vec as identity permutation: p[i] = i;

  permutation_vec.resize(nb_qd_pts);
  for (Uint i = 0; i < permutation_vec.size(); ++i)
  {
    permutation_vec[i] = i;
  }

  // If this permutation is identity, we are done
  if (permutation_code.is_identity(ElemShape::Line))
  {
    return;
  }

  // Temporary vector
  std::vector<Uint> old_permutation(nb_qd_pts);
  std::vector<Uint> new_permutation(nb_qd_pts);

  // Apply all flips first
  for (Uint i = 0; i < permutation_code.nb_flips(); ++i)
  {
    old_permutation.swap(permutation_vec);

    // Fill flip permutation:
    new_permutation.resize(nb_qd_pts);
    for (Uint q = 0; q < nb_qd_pts; ++q)
    {
      new_permutation[q] = nb_qd_pts - q - 1;
    }

    for (Uint j = 0; j < new_permutation.size(); ++j)
    {
      permutation_vec[j] = old_permutation[new_permutation[j]];
    }
  }

  // Then apply all rotations
  for (Uint i = 0; i < permutation_code.nb_rotations(); ++i)
  {
    old_permutation.swap(permutation_vec);

    // Fill rotation permutation:
    new_permutation.resize(nb_qd_pts);
    for (Uint q = 0; q < nb_qd_pts; ++q)
    {
      new_permutation[q] = nb_qd_pts - q - 1;
    }

    for (Uint j = 0; j < new_permutation.size(); ++j)
    {
      permutation_vec[j] = old_permutation[new_permutation[j]];
    }
  }
}

// ----------------------------------------------------------------------------

void compute_triag_barycentric_coords(const math::DenseSVec<Real, 2> &V0,
                                      const math::DenseSVec<Real, 2> &V1,
                                      const math::DenseSVec<Real, 2> &V2,
                                      const math::DenseDMat<Real> &cart_coords,
                                      math::DenseDMat<Real> &bar_coords)
{
  // We suppose that each row of 'cart_coords' contains coordinates of one
  // quadrature point
  const Uint nb_qd_pts = cart_coords.rows();
  bar_coords.resize(nb_qd_pts, 3);

  // transf_mat is a matrix that transforms from cartesian coordinates
  // into barycentric coordinates
  math::DenseSMat<Real, 2, 2> tmp, transf_mat;
  tmp(0, 0) = V1[X0] - V0[X0];
  tmp(0, 1) = V2[X0] - V0[X0];
  tmp(1, 0) = V1[X1] - V0[X1];
  tmp(1, 1) = V2[X1] - V0[X1];

  tmp.inv(transf_mat);

  math::DenseSVec<Real, 2> diff, result;

  for (Uint qd_pt = 0; qd_pt < nb_qd_pts; ++qd_pt)
  {
    diff[X0] = cart_coords(qd_pt, X0) - V0[X0];
    diff[X1] = cart_coords(qd_pt, X1) - V0[X1];

    // Result contains lambda_1 and lambda_2, lambda_0 can be obtained as
    // lambda_0 = 1 - lambda_1 - lambda_2
    result = transf_mat * diff;

    bar_coords(qd_pt, 0) = 1. - result[0] - result[1];
    bar_coords(qd_pt, 1) = result[0];
    bar_coords(qd_pt, 2) = result[1];
  }
}

// ----------------------------------------------------------------------------

void canonical_triag_quadrature_rotation_permutation(const math::DenseDMat<Real> &cart_coords,
                                                     std::vector<Uint> &pvec)
{
  // Algorithm idea: express original quadrature points in barycentric
  // coordinates Then rotate the quadrature coordinates by performing cyclic
  // permutation of the corresponding barycentric coords UNDER THE ASSUMPTION
  // THAT THE QUADRATURE POINTS HAVE SYMMETRIC DISTRIBUTION, the new points
  // should map on the old points, therefore search for the correspondence
  // between the indexes of the original and rotated points - this is the
  // sought permutation

  // --------------------------------------------------------------------------
  // Step 1: compute the barycentric coordinates in the original quadrature
  // --------------------------------------------------------------------------

  math::DenseDMat<Real> barycentric_coords;

  // We suppose that the vertices of the reference triangle are
  // (-1,-1), (1,-1), (-1,1)
  math::DenseSVec<Real, 2> V0_ref, V1_ref, V2_ref;
  V0_ref[X0] = -1.;
  V0_ref[X1] = -1.;

  V1_ref[X0] = 1.;
  V1_ref[X1] = -1.;

  V2_ref[X0] = -1.;
  V2_ref[X1] = 1.;

  compute_triag_barycentric_coords(V0_ref, V1_ref, V2_ref, cart_coords, barycentric_coords);

  // --------------------------------------------------------------------------
  // Step 2: express the quadrature point coordinates in equilateral triangle
  //        - build the rotated coordinates in equilateral triangle
  // --------------------------------------------------------------------------

  math::DenseSVec<Real, 2> V0_equ, V1_equ, V2_equ;
  V0_equ[X0] = -1.;
  V0_equ[X1] = -1.;

  V1_equ[X0] = 1.;
  V1_equ[X1] = -1.;

  V2_equ[X0] = 0.;
  V2_equ[X1] = -1. + std::sqrt(3.);

  math::DenseDMat<Real> cart_coords_equidist;
  cart_coords_equidist.resize(cart_coords.rows(), cart_coords.cols());

  math::DenseDMat<Real> cart_coords_rot;
  cart_coords_rot.resize(cart_coords.rows(), cart_coords.cols());

  const Uint nb_qd_pts = cart_coords.rows();

  for (Uint qd_pt = 0; qd_pt < nb_qd_pts; ++qd_pt)
  {
    cart_coords_equidist(qd_pt, X0) = V0_equ[X0] * barycentric_coords(qd_pt, 0) +
                                      V1_equ[X0] * barycentric_coords(qd_pt, 1) +
                                      V2_equ[X0] * barycentric_coords(qd_pt, 2);
    cart_coords_equidist(qd_pt, X1) = V0_equ[X1] * barycentric_coords(qd_pt, 0) +
                                      V1_equ[X1] * barycentric_coords(qd_pt, 1) +
                                      V2_equ[X1] * barycentric_coords(qd_pt, 2);

    // The 'rotated' coordinates take the linear combination of barycentric
    // coordinates with index permutation 0->1, 1->2, 2->0
    cart_coords_rot(qd_pt, X0) = V0_equ[X0] * barycentric_coords(qd_pt, 1) +
                                 V1_equ[X0] * barycentric_coords(qd_pt, 2) +
                                 V2_equ[X0] * barycentric_coords(qd_pt, 0);
    cart_coords_rot(qd_pt, X1) = V0_equ[X1] * barycentric_coords(qd_pt, 1) +
                                 V1_equ[X1] * barycentric_coords(qd_pt, 2) +
                                 V2_equ[X1] * barycentric_coords(qd_pt, 0);
  }

  // --------------------------------------------------------------------------
  // Build the correspondence vector between 'old' and rotated coords - this
  // is the permutation we are looking for
  // --------------------------------------------------------------------------
  pvec.resize(nb_qd_pts, 0);

  math::DenseSVec<Real, 2> diff;
  Real min_point_dist = 0.0;

  for (Uint ipt = 0; ipt < nb_qd_pts; ++ipt)
  {
    min_point_dist                  = 1.e6;
    bool failed_to_find_permuted_pt = true;

    for (Uint jpt = 0; jpt < nb_qd_pts; ++jpt)
    {
      // If the points match, we have a new entry for the permutation
      // vector: diff = cart_coords_equidist.const_row(ipt) -
      // cart_coords_equidist.const_row(jpt);
      diff[X0] = cart_coords_equidist(ipt, X0) - cart_coords_rot(jpt, X0);
      diff[X1] = cart_coords_equidist(ipt, X1) - cart_coords_rot(jpt, X1);

      const Real point_dist = math::norm_max(diff);
      min_point_dist        = std::min(point_dist, min_point_dist);

      if (point_dist < 1.e-11)
      {
        // Now the point ipt (original index) is on position jpt (index
        // after rotation)
        pvec[jpt] = ipt;

        failed_to_find_permuted_pt = false;
        break;
      }
    } // Loop over jpt

    if (failed_to_find_permuted_pt == true)
    {
      std::cerr << "compute_quadrature_rotation_permutation: failed to find "
                   "permuted point \n"
                << "for quadrature point nr. " << ipt << " (nearest match " << min_point_dist << ")"
                << std::endl;
    }
  } // Loop over ipt
}

// ----------------------------------------------------------------------------

void canonical_triag_quadrature_flip_permutation(const math::DenseDMat<Real> &cart_coords,
                                                 std::vector<Uint> &pvec)
{
  // Algorithm idea: express original quadrature points in barycentric
  // coordinates Then flip the quadrature coordinates by swapping the
  // barycentric coordinates corresponding to vertices 0 and 1 UNDER THE
  // ASSUMPTION THAT THE QUADRATURE POINTS HAVE SYMMETRIC DISTRIBUTION, the
  // new points should map on the old points, therefore search for the
  // correspondence between the indexes of the original and rotated points -
  // this is the sought permutation

  // --------------------------------------------------------------------------
  // Step 1: compute the barycentric coordinates in the original quadrature
  // --------------------------------------------------------------------------

  math::DenseDMat<Real> barycentric_coords;

  // We suppose that the vertices of the reference triangle are
  // (-1,-1), (1,-1), (-1,1)
  math::DenseSVec<Real, 2> V0_ref, V1_ref, V2_ref;
  V0_ref[X0] = -1.;
  V0_ref[X1] = -1.;

  V1_ref[X0] = 1.;
  V1_ref[X1] = -1.;

  V2_ref[X0] = -1.;
  V2_ref[X1] = 1.;

  compute_triag_barycentric_coords(V0_ref, V1_ref, V2_ref, cart_coords, barycentric_coords);

  // --------------------------------------------------------------------------
  // Step 2: express the quadrature point coordinates in equilateral triangle
  //        - build the flipped coordinates in equilateral triangle
  // --------------------------------------------------------------------------

  math::DenseSVec<Real, 2> V0_equ, V1_equ, V2_equ;
  V0_equ[X0] = -1.;
  V0_equ[X1] = -1.;

  V1_equ[X0] = 1.;
  V1_equ[X1] = -1.;

  V2_equ[X0] = 0.;
  V2_equ[X1] = -1. + std::sqrt(3.);

  math::DenseDMat<Real> cart_coords_equidist;
  cart_coords_equidist.resize(cart_coords.rows(), cart_coords.cols());

  math::DenseDMat<Real> cart_coords_flip;
  cart_coords_flip.resize(cart_coords.rows(), cart_coords.cols());

  const Uint nb_qd_pts = cart_coords.rows();

  for (Uint qd_pt = 0; qd_pt < nb_qd_pts; ++qd_pt)
  {
    cart_coords_equidist(qd_pt, X0) = V0_equ[X0] * barycentric_coords(qd_pt, 0) +
                                      V1_equ[X0] * barycentric_coords(qd_pt, 1) +
                                      V2_equ[X0] * barycentric_coords(qd_pt, 2);
    cart_coords_equidist(qd_pt, X1) = V0_equ[X1] * barycentric_coords(qd_pt, 0) +
                                      V1_equ[X1] * barycentric_coords(qd_pt, 1) +
                                      V2_equ[X1] * barycentric_coords(qd_pt, 2);

    // The 'rotated' coordinates take the linear combination of barycentric
    // coordinates with index permutation 0->1, 1->2, 2->0
    cart_coords_flip(qd_pt, X0) = V0_equ[X0] * barycentric_coords(qd_pt, 1) +
                                  V1_equ[X0] * barycentric_coords(qd_pt, 0) +
                                  V2_equ[X0] * barycentric_coords(qd_pt, 2);
    cart_coords_flip(qd_pt, X1) = V0_equ[X1] * barycentric_coords(qd_pt, 1) +
                                  V1_equ[X1] * barycentric_coords(qd_pt, 0) +
                                  V2_equ[X1] * barycentric_coords(qd_pt, 2);
  }

  // --------------------------------------------------------------------------
  // Build the correspondence vector between 'old' and rotated coords - this
  // is the permutation we are looking for
  // --------------------------------------------------------------------------
  pvec.resize(nb_qd_pts, 0);

  math::DenseSVec<Real, 2> diff;
  Real min_point_dist = 0.0;

  for (Uint ipt = 0; ipt < nb_qd_pts; ++ipt)
  {
    min_point_dist                  = 1.e6;
    bool failed_to_find_permuted_pt = true;

    for (Uint jpt = 0; jpt < nb_qd_pts; ++jpt)
    {
      // If the points match, we have a new entry for the permutation
      // vector: diff = cart_coords_equidist.const_row(ipt) -
      // cart_coords_equidist.const_row(jpt);
      diff[X0] = cart_coords_equidist(ipt, X0) - cart_coords_flip(jpt, X0);
      diff[X1] = cart_coords_equidist(ipt, X1) - cart_coords_flip(jpt, X1);

      const Real point_dist = math::norm_max(diff);
      min_point_dist        = std::min(point_dist, min_point_dist);

      if (point_dist < 1.e-11)
      {
        // Now the point ipt (original index) is on position jpt (index
        // after flip)
        pvec[jpt]                  = ipt;
        failed_to_find_permuted_pt = false;
        break;
      }
    } // Loop over jpt

    if (failed_to_find_permuted_pt == true)
    {
      std::cerr << "compute_quadrature_flip_permutation: failed to find "
                   "permuted point \n"
                << "for quadrature point nr. " << ipt << " (nearest match " << min_point_dist << ")"
                << std::endl;
    }
  } // Loop over ipt
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace mesh

} // namespace pdekit
