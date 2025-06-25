#ifndef PDEKIT_Mesh_Quadrature_Transform_Utils_hpp
#define PDEKIT_Mesh_Quadrature_Transform_Utils_hpp

#include "math/DenseDMat.hpp"
#include "math/DenseSVec.hpp"
#include "mesh/EntityRealignCode.hpp"

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
                                      std::vector<Uint> &permutation_vec);

/// @param V0, V1, V2  ... vertices of triangle in which 'cart_coords' -
/// Cartesian
///                        coordinates of the quadrature points are defined
/// @param cart_coords ... Cartesian coordinates of the barycentric coordinates
///                        The matrix should have as many rows as there are
///                        quadrature points, and 2 colums
/// @param bar_coords  ... resulting barycentric coordinates

void compute_triag_barycentric_coords(const math::DenseSVec<Real, 2> &V0,
                                      const math::DenseSVec<Real, 2> &V1,
                                      const math::DenseSVec<Real, 2> &V2,
                                      const math::DenseDMat<Real> &cart_coords,
                                      math::DenseDMat<Real> &bar_coords);

// ----------------------------------------------------------------------------

/// Takes coordinates of quadrature and computes a permutation vector which says
/// how to reorder the quadrature points so that they correspond to the triangle
/// rotated clockwise so that the triangle vertex 0 becomes vertex 1, 1 becomes
/// 2 and 2 becomes 0 (i.e. cyclic permutation 0->1, 1->2, 2->0 which
/// corresponds to rotation by 120 degrees)

void canonical_triag_quadrature_rotation_permutation(const math::DenseDMat<Real> &cart_coords,
                                                     std::vector<Uint> &pvec);

// ----------------------------------------------------------------------------

/// Takes coordinates of quadrature and computes a permutation vector which says
/// how to reorder the quadrature points so that they correspond to a 'flipped
/// triangle' where it's vertices 0 and 1 have been swapped

void canonical_triag_quadrature_flip_permutation(const math::DenseDMat<Real> &cart_coords,
                                                 std::vector<Uint> &pvec);

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace mesh

} // namespace pdekit

#endif
