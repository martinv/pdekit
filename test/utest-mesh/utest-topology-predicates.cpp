/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE topology_predicates_test
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iostream>

/// PDEKIT headers
#include "mesh/TopologyPredicates.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

typedef TopologyPredicates TP;

BOOST_AUTO_TEST_CASE(topology_predicates_utest)
{

  /// --------------------------------------------------------------------------
  /// Compute triangle permutation
  /// a) rotation
  /// --------------------------------------------------------------------------

  Uint p3tri_perm_rot[10] = {2, 0, 1, 7, 8, 3, 4, 5, 6, 9};
  Uint p4tri_perm_rot[15] = {2, 0, 1, 9, 10, 11, 3, 4, 5, 6, 7, 8, 14, 12, 13};
  Uint p5tri_perm_rot[21] = {2, 0, 1, 11, 12, 13, 14, 3,  4,  5, 6,
                             7, 8, 9, 10, 17, 15, 16, 20, 18, 19};
  Uint p6tri_perm_rot[28] = {2, 0,  1,  13, 14, 15, 16, 17, 3,  4,  5,  6,  7,  8,
                             9, 10, 11, 12, 20, 18, 19, 25, 26, 21, 22, 23, 24, 27};
  Uint p7tri_perm_rot[36] = {2,  0,  1,  15, 16, 17, 18, 19, 20, 3,  4,  5,
                             6,  7,  8,  9,  10, 11, 12, 13, 14, 23, 21, 22,
                             30, 31, 32, 24, 25, 26, 27, 28, 29, 35, 33, 34};

  const EntityRealignCode rotation_perm_code = EntityRealignCode::single_rotation(ElemShape::Triag);
  const EntityRealignCode flip_perm_code     = EntityRealignCode::single_flip(ElemShape::Triag);

  EntityDofRealign p(
      std::make_pair(PointSetTag::string_to_tag("Triag-P3-Equidist"), rotation_perm_code));
  std::cout << "Rotation permutation for p3 triangle = " << std::endl;
  p.get().print();
  std::cout << std::endl;

  // Check rotation permutation for p3 triangle
  BOOST_CHECK(p.get().size() == 10);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p3tri_perm_rot[i]);
  }

  // Check rotation permutation for p4 triangle
  p.change_type(PointSetTag::string_to_tag("Triag-P4-Equidist"), rotation_perm_code);
  BOOST_CHECK(p.get().size() == 15);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p4tri_perm_rot[i]);
  }

  // Check rotation permutation for p5 triangle
  p.change_type(PointSetTag::string_to_tag("Triag-P5-Equidist"), rotation_perm_code);
  BOOST_CHECK(p.get().size() == 21);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p5tri_perm_rot[i]);
  }

  // Check rotation permutation for p6 triangle
  p.change_type(PointSetTag::string_to_tag("Triag-P6-Equidist"), rotation_perm_code);
  BOOST_CHECK(p.get().size() == 28);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p6tri_perm_rot[i]);
  }

  // Check rotation permutation for p4 triangle
  p.change_type(PointSetTag::string_to_tag("Triag-P7-Equidist"), rotation_perm_code);
  BOOST_CHECK(p.get().size() == 36);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p7tri_perm_rot[i]);
  }

  /// --------------------------------------------------------------------------
  /// Compute triangle permutation
  /// b) flip
  /// --------------------------------------------------------------------------

  Uint p3tri_perm_flip[10] = {1, 0, 2, 4, 3, 8, 7, 6, 5, 9};
  Uint p4tri_perm_flip[15] = {1, 0, 2, 5, 4, 3, 11, 10, 9, 8, 7, 6, 13, 12, 14};
  Uint p5tri_perm_flip[21] = {1,  0, 2, 6, 5,  4,  3,  14, 13, 12, 11,
                              10, 9, 8, 7, 16, 15, 17, 18, 20, 19};
  Uint p6tri_perm_flip[28] = {1,  0,  2, 7, 6,  5,  4,  3,  17, 16, 15, 14, 13, 12,
                              11, 10, 9, 8, 19, 18, 20, 22, 21, 26, 25, 24, 23, 27};
  Uint p7tri_perm_flip[36] = {1,  0,  2,  8,  7,  6,  5,  4,  3,  20, 19, 18,
                              17, 16, 15, 14, 13, 12, 11, 10, 9,  22, 21, 23,
                              26, 25, 24, 32, 31, 30, 29, 28, 27, 34, 33, 35};

  // Check flip permutation for p3 triangle
  p.change_type(PointSetTag::string_to_tag("Triag-P3-Equidist"), flip_perm_code);
  BOOST_CHECK(p.get().size() == 10);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p3tri_perm_flip[i]);
  }

  // Check flip permutation for p4 triangle
  p.change_type(PointSetTag::string_to_tag("Triag-P4-Equidist"), flip_perm_code);
  BOOST_CHECK(p.get().size() == 15);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p4tri_perm_flip[i]);
  }

  // Check flip permutation for p5 triangle
  p.change_type(PointSetTag::string_to_tag("Triag-P5-Equidist"), flip_perm_code);
  BOOST_CHECK(p.get().size() == 21);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p5tri_perm_flip[i]);
  }
  // Check flip permutation for p6 triangle
  p.change_type(PointSetTag::string_to_tag("Triag-P6-Equidist"), flip_perm_code);
  BOOST_CHECK(p.get().size() == 28);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p6tri_perm_flip[i]);
  }

  // Check flip permutation for p7 triangle
  p.change_type(PointSetTag::string_to_tag("Triag-P7-Equidist"), flip_perm_code);
  BOOST_CHECK(p.get().size() == 36);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p7tri_perm_flip[i]);
  }

  /// --------------------------------------------------------------------------
  /// Compute quadrilateral permutation
  /// a) rotation
  /// --------------------------------------------------------------------------

  Uint p3quad_perm_rot[16] = {3, 0, 1, 2, 10, 11, 4, 5, 6, 7, 8, 9, 15, 12, 13, 14};
  Uint p4quad_perm_rot[25] = {3,  0,  1,  2,  13, 14, 15, 4,  5,  6,  7,  8, 9,
                              10, 11, 12, 19, 16, 17, 18, 23, 20, 21, 22, 24};
  Uint p5quad_perm_rot[36] = {3,  0,  1,  2,  16, 17, 18, 19, 4,  5,  6,  7,
                              8,  9,  10, 11, 12, 13, 14, 15, 23, 20, 21, 22,
                              30, 31, 24, 25, 26, 27, 28, 29, 35, 32, 33, 34};
  Uint p6quad_perm_rot[49] = {3,  0,  1,  2,  19, 20, 21, 22, 23, 4,  5,  6,  7,  8,  9,  10, 11,
                              12, 13, 14, 15, 16, 17, 18, 27, 24, 25, 26, 37, 38, 39, 28, 29, 30,
                              31, 32, 33, 34, 35, 36, 43, 40, 41, 42, 47, 44, 45, 46, 48};

  p.change_type(PointSetTag::string_to_tag("Quad-P3-Equidist"), rotation_perm_code);

  // Check rotation permutation for p3 quad
  BOOST_CHECK(p.get().size() == 16);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p3quad_perm_rot[i]);
  }

  p.change_type(PointSetTag::string_to_tag("Quad-P4-Equidist"), rotation_perm_code);

  // Check rotation permutation for p4 quad
  BOOST_CHECK(p.get().size() == 25);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p4quad_perm_rot[i]);
  }

  p.change_type(PointSetTag::string_to_tag("Quad-P5-Equidist"), rotation_perm_code);

  // Check rotation permutation for p5 quad
  BOOST_CHECK(p.get().size() == 36);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p5quad_perm_rot[i]);
  }

  p.change_type(PointSetTag::string_to_tag("Quad-P6-Equidist"), rotation_perm_code);

  // Check rotation permutation for p6 quad
  BOOST_CHECK(p.get().size() == 49);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p6quad_perm_rot[i]);
  }

  /// --------------------------------------------------------------------------
  /// Compute quadrilateral permutation
  /// b) flip
  /// --------------------------------------------------------------------------

  Uint p3quad_perm_flip[16] = {1, 0, 3, 2, 5, 4, 11, 10, 9, 8, 7, 6, 13, 12, 15, 14};
  Uint p4quad_perm_flip[25] = {1, 0, 3, 2,  6,  5,  4,  15, 14, 13, 12, 11, 10,
                               9, 8, 7, 17, 16, 19, 18, 20, 23, 22, 21, 24};
  Uint p5quad_perm_flip[36] = {1,  0,  3,  2,  7,  6,  5,  4,  19, 18, 17, 16,
                               15, 14, 13, 12, 11, 10, 9,  8,  21, 20, 23, 22,
                               25, 24, 31, 30, 29, 28, 27, 26, 33, 32, 35, 34};
  Uint p6quad_perm_flip[49] = {1,  0,  3,  2,  8,  7,  6,  5,  4,  23, 22, 21, 20, 19, 18, 17, 16,
                               15, 14, 13, 12, 11, 10, 9,  25, 24, 27, 26, 30, 29, 28, 39, 38, 37,
                               36, 35, 34, 33, 32, 31, 41, 40, 43, 42, 44, 47, 46, 45, 48};

  p.change_type(PointSetTag::string_to_tag("Quad-P3-Equidist"), flip_perm_code);

  // Check rotation permutation for p3 quad
  BOOST_CHECK(p.get().size() == 16);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p3quad_perm_flip[i]);
  }

  p.change_type(PointSetTag::string_to_tag("Quad-P4-Equidist"), flip_perm_code);

  // Check rotation permutation for p4 quad
  BOOST_CHECK(p.get().size() == 25);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p4quad_perm_flip[i]);
  }

  p.change_type(PointSetTag::string_to_tag("Quad-P5-Equidist"), flip_perm_code);
  p.get().print();
  std::cout << std::endl;

  // Check rotation permutation for p5 quad
  BOOST_CHECK(p.get().size() == 36);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p5quad_perm_flip[i]);
  }

  p.change_type(PointSetTag::string_to_tag("Quad-P6-Equidist"), flip_perm_code);

  // Check rotation permutation for p6 quad
  BOOST_CHECK(p.get().size() == 49);

  for (Uint i = 0; i < p.get().size(); ++i)
  {
    BOOST_CHECK(p.get().vertex(i) == p6quad_perm_flip[i]);
  }

  /// --------------------------------------------------------------------------
  /// Match two edges (1d entities)
  /// --------------------------------------------------------------------------

  Uint p5line_array[6]     = {45, 85, 12, 11, 8, 19};
  Uint p5line_array_rot[6] = {85, 45, 19, 8, 11, 12};

  MeshEntity p5line(common::ArrayView<const Uint, _1D, Uint>(&p5line_array[0], 6), 0,
                    PointSetTag::string_to_tag("Line-P5-Equidist"));
  MeshEntity p5line_rot(common::ArrayView<const Uint, _1D, Uint>(&p5line_array_rot[0], 6), 0,
                        PointSetTag::string_to_tag("Line-P5-Equidist"));

  EntityDofRealign permutation;

  const bool p5lines_match = TP::entities_match(p5line, p5line, permutation);
  BOOST_CHECK(p5lines_match == true);

  const bool p5lines_match_rot = TP::entities_match_reverse(p5line, p5line_rot, permutation);
  BOOST_CHECK(p5lines_match_rot == true);

  /// --------------------------------------------------------------------------
  /// Match two p5 triangles without having to flip them
  /// --------------------------------------------------------------------------

  Uint p3triag_array[10]         = {100, 101, 102, 20, 21, 30, 31, 40, 41, 1};
  Uint p3triag_array_rot[10]     = {101, 102, 100, 30, 31, 40, 41, 20, 21, 1};
  Uint p3triag_array_rot_rot[10] = {102, 100, 101, 40, 41, 20, 21, 30, 31, 1};
  Uint p3triag_array_nomatch[10] = {102, 100, 101, 40, 51, 20, 21, 30, 31, 1};

  MeshEntity p3tr(common::ArrayView<const Uint, _1D, Uint>(&p3triag_array[0], 10), 0,
                  PointSetTag::string_to_tag("Triag-P3-Equidist"));
  MeshEntity p3tr_rot(common::ArrayView<const Uint, _1D, Uint>(&p3triag_array_rot[0], 10), 0,
                      PointSetTag::string_to_tag("Triag-P3-Equidist"));
  MeshEntity p3tr_rot_rot(common::ArrayView<const Uint, _1D, Uint>(&p3triag_array_rot_rot[0], 10),
                          0, PointSetTag::string_to_tag("Triag-P3-Equidist"));
  MeshEntity p3tr_nomatch(common::ArrayView<const Uint, _1D, Uint>(&p3triag_array_nomatch[0], 10),
                          0, PointSetTag::string_to_tag("Triag-P3-Equidist"));

  const bool p3triags_match = TP::entities_match(p3tr, p3tr, permutation);
  BOOST_CHECK(p3triags_match == true);

  const bool p3triags_match_rot = TP::entities_match(p3tr, p3tr_rot, permutation);
  BOOST_CHECK(p3triags_match_rot == true);

  const bool p3triags_match_rot_rot = TP::entities_match(p3tr, p3tr_rot_rot, permutation);
  BOOST_CHECK(p3triags_match_rot_rot == true);

  const bool p3triags_nomatch = TP::entities_match(p3tr, p3tr_nomatch, permutation);
  BOOST_CHECK(p3triags_nomatch == false);

  /// --------------------------------------------------------------------------
  /// Match two p5 triangles when a flip is required
  /// --------------------------------------------------------------------------

  Uint p3triag_array_f[10]         = {101, 100, 102, 21, 20, 41, 40, 31, 30, 1};
  Uint p3triag_array_f_rot[10]     = {102, 101, 100, 31, 30, 21, 20, 41, 40, 1};
  Uint p3triag_array_f_rot_rot[10] = {100, 102, 101, 41, 40, 31, 30, 21, 20, 1};

  MeshEntity p3tr_f(common::ArrayView<const Uint, _1D, Uint>(&p3triag_array_f[0], 10), 0,
                    PointSetTag::string_to_tag("Triag-P3-Equidist"));
  MeshEntity p3tr_f_rot(common::ArrayView<const Uint, _1D, Uint>(&p3triag_array_f_rot[0], 10), 0,
                        PointSetTag::string_to_tag("Triag-P3-Equidist"));
  MeshEntity p3tr_f_rot_rot(
      common::ArrayView<const Uint, _1D, Uint>(&p3triag_array_f_rot_rot[0], 10), 0,
      PointSetTag::string_to_tag("Triag-P3-Equidist"));

  // When using the 'reverse matching', an element should not be considered
  // identical to itself
  const bool p3triags_match_false_positive = TP::entities_match_reverse(p3tr, p3tr, permutation);
  BOOST_CHECK(p3triags_match_false_positive == false);

  const bool p3triags_match_f = TP::entities_match_reverse(p3tr, p3tr_f, permutation);
  BOOST_CHECK(p3triags_match_f == true);

  const bool p3triags_match_f_rot = TP::entities_match_reverse(p3tr, p3tr_f_rot, permutation);
  BOOST_CHECK(p3triags_match_f_rot == true);

  const bool p3triags_match_f_rot_rot =
      TP::entities_match_reverse(p3tr, p3tr_f_rot_rot, permutation);
  BOOST_CHECK(p3triags_match_f_rot_rot == true);

  /// --------------------------------------------------------------------------
  /// Match two p3 triangles which are faces of neighboring p3 tetrahedra
  /// --------------------------------------------------------------------------

  Uint p3tet_array[20] = {0,  1,  2,  3,  14, 15, 16,  17,  18,  19,
                          20, 21, 22, 23, 24, 25, 100, 101, 102, 103};
  MeshEntity p3tet(common::ArrayView<const Uint, _1D, Uint>(&p3tet_array[0], 20), 0,
                   PointSetTag::string_to_tag("Tetra-P3-Equidist"));

  std::cout << "P3 tetra:" << p3tet << std::endl;

  MeshEntity face102 = p3tet.sub_entity(_2D, 0);
  std::cout << "Face 102 = " << face102 << std::endl;

  MeshEntity edge10 = face102.sub_entity(_1D, 0);
  std::cout << "Edge 10 = " << edge10 << std::endl;
}
