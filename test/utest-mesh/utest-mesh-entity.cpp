/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mesh_entity_test
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iostream>

/// PDEKIT headers
#include "mesh/MeshEntity.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(triangular_mesh_entity_utest)
{
  const Uint p1_triag_data[3] = {10, 11, 12};
  MeshEntity p1_tri(common::ArrayView<const Uint, _1D, Uint>(&p1_triag_data[0], 3), 0,
                    PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist));
  p1_tri.print_reference_topology();

  BOOST_CHECK_EQUAL(p1_tri.nb_vert(), 3u);

  for (Uint v = 0; v < p1_tri.nb_vert(); ++v)
  {
    BOOST_CHECK_EQUAL(p1_tri.vertex(v), p1_triag_data[v]);
  }

  // Triangle has 3 edges
  BOOST_CHECK_EQUAL(p1_tri.nb_sub_elements(_1D), 3u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(tetrahedral_mesh_entity_utest)
{
  const Uint p1_tet_data[4] = {10, 11, 12, 13};
  MeshEntity p1_tet(common::ArrayView<const Uint, _1D, Uint>(&p1_tet_data[0], 4), 0,
                    PointSetTag(ElemShape::Tetra, P1, PointSetID::Equidist));
  p1_tet.print_reference_topology();

  BOOST_CHECK_EQUAL(p1_tet.nb_vert(), 4u);

  for (Uint v = 0; v < p1_tet.nb_vert(); ++v)
  {
    BOOST_CHECK_EQUAL(p1_tet.vertex(v), p1_tet_data[v]);
  }

  // Check faces - tetrahedra has 4 faces
  BOOST_CHECK_EQUAL(p1_tet.nb_sub_elements(_2D), 4u);
  MeshEntity tmp = p1_tet;
  tmp.local_transform(_2D, 0);
  BOOST_CHECK_EQUAL(tmp.vertex(0), 10u);
  BOOST_CHECK_EQUAL(tmp.vertex(1), 12u);
  BOOST_CHECK_EQUAL(tmp.vertex(2), 11u);

  tmp = p1_tet;
  tmp.local_transform(_2D, 1);
  BOOST_CHECK_EQUAL(tmp.vertex(0), 10u);
  BOOST_CHECK_EQUAL(tmp.vertex(1), 11u);
  BOOST_CHECK_EQUAL(tmp.vertex(2), 13u);

  tmp = p1_tet;
  tmp.local_transform(_2D, 2);
  BOOST_CHECK_EQUAL(tmp.vertex(0), 10u);
  BOOST_CHECK_EQUAL(tmp.vertex(1), 13u);
  BOOST_CHECK_EQUAL(tmp.vertex(2), 12u);

  tmp = p1_tet;
  tmp.local_transform(_2D, 3);
  BOOST_CHECK_EQUAL(tmp.vertex(0), 13u);
  BOOST_CHECK_EQUAL(tmp.vertex(1), 11u);
  BOOST_CHECK_EQUAL(tmp.vertex(2), 12u);

  // Now take the last face (tmp) and cast it to edges
  std::cout << "Last face = " << tmp << std::endl;
  MeshEntity tmp2 = tmp;
  tmp2.local_transform(_1D, 0);
  std::cout << tmp2 << std::endl;

  // Tetrahedra has 6 edges
  BOOST_CHECK_EQUAL(p1_tet.nb_sub_elements(_1D), 6u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mesh_entity_utest)
{
  /// --------------------------------------------------------------------------
  /// Create 2 triangles (2d entities)
  /// --------------------------------------------------------------------------

  Uint p2triag_data_a[6] = {10, 11, 12, 100, 101, 102};
  Uint p2triag_data_b[6] = {13, 14, 15, 200, 201, 202};

  MeshEntity p2tri_a(common::ArrayView<const Uint, _1D, Uint>(&p2triag_data_a[0], 6), 0,
                     PointSetTag::string_to_tag("Triag-P2-Equidist"));
  MeshEntity p2tri_b(common::ArrayView<const Uint, _1D, Uint>(&p2triag_data_b[0], 6), 1,
                     PointSetTag::string_to_tag("Triag-P2-Equidist"));

  std::cout << "First p2 triangle = " << p2tri_a << std::endl;
  for (Uint i = 0; i < p2tri_a.nb_vert(); ++i)
  {
    std::cout << "  Vertex " << p2tri_a.vertex(i) << " is p1 vertex " << p2tri_a.vert_is_p1(i)
              << std::endl;
  }

  std::cout << "Second p2 triangle = " << p2tri_b << std::endl;
  for (Uint i = 0; i < p2tri_b.nb_vert(); ++i)
  {
    std::cout << "  Vertex " << p2tri_b.vertex(i) << " is p1 vertex " << p2tri_b.vert_is_p1(i)
              << std::endl;
  }

  MeshEntity p2line_a = p2tri_a;
  p2line_a.local_transform(_1D, 0);

  std::cout << "Edges incident to first edge of p2 triag:" << std::endl;
  std::cout << p2line_a.local_incident_entities(_1D) << std::endl;

  std::cout << "First edge of the first p2 triangle:" << std::endl;
  std::cout << p2line_a << std::endl;
  for (Uint i = 0; i < p2line_a.nb_vert(); ++i)
  {
    std::cout << "  Vertex " << p2line_a.vertex(i) << " is p1 vertex " << p2line_a.vert_is_p1(i)
              << std::endl;
  }

  /// --------------------------------------------------------------------------
  /// Generate a sample p3 triangle
  /// --------------------------------------------------------------------------

  Uint p3triag_array[10] = {100, 101, 102, 20, 21, 30, 31, 40, 41, 1};
  MeshEntity p3triag(common::ArrayView<const Uint, _1D, Uint>(&p3triag_array[0], 10), 0,
                     PointSetTag::string_to_tag("Triag-P3-Equidist"));

  p3triag.print_reference_topology();

  /// --------------------------------------------------------------------------
  /// Generate a sample p5 tetrahedra
  /// --------------------------------------------------------------------------

  Uint p5tet_array[56] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
                          14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                          28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
                          42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55};

  Uint face_021_check[21] = {0, 2, 1, 15, 14, 13, 12, 11, 10, 9, 8,
                             7, 6, 5, 4,  28, 29, 30, 31, 32, 33};
  Uint face_013_check[21] = {0,  1,  3,  4,  5,  6,  7,  27, 26, 25, 24,
                             16, 17, 18, 19, 34, 35, 36, 37, 38, 39};
  Uint face_032_check[21] = {0,  3,  2,  19, 18, 17, 16, 20, 21, 22, 23,
                             12, 13, 14, 15, 40, 41, 42, 43, 44, 45};
  Uint face_312_check[21] = {3,  1,  2,  24, 25, 26, 27, 8,  9,  10, 11,
                             23, 22, 21, 20, 46, 47, 48, 49, 50, 51};

  MeshEntity p5tetra(common::ArrayView<const Uint, _1D, Uint>(&p5tet_array[0], 56), 0,
                     PointSetTag::string_to_tag("Tetra-P5-Equidist"));

  // Get the first face (1-0-2) of p5 tetra, check that its vertices are
  // correct
  MeshEntity tet_face_021 = p5tetra.sub_entity(_2D, 0);
  for (Uint i = 0; i < 21; ++i)
  {
    BOOST_CHECK(tet_face_021.vertex(i) == face_021_check[i]);
  }

  // Get the second face (0-1-3) of p5 tetra, check that its vertices are
  // correct
  MeshEntity tet_face_013 = p5tetra.sub_entity(_2D, 1);
  for (Uint i = 0; i < 21; ++i)
  {
    BOOST_CHECK(tet_face_013.vertex(i) == face_013_check[i]);
  }
  std::cout << std::endl;

  // Get the second face (0-3-2) of p5 tetra, check that its vertices are
  // correct
  MeshEntity tet_face_032 = p5tetra.sub_entity(_2D, 2);
  for (Uint i = 0; i < 21; ++i)
  {
    BOOST_CHECK(tet_face_032.vertex(i) == face_032_check[i]);
  }

  // Get the second face (3-1-2) of p5 tetra, check that its vertices are
  // correct
  MeshEntity tet_face_312 = p5tetra.sub_entity(_2D, 3);
  for (Uint i = 0; i < 21; ++i)
  {
    BOOST_CHECK(tet_face_312.vertex(i) == face_312_check[i]);
  }

  // Check for edges. Not that it should be possible to obtain the edges just
  // by transforming the face entities into their own sub-entities. We are
  // going to test it here.
  MeshEntity edge02_of_face021 = tet_face_021.sub_entity(_1D, 0);
  std::cout << "Edge 1-0 of face 1-0-2:" << std::endl;
  std::cout << edge02_of_face021 << std::endl;

  p5tetra.print_reference_topology();
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(standard_region_utest)
{
  mesh::StdRegion std_reg;
  std_reg.change_type(mesh::PointSetTag(ElemShape::Triag, P2, PointSetID::Equidist));

  std::cout << "Triangle normals: " << std::endl;
  for (Uint f = 0; f < std_reg.get().nb_entities(_1D); ++f)
  {

    std::cout << std_reg.get().facet_normal(f) << std::endl;
  }

  std_reg.change_type(mesh::PointSetTag(ElemShape::Tetra, P5, PointSetID::Equidist));

  std::cout << "Tetra normals: " << std::endl;
  for (Uint f = 0; f < std_reg.get().nb_entities(_2D); ++f)
  {

    std::cout << std_reg.get().facet_normal(f) << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mesh_entity_vertex_iterator_utest)
{
  const Uint p1_tet_data[4] = {10, 11, 12, 13};
  MeshEntity p1_tet(common::ArrayView<const Uint, _1D, Uint>(&p1_tet_data[0], 4), 0,
                    PointSetTag(ElemShape::Tetra, P1, PointSetID::Equidist));
  p1_tet.print_reference_topology();

  BOOST_CHECK_EQUAL(p1_tet.nb_vert(), 4u);

  for (Uint v = 0; v < p1_tet.nb_vert(); ++v)
  {
    BOOST_CHECK_EQUAL(p1_tet.vertex(v), p1_tet_data[v]);
  }

  Uint v = 0;
  for (MeshEntity::const_iterator vert_it = p1_tet.cbegin(); vert_it != p1_tet.cend();
       ++vert_it, ++v)
  {
    BOOST_CHECK_EQUAL(p1_tet.vertex(v), p1_tet_data[v]);
    BOOST_CHECK_EQUAL(p1_tet.vertex(v), p1_tet_data[v]);
  }

  // Check faces - tetrahedra has 4 faces
  BOOST_CHECK_EQUAL(p1_tet.nb_sub_elements(_2D), 4u);
  MeshEntity tmp = p1_tet;
  tmp.local_transform(_2D, 0);
  const Uint face0_check[3] = {10u, 12u, 11u};

  v = 0;
  for (MeshEntity::const_iterator vert_it = tmp.cbegin(); vert_it != tmp.cend(); ++vert_it, ++v)
  {
    BOOST_CHECK_EQUAL(tmp.vertex(v), face0_check[v]);
    BOOST_CHECK_EQUAL(tmp.vertex(v), face0_check[v]);
  }

  tmp = p1_tet;
  tmp.local_transform(_2D, 1);
  const Uint face1_check[3] = {10u, 11u, 13u};

  v = 0;
  for (MeshEntity::const_iterator vert_it = tmp.cbegin(); vert_it != tmp.cend(); ++vert_it, ++v)
  {
    BOOST_CHECK_EQUAL(tmp.vertex(v), face1_check[v]);
    BOOST_CHECK_EQUAL(tmp.vertex(v), face1_check[v]);
  }

  tmp = p1_tet;
  tmp.local_transform(_2D, 2);
  const Uint face2_check[3] = {10u, 13u, 12u};

  v = 0;
  for (MeshEntity::const_iterator vert_it = tmp.cbegin(); vert_it != tmp.cend(); ++vert_it, ++v)
  {
    BOOST_CHECK_EQUAL(tmp.vertex(v), face2_check[v]);
    BOOST_CHECK_EQUAL(tmp.vertex(v), face2_check[v]);
  }

  tmp = p1_tet;
  tmp.local_transform(_2D, 3);
  const Uint face3_check[3] = {13u, 11u, 12u};

  v = 0;
  for (MeshEntity::const_iterator vert_it = tmp.cbegin(); vert_it != tmp.cend(); ++vert_it, ++v)
  {
    BOOST_CHECK_EQUAL(tmp.vertex(v), face3_check[v]);
    BOOST_CHECK_EQUAL(tmp.vertex(v), face3_check[v]);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mesh_entity_vertex_hash)
{
  const Uint p1_tet_data[4] = {10, 11, 12, 13};
  MeshEntity p1_tet(common::ArrayView<const Uint, _1D, Uint>(&p1_tet_data[0], 4), 0,
                    PointSetTag(ElemShape::Tetra, P1, PointSetID::Equidist));

  const Uint vert_hash = MeshEntityVertHasher::hash(p1_tet);
  BOOST_CHECK_EQUAL(vert_hash, 165901833u);
}

// ----------------------------------------------------------------------------
