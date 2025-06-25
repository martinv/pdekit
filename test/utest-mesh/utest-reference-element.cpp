/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE reference_element_test
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iostream>

/// PDEKIT headers
#include "mesh/std_region/StdRegion.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

BOOST_AUTO_TEST_CASE(reference_element_utest)
{
  StdRegion elem;

  // Switch to undefined reference element and check that its a valid element
  elem.change_type(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined));
  // elem.change_type(Undefined,_0D,Undefined); // This should also
  // compile - change_type
  // accepts template parameter pack as argument

  BOOST_CHECK_EQUAL(elem.get().nb_nodes(), 0U);
  BOOST_CHECK_EQUAL(elem.get().topo_dim(), _0D);

  // Switch to a 'normal' element

  PointSetTag elem_type_tag(ElemShape::Quad, _2D, PointSetID::Warpblend);
  elem.change_type(elem_type_tag);

  const StdRegion::value_type::coordinates_type &coordinates = elem.get().coordinates();

  std::cout << "Coordinates of 2D Warpblend Quadrilateral = " << std::endl
            << coordinates << std::endl;

  for (Uint f = 0; f < elem.get().nb_entities(_1D); ++f)
  {
    std::cout << "Face " << f << ":" << std::endl;
    std::cout << *elem.get().elem_entity(_1D, f) << std::endl;
  }

  // Perform some checks
  BOOST_CHECK_EQUAL(elem.get().pt_set_id(), elem_type_tag);
  BOOST_CHECK_EQUAL(elem.get().nb_nodes(), 9U);
  BOOST_CHECK_EQUAL(elem.get().topo_dim(), _2D);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(reference_line_topology_utest)
{
  StdRegion ref_line;

  // ----------------------------
  // REFERENCE LINE, ORDER 1 TO 5
  // ----------------------------

  for (Uint order = 1; order <= 5; ++order)
  {
    ref_line.change_type(PointSetTag(ElemShape::Line, order, PointSetID::Equidist));

    std::shared_ptr<StdRegionEntity const> line_as_entity = ref_line.get().elem_entity(_1D, 0);
    BOOST_CHECK_EQUAL((*line_as_entity).nb_vert(), order + 1);

    for (Uint v = 0; v < (*line_as_entity).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*line_as_entity).vertex(v), v);
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(reference_triag_topology_utest)
{
  StdRegion ref_triag;

  // --------------------------------
  // REFERENCE TRIANGLE, ORDER 1 TO 5
  // --------------------------------

  for (Uint order = 1; order <= 5; ++order)
  {
    ref_triag.change_type(PointSetTag(ElemShape::Triag, order, PointSetID::Equidist));

    BOOST_CHECK_EQUAL(ref_triag.get().nb_entities(_1D), 3u);

    Uint internal_vertex_in_edge = 3;

    for (Uint e = 0; e < ref_triag.get().nb_entities(_1D); ++e)
    {
      std::shared_ptr<StdRegionEntity const> tri_edge = ref_triag.get().elem_entity(_1D, e);
      BOOST_CHECK_EQUAL((*tri_edge).nb_vert(), order + 1);

      // Check first vertex of edge (first extremity)
      BOOST_CHECK_EQUAL((*tri_edge).vertex(0), e);
      // Check second vertex of edge (second extremity)
      BOOST_CHECK_EQUAL((*tri_edge).vertex(1), (e + 1) % 3);

      for (Uint v = 2; v < (*tri_edge).nb_vert(); ++v)
      {
        BOOST_CHECK_EQUAL((*tri_edge).vertex(v), internal_vertex_in_edge);
        internal_vertex_in_edge++;
      }
    }
  } // Loop over polynomial orders
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(reference_quad_topology_utest)
{
  StdRegion ref_quad;

  // -------------------------------------
  // REFERENCE QUADRILATERAL, ORDER 1 TO 5
  // -------------------------------------

  for (Uint order = 1; order <= 5; ++order)
  {
    ref_quad.change_type(PointSetTag(ElemShape::Quad, order, PointSetID::Equidist));

    BOOST_CHECK_EQUAL(ref_quad.get().nb_entities(_1D), 4u);

    Uint internal_vertex_in_edge = 4;

    for (Uint e = 0; e < ref_quad.get().nb_entities(_1D); ++e)
    {
      std::shared_ptr<StdRegionEntity const> quad_edge = ref_quad.get().elem_entity(_1D, e);
      BOOST_CHECK_EQUAL((*quad_edge).nb_vert(), order + 1);

      // Check first vertex of edge (first extremity)
      BOOST_CHECK_EQUAL((*quad_edge).vertex(0), e);
      // Check second vertex of edge (second extremity)
      BOOST_CHECK_EQUAL((*quad_edge).vertex(1), (e + 1) % 4);

      for (Uint v = 2; v < (*quad_edge).nb_vert(); ++v)
      {
        BOOST_CHECK_EQUAL((*quad_edge).vertex(v), internal_vertex_in_edge);
        internal_vertex_in_edge++;
      }
    }
  } // Loop over polynomial orders
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(reference_tetra_topology_utest)
{
  StdRegion ref_tetra;

  // -----------------
  // P1 EQUIDIST TETRA
  // -----------------
  ref_tetra.change_type(PointSetTag(ElemShape::Tetra, P1, PointSetID::Equidist));

  // CHECK ELEMENT FACES
  const std::vector<std::vector<Uint>> p1_tet_faces = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {3, 1, 2}};

  BOOST_CHECK_EQUAL(ref_tetra.get().nb_entities(_2D), 4u);

  for (Uint f = 0; f < ref_tetra.get().nb_entities(_2D); ++f)
  {
    std::shared_ptr<StdRegionEntity const> p1_tet_face = ref_tetra.get().elem_entity(_2D, f);
    BOOST_CHECK_EQUAL((*p1_tet_face).nb_vert(), 3u);

    for (Uint v = 0; v < (*p1_tet_face).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p1_tet_face).vertex(v), p1_tet_faces[f][v]);
    }
  }

  // CHECK ELEMENT EDGES
  const std::vector<std::vector<Uint>> p1_tet_edges = {{0, 1}, {1, 2}, {2, 0},
                                                       {3, 0}, {3, 2}, {3, 1}};

  BOOST_CHECK_EQUAL(ref_tetra.get().nb_entities(_1D), 6u);

  for (Uint e = 0; e < ref_tetra.get().nb_entities(_1D); ++e)
  {
    std::shared_ptr<StdRegionEntity const> p1_tet_edge = ref_tetra.get().elem_entity(_1D, e);
    BOOST_CHECK_EQUAL((*p1_tet_edge).nb_vert(), 2u);

    for (Uint v = 0; v < (*p1_tet_edge).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p1_tet_edge).vertex(v), p1_tet_edges[e][v]);
    }
  }

  // -----------------
  // P2 EQUIDIST TETRA
  // -----------------
  ref_tetra.change_type(PointSetTag(ElemShape::Tetra, P2, PointSetID::Equidist));

  // CHECK ELEMENT FACES
  const std::vector<std::vector<Uint>> p2_tet_faces = {
      {0, 2, 1, 6, 5, 4}, {0, 1, 3, 4, 9, 7}, {0, 3, 2, 7, 8, 6}, {3, 1, 2, 9, 5, 8}};

  BOOST_CHECK_EQUAL(ref_tetra.get().nb_entities(_2D), 4u);

  for (Uint f = 0; f < ref_tetra.get().nb_entities(_2D); ++f)
  {
    std::shared_ptr<StdRegionEntity const> p2_tet_face = ref_tetra.get().elem_entity(_2D, f);
    BOOST_CHECK_EQUAL((*p2_tet_face).nb_vert(), 6u);

    for (Uint v = 0; v < (*p2_tet_face).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p2_tet_face).vertex(v), p2_tet_faces[f][v]);
    }
  }

  // CHECK ELEMENT EDGES
  const std::vector<std::vector<Uint>> p2_tet_edges = {{0, 1, 4}, {1, 2, 5}, {2, 0, 6},
                                                       {3, 0, 7}, {3, 2, 8}, {3, 1, 9}};

  BOOST_CHECK_EQUAL(ref_tetra.get().nb_entities(_1D), 6u);

  for (Uint e = 0; e < ref_tetra.get().nb_entities(_1D); ++e)
  {
    std::shared_ptr<StdRegionEntity const> p2_tet_edge = ref_tetra.get().elem_entity(_1D, e);
    BOOST_CHECK_EQUAL((*p2_tet_edge).nb_vert(), 3u);

    for (Uint v = 0; v < (*p2_tet_edge).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p2_tet_edge).vertex(v), p2_tet_edges[e][v]);
    }
  }

  // -----------------
  // P3 EQUIDIST TETRA
  // -----------------
  ref_tetra.change_type(PointSetTag(ElemShape::Tetra, P3, PointSetID::Equidist));

  // CHECK ELEMENT FACES
  const std::vector<std::vector<Uint>> p3_tet_faces = {{0, 2, 1, 9, 8, 7, 6, 5, 4, 16},
                                                       {0, 1, 3, 4, 5, 15, 14, 10, 11, 17},
                                                       {0, 3, 2, 11, 10, 12, 13, 8, 9, 18},
                                                       {3, 1, 2, 14, 15, 6, 7, 13, 12, 19}};

  BOOST_CHECK_EQUAL(ref_tetra.get().nb_entities(_2D), 4u);

  for (Uint f = 0; f < ref_tetra.get().nb_entities(_2D); ++f)
  {
    std::shared_ptr<StdRegionEntity const> p3_tet_face = ref_tetra.get().elem_entity(_2D, f);
    BOOST_CHECK_EQUAL((*p3_tet_face).nb_vert(), 10u);

    for (Uint v = 0; v < (*p3_tet_face).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p3_tet_face).vertex(v), p3_tet_faces[f][v]);
    }
  }

  // CHECK ELEMENT EDGES
  const std::vector<std::vector<Uint>> p3_tet_edges = {
      {0, 1, 4, 5}, {1, 2, 6, 7}, {2, 0, 8, 9}, {3, 0, 10, 11}, {3, 2, 12, 13}, {3, 1, 14, 15}};

  BOOST_CHECK_EQUAL(ref_tetra.get().nb_entities(_1D), 6u);

  for (Uint e = 0; e < ref_tetra.get().nb_entities(_1D); ++e)
  {
    std::shared_ptr<StdRegionEntity const> p3_tet_edge = ref_tetra.get().elem_entity(_1D, e);
    BOOST_CHECK_EQUAL((*p3_tet_edge).nb_vert(), 4u);

    for (Uint v = 0; v < (*p3_tet_edge).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p3_tet_edge).vertex(v), p3_tet_edges[e][v]);
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(reference_hexa_topology_utest)
{
  std::cout << "********************************************" << std::endl;
  std::cout << "            HEXA TOPOLOGY UTEST             " << std::endl;
  std::cout << "********************************************" << std::endl;
  StdRegion ref_hexa;

  // -----------------
  // P1 EQUIDIST HEXA
  // -----------------
  ref_hexa.change_type(PointSetTag(ElemShape::Hexa, P1, PointSetID::Equidist));

  // CHECK ELEMENT FACES
  const std::vector<std::vector<Uint>> p1_hex_faces = {{0, 3, 2, 1}, {0, 1, 5, 4}, {0, 4, 7, 3},
                                                       {1, 2, 6, 5}, {2, 3, 7, 6}, {4, 5, 6, 7}};

  BOOST_CHECK_EQUAL(ref_hexa.get().nb_entities(_2D), 6u);
  for (Uint f = 0; f < ref_hexa.get().nb_entities(_2D); ++f)
  {
    std::shared_ptr<StdRegionEntity const> p1_hex_face = ref_hexa.get().elem_entity(_2D, f);
    BOOST_CHECK_EQUAL((*p1_hex_face).nb_vert(), 4u);

    for (Uint v = 0; v < (*p1_hex_face).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p1_hex_face).vertex(v), p1_hex_faces[f][v]);
    }
  }

  // CHECK ELEMENT EDGES
  const std::vector<std::vector<Uint>> p1_hex_edges = {{0, 1}, {0, 3}, {0, 4}, {1, 2},
                                                       {1, 5}, {2, 3}, {2, 6}, {3, 7},
                                                       {4, 5}, {4, 7}, {5, 6}, {6, 7}};

  BOOST_CHECK_EQUAL(ref_hexa.get().nb_entities(_1D), 12u);

  for (Uint e = 0; e < ref_hexa.get().nb_entities(_1D); ++e)
  {
    std::shared_ptr<StdRegionEntity const> p1_hex_edge = ref_hexa.get().elem_entity(_1D, e);
    BOOST_CHECK_EQUAL((*p1_hex_edge).nb_vert(), 2u);

    for (Uint v = 0; v < (*p1_hex_edge).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p1_hex_edge).vertex(v), p1_hex_edges[e][v]);
    }
  }

  // -----------------
  // P2 EQUIDIST HEXA
  // -----------------
  ref_hexa.change_type(PointSetTag(ElemShape::Hexa, P2, PointSetID::Equidist));

  // CHECK ELEMENT FACES
  const std::vector<std::vector<Uint>> p2_hex_faces = {
      {0, 3, 2, 1, 9, 13, 11, 8, 20},   {0, 1, 5, 4, 8, 12, 16, 10, 21},
      {0, 4, 7, 3, 10, 17, 15, 9, 22},  {1, 2, 6, 5, 11, 14, 18, 12, 23},
      {2, 3, 7, 6, 13, 15, 19, 14, 24}, {4, 5, 6, 7, 16, 18, 19, 17, 25}};

  BOOST_CHECK_EQUAL(ref_hexa.get().nb_entities(_2D), 6u);
  for (Uint f = 0; f < ref_hexa.get().nb_entities(_2D); ++f)
  {
    std::shared_ptr<StdRegionEntity const> p2_hex_face = ref_hexa.get().elem_entity(_2D, f);
    BOOST_CHECK_EQUAL((*p2_hex_face).nb_vert(), 9u);

    for (Uint v = 0; v < (*p2_hex_face).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p2_hex_face).vertex(v), p2_hex_faces[f][v]);
    }
  }

  // CHECK ELEMENT EDGES
  const std::vector<std::vector<Uint>> p2_hex_edges = {
      {0, 1, 8},  {0, 3, 9},  {0, 4, 10}, {1, 2, 11}, {1, 5, 12}, {2, 3, 13},
      {2, 6, 14}, {3, 7, 15}, {4, 5, 16}, {4, 7, 17}, {5, 6, 18}, {6, 7, 19}};

  BOOST_CHECK_EQUAL(ref_hexa.get().nb_entities(_1D), 12u);

  for (Uint e = 0; e < ref_hexa.get().nb_entities(_1D); ++e)
  {
    std::shared_ptr<StdRegionEntity const> p2_hex_edge = ref_hexa.get().elem_entity(_1D, e);
    BOOST_CHECK_EQUAL((*p2_hex_edge).nb_vert(), 3u);

    for (Uint v = 0; v < (*p2_hex_edge).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p2_hex_edge).vertex(v), p2_hex_edges[e][v]);
    }
  }

  // -----------------
  // P3 EQUIDIST HEXA
  // -----------------
  ref_hexa.change_type(PointSetTag(ElemShape::Hexa, P3, PointSetID::Equidist));

  // CHECK ELEMENT FACES
  const std::vector<std::vector<Uint>> p3_hex_faces = {
      {0, 3, 2, 1, 10, 11, 19, 18, 15, 14, 9, 8, 32, 33, 34, 35},
      {0, 1, 5, 4, 8, 9, 16, 17, 25, 24, 13, 12, 36, 37, 38, 39},
      {0, 4, 7, 3, 12, 13, 26, 27, 23, 22, 11, 10, 40, 41, 42, 43},
      {1, 2, 6, 5, 14, 15, 20, 21, 29, 28, 17, 16, 44, 45, 46, 47},
      {2, 3, 7, 6, 18, 19, 22, 23, 31, 30, 21, 20, 48, 49, 50, 51},
      {4, 5, 6, 7, 24, 25, 28, 29, 30, 31, 27, 26, 52, 53, 54, 55}};

  BOOST_CHECK_EQUAL(ref_hexa.get().nb_entities(_2D), 6u);
  for (Uint f = 0; f < ref_hexa.get().nb_entities(_2D); ++f)
  {
    std::shared_ptr<StdRegionEntity const> p3_hex_face = ref_hexa.get().elem_entity(_2D, f);
    BOOST_CHECK_EQUAL((*p3_hex_face).nb_vert(), 16u);

    for (Uint v = 0; v < (*p3_hex_face).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p3_hex_face).vertex(v), p3_hex_faces[f][v]);
    }
  }

  // CHECK ELEMENT EDGES
  const std::vector<std::vector<Uint>> p3_hex_edges = {
      {0, 1, 8, 9},   {0, 3, 10, 11}, {0, 4, 12, 13}, {1, 2, 14, 15},
      {1, 5, 16, 17}, {2, 3, 18, 19}, {2, 6, 20, 21}, {3, 7, 22, 23},
      {4, 5, 24, 25}, {4, 7, 26, 27}, {5, 6, 28, 29}, {6, 7, 30, 31}};

  BOOST_CHECK_EQUAL(ref_hexa.get().nb_entities(_1D), 12u);

  for (Uint e = 0; e < ref_hexa.get().nb_entities(_1D); ++e)
  {
    std::shared_ptr<StdRegionEntity const> p3_hex_edge = ref_hexa.get().elem_entity(_1D, e);
    BOOST_CHECK_EQUAL((*p3_hex_edge).nb_vert(), 4u);

    for (Uint v = 0; v < (*p3_hex_edge).nb_vert(); ++v)
    {
      BOOST_CHECK_EQUAL((*p3_hex_edge).vertex(v), p3_hex_edges[e][v]);
    }
  }

  // std::cout << "Coordinates:" << std::endl;
  // std::cout << *ref_hexa.get().coordinates() << std::endl;
}
