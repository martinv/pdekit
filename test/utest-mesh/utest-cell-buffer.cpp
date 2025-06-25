/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cell_buffer_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers
#include "mesh/CellBuffer.hpp"

using namespace pdekit;

BOOST_AUTO_TEST_CASE(cell_buffer_utest)
{
  std::vector<Real> one_cell_coords;
  mesh::CellBuffer<_2D, _2D> buffer;
  buffer.reserve(22, 6);

  std::vector<Uint> cell_vertices = {1, 2, 3};
  one_cell_coords.resize(cell_vertices.size() * _2D);
  buffer.push_back_cell(0, mesh::PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist),
                        cell_vertices, 0u, one_cell_coords);

  cell_vertices = {4, 5, 6};
  one_cell_coords.resize(cell_vertices.size() * _2D);
  buffer.push_back_cell(1, mesh::PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist),
                        cell_vertices, 0u, one_cell_coords);

  cell_vertices = {7, 8, 9};
  one_cell_coords.resize(cell_vertices.size() * _2D);
  buffer.push_back_cell(2, mesh::PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist),
                        cell_vertices, 0u, one_cell_coords);

  cell_vertices = {10, 11, 12};
  one_cell_coords.resize(cell_vertices.size() * _2D);
  buffer.push_back_cell(3, mesh::PointSetTag(ElemShape::Triag, P1, PointSetID::Warpblend),
                        cell_vertices, 0u, one_cell_coords);

  cell_vertices = {13, 14, 15, 16};
  one_cell_coords.resize(cell_vertices.size() * _2D);
  buffer.push_back_cell(4, mesh::PointSetTag(ElemShape::Quad, P1, PointSetID::Equidist),
                        cell_vertices, 0u, one_cell_coords);

  cell_vertices = {17, 18, 19, 20, 21, 22};
  one_cell_coords.resize(cell_vertices.size() * _2D);
  buffer.push_back_cell(5, mesh::PointSetTag(ElemShape::Triag, P2, PointSetID::Warpblend),
                        cell_vertices, 0u, one_cell_coords);

  BOOST_CHECK_EQUAL(buffer.nb_active_cells(), 6U);

  const mesh::MeshEntity triag = buffer.active_cell(mesh::ActiveIdx(1));
  std::cout << "triag = " << triag << std::endl;

  BOOST_CHECK_EQUAL(triag.nb_vert(), 3U);
  BOOST_CHECK_EQUAL(triag.vertex(0), 4U);
  BOOST_CHECK_EQUAL(triag.vertex(1), 5U);
  BOOST_CHECK_EQUAL(triag.vertex(2), 6U);

  mesh::MeshEntity triag2;
  buffer.fill_cell(mesh::ActiveIdx(1), triag2);
  BOOST_CHECK_EQUAL(triag2.nb_vert(), 3U);
  BOOST_CHECK_EQUAL(triag2.vertex(0), 4U);
  BOOST_CHECK_EQUAL(triag2.vertex(1), 5U);
  BOOST_CHECK_EQUAL(triag2.vertex(2), 6U);

  const mesh::MeshEntity quad = buffer.active_cell(mesh::ActiveIdx(4));
  std::cout << "quad = " << quad << std::endl;
  BOOST_CHECK_EQUAL(quad.nb_vert(), 4U);
  BOOST_CHECK_EQUAL(quad.vertex(0), 13U);
  BOOST_CHECK_EQUAL(quad.vertex(1), 14U);
  BOOST_CHECK_EQUAL(quad.vertex(2), 15U);
  BOOST_CHECK_EQUAL(quad.vertex(3), 16U);

  mesh::MeshEntity quad2;
  buffer.fill_cell(mesh::ActiveIdx(4), quad2);
  BOOST_CHECK_EQUAL(quad2.nb_vert(), 4U);
  BOOST_CHECK_EQUAL(quad2.vertex(0), 13U);
  BOOST_CHECK_EQUAL(quad2.vertex(1), 14U);
  BOOST_CHECK_EQUAL(quad2.vertex(2), 15U);
  BOOST_CHECK_EQUAL(quad2.vertex(3), 16U);

  std::cout << "Printing buffer " << std::endl;
  buffer.print(true);
  std::cout << "Buffer printed" << std::endl;

  for (Uint c = 0; c < buffer.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = buffer.active_cell(mesh::ActiveIdx(c));
    std::cout << "--> " << cell << std::endl;
  }

  std::cout << std::endl;

  for (auto it = std::begin(buffer); it != std::end(buffer); ++it)
  {
    const mesh::MeshEntity cell = it->mesh_entity();
    std::cout << "==> " << cell << std::endl;
  }

  std::cout << std::endl;

  for (const auto &cell : buffer)
  {
    std::cout << "===> " << cell.mesh_entity() << std::endl;
  }

  mesh::CellBuffer<_2D, _2D> buffer2;
  buffer2 = buffer;
  buffer2.print(true);
}

// ----------------------------------------------------------------------------
