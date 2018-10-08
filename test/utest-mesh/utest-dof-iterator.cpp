/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE dof_iterator_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers
#include "common/Constants.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

typedef Tria<Cart2D> MeshType;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dof_iterator_validity_utest)
{
  MeshType mesh2d("mesh2d");

  const std::string infilename = "rectangle_mixed_elem_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, mesh2d, "geo_dofs");

  const MeshType::dof_storage_type &geo_cell_dofs = *mesh2d.dof_storage("geo_dofs");

  typedef typename result_of::dof_map_t<Cart2D>::const_dof_iterator dof_iterator_type;

  dof_iterator_type it1 = geo_cell_dofs.begin();
  dof_iterator_type it2 = geo_cell_dofs.begin();

  bool result = (it1 == it2);
  // std::cout << "Result = " << result << " [should be true]" << std::endl;
  BOOST_CHECK_EQUAL(result, true);

  ++it2;
  result = (it1 == it2);
  BOOST_CHECK_EQUAL(result, false);
  // std::cout << "Result = " << result << " [should be false]" << std::endl;

  BOOST_CHECK_EQUAL(dof_iterator_type::distance(it1, it2), 1);

  --it2;

  const MeshEntity cell1 = it1->mesh_entity();
  const MeshEntity cell2 = (*(it2.operator->())).mesh_entity();

  std::cout << "Cell 1 = " << cell1 << std::endl;
  std::cout << "Cell 2 = " << cell2 << std::endl;
  std::cout << "Element types:" << std::endl;

  for (Uint i = 0; i < geo_cell_dofs.all_active_dof_groups().size(); ++i)
  {
    const mesh::MeshEntity cell = (*geo_cell_dofs.all_active_dof_groups()[i].begin()).mesh_entity();

    std::cout << cell.pt_set_id().as_string() << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dof_iterator_loop_utest)
{
  // std::cout << std::endl << "************ RUNNING CELL ITERATOR TEST
  // ************" << std::endl;

  MeshType mesh2d("mesh2d");

  // Mesh root("root");

  const std::string infilename = "rectangle_mixed_elem_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, mesh2d, "geo_dofs");

  const MeshType::dof_storage_type &geo_cell_dofs = *mesh2d.dof_storage("geo_dofs");

  std::cout << std::endl;
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Iterating over the cells group by group:" << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  for (MeshType::dof_storage_type::const_dof_range_typed cells :
       geo_cell_dofs.all_active_dof_groups())
  {
    for (MeshType::dof_storage_type::const_dof_iterator_typed dof_iter = cells.begin();
         dof_iter != cells.end(); ++dof_iter)
    {
      // std::cout << *dof_iter << std::endl;

      const MeshEntity cell = dof_iter->mesh_entity();
      std::cout << "[" << cell.idx() << "]  " << cell << std::endl;
      for (Uint i = 0; i < cell.nb_sub_elements(_1D); ++i)
      {
        std::cout << " " << cell.sub_entity(_1D, i) << std::endl;
      }
      std::cout << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Iteration over all cells of all groups in one loop:" << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  for (MeshType::dof_storage_type::const_dof_iterator dof_iter = geo_cell_dofs.begin();
       dof_iter != geo_cell_dofs.end(); ++dof_iter)
  {
    const MeshEntity cell = dof_iter->mesh_entity();
    std::cout << "[" << cell.idx() << "]  " << cell << std::endl;
  }

  std::cout << std::endl;
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Testing direct access to cells:" << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  for (Uint i = 0; i < geo_cell_dofs.nb_active_cells(); ++i)
  {
    const MeshEntity cell = geo_cell_dofs.active_cell(ActiveIdx(i));
    std::cout << "[" << cell.idx() << "]  " << cell << std::endl;
  }

  std::cout << std::endl;

  // std::cout << std::endl << "************ FINISHED CELL ITERATOR TEST
  // ************" << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dof_iterator_synchronization_utest)
{
  MeshType mesh2d("mesh2d");

  const std::string infilename = "rectangle_mixed_elem_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, mesh2d, "geo_dofs_p1");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs_p1 = mesh2d.dof_storage("geo_dofs_p1");
  common::PtrHandle<MeshType::dof_storage_type> geo_dofs_p2 =
      mesh2d.create_dof_storage("geo_dofs_p2");

  MeshType::dof_storage_type::clone_continuous(mesh2d, *geo_dofs_p1, *geo_dofs_p2, P2,
                                               PointSetID::Equidist);

  const MeshType::dof_storage_type &cell_dofs_p1 = *geo_dofs_p1;
  const MeshType::dof_storage_type &cell_dofs_p2 = *geo_dofs_p2;

  MeshType::dof_storage_type::const_dof_iterator it2 = cell_dofs_p2.cbegin();

  for (MeshType::dof_storage_type::const_dof_range_typed const &dof_range :
       cell_dofs_p1.all_active_dof_groups())
  {
    for (MeshType::dof_storage_type::const_dof_iterator_typed it1 = dof_range.begin();
         it1 != dof_range.end(); ++it1)
    {
      synchronize_dof_iterators(it1, it2);
      std::cout << it1->mesh_entity() << " [synchronized] " << it2->mesh_entity() << std::endl;

      const MeshEntity cell_p1 = it1->mesh_entity();
      const MeshEntity cell_p2 = it2->mesh_entity();

      BOOST_CHECK_LE(cell_p1.nb_vert(), cell_p2.nb_vert());

      for (Uint v = 0; v < cell_p1.nb_vert(); ++v)
      {
        BOOST_CHECK_EQUAL(cell_p1.vertex(v), cell_p2.vertex(v));
      }
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(bdry_dof_iterator_loop_utest_2D)
{
  MeshType mesh2d("mesh2d");

  const std::string infilename = "rectangle_mixed_elem_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, mesh2d, "geo_dofs");

  const MeshType::dof_storage_type &geo_cell_dofs = *mesh2d.dof_storage("geo_dofs");
  const MeshBoundarySet<Cart2D> mesh_boundary     = mesh2d.all_boundaries();

  typedef BoundaryFacets<Cart2D, Cart2D::TDIM - 1>::const_dof_iterator dof_iterator_type;

  const std::shared_ptr<BoundaryFacets<Cart2D, Cart2D::TDIM - 1>> bottom_bdry =
      mesh_boundary.domain("Bottom");

  std::cout << "Bottom boundary" << std::endl;

  const std::vector<std::pair<Uint, Uint>> bottom_bdry_check = {{0, 6}, {6, 7},  {7, 8},   {8, 1},
                                                                {1, 9}, {9, 10}, {10, 11}, {11, 2}};

  Uint i = 0;
  for (dof_iterator_type it = bottom_bdry->cbegin(geo_cell_dofs);
       it != bottom_bdry->cend(geo_cell_dofs); ++it)
  {
    const MeshEntity bcell = it->mesh_entity();
    BOOST_CHECK_EQUAL(bcell.vertex(0), bottom_bdry_check[i].first);
    BOOST_CHECK_EQUAL(bcell.vertex(1), bottom_bdry_check[i].second);
    i++;
    std::cout << bcell << std::endl;
  }

  // Loop over the same cells again, this time directly extracting each cell
  for (Uint ac = 0; ac < bottom_bdry->nb_active_cells(); ++ac)
  {
    const MeshEntity bcell = bottom_bdry->active_cell(geo_cell_dofs, ActiveIdx(ac));
    BOOST_CHECK_EQUAL(bcell.vertex(0), bottom_bdry_check[ac].first);
    BOOST_CHECK_EQUAL(bcell.vertex(1), bottom_bdry_check[ac].second);
  }

  const std::shared_ptr<BoundaryFacets<Cart2D, Cart2D::TDIM - 1>> right_bdry =
      mesh_boundary.domain("Right");

  std::cout << "Right boundary" << std::endl;

  const std::vector<std::pair<Uint, Uint>> right_bdry_check = {
      {2, 12}, {12, 13}, {13, 14}, {14, 3}};

  i = 0;
  for (dof_iterator_type it = right_bdry->cbegin(geo_cell_dofs);
       it != right_bdry->cend(geo_cell_dofs); ++it)
  {
    const MeshEntity bcell = it->mesh_entity();
    BOOST_CHECK_EQUAL(bcell.vertex(0), right_bdry_check[i].first);
    BOOST_CHECK_EQUAL(bcell.vertex(1), right_bdry_check[i].second);
    i++;
    std::cout << bcell << std::endl;
  }

  // Loop over the same cells again, this time directly extracting each cell
  for (Uint ac = 0; ac < right_bdry->nb_active_cells(); ++ac)
  {
    const MeshEntity bcell = right_bdry->active_cell(geo_cell_dofs, ActiveIdx(ac));
    BOOST_CHECK_EQUAL(bcell.vertex(0), right_bdry_check[ac].first);
    BOOST_CHECK_EQUAL(bcell.vertex(1), right_bdry_check[ac].second);
  }

  const std::shared_ptr<BoundaryFacets<Cart2D, Cart2D::TDIM - 1>> top_bdry =
      mesh_boundary.domain("Top");

  std::cout << "Top boundary" << std::endl;

  const std::vector<std::pair<Uint, Uint>> top_bdry_check = {{3, 15}, {15, 16}, {16, 17}, {17, 4},
                                                             {4, 18}, {18, 19}, {19, 20}, {20, 5}};

  i = 0;
  for (dof_iterator_type it = top_bdry->cbegin(geo_cell_dofs); it != top_bdry->cend(geo_cell_dofs);
       ++it)
  {
    const MeshEntity bcell = it->mesh_entity();
    BOOST_CHECK_EQUAL(bcell.vertex(0), top_bdry_check[i].first);
    BOOST_CHECK_EQUAL(bcell.vertex(1), top_bdry_check[i].second);
    i++;
    std::cout << bcell << std::endl;
  }

  // Loop over the same cells again, this time directly extracting each cell
  for (Uint ac = 0; ac < top_bdry->nb_active_cells(); ++ac)
  {
    const MeshEntity bcell = top_bdry->active_cell(geo_cell_dofs, ActiveIdx(ac));
    BOOST_CHECK_EQUAL(bcell.vertex(0), top_bdry_check[ac].first);
    BOOST_CHECK_EQUAL(bcell.vertex(1), top_bdry_check[ac].second);
  }

  const std::shared_ptr<BoundaryFacets<Cart2D, Cart2D::TDIM - 1>> left_bdry =
      mesh_boundary.domain("Left");

  std::cout << "Left boundary" << std::endl;

  const std::vector<std::pair<Uint, Uint>> left_bdry_check = {{5, 21}, {21, 22}, {22, 23}, {23, 0}};

  i = 0;
  for (dof_iterator_type it = left_bdry->cbegin(geo_cell_dofs);
       it != left_bdry->cend(geo_cell_dofs); ++it)
  {
    const MeshEntity bcell = it->mesh_entity();
    BOOST_CHECK_EQUAL(bcell.vertex(0), left_bdry_check[i].first);
    BOOST_CHECK_EQUAL(bcell.vertex(1), left_bdry_check[i].second);
    i++;
    std::cout << bcell << std::endl;
  }

  // Loop over the same cells again, this time directly extracting each cell
  for (Uint ac = 0; ac < left_bdry->nb_active_cells(); ++ac)
  {
    const MeshEntity bcell = left_bdry->active_cell(geo_cell_dofs, ActiveIdx(ac));
    BOOST_CHECK_EQUAL(bcell.vertex(0), left_bdry_check[ac].first);
    BOOST_CHECK_EQUAL(bcell.vertex(1), left_bdry_check[ac].second);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(bdry_dof_iterator_loop_utest_3D)
{
  Tria<Cart3D> mesh3d("mesh3d");

  const std::string infilename = "unit_cube_mixed_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, mesh3d, "geo_dofs");

  const Tria<Cart3D>::dof_storage_type &geo_cell_dofs = *mesh3d.dof_storage("geo_dofs");
  const MeshBoundarySet<Cart3D> mesh_boundary         = mesh3d.all_boundaries();

  typedef BoundaryFacets<Cart3D, Cart3D::TDIM - 1>::const_dof_iterator_typed dof_iterator_type;

  const std::shared_ptr<BoundaryFacets<Cart3D, Cart3D::TDIM - 1>> top_bdry =
      mesh_boundary.domain("top");

  std::cout << "Top boundary (iteration by cell types)" << std::endl;

  std::vector<common::IteratorRange<dof_iterator_type>> iter_ranges;
  top_bdry->all_bdry_dof_ranges(geo_cell_dofs, iter_ranges);

  for (Uint r = 0; r < iter_ranges.size(); ++r)
  {
    std::cout << "Range " << r << std::endl;
    for (dof_iterator_type it = iter_ranges[r].begin(); it != iter_ranges[r].end(); ++it)
    {
      std::cout << it->mesh_entity() << " [" << it->mesh_entity().pt_set_id().as_string() << "]"
                << std::endl;
    }
  }
}

// ================================ TIMINGS ===================================

BOOST_AUTO_TEST_CASE(dof_iterator_timings_utest)
{
  MeshType mesh2d("mesh2d");

  mesh::MeshCreator::make_unit_quad(mesh2d, "geo_dofs", 1001,
                                    false); // 1 million quads

  typedef typename result_of::dof_map_t<Cart2D> cell_storage_type;
  cell_storage_type const &geo_cell_dofs = *mesh2d.dof_storage("geo_dofs");

  std::cout << "Number of cells in mesh = " << geo_cell_dofs.nb_active_cells() << std::endl;

  clock_t start, end;
  Real elapsed;

  // Loop over cells based on their types

  start = clock();

  for (MeshType::dof_storage_type::const_dof_range_typed cells :
       geo_cell_dofs.all_active_dof_groups())
  {
    for (MeshType::dof_storage_type::const_dof_iterator_typed dof_iter = cells.begin();
         dof_iter != cells.end(); ++dof_iter)
    {
      const MeshEntity cell = dof_iter->mesh_entity();
      // BOOST_CHECK_EQUAL(cell.nb_vert(), 4u);
      // std::cout << "[" << cell.idx() << "]  " << cell << std::endl;
    }
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "1 loop over all cells going type by type: " << elapsed << " sec" << std::endl;

  // ----------------------------
  // Loop over cells, no distinction of cell groups

  start = clock();

  for (typename cell_storage_type::const_dof_iterator dof_iter = geo_cell_dofs.cbegin();
       dof_iter != geo_cell_dofs.cend(); ++dof_iter)
  {
    const MeshEntity cell = dof_iter->mesh_entity();
    // BOOST_CHECK_EQUAL(cell.nb_vert(), 4u);
    // std::cout << "[" << cell.idx() << "]  " << cell << std::endl;
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "1 loop over all cells, types vary as we go from cell to cell: " << elapsed << " sec"
            << std::endl;

  // ----------------------------
  // Loop over cells, directly access n-th cell

  start = clock();

  for (Uint c = 0; c < geo_cell_dofs.nb_active_cells(); ++c)
  {
    const MeshEntity cell = geo_cell_dofs.active_cell(ActiveIdx(c));
    // BOOST_CHECK_EQUAL(cell.nb_vert(), 4u);
    // std::cout << "[" << cell.idx() << "]  " << cell << std::endl;
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "1 loop over all cells directly extracting n-th cell: " << elapsed << " sec"
            << std::endl;
}

// ----------------------------------------------------------------------------
