/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE scalar_mesh_function_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

/// PDEKIT headers
#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;

struct ScalarMeshFunctionFixture
{
  /// common setup for each test case
  ScalarMeshFunctionFixture()
  {
  }

  ~ScalarMeshFunctionFixture()
  {
  }

  // MeshCreator m_creator;
  mesh::gmsh::GmshWriter m_gmsh_writer;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(ScalarMeshFunction_TestSuite, ScalarMeshFunctionFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(scalar_mesh_function_adapt_utest)
{
  mesh::Tria<mesh::Cart2D> unit_triag_mesh("unit_triag_mesh");
  mesh::MeshCreator::make_unit_triangle(unit_triag_mesh, "geo_dofs", 30);

  typedef result_of::dof_map_t<mesh::Cart2D> cell_dofs_type;
  common::PtrHandle<cell_dofs_type> geo_cell_dofs = unit_triag_mesh.dof_storage("geo_dofs");
  common::PtrHandle<cell_dofs_type> sol_cell_dofs = unit_triag_mesh.create_dof_storage("sol_dofs");

  // Create discontinuous dofs for solution
  cell_dofs_type::clone_discontinuous(unit_triag_mesh, *geo_cell_dofs, *sol_cell_dofs, P2,
                                      PointSetID::Equidist);

  const Uint nb_cells = (*sol_cell_dofs).nb_active_cells();
  // const Uint nb_facets = sol_facets.nb_cells();

  std::vector<Uint> cell_p_order(nb_cells, P1);
  for (Uint c = 0; c < nb_cells / 2; ++c)
  {
    cell_p_order[c] = P2;
    // cell_p_order[9] = P3;
  }

  // Set up mesh adapter
  mesh::MeshAdaptSequence<mesh::Cart2D> adapter;
  adapter.define_p_adapt_ops(cell_p_order);

  interpolation::ScalarMeshFunction<Real> f1("space", "f1");
  f1.resize((*sol_cell_dofs).nb_nodes());
  f1.fill(0.0);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < (*sol_cell_dofs).nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<mesh::Cart2D> tcell_view =
        (*sol_cell_dofs).tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity active_cell = (*sol_cell_dofs).active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> active_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), active_cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < active_cell.nb_vert(); ++n)
    {
      const math::DenseConstVecView<Real> node_coord = active_coords.row_transpose(n);
      f1[active_cell.vertex(n)]                      = node_coord[X0] + node_coord[X1];
    }
  }

  f1.adapt(*sol_cell_dofs, adapter);
  (*sol_cell_dofs).adapt(adapter);

  m_gmsh_writer.write_mesh_to_file(unit_triag_mesh, "sol_dofs",
                                   "p_adapted_triag_mesh_w_scalar_func.msh");
  m_gmsh_writer.append_nodal_function_to_file(
      unit_triag_mesh, "p_adapted_triag_mesh_w_scalar_func.msh", f1, "u_p_adapted");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
