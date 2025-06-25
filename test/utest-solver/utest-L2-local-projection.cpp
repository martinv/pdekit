/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Local L2 - projection_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iomanip>
#include <iostream>

#include "interpolation/L2ProjectionLocal.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/MathConstants.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

struct L2Local_Fixture
{
  /// common setup for each test case
  L2Local_Fixture()
  {
    /// arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~L2Local_Fixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;

  gmsh::GmshWriter m_gmsh_writer;

  interpolation::L2ProjectionLocal m_projection;

  struct PolyRule
  {
    // Return: a tuple consisting of
    // 1) Shape function tag of source shape function (which we're
    // projecting FROM) 2) Shape function tag of target polynomial (i.e
    // space we're projecting TO)
    std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> sf_tags(const mesh::PointSetTag &geo_tag,
                                                         const mesh::PointSetTag &sol_tag_src,
                                                         const mesh::PointSetTag &sol_tag_tgt) const
    {
      const mesh::sf::SFTag src_poly_space_tag(sol_tag_src.elem_shape(), SFunc::Lagrange,
                                               sol_tag_src.poly_order(), ModalBasis::Modal);

      const mesh::sf::SFTag tgt_poly_space_tag(sol_tag_tgt.elem_shape(), SFunc::Lagrange,
                                               sol_tag_tgt.poly_order(), ModalBasis::Modal);

      const std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> result =
          std::make_tuple(src_poly_space_tag, tgt_poly_space_tag);
      return result;
    }

    // Return: a tuple consisting of
    // 1) Quadrature order to be used in order to construct LHS matrix
    // operator 2) Quadrature order to be used in order to construct RHS
    std::tuple<Uint, Uint> quadrature_orders(const mesh::PointSetTag &geo_tag,
                                             const mesh::PointSetTag &sol_tag_src,
                                             const mesh::PointSetTag &sol_tag_tgt) const
    {
      const Uint src_quad_order =
          std::max(sol_tag_src.poly_order() + sol_tag_tgt.poly_order(), geo_tag.poly_order());
      const Uint tgt_quad_order = std::max(2 * sol_tag_tgt.poly_order(), geo_tag.poly_order());

      const std::tuple<Uint, Uint> result = std::make_tuple(src_quad_order, tgt_quad_order);
      return result;
    }
  };
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(L2ProjectionLocal_TestSuite, L2Local_Fixture)

// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(check_L2_projection_Local_2D)
{
  typedef Tria<Cart2D> MeshType;

  // --------------------------------------------------------------------------
  // 1) CREATE THE FINE LEVEL ('SOURCE') MESH FROM WHICH WE ARE PROJECTING
  // --------------------------------------------------------------------------

  std::shared_ptr<MeshType> mesh2D(new MeshType("mesh2D"));

  // Create a quad domain [0,1] x [0,1] meshed with triangles
  MeshCreator::make_unit_quad(*mesh2D, "source_dofs", 30, true);

  // Increase the polynomial order
  common::PtrHandle<MeshType::dof_storage_type> source_dofs = mesh2D->dof_storage("source_dofs");
  (*source_dofs).upgrade(*mesh2D, P2, PointSetID::Warpblend);

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> solution_src =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "u_source");

  /// Configure the function space for solution

  const result_of::dof_map_t<Cart2D> &source_cell_dofs = *source_dofs;
  const Uint nb_nodes_src                              = source_cell_dofs.nb_nodes();

  solution_src->resize(3, nb_nodes_src);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < source_cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<Cart2D> tcell_view = source_cell_dofs.tcell(mesh::ActiveIdx(c));
    const MeshEntity cell = source_cell_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      const math::DenseConstVecView<Real> point = cell_coords.row_transpose(n);
      typename interpolation::VectorMeshFunction<Real>::entry_type nodal_value =
          solution_src->value(cell.vertex(n));

      nodal_value[0] = std::sin(math::pi * point[X]) * std::sin(math::pi * point[Y]);
      nodal_value[1] = std::sin(2. * math::pi * point[X0]) * std::sin(3. * math::pi * point[X1]);
      nodal_value[2] = std::sin(4. * math::pi * point[X0]) * std::sin(6. * math::pi * point[X1]);
    }
  }

  // --------------------------------------------------------------------------
  // 2) CREATE THE COARSE LEVEL ('TARGET') MESH TO WHICH WE ARE PROJECTING
  // --------------------------------------------------------------------------

  common::PtrHandle<MeshType::dof_storage_type> target_dofs =
      mesh2D->create_dof_storage("target_dofs");
  MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *source_dofs, *target_dofs, P1,
                                                  PointSetID::Equidist);

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> solution_tgt =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "u_target");

  /// Configure the function space for solution
  const Uint nb_nodes_tgt = (*target_dofs).nb_nodes();

  solution_tgt->resize(3, nb_nodes_tgt);
  solution_tgt->fill(0.0);

  /*
  for(Uint n = 0; n < nb_nodes_tgt; ++n)
  {
    const Real x = coordinates_tgt(n,X);
    const Real y = coordinates_tgt(n,Y);

    (*solution_tgt)(0,n) = std::sin(2.*PI*x) * std::sin(3.*PI*y);
  }
  */

  // --------------------------------------------------------------------------
  // 3) PROJECT THE SOLUTION FROM SOURCE TO TARGET
  // --------------------------------------------------------------------------

  m_projection.build_projection_operator<Cart2D>(*mesh2D, *source_dofs, *target_dofs, PolyRule{});
  m_projection.project<Cart2D>(*mesh2D, *source_dofs, *target_dofs, PolyRule{}, *solution_src,
                               *solution_tgt);

  std::cout << "Projection finished" << std::endl;

  // --------------------------------------------------------------------------
  // 4) WRITE THE DATA TO GMSH FILES
  // --------------------------------------------------------------------------

  m_gmsh_writer.write_mesh_to_file(*mesh2D, "source_dofs", "source_mesh_L2_local.msh");
  m_gmsh_writer.append_nodal_function_to_file(*mesh2D, "source_mesh_L2_local.msh", *solution_src,
                                              "u_source");

  m_gmsh_writer.write_mesh_to_file(*mesh2D, "target_dofs", "target_mesh_L2_local.msh");
  m_gmsh_writer.append_nodal_function_to_file(*mesh2D, "target_mesh_L2_local.msh", *solution_tgt,
                                              "u_target");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
