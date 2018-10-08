/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Local_residual_restriction_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iomanip>
#include <iostream>

#include "common/MPI/MPIEnv.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/ResidualRestrictionLocal.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/MathConstants.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

struct LocalResidualRestriction_Fixture
{
  /// common setup for each test case
  LocalResidualRestriction_Fixture()
  {
    /// arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~LocalResidualRestriction_Fixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;

  gmsh::GmshWriter m_gmsh_writer;

  interpolation::ResidualRestrictionLocal m_projection;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(ResidualRestrictionLocal_TestSuite, LocalResidualRestriction_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(check_residual_restriction_local_2D)
{
  typedef Tria<Cart2D> MeshType;

  // --------------------------------------------------------------------------
  // 1) CREATE THE FINE LEVEL ('SOURCE') MESH FROM WHICH WE ARE PROJECTING
  // --------------------------------------------------------------------------

  std::shared_ptr<MeshType> mesh2D(new MeshType("source_mesh"));

  // Create a quad domain [0,1] x [0,1] meshed with triangles
  MeshCreator::make_unit_quad(*mesh2D, "source_dofs", 30, true);

  // Increase the polynomial order
  common::PtrHandle<MeshType::dof_storage_type> source_dofs = mesh2D->dof_storage("source_dofs");
  (*source_dofs).upgrade(*mesh2D, P4, PointSetID::Warpblend);

  common::PtrHandle<MeshType::dof_storage_type> target_dofs_discontinuous =
      mesh2D->create_dof_storage("target_dofs_discontinuous");
  MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *source_dofs, *target_dofs_discontinuous,
                                                  P3, PointSetID::Equidist);

  common::PtrHandle<MeshType::dof_storage_type> target_dofs_continuous =
      mesh2D->create_dof_storage("target_dofs_continuous");
  MeshType::dof_storage_type::clone_continuous(*mesh2D, *source_dofs, *target_dofs_continuous, P3,
                                               PointSetID::Equidist);

  interpolation::FunctionSpace<Cart2D>::ptr fs_residual_src =
      std::make_shared<interpolation::FunctionSpace<Cart2D>>();

  // Here the quadrature order has to be the square of what's the element
  // polynomial order: we're going to
  // integrate the mass matrix entries \varphi_i \cdot \varphi_j
  fs_residual_src->set_reference_fe_values((*source_dofs).as_range(),
                                           [=](const ElemShape shape, const Uint order) {
                                             return mesh::sf::SFTag(shape, SFunc::Lagrange, order,
                                                                    ModalBasis::Modal);
                                           },
                                           [=](const ElemShape shape, const Uint order) {
                                             return mesh::PointSetTag(shape, P4, PointSetID::Gauss);
                                           });

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> residual_src =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "res_source");

  /// Configure the function space for solution

  const Uint nb_nodes_src = (*source_dofs).nb_nodes();

  residual_src->resize(3, nb_nodes_src);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < (*source_dofs).nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<Cart2D> topo_view = (*source_dofs).tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell                    = (*source_dofs).active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        topo_view.pt_set_id(), cell.pt_set_id(), topo_view.coordinates());

    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      const math::DenseConstVecView<Real> point = cell_coords.row_transpose(n);

      typename interpolation::VectorMeshFunction<Real>::entry_type nodal_value =
          residual_src->value(cell.vertex(n));

      nodal_value[0] = std::sin(math::pi * point[X0]) * std::sin(math::pi * point[X1]);
      nodal_value[1] = std::sin(2. * math::pi * point[X0]) * std::sin(3. * math::pi * point[X1]);
      nodal_value[2] = std::sin(4. * math::pi * point[X0]) * std::sin(6. * math::pi * point[X1]);
    }
  }

  // --------------------------------------------------------------------------
  // 2) CREATE THE COARSE LEVEL ('TARGET') MESH TO WHICH WE ARE PROJECTING
  //    In this case, we consider DISCONTINUOUS mesh
  // --------------------------------------------------------------------------

  std::shared_ptr<MeshType> target_mesh_discontinuous(new MeshType("target_mesh_discontinuous"));
  // Create a quad domain [0,1] x [0,1] meshed with triangles

  interpolation::FunctionSpace<Cart2D>::ptr fs_residual_tgt =
      std::make_shared<interpolation::FunctionSpace<Cart2D>>();
  fs_residual_tgt->set_reference_fe_values((*target_dofs_discontinuous).as_range(),
                                           [=](const ElemShape shape, const Uint order) {
                                             return mesh::sf::SFTag(shape, SFunc::Lagrange, order,
                                                                    ModalBasis::Modal);
                                           },
                                           [=](const ElemShape shape, const Uint order) {
                                             return mesh::PointSetTag(shape, P1, PointSetID::Gauss);
                                           });

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> residual_tgt =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "res_target");

  const Uint nb_nodes_tgt_discontinuous = (*target_dofs_discontinuous).nb_nodes();
  residual_tgt->resize(3, nb_nodes_tgt_discontinuous);
  residual_tgt->fill(0.0);

  // --------------------------------------------------------------------------
  // 3) PROJECT THE SOLUTION FROM SOURCE TO TARGET
  // --------------------------------------------------------------------------

  m_projection.build_projection_operator<Cart2D>(*mesh2D, *source_dofs, *target_dofs_discontinuous);
  m_projection.project<Cart2D>(*mesh2D, *source_dofs, *target_dofs_discontinuous, *residual_src,
                               *residual_tgt);

  std::cout << "Local residual restriction [discontinuous] finished" << std::endl;

  // --------------------------------------------------------------------------
  // 4) WRITE THE DATA TO GMSH FILES
  // --------------------------------------------------------------------------

  m_gmsh_writer.write_mesh_to_file(*mesh2D, "source_dofs",
                                   "source_mesh_residual_restriction_local.msh");
  m_gmsh_writer.append_nodal_function_to_file(*mesh2D, "source_mesh_residual_restriction_local.msh",
                                              *residual_src, "r_source");

  m_gmsh_writer.write_mesh_to_file(*mesh2D, "target_dofs_discontinuous",
                                   "target_mesh_discontinuous_residual_restriction_local.msh");
  m_gmsh_writer.append_nodal_function_to_file(
      *mesh2D, "target_mesh_discontinuous_residual_restriction_local.msh", *residual_tgt,
      "r_target");

  // --------------------------------------------------------------------------
  // 5) CREATE THE COARSE LEVEL ('TARGET') MESH TO WHICH WE ARE PROJECTING
  //    In this case, we consider CONTINUOUS mesh
  // --------------------------------------------------------------------------

  const Uint nb_nodes_tgt_continuous = (*target_dofs_continuous).nb_nodes();
  residual_tgt->resize(3, nb_nodes_tgt_continuous);
  residual_tgt->fill(0.0);

  // --------------------------------------------------------------------------
  // 3) PROJECT THE SOLUTION FROM SOURCE TO TARGET
  // --------------------------------------------------------------------------

  m_projection.build_projection_operator<Cart2D>(*mesh2D, *source_dofs, *target_dofs_continuous);
  m_projection.project<Cart2D>(*mesh2D, *source_dofs, *target_dofs_continuous, *residual_src,
                               *residual_tgt);

  std::cout << "Local residual restriction [continuous] finished" << std::endl;

  // --------------------------------------------------------------------------
  // 4) WRITE THE DATA TO GMSH FILES
  // --------------------------------------------------------------------------

  m_gmsh_writer.write_mesh_to_file(*mesh2D, "target_dofs_continuous",
                                   "target_mesh_continuous_residual_restriction_local.msh");
  m_gmsh_writer.append_nodal_function_to_file(
      *mesh2D, "target_mesh_continuous_residual_restriction_local.msh", *residual_tgt, "r_target");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
