/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE operator_residual_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <ctime>
#include <iomanip>
#include <iostream>

#include "interpolation/OperatorResidualCellwise.hpp"
#include "interpolation/OperatorResidualNodal.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

struct PhysicsMockup2D
{
  enum
  {
    NEQ = 4
  };

  enum
  {
    DIM = 2
  };

  struct Properties
  {
    math::DenseSVec<Real, _2D> coords;
    math::DenseSVec<Real, 4> vars;
    math::DenseSMat<Real, 4, _2D> grad_vars;
  };

  template <typename V1, typename V2, typename M1>
  static void compute_properties(const V1 &coord, const V2 &sol, const M1 &grad_vars, Properties &p)
  {
    p.coords    = coord;     // cache the coordinates locally
    p.vars      = sol;       // cache the variables locally
    p.grad_vars = grad_vars; // cache the gradient of variables locally
  }

  template <typename M1>
  static void flux(const Properties &p, M1 &flux)
  {
    flux(0, X) = p.vars[X] * p.vars[X]; // x^2
    flux(1, X) = p.vars[X] * p.vars[X];
    flux(2, X) = p.vars[X] * p.vars[X];
    flux(3, X) = p.vars[X] * p.vars[X];

    flux(0, Y) = p.vars[Y] * p.vars[Y]; // y^2
    flux(1, Y) = p.vars[Y] * p.vars[Y];
    flux(2, Y) = p.vars[Y] * p.vars[Y];
    flux(3, Y) = p.vars[Y] * p.vars[Y];
  }

  template <typename V1>
  static void flux_divergence(const Properties &p, V1 &flux_div)
  {
    flux_div[0] = 2. * (p.vars[0] + p.vars[1]); // 2 * (u[0] + u[1])
    flux_div[1] = 2. * (p.vars[0] + p.vars[1]); // 2 * (u[0] + u[1])
    flux_div[2] = 2. * (p.vars[0] + p.vars[1]); // 2 * (u[0] + u[1])
    flux_div[3] = 2. * (p.vars[0] + p.vars[1]); // 2 * (u[0] + u[1])
  }
};

// ----------------------------------------------------------------------------

struct PhysicsMockup3D
{
  enum
  {
    NEQ = 5
  };

  enum
  {
    DIM = 3
  };

  struct Properties
  {
    math::DenseSVec<Real, _3D> coords;
    math::DenseSVec<Real, 5> vars;
    math::DenseSMat<Real, 5, _3D> grad_vars;
  };

  template <typename V1, typename V2, typename M1>
  static void compute_properties(const V1 &coord, const V2 &sol, const M1 &grad_vars, Properties &p)
  {
    p.coords    = coord;     // cache the coordinates locally
    p.vars      = sol;       // cache the variables locally
    p.grad_vars = grad_vars; // cache the gradient of variables locally
  }

  template <typename M1>
  static void flux(const Properties &p, M1 &flux)
  {
    flux(0, X) = p.vars[X] * p.vars[X]; // x^2
    flux(1, X) = p.vars[X] * p.vars[X];
    flux(2, X) = p.vars[X] * p.vars[X];
    flux(3, X) = p.vars[X] * p.vars[X];
    flux(4, X) = p.vars[X] * p.vars[X];

    flux(0, Y) = p.vars[Y] * p.vars[Y]; // y^2
    flux(1, Y) = p.vars[Y] * p.vars[Y];
    flux(2, Y) = p.vars[Y] * p.vars[Y];
    flux(3, Y) = p.vars[Y] * p.vars[Y];
    flux(4, Y) = p.vars[Y] * p.vars[Y];

    flux(0, Z) = p.vars[Z] * p.vars[Z]; // z^2
    flux(1, Z) = p.vars[Z] * p.vars[Z];
    flux(2, Z) = p.vars[Z] * p.vars[Z];
    flux(3, Z) = p.vars[Z] * p.vars[Z];
    flux(4, Z) = p.vars[Z] * p.vars[Z];
  }

  template <typename V1>
  static void flux_divergence(const Properties &p, V1 &flux_div)
  {
    flux_div[0] = 2. * (p.vars[0] + p.vars[1] + p.vars[2]); // 2 * (u[0] + u[1] + u[2])
    flux_div[1] = 2. * (p.vars[0] + p.vars[1] + p.vars[2]); // 2 * (u[0] + u[1] + u[2])
    flux_div[2] = 2. * (p.vars[0] + p.vars[1] + p.vars[2]); // 2 * (u[0] + u[1] + u[2])
    flux_div[3] = 2. * (p.vars[0] + p.vars[1] + p.vars[2]); // 2 * (u[0] + u[1] + u[2])
    flux_div[4] = 2. * (p.vars[0] + p.vars[1] + p.vars[2]); // 2 * (u[0] + u[1] + u[2])
  }
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
void test_operator_residual(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                            interpolation::VectorMeshFunction<Real> &operator_residual_nodal,
                            interpolation::VectorMeshFunction<Real> &operator_residual_cellwise)
{
  interpolation::OperatorResidualNodal<MeshConfig, Physics> operator_res_nodal;
  interpolation::OperatorResidualCellwise<MeshConfig, Physics> operator_res_cellwise;
  operator_res_nodal.setup(sol_dofs);
  operator_res_cellwise.setup(sol_dofs);

  interpolation::VectorMeshFunction<Real> solution("", "solution");
  solution.resize(Physics::NEQ, sol_dofs.nb_nodes());

  mesh::adapt::LocalInterpolator loc_interpolator;

  // Fill the solution field
  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell                         = sol_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      const math::DenseConstVecView<Real> node_coord = cell_coords.row_transpose(n);
      interpolation::VectorMeshFunction<Real>::entry_type function_nodal_value =
          solution.value(cell.vertex(n));

      // node_coord.size() should be equal to Physics::DIM
      for (Uint i = 0; i < node_coord.size(); ++i)
      {
        function_nodal_value[i] = node_coord[i];
      }

      for (Uint i = node_coord.size(); i < Physics::NEQ; ++i)
      {
        function_nodal_value[i] = 1.0;
      }
    }
  }
  operator_res_nodal.evaluate(sol_dofs, solution, operator_residual_nodal);
  operator_res_cellwise.evaluate(sol_dofs, solution, operator_residual_cellwise);

  typename Physics::Properties phys_props;
  math::DenseDVec<Real> flux_div(Physics::NEQ);

  // Main loop over cells: compare computed and expected values of flux
  // divergence
  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell                     = sol_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), sol_cell.pt_set_id(), tcell_view.coordinates());

    for (Uint v = 0; v < sol_cell.nb_vert(); ++v)
    {
      const math::DenseConstVecView<Real> node_coord = cell_coords.row_transpose(v);
      phys_props.coords                              = node_coord;
      // std::cout << "Coordinates = " << phys_props.coords << std::endl;
      phys_props.vars = solution.const_value(sol_cell.vertex(v));
      // phys_props.grad_vars = ...

      Physics::flux_divergence(phys_props, flux_div);

      interpolation::VectorMeshFunction<Real>::const_entry_type computed_flux_div =
          operator_residual_nodal.const_value(sol_cell.vertex(v));

      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        BOOST_CHECK_LE(std::abs(flux_div[eq] - computed_flux_div[eq]), 1.e-10);
      }

      /*
      std::cout << "Reference = " << flux_div << std::endl;
      std::cout << "Computed  = " <<
      operator_residual_nodal.const_value(sol_cell.vertex(v))
                << std::endl;
      std::cout << std::endl;
      */

    } // Loop over vertices of solution cell

    // std::cout <<
    // "---------------------------------------------------------------" <<
    // std::endl;

  } // Loop over cells in dof handler
} // test_operator_residual

typedef mesh::Cart2D MeshConfig2D;
typedef mesh::Cart3D MeshConfig3D;

typedef mesh::Tria<MeshConfig2D> MeshType2D;
typedef mesh::Tria<MeshConfig3D> MeshType3D;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(operator_residual_2D_utest)
{
  MeshType2D mesh2d("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_tri_p2.msh", mesh2d, "geo_dofs");
  common::PtrHandle<MeshType2D::dof_storage_type> geo_dofs = mesh2d.dof_storage("geo_dofs");
  common::PtrHandle<MeshType2D::dof_storage_type> sol_dofs = mesh2d.create_dof_storage("sol_dofs");

  MeshType2D::dof_storage_type::clone_discontinuous(mesh2d, *geo_dofs, *sol_dofs, P4,
                                                    PointSetID::Warpblend);

  interpolation::VectorMeshFunction<Real> operator_residual_nodal("", "operator_res_nodal");
  interpolation::VectorMeshFunction<Real> operator_residual_cellwise("", "operator_res_cellwise");

  test_operator_residual<MeshConfig2D, PhysicsMockup2D, _2D>(*sol_dofs, operator_residual_nodal,
                                                             operator_residual_cellwise);

  mesh::gmsh::GmshWriter mesh_writer;
  mesh_writer.write_mesh_to_file(mesh2d, "sol_dofs", "operator_res_2D.msh");
  mesh_writer.append_nodal_function_to_file(mesh2d, "operator_res_2D.msh", operator_residual_nodal,
                                            "op_res_nodal");
  mesh_writer.append_cell_function_to_file(mesh2d, "operator_res_2D.msh",
                                           operator_residual_cellwise, "op_res_cellwise");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(operator_residual_3D_utest)
{
  MeshType3D mesh3d("cube3D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_cube_tet_p3.msh", mesh3d, "geo_dofs");
  common::PtrHandle<MeshType3D::dof_storage_type> geo_dofs = mesh3d.dof_storage("geo_dofs");
  common::PtrHandle<MeshType3D::dof_storage_type> sol_dofs = mesh3d.create_dof_storage("sol_dofs");

  MeshType3D::dof_storage_type::clone_continuous(mesh3d, *geo_dofs, *sol_dofs, P4,
                                                 PointSetID::Warpblend);

  interpolation::VectorMeshFunction<Real> operator_residual_nodal("", "operator_res_nodal");
  interpolation::VectorMeshFunction<Real> operator_residual_cellwise("", "operator_res_cellwise");

  test_operator_residual<MeshConfig3D, PhysicsMockup3D, _3D>(*sol_dofs, operator_residual_nodal,
                                                             operator_residual_cellwise);

  mesh::gmsh::GmshWriter mesh_writer;
  mesh_writer.write_mesh_to_file(mesh3d, "sol_dofs", "operator_res_3D.msh");
  mesh_writer.append_nodal_function_to_file(mesh3d, "operator_res_3D.msh", operator_residual_nodal,
                                            "op_res_nodal");
  mesh_writer.append_cell_function_to_file(mesh3d, "operator_res_3D.msh",
                                           operator_residual_cellwise, "op_res_cellwise");
}

// ----------------------------------------------------------------------------
