#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/scalar/Burgers2D.hpp"
#include "solver/SolverIO.hpp"
#include "solver/art_visc/ArtificialViscosity.hpp"
#include "solver/rdm/DRDMethodExplicit.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/blending_coeff/RDBlendingCoeff.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::Burgers2D phys_model;

// ----------------------------------------------------------------------------

Real left_bc_expr(const math::DenseConstVecView<Real> &point_coord,
                  const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                  const Uint component)
{
  return 1.5;
}

Real right_bc_expr(const math::DenseConstVecView<Real> &point_coord,
                   const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                   const Uint component)
{
  return -0.5;
}

Real top_bc_expr(const math::DenseConstVecView<Real> &point_coord,
                 const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                 const Uint component)
{
  return 1.5;
}

// ----------------------------------------------------------------------------

int main()
{
  // ------------------------------------------------------
  // SELECT SOLVER TYPE
  // ------------------------------------------------------

  typedef rdm::PGRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;
  // typedef rdm::DRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA,
  // rdm::FacetDG> scheme_type;
  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme = std::make_shared<scheme_type>();

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  // const std::string infilename = "square_quad_p3.msh";
  const std::string infilename = "square_tri_p3.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  if (scheme->is_continuous())
  {
    MeshType::dof_storage_type::clone_continuous(*mesh2D, *geo_dofs, *sol_dofs, P3,
                                                 PointSetID::Warpblend);
  }
  else
  {
    MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P3,
                                                    PointSetID::Equidist);
  }

  std::vector<Int> reordering;
  if (scheme->is_continuous())
  {
    scheme->compute_node_reordering(*mesh2D, *sol_dofs, reordering);
    (*sol_dofs).renumber_dofs(reordering);
  }
  else
  {
    scheme->compute_cell_reordering(*mesh2D, *sol_dofs, reordering);
    (*sol_dofs).renumber_dofs_blockwise(reordering);
  }

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P6;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  vector_function::ptr solution = std::make_shared<vector_function>("", "solution");
  vector_function::ptr residual = std::make_shared<vector_function>("", "residual");

  const Uint nb_nodes = (*sol_dofs).nb_nodes();

  solution->resize(phys_model::NEQ, nb_nodes);
  residual->resize(phys_model::NEQ, nb_nodes);

  // ------------------------------------------------------
  // DISCONTINUITY SENSOR SETUP
  // ------------------------------------------------------

  solver::ArtificialViscosity<MeshConfig> art_visc_operator;
  art_visc_operator.prepare_indicator_operators(*sol_dofs);

  scalar_function::ptr indicator = std::make_shared<scalar_function>("", "smoothness_indicator");
  indicator->resize((*geo_dofs).nb_active_cells());

  // scalar_function::ptr art_visc = std::make_shared<scalar_function>("",
  // "artificial_viscosity"); art_visc->resize((*geo_dofs).nb_active_cells());

  scalar_function::ptr art_visc_nodal =
      std::make_shared<scalar_function>("", "artificial_viscosity_nodal");
  art_visc_nodal->resize((*sol_dofs).nb_nodes());

  // ------------------------------------------------------
  // BLENDING COEFFICIENT FOR NONLINEAR SCHEMES
  // ------------------------------------------------------

  rdm::RDBlendingCoeff<MeshConfig> rd_blending_coeff_operator;

  // Types of viscosity coefficient: "LF_Blend", "ViscosityIndicator",
  // "WongJansen"
  const std::string blend_coeff_type = "LF_Blend";
  rd_blending_coeff_operator.setup(*sol_dofs, blend_coeff_type);

  // This is reference value for Wong-Jansen blending coefficient
  rd_blending_coeff_operator.set_param("u_inf", 1.0);

  scalar_function::ptr blending_coeff =
      std::make_shared<scalar_function>("", "blending_coefficient");
  blending_coeff->resize((*geo_dofs).nb_nodes());

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residual);
  scheme->set_artificial_viscosity(art_visc_nodal);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  std::shared_ptr<scheme_type::bc_base_type> bottom =
      scheme->add_boundary_condition("WeakDirichlet", "bottom_bc", "bottom");
  bottom->set_expression(&BurgersInlet2D::value);

  /*
  std::shared_ptr<scheme_type::bc_base_type> top =
      scheme->add_boundary_condition("StrongDirichlet", "top_bc", "top");
  top->set_expression(&top_bc_expr);
  */

  std::shared_ptr<scheme_type::bc_base_type> left =
      scheme->add_boundary_condition("WeakDirichlet", "left_bc", "left");
  left->set_expression(&left_bc_expr);

  std::shared_ptr<scheme_type::bc_base_type> right =
      scheme->add_boundary_condition("WeakDirichlet", "right_bc", "right");
  right->set_expression(&right_bc_expr);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  typedef vector_function::entry_type node_value_type;
  typedef vector_function::const_entry_type const_node_value_type;

  math::DenseDVec<Real> res_L2_norm(1);

  const std::clock_t cpu_start = std::clock();
  auto wall_start              = std::chrono::high_resolution_clock::now();

  for (Uint iter = 0; iter < 4000; ++iter) // 4000
  {
    const Real CFL = 0.5;
    scheme->assemble_lhs_and_rhs(CFL);
    scheme->solve({});

    /*
    rd_blending_coeff_operator.calculate(geo_mesh->topology().dof_storage(),
                                         sol_mesh->topology().dof_storage(),
    *solution, *blending_coeff);
    */

    if ((iter > 0) && ((iter % 5) == 0))
    {
      // art_visc_operator.compute_artificial_viscosity(*geo_dofs,
      // *sol_dofs, *solution, *art_visc);
      art_visc_operator.compute_artificial_viscosity_nodal(*mesh2D, *sol_dofs, *solution,
                                                           *art_visc_nodal);
    }

    const std::clock_t cpu_end = std::clock();
    const double cpu_duration  = (cpu_end - cpu_start) / (double)CLOCKS_PER_SEC;

    if (iter % 10 == 0)
    {
      std::cout.precision(15);
      norm_L2(*residual, res_L2_norm);

      solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);
    }
  }

  const std::clock_t cpu_end = std::clock();
  auto wall_end              = std::chrono::high_resolution_clock::now();

  const double cpu_duration = (cpu_end - cpu_start) / (double)CLOCKS_PER_SEC;
  auto wclock_duration = std::chrono::duration<double, std::milli>(wall_end - wall_start).count();

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop)        = " << cpu_duration << " s" << std::endl;
  std::cout << "Wall clock time (main iteration loop) = " << wclock_duration * 1.e-3 << " s"
            << std::endl;

  // Compute the smoothness indicator and artificial viscosity once more to
  // store for visualization

  art_visc_operator.compute_indicator_normalized(*sol_dofs, *solution, *indicator);
  // art_visc_operator.compute_artificial_viscosity(*geo_dofs, *sol_dofs,
  // *solution, *art_visc);

  art_visc_operator.compute_artificial_viscosity_nodal(*mesh2D, *sol_dofs, *solution,
                                                       *art_visc_nodal);

  // Compute the blending coefficient and store its values for visualization
  rd_blending_coeff_operator.calculate(*mesh2D, *sol_dofs, *solution, *blending_coeff);

  /// Write the output to file

  const std::string outfilename = "output_Burgers_2D.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "solution");
  meshwriter.append_cell_function_to_file(*mesh2D, outfilename, *indicator, "smoothness_indicator");
  // meshwriter.append_cell_function_to_file(*mesh2D, outfilename, *art_visc,
  // "artificial_viscosity");
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *art_visc_nodal,
                                           "artificial_viscosity_nodal");
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *blending_coeff,
                                           "blending_coefficient");

  /*
  vtk::VtkWriter vtkwriter;
  vtkwriter.append_nodal_function_to_file(*geo_mesh, "output.vtu", *solution,
  "solution");
  */

  return 0;
}
