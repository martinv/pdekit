#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>

// #include <sys/resource.h>

#include "common/PDEKit.hpp"
#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/ErrorEstimator.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/scalar/RotationAdvection2D.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodExplicit.hpp"
#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/bc/WeakDirichletBC.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"
#include "tools/MeshInspector.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::RotationAdvection2D<physics::AroundOrigin2D> phys_model;

// ----------------------------------------------------------------------------

int main()
{
  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  // struct rusage r_usage;

  std::shared_ptr<common::Component> mesh_component = common::make_component<MeshType>("mesh2D");
  std::shared_ptr<MeshType> mesh2D                  = mesh_component->handle<MeshType>();

  const std::string infilename = "advection-tri.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  mesh2D->change_std_region_types([](const PointSetTag &tag_old) -> PointSetTag {
    return PointSetTag(tag_old.elem_shape(), P2, tag_old.ref_topology());
  });

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  // Create a discontinuous connectivity for solution
  MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P2,
                                                  PointSetID::Equidist);

  // First adaptation pass
  MeshAdaptSequence<Cart2D> adapt_schedule;

  const Uint nb_old_cells = (*geo_dofs).nb_active_cells();

  std::vector<CellTransform> cell_adapt_ops(nb_old_cells);
  cell_adapt_ops.assign(nb_old_cells, CellTransform::NO_TRANS);

  // Refine cells around [0.7;0.6]

  for (Uint ac = 0; ac < (*sol_dofs).nb_active_cells(); ++ac)
  {
    const mesh::CellTopologyView<MeshConfig> tcell         = mesh2D->active_cell(ActiveIdx(ac));
    const mesh::CellGeometry<MeshConfig::GDIM> cell_coords = tcell.coordinates();
    Real cx0                                               = 0.0;
    Real cx1                                               = 0.0;

    const Uint nb_vert = cell_coords.size();

    for (Uint v = 0; v < nb_vert; ++v)
    {
      const math::DenseVecView<const Real> node_coords = cell_coords.const_node_view(v);
      cx0 += node_coords[X0];
      cx1 += node_coords[X1];
    }
    cx0 /= nb_vert;
    cx1 /= nb_vert;

    const Real dist = std::sqrt((cx0 - 0.7) * (cx0 - 0.7) + (cx1 - 0.6) * (cx1 - 0.6));

    if (dist < 0.2)
    {
      cell_adapt_ops[ac] = CellTransform::UNIFORM_REFINE;
    }
    else
    {
      cell_adapt_ops[ac] = CellTransform::NO_TRANS;
    }
  }

  // cell_adapt_ops[40] = CellAdaptOpID::UNIFORM_REFINE;
  // cell_adapt_ops[41] = CellAdaptOpID::UNIFORM_REFINE;

  adapt_schedule.define_h_adapt_ops(*mesh2D, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);
  mesh2D->adapt(adapt_schedule);
  (*geo_dofs).adapt(adapt_schedule);
  (*sol_dofs).adapt(adapt_schedule);

  /*
  tools::MeshInspector<MeshConfig> mesh_inspector;

  const Real mesh_volume = mesh_inspector.check_mesh_consistency(*sol_mesh,
  4); std::cout << "Mesh inspector: mesh volume = " << mesh_volume <<
  std::endl; mesh_inspector.print_active_facets(*sol_mesh);

  mesh_inspector.print_active_topo_cell_info(*sol_mesh, 544);

  interpolation::VectorMeshFunction<Real> cell_colors("", "cell_colors");
  mesh_inspector.color_mesh(*sol_mesh, cell_colors);

  sol_mesh->topology().cells().print_complete_skeleton(MeshConfig::TDIM-1);
  */

  /*
  mesh::MeshBoundarySet<MeshConfig> const &geo_mesh_boundaries =
      geo_mesh->topology().all_boundaries();
  mesh::MeshBoundarySet<MeshConfig> const &sol_mesh_boundaries =
      sol_mesh->topology().all_boundaries();
  */

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P4;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  vector_function::ptr solution = std::make_shared<vector_function>("", "solution");
  vector_function::ptr residual = std::make_shared<vector_function>("", "residual");

  // Choose scheme, configure the function spaces (solution, residual,
  // update coefficients

  const Uint nb_nodes = (*sol_dofs).nb_nodes();

  solution->resize(phys_model::NEQ, nb_nodes);
  residual->resize(phys_model::NEQ, nb_nodes);

  //  interpolation::FunctionSpace::Function::Ptr f_ptr =
  // fs_solution->function("residual");
  //  std::cout << "Got function " << f_ptr->name() << std::endl;

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  typedef rdm::DRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA, rdm::FacetDG> scheme_type;
  scheme_type scheme;

  /*
  scheme.configure_spaces(fs_geometry_cells, fs_solution_cells,
  fs_geometry_cell_contours, fs_solution_cell_contours, fs_geometry_facets,
  fs_solution_facets);
  */

  scheme.configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme.initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);

  scheme.set_vec_function(SolverVecFn::solution, solution);
  scheme.set_vec_function(SolverVecFn::residuals, residual);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  /*
  // Initial condition
  solver::InitialCondition<MeshConfig> ic("hat_initial_condition", sol_mesh);

  ic.set_domain(_2D, "fluid");
  ic.set_expression(&CosineHat2D::value);
  ic.apply_function(*solution);
  */

  /*
  result_of::geometry<MeshConfig>::type const &sol_geometry =
  sol_mesh->geometry();

  for (Uint c = 0; c < sol_cells.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = sol_cells.cell(c);

    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      const math::ConstVectorBlock<Real> node_coords =
  sol_geometry.node(cell.vertex(v));
      // const Real init_value = std::sin(node_coords[X0] * node_coords[X1]);
      const Real init_value = node_coords[X0] * node_coords[X0] +
  node_coords[X1] * node_coords[X1];
      (*solution)(0, cell.vertex(v)) = init_value;
    }
  }
  */

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  std::shared_ptr<scheme_type::bc_base_type> farfield =
      scheme.add_boundary_condition("StrongDirichlet", "strong_farfield", "farfield");
  farfield->set_expression(&Zero::value);

  std::shared_ptr<scheme_type::bc_base_type> strong_inlet =
      scheme.add_boundary_condition("StrongDirichlet", "strong_inlet", "inlet");
  strong_inlet->set_expression(&CosineHat2D::value);

  /*
  std::shared_ptr<scheme_type::bc_base_type> weak_inlet =
      scheme.add_boundary_condition("WeakDirichlet", "weak_inlet", "inlet");
  weak_inlet->set_expression(&CosineHat2D::value);
  */

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::cout.precision(15);

  math::DenseDVec<Real> res_L1_norm(1);

  start = std::chrono::system_clock::now();

  for (Uint iter = 0; iter < 950; ++iter) // 950
  {
    const Real CFL = 0.5;

    scheme.assemble_lhs_and_rhs(CFL);
    scheme.solve({});

    norm_L1(*residual, res_L1_norm);

    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L1_norm);
  }

  end                                 = std::chrono::system_clock::now();
  std::chrono::duration<Real> elapsed = end - start;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed.count() << " s" << std::endl;
  // getrusage(RUSAGE_SELF, &r_usage);
  // std::cout << "Memory usage: " << r_usage.ru_maxrss << " Kb" << std::endl;

  /// Compute the error in the solution
  Rotation2DSolution exact_solution;
  interpolation::ErrorEstimator<MeshConfig> error_estim(mesh2D);
  error_estim.compute_pointwise_L2_error(exact_solution, *sol_dofs, *solution);
  error_estim.compute_L1_error(exact_solution, *mesh2D, *sol_dofs, *solution);
  error_estim.compute_L2_error(exact_solution, *mesh2D, *sol_dofs, *solution);
  error_estim.compute_infty_error(exact_solution, *mesh2D, *sol_dofs, *solution);

  /// Write the output to file
  const std::string outfilename = "output_advection_rds_2D_h_adapt.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "solution");
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *residual, "residual");
  // meshwriter.append_cell_function_to_file(*sol_mesh, outfilename,
  // cell_colors, "cell_color");
  // meshwriter.append_nodal_function_to_file(*sol_mesh, outfilename,
  // *update_coeff, "update_coeff");

  /*
  vtk::VtkWriter vtkwriter;
  vtkwriter.append_nodal_function_to_file(*sol_mesh,
  "output_advection_rds_2D_h_adapt.vtu", *solution, "solution");
  */

  return 0;
}
