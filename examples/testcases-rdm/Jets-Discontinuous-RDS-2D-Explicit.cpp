#include <ctime>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodExplicit.hpp"
#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<Cart2D> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::Euler2DCons phys_model;

// ----------------------------------------------------------------------------

int main()
{
  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  const std::string infilename = "trapezoid-tri-p1.msh";
  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                                  PointSetID::Equidist);

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P4;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  vector_function::ptr solution = std::make_shared<vector_function>("", "solution");
  vector_function::ptr residual = std::make_shared<vector_function>("", "residual");

  /*
  Uint nb_nodes = 0;

  typedef result_of::dof_handler<MeshConfig>::type cells_type;

  for (cells_type::const_dof_iterator dof_it = sol_cells.begin(); dof_it !=
  sol_cells.end();
       ++dof_it)
  {
    const MeshEntity cell = *dof_it;
    nb_nodes += cell.nb_vert();
  }
  */

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
  scheme.configure_spaces(fs_solution_cells, fs_solution_cells,
  fs_solution_cell_contours, fs_solution_cell_contours, fs_solution_facets,
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

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh2D);
  ic.set_domain(_2D, "domain");
  ic.set_expression(&RiemannFansInlet2D::value);
  ic.apply_function("sol_dofs",
                    *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  /*
  BoundaryFacets<MeshConfig> const &bottom_geo_bdry =
  *(geo_mesh_boundaries.domain("bottom"));
  interpolation::FunctionSpace<MeshConfig, _1D>::ptr bottom_geo_space =
      std::make_shared<interpolation::FunctionSpace<MeshConfig, _1D>>();
  bottom_geo_space->set_reference_fe_values(bottom_geo_bdry,
  geo_mesh->topology().dof_storage(), SFunc::Lagrange, quadrature_order,
  PointSetID::Gauss);

  BoundaryFacets<MeshConfig> const &bottom_sol_bdry =
  *(sol_mesh_boundaries.domain("bottom"));
  interpolation::FunctionSpace<MeshConfig, _1D>::ptr bottom_sol_space =
      std::make_shared<interpolation::FunctionSpace<MeshConfig, _1D>>();
  bottom_sol_space->set_reference_fe_values(bottom_sol_bdry,
  sol_mesh->topology().dof_storage(), SFunc::Lagrange, quadrature_order,
  PointSetID::Gauss);

  rdm::StrongDirichletBC<MeshConfig, phys_model, _1D> inlet("strong_inlet");
  inlet.configure_mesh_data(geo_topology.cells_ptr(),
  geo_topology.all_boundaries_ptr(), geo_mesh->geometry_ptr(),
  geo_topology.dof_storage_ptr(), sol_mesh->geometry_ptr(),
  sol_topology.dof_storage_ptr(), "bottom");
  inlet.configure_spaces(bottom_geo_space, bottom_sol_space, solution);
  inlet.set_residuals(nodal_residual);
  inlet.set_update_coeffs(update_coeff);
  inlet.set_expression(&RiemannFansInlet2D::value);
  */

  std::shared_ptr<scheme_type::bc_base_type> inlet =
      scheme.add_boundary_condition("StrongDirichlet", "strong_inlet", "bottom");
  inlet->set_expression(&RiemannFansInlet2D::value);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;

  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  start = clock();
  for (Uint iter = 0; iter < 250; ++iter)
  {
    const Real CFL = 0.3;

    scheme.assemble_lhs_and_rhs(CFL);
    scheme.solve({});
    scheme.compute_residual_norm(res_L2_norm);

    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Write the output to file

  const std::string outfilename = "output_jets2d_discontinuous.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "u");

  return 0;
}
