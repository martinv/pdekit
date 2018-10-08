#include <ctime>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
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

  MeshType::dof_storage_type::clone_continuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                               PointSetID::Equidist);

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P4;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  vector_function::ptr solution = std::make_shared<vector_function>("", "solution");
  vector_function::ptr residual = std::make_shared<vector_function>("", "residual");

  const Uint nb_nodes = (*sol_dofs).nb_nodes();

  solution->resize(phys_model::NEQ, nb_nodes);
  residual->resize(phys_model::NEQ, nb_nodes);

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  typedef rdm::PGRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;
  scheme_type scheme;

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
  BoundaryFacets<MeshConfig> const &bottom_bdry =
  *(boundaries.domain("bottom")); interpolation::FunctionSpace<MeshConfig,
  _1D>::ptr bottom_bdry_space =
      std::make_shared<interpolation::FunctionSpace<MeshConfig, _1D>>();
  bottom_bdry_space->set_reference_fe_values(bottom_bdry,
  mesh->topology().dof_storage(), SFunc::Lagrange, quadrature_order,
  PointSetID::Gauss);

  rdm::StrongDirichletBC<MeshConfig, physics::Euler2DCons, _1D>
  inlet("strong_inlet"); inlet.configure_mesh_data(topology.cells_ptr(),
  topology.all_boundaries_ptr(), mesh->geometry_ptr(),
  topology.dof_storage_ptr(), mesh->geometry_ptr(),
  topology.dof_storage_ptr(), "bottom");
  inlet.configure_spaces(bottom_bdry_space, bottom_bdry_space, solution);
  inlet.set_residuals(residual); inlet.set_update_coeffs(update_coeff);
  inlet.set_expression(&RiemannFansInlet2D::value);
  */

  /*
  // Impose the inlet BC strongly:
  std::shared_ptr<scheme_type::bc_base_type> inlet =
      scheme.add_boundary_condition("StrongDirichlet", "strong_inlet",
  "bottom"); inlet->set_expression(&RiemannFansInlet2D::value);
  */

  // Impose the inlet BC weakly:
  std::shared_ptr<scheme_type::bc_base_type> inlet =
      scheme.add_boundary_condition("WeakDirichlet", "weak_inlet", "bottom");
  inlet->set_expression(&RiemannFansInlet2D::value);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;

  math::DenseDVec<Real> res_L2_norm;

  start = clock();
  for (Uint iter = 0; iter < 100; ++iter)
  {
    const Real CFL = 0.5;

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

  const std::string outfilename = "output_jets2d_continuous.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "u");

  //  std::cout << "Solution: " << std::endl;
  //  std::cout << *solution << std::endl;

  //  std::cout << "Residuals: " << std::endl;
  //  std::cout << *residual << std::endl;

  //  std::cout << "Update coefficient: " << std::endl;
  //  std::cout << *update_coeff << std::endl;

  return 0;
}
