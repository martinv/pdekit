#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "common/PDEKit.hpp"
#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/ErrorEstimator.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/scalar/RotationAdvection3D.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodImplicit.hpp"
#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart3D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::RotationAdvection3D phys_model;

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // ------------------------------------------------------
  // INITIALIZE ENVIRONMENT
  // ------------------------------------------------------

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh3D = std::make_shared<MeshType>("mesh3d");

  const std::string infilename = "advection-tet.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh3D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh3D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh3D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_discontinuous(*mesh3D, *geo_dofs, *sol_dofs, P1,
                                                  PointSetID::Warpblend);

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P2;

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

  typedef rdm::DRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA, rdm::FacetDG> scheme_type;
  scheme_type scheme;

  scheme.configure_mesh_data(mesh3D, geo_dofs, sol_dofs);
  scheme.initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);

  scheme.set_vec_function(SolverVecFn::solution, solution);
  scheme.set_vec_function(SolverVecFn::residuals, residual);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  std::shared_ptr<scheme_type::bc_base_type> farfield =
      scheme.add_boundary_condition("StrongDirichlet", "strong_farfield", "farfield");
  farfield->set_expression(&Zero::value);

  std::shared_ptr<scheme_type::bc_base_type> strong_inlet =
      scheme.add_boundary_condition("StrongDirichlet", "strong_inlet", "inlet");
  strong_inlet->set_expression(&CosineHat3D::value);

  /*
  std::shared_ptr<scheme_type::bc_base_type> weak_inlet =
      scheme.add_boundary_condition("WeakDirichlet", "weak_inlet", "inlet");
  weak_inlet->set_expression(&CosineHat3D::value);
  */

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::cout.precision(15);

  math::DenseDVec<Real> res_L2_norm(1);

  start = std::chrono::system_clock::now();

  for (Uint iter = 0; iter < 50; ++iter)
  {
    const Real CFL = iter < 5 ? 3.0 : std::min(25.0, 4.0 + std::pow(1.30, iter - 5));

    scheme.set_cfl(CFL);

    if (iter % 5 == 0)
    {
      scheme.assemble_lhs_and_rhs(CFL);
      scheme.solve({SolverOption::RecomputePreconditioner});
    }
    else
    {
      scheme.assemble_rhs();
      scheme.solve({});
    }

    scheme.compute_residual_norm(res_L2_norm);

    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);
  }

  end                                 = std::chrono::system_clock::now();
  std::chrono::duration<Real> elapsed = end - start;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed.count() << " s" << std::endl;

  /// Compute the error in the solution

  /*
  Rotation3DSolution exact_solution;
  interpolation::ErrorEstimator error_estim(mesh);
  error_estim.compute_pointwise_error<Rotation3DSolution,MeshType>(exact_solution,*solution);
  */

  /// Write the output to file

  const std::string outfilename = "output_advection_p1_3D_discontinuous.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh3D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh3D, outfilename, *solution, "solution");
  meshwriter.append_nodal_function_to_file(
      *mesh3D, outfilename, *scheme.scal_function(SolverScalFn::time_step), "time_step");

  mpi_env.finalize();

  return 0;
}
