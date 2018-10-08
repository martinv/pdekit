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
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/time/ExplicitTimeStepper.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart3D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::RotationAdvection3D phys_model;

// ============================================================================

int main()
{
  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh3D = std::make_shared<MeshType>("mesh3D");

  const std::string infilename = "advection-tet.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh3D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh3D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh3D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_continuous(*mesh3D, *geo_dofs, *sol_dofs, P1,
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

  typedef rdm::PGRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;
  scheme_type scheme;

  scheme.configure_mesh_data(mesh3D, geo_dofs, sol_dofs);
  scheme.initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme.set_vec_function(SolverVecFn::solution, solution);
  scheme.set_vec_function(SolverVecFn::residuals, residual);

  time::ExplicitEuler time_stepper;
  time_stepper.setup(phys_model::NEQ, nb_nodes);

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

  clock_t start, end;
  Real elapsed;
  std::cout.precision(15);

  math::DenseDVec<Real> res_L1_norm(1);

  start = clock();

  for (Uint iter = 0; iter < 500; ++iter)
  {
    /*
    residual->fill(0.0);
    time_update->reset_wave_speeds();

    scheme.iterate();

    for (Uint n = 0; n < solution->nb_entries(); ++n)
    {
      node_value_type sol_value = solution->value(n);
      const_node_value_type nodal_residual = residual->const_value(n);
      const Real nodal_wave_speed = (*time_update).wave_speed(n);

      for (Uint eq = 0; eq < physics::RotationAdvection3D::NEQ; ++eq)
      {
        sol_value[eq] -= 0.5 / nodal_wave_speed * nodal_residual[eq];
      }
    }
    */

    time_stepper.advance_in_time(scheme, 0.5);

    if (iter % 10 == 0)
    {
      norm_L1(*residual, res_L1_norm);

      std::cout << "Iter = " << std::setw(4) << iter << ", res = " << std::setw(20)
                << res_L1_norm[0] << " , log(res) = " << std::setw(20) << std::log(res_L1_norm[0])
                << std::endl;
    }
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Compute the error in the solution

  /*
  Rotation3DSolution exact_solution;
  interpolation::ErrorEstimator error_estim(mesh);
  error_estim.compute_pointwise_error<Rotation3DSolution,MeshType>(exact_solution,*solution);
  */

  /// Write the output to file

  const std::string outfilename = "output_advection_p1_3D.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh3D, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*mesh3D, outfilename, *solution, "solution");

  return 0;
}
