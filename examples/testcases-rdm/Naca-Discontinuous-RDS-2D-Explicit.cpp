#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodExplicit.hpp"
#include "solver/rdm/bc/WeakFarfield.hpp"
#include "solver/rdm/bc/WeakWall.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"

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

  const std::string infilename = "naca0012_p2_tri.msh";
  // const std::string infilename = "naca0012_p2_mixed_elem.msh";
  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P2,
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

  typedef rdm::DRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA, rdm::FacetDG> scheme_type;
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
  ic.set_domain(_2D, "InnerCells");
  math::DenseDVec<Real> ic_state(phys_model::NEQ);

  ic_state[0] = 1.22503;
  ic_state[1] = 208.30559;
  ic_state[2] = 7.27419;
  ic_state[3] = 271044.38;

  ic.apply_values("sol_dofs", ic_state,
                  *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  scheme.add_boundary_condition("WeakWall", "wall_bc", "Wall");

  // std::shared_ptr<bc_base_type> farfield_bc =
  auto farfield_bc = scheme.add_boundary_condition("WeakFarfield", "farfield_bc", "Farfield");
  farfield_bc->set_reference_state(ic_state);
  farfield_bc->print_bc_parameters();

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  math::DenseDVec<Real> residual_norm(phys_model::NEQ);

  clock_t start, end;
  Real elapsed;

  start = clock();

  for (Uint iter = 0; iter <= 10000; ++iter)
  {
    const Real CFL = 0.5;
    scheme.assemble_lhs_and_rhs(CFL);
    scheme.solve({});

    if (iter % 10 == 0)
    {
      interpolation::norm_L2(*residual, residual_norm);
      solver::SolverIO::print_iter_and_res_norm(iter, CFL, residual_norm);
    }
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  // ------------------------------------------------------
  // WRITE THE OUTPUT TO FILE
  // ------------------------------------------------------

  const std::string outfilename = "output_naca_discontinuous_tri_p1.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "u");

  // Write Mach number values to file
  scalar_function::ptr Mach_number = std::make_shared<scalar_function>("", "Ma");
  Mach_number->resize(nb_nodes);

  physics::Euler2DCons physical_model;
  physics::Euler2DProperties properties;
  physics::Euler2DProperties::SolGradM gradient_matrix;

  /*
  for (Uint c = 0; c < (*sol_dofs).nb_active_cells(); ++c)
  {
    const mesh::MeshEntity active_cell = (*sol_dofs).active_cell(mesh::ActiveIdx(c));
    const mesh::DofCoordinates<MeshConfig::GDIM> active_coords =
        (*sol_dofs).active_cell_coords(active_cell);

    for (Uint n = 0; n < active_cell.nb_vert(); ++n)
    {
      physical_model.compute_properties(active_coords.c(n),
                                        solution->const_value(active_cell.vertex(n)),
                                        gradient_matrix, properties);
      (*Mach_number)[active_cell.vertex(n)] = properties.Ma;
    }
  }
  */

  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution, [](const phys_model::Properties &props) { return props.Ma; },
      *Mach_number);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *Mach_number, "Ma");
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *residual, "Res");

  return 0;
}
