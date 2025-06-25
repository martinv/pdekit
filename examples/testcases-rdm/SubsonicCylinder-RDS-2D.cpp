#include <ctime>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodImplicit.hpp"
#include "solver/rdm/PGRDMethodImplicit.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters/FacetSplitters.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::interpolation;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::Euler2DCons phys_model;

// ----------------------------------------------------------------------------

class CylinderFarfieldBC
{
  public:
  // Constructor
  CylinderFarfieldBC();

  // Destructor
  ~CylinderFarfieldBC();

  static Real Dirichlet_bc_value(const math::DenseConstVecView<Real> &node,
                                 const vector_function::const_entry_type &solution,
                                 const Uint component);
};

CylinderFarfieldBC::CylinderFarfieldBC()
{
}

CylinderFarfieldBC::~CylinderFarfieldBC()
{
}

Real CylinderFarfieldBC::Dirichlet_bc_value(const math::DenseConstVecView<Real> &node,
                                            const vector_function::const_entry_type &solution,
                                            const Uint component)
{
  Real result = 0.0;

  switch (component)
  {

    case 0:
      result = 1.2250;
      break;

    case 1:
      result = 158.407;
      break;

    case 2:
      result = 0.0;
      break;

    case 3:
      result = 263554.0;
      break;
  };

  return result;
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // ------------------------------------------------------
  // INITIALIZE ENVIRONMENT
  // ------------------------------------------------------

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // DEFINE METHOD TYPE
  // ------------------------------------------------------

  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;
  // typedef rdm::DRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA,
  // rdm::FacetDG> scheme_type;

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme = std::make_shared<scheme_type>();

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  // Remark: for p3 mesh, use CFL = 10.0 at least for the first 3000
  // iterations (continuous solver)

  // const std::string infilename = "cylinder2_regular_p1.msh";
  const std::string infilename = "cylinder2_regular_p2p2.msh";
  // const std::string infilename = "cylinder2_regular_p3p3.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  if (scheme->is_continuous())
  {
    MeshType::dof_storage_type::clone_continuous(*mesh2D, *geo_dofs, *sol_dofs, P2,
                                                 PointSetID::Warpblend);
  }
  else
  {
    MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P2,
                                                    PointSetID::Warpblend);
  }

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

  scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residual);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh2D);
  ic.set_domain(_2D, "InnerCells");
  math::DenseDVec<Real> ic_state(physics::Euler2DCons::NEQ);

  ic_state[0] = 1.2250;
  ic_state[1] = 158.407;
  ic_state[2] = 0.0;
  ic_state[3] = 263554.0;

  std::cout << "Inlet state = " << ic_state << std::endl;

  ic.apply_values("sol_dofs", ic_state,
                  *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  scheme->add_boundary_condition("WeakWall", "bottom_wall_bc", "LowerCylinder");
  scheme->add_boundary_condition("WeakWall", "top_wall_bc", "UpperCylinder");

  // std::shared_ptr<bc_base_type> farfield_bc =
  /*
  auto farfield_bc = scheme->add_boundary_condition("WeakFarfield",
  "farfield_bc", "Farfield"); farfield_bc->set_reference_state(ic_state);
  farfield_bc->print_bc_parameters();
  */

  /*
  std::shared_ptr<scheme_type::bc_base_type> weak_Dirichlet_farfield =
      scheme->add_boundary_condition("WeakDirichlet", "farfield_bc",
  "Farfield");
  weak_Dirichlet_farfield->set_expression(&CylinderFarfieldBC::Dirichlet_bc_value);
  */

  std::shared_ptr<scheme_type::bc_base_type> strong_Dirichlet_farfield =
      scheme->add_boundary_condition("StrongDirichlet", "farfield_bc", "Farfield");
  strong_Dirichlet_farfield->set_expression(&CylinderFarfieldBC::Dirichlet_bc_value);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;

  std::cout.precision(15);
  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  const Real CFLMax = scheme->is_continuous() ? 25.0 : 15.0;

  start = clock();
  for (Uint iter = 0; iter < 150; ++iter)
  {
    const Real CFL = iter < 10 ? 15.0 : std::min(CFLMax, 10.0 + std::pow(1.50, iter - 10));
    // const Real CFL = 10.0;

    scheme->set_cfl(CFL);

    if (iter % 5 == 0)
    {
      scheme->assemble_lhs_and_rhs(CFL);
      scheme->solve({SolverOption::RecomputePreconditioner});
    }
    else
    {
      scheme->assemble_rhs();
      scheme->solve({});
    }

    scheme->compute_residual_norm(res_L2_norm);

    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L2_norm);
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Write the output to file

  const std::string outfilename = "output_subsonic_cylinder_rds_2D.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "u");

  /// Write Mach number values to file
  scalar_function::ptr Mach_number = std::make_shared<scalar_function>("", "Ma");
  Mach_number->resize(nb_nodes);

  solver::PostprocessingUtils<MeshConfig>::compute_Euler_quantity<phys_model>(
      *sol_dofs, *solution, [](const phys_model::Properties &props) { return props.Ma; },
      *Mach_number);

  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *Mach_number, "Ma");
  meshwriter.append_nodal_function_to_file(
      *mesh2D, outfilename, *scheme->vec_function(SolverVecFn::residuals), "residuals");

  mpi_env.finalize();

  return 0;
}
