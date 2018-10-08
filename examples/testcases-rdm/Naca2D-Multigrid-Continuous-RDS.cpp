#include <ctime>
#include <iostream>
#include <memory>

#include "common/MPI/MPIEnv.hpp"
#include "common/PDEKit.hpp"
#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/multigrid/PMultigridRDM.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::Euler2DCons phys_model;

// ----------------------------------------------------------------------------

void print_norm(const Uint level, const math::DenseDVec<Real> &norm)
{
  std::cout << "   Res [P" << level << "] =";

  for (Uint i = 0; i < norm.size(); ++i)
  {
    std::cout << " " << std::setw(15) << norm[i];
  }

  std::cout << " , log(res) =";

  for (Uint i = 0; i < norm.size(); ++i)
  {
    std::cout << " " << std::setw(15) << std::log(norm[i]);
  }
  std::cout << std::endl;
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // ------------------------------------------------------
  // INITIALIZE ENVIRONMENT
  // ------------------------------------------------------

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // READ MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  const std::string infilename = "naca0012_p2_tri.msh";
  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  typedef rdm::multigrid::PMultigridRDM<MeshConfig, phys_model, rdm::PGLDA> multigrid_type;
  multigrid_type multigrid;

  const std::vector<Uint> levels = {P0, P1, P2};
  multigrid.build_mesh_p_sequence(mesh2D, "geo_dofs", "P2-P1-P0-P1-P2");
  // multigrid.define_nb_smoothing_iters_per_level({ 10, 10, 50, 10, 10 });
  // multigrid.define_nb_smoothing_iters_per_level({ 3, 3, 10, 3, 3 });
  multigrid.define_nb_smoothing_iters_per_level({10, 10, 10, 10, 10});

  math::DenseDVec<Real> ic_state(phys_model::NEQ);

  ic_state[0] = 1.22503;
  ic_state[1] = 208.30559;
  ic_state[2] = 7.27419;
  ic_state[3] = 271044.38;

  for (Uint l = 0; l < multigrid.nb_levels(); ++l)
  {
    // const std::shared_ptr<multigrid_type::bc_base_type> bc_farfield =

    multigrid.add_boundary_condition("WeakWall", "wall_bc", "Wall", levels[l]);
    std::shared_ptr<multigrid_type::bc_base_type> bc_farfield =
        multigrid.add_boundary_condition("WeakFarfield", "farfield_bc", "Farfield", levels[l]);
    bc_farfield->set_reference_state(ic_state);
  }

  multigrid.apply_initial_condition("InnerCells", ic_state);

  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);
  std::cout.precision(10);

  clock_t start, end;
  Real elapsed;
  start = clock();

  // multigrid.mg_startup();

  for (Uint i = 0; i < 20; ++i)
  {
    std::cout << "Cycle = " << i << std::endl;
    multigrid.make_cycle();

    for (Uint l = 0; l < levels.size(); ++l)
    {
      multigrid.compute_residual_norm(levels[l], res_L2_norm);
      print_norm(levels[l], res_L2_norm);
    }
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  multigrid.write_sequence_to_files();

  mpi_env.finalize();

  return 0;
}
