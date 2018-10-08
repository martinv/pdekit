#include <chrono>
#include <iostream>
#include <memory>

#include "common/MPI/MPIEnv.hpp"
#include "common/PDEKit.hpp"
#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "physics/scalar/RotationAdvection2D.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/multigrid/PMultigridRDM.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;

int main(int argc, char *argv[])
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh_support");

  const std::string infilename = "advection-tri.msh";
  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  typedef physics::RotationAdvection2D<physics::AroundOrigin2D> phys_model;
  typedef rdm::multigrid::PMultigridRDM<MeshConfig, phys_model, rdm::PGLDA> multigrid_type;
  multigrid_type multigrid;

  // multigrid.build_mesh_p_sequence(mesh2D, levels);
  multigrid.build_mesh_p_sequence(mesh2D, "geo_dofs", "P2-P1-P0-P1-P2");
  multigrid.define_nb_smoothing_iters_per_level({10, 10, 50, 10, 10});

  const std::vector<Uint> levels = {P0, P1, P2};

  for (Uint l = 0; l < multigrid.nb_levels(); ++l)
  {
    std::shared_ptr<multigrid_type::bc_base_type> farfield_bc = multigrid.add_boundary_condition(
        "StrongDirichlet", "strong_farfield", "farfield", levels[l]);
    farfield_bc->set_expression(&Zero::value);

    const std::shared_ptr<multigrid_type::bc_base_type> inlet_bc =
        multigrid.add_boundary_condition("StrongDirichlet", "inlet_bc", "inlet", levels[l]);
    inlet_bc->set_expression(&CosineHat2D::value);
  }

  multigrid.apply_initial_condition("fluid", &CosineHat2D::value);

  math::DenseDVec<Real> res_L1_norm(1);
  std::cout.precision(15);

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::milliseconds elapsed;
  start = std::chrono::high_resolution_clock::now();

  // multigrid.mg_startup();

  for (Uint i = 0; i < 20; ++i)
  {
    std::cout << "Cycle = " << i << std::endl;
    multigrid.make_cycle();

    for (Uint l = 0; l < levels.size(); ++l)
    {
      multigrid.compute_residual_norm(levels[l], res_L1_norm);
      std::cout << "   Res [P" << levels[l] << "] = " << std::setw(20) << res_L1_norm[0]
                << " , log(res) = " << std::setw(20) << std::log(res_L1_norm[0]) << std::endl;
    }
  }

  /*
  std::cout << "*****************************************" << std::endl;
  multigrid.smooth_single_level(P0, 1000);

  for (Uint l = 0; l < levels.size(); ++l)
  {
    multigrid.compute_residual_norm(levels[l], res_L1_norm);
    std::cout << "   Res [P" << levels[l] << "] = " << std::setw(20) <<
  res_L1_norm[0]
              << " , log(res) = " << std::setw(20) << std::log(res_L1_norm[0])
  << std::endl;
  }
  std::cout << "*****************************************" << std::endl;
  */

  end     = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  std::cout.setf(std::ios::fixed);
  std::cout.precision(7);
  std::cout << "CPU time (main iteration loop) = " << elapsed.count() << " ms" << std::endl;

  multigrid.write_sequence_to_files();

  mpi_env.finalize();

  return 0;
}
