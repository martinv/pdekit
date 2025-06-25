/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE coarse_scale_correction_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>
#include <map>
#include <memory>

/// PDEKIT headers
#include "common/MPI/MPIEnv.hpp"
#include "graph/GraphReordering.hpp"
#include "interpolation/CoarseScaleCorrectionOpBuilder.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "linear_system/CoarseScaleCorrection.hpp"
#include "linear_system/LSTpetra.hpp"
#include "linear_system/TpetraDofMap.hpp"
#include "math/MathConstants.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

struct MPIFixture
{
  MPIFixture()
  {
    int argc    = boost::unit_test::framework::master_test_suite().argc;
    char **argv = boost::unit_test::framework::master_test_suite().argv;

    common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);
    is_initialized                              = mpi_env.is_initialized();
    is_finalized                                = mpi_env.is_finalized();
  }

  ~MPIFixture()
  {
    common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance();
    mpi_env.finalize();

    std::cout << "MPI environment is finalized (should be true): " << mpi_env.is_finalized()
              << std::endl;
  }

  static bool is_initialized;
  static bool is_finalized;
};

bool MPIFixture::is_initialized = false;
bool MPIFixture::is_finalized   = false;

// ----------------------------------------------------------------------------

struct CoarseScaleCorrection_Fixture
{
  /// common setup for each test case
  CoarseScaleCorrection_Fixture()
  {
  }

  /// common tear-down for each test case
  ~CoarseScaleCorrection_Fixture()
  {
  }

  template <typename MeshConfig>
  void compute_cell_reordering(const mesh::Tria<MeshConfig> &cell_topology,
                               typename pdekit::result_of::dof_map_t<MeshConfig> const &dofs,
                               std::vector<Uint> &cell_reordering,
                               std::vector<Uint> &dof_reordering)
  {
    common::BlockArray<Uint, Uint> mesh_dual_graph_crs;

    cell_topology.build_dual_graph_undirected(mesh_dual_graph_crs);

    // For each old active cell, 'm_cell_reordering' stores the new active
    // index on position i Note that we're actually NOT going to renumber
    // the cells, only change the order in which the cells are assigned
    // their dof numbers

    // Inverse map to the above: for each new active index on position i,
    // store what the old active index was
    std::vector<Uint> cell_reordering_new_to_old;
    graph::GraphReordering::compute_reverse_cuthill_mckee(mesh_dual_graph_crs, cell_reordering);

    // Generate the inverse map
    cell_reordering_new_to_old.resize(cell_reordering.size());

    for (Uint i = 0; i < cell_reordering.size(); ++i)
    {
      cell_reordering_new_to_old[cell_reordering[i]] = i;
    }

    dof_reordering.resize(dofs.nb_nodes());

    Uint dof_idx = 0;
    for (Uint ac = 0; ac < cell_reordering_new_to_old.size(); ++ac)
    {
      const Uint old_ac_id        = cell_reordering_new_to_old[ac];
      const mesh::MeshEntity cell = dofs.active_cell(mesh::ActiveIdx(old_ac_id));
      for (Uint v = 0; v < cell.nb_vert(); ++v)
      {
        dof_reordering[cell.vertex(v)] = dof_idx++;
      }
    }
  }

  void vec_function_to_trilinos_vector(const interpolation::VectorMeshFunction<Real> &vec_function,
                                       ls::TpetraMultiVector<Real> &tpetra_vec)
  {
    using const_function_entry = interpolation::VectorMeshFunction<Real>::const_entry_type;

    const Uint nb_fields = vec_function.nb_fields();

    for (Uint i = 0; i < vec_function.nb_entries(); ++i)
    {
      const const_function_entry entry = vec_function.const_value(i);

      for (Uint f = 0; f < nb_fields; ++f)
      {
        tpetra_vec.insert_value(i * nb_fields + f, entry[f], 0);
      }
    }
  }

  void trilinos_vector_to_vec_function(const ls::TpetraMultiVector<Real> &tpetra_vec,
                                       const Uint nb_fields,
                                       interpolation::VectorMeshFunction<Real> &vec_function)
  {
    using function_entry = interpolation::VectorMeshFunction<Real>::entry_type;

    const Uint nb_entries = tpetra_vec.size() / nb_fields;
    vec_function.resize(nb_fields, nb_entries);

    for (Uint i = 0; i < nb_entries; ++i)
    {
      function_entry entry = vec_function.value(i);

      for (Uint f = 0; f < nb_fields; ++f)
      {
        entry[f] = tpetra_vec.value(i * nb_fields + f, 0);
      }
    }
  }
};

// ----------------------------------------------------------------------------

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_FIXTURE_TEST_SUITE(CoarseScaleCorrection_TestSuite, CoarseScaleCorrection_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(coarse_scale_correction_utest_2D)
{
  BOOST_CHECK_EQUAL(MPIFixture::is_initialized, true);
  BOOST_CHECK_EQUAL(MPIFixture::is_finalized, false);

  const std::string infilename = "unit_square_tri_p2.msh";

  mesh::gmsh::GmshReader gmsh_reader;
  mesh::Tria<Cart2D>::shared_ptr input_mesh = std::make_shared<mesh::Tria<Cart2D>>("mesh");
  gmsh_reader.read_mesh_from_file(infilename, *input_mesh, "geo_dofs");

  common::PtrHandle<mesh::DofMap<Cart2D>> geo_dofs = input_mesh->dof_storage("geo_dofs");

  common::PtrHandle<mesh::DofMap<Cart2D>> sol_dofs_coarse =
      input_mesh->create_dof_storage("sol_dofs_coarse");
  mesh::DofMap<Cart2D>::clone_discontinuous(*input_mesh, *geo_dofs, *sol_dofs_coarse, P1,
                                            PointSetID::Warpblend);

  common::PtrHandle<mesh::DofMap<Cart2D>> sol_dofs_fine =
      input_mesh->create_dof_storage("sol_dofs_fine");
  mesh::DofMap<Cart2D>::clone_discontinuous(*input_mesh, *geo_dofs, *sol_dofs_fine, P3,
                                            PointSetID::Warpblend);

  std::vector<Uint> cell_reordering;
  std::vector<Uint> dof_reordering;

  compute_cell_reordering(*input_mesh, *sol_dofs_fine, cell_reordering, dof_reordering);

  (*sol_dofs_coarse).renumber_dofs_blockwise(cell_reordering);
  (*sol_dofs_fine).renumber_dofs_blockwise(cell_reordering);

  // ---------------------
  // Create Tpetra vectors
  // ---------------------

  const Uint NEQ = 4;

  const Uint nb_sol_dofs_coarse = (*sol_dofs_coarse).nb_nodes();
  const Uint nb_sol_dofs_fine   = (*sol_dofs_fine).nb_nodes();

  std::shared_ptr<ls::TpetraComm<Int>> tpetra_comm =
      std::make_shared<ls::TpetraComm<Int>>(common::mpi::MPIEnv::instance().comm());

  std::shared_ptr<ls::TpetraDofMap<Int>> vec_map_coarse =
      std::make_shared<ls::TpetraDofMap<Int>>(nb_sol_dofs_coarse * NEQ, 0, *tpetra_comm);

  std::shared_ptr<ls::TpetraDofMap<Int>> vec_map_fine =
      std::make_shared<ls::TpetraDofMap<Int>>(nb_sol_dofs_fine * NEQ, 0, *tpetra_comm);

  ls::TpetraMultiVector<Real> tpetra_vec_coarse(*vec_map_coarse, 1);
  ls::TpetraMultiVector<Real> tpetra_vec_fine(*vec_map_fine, 1);

  // --------------------------------------
  // Prepare coarse scale correction object
  // --------------------------------------

  ls::CoarseScaleCorrection correction;
  common::BlockArray<Uint, Uint> mesh_dual_graph_crs;

  input_mesh->build_dual_graph_undirected(mesh_dual_graph_crs);

  std::vector<Uint> block_sizes_coarse((*sol_dofs_coarse).nb_active_cells());
  std::vector<Uint> block_sizes_fine((*sol_dofs_fine).nb_active_cells());

  for (Uint ac = 0; ac < (*sol_dofs_coarse).nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity cell_coarse = (*sol_dofs_coarse).active_cell(mesh::ActiveIdx(ac));
    const mesh::MeshEntity cell_fine   = (*sol_dofs_fine).active_cell(mesh::ActiveIdx(ac));

    block_sizes_coarse[ac] = NEQ * cell_coarse.nb_vert();
    block_sizes_fine[ac]   = NEQ * cell_fine.nb_vert();
  }

  std::unique_ptr<ls::LocalTransferOps<Real>> local_ops(new ls::LocalTransferOps<Real>());
  interpolation::CoarseScaleCorrectionOpBuilder::build_restriction_ops<Cart2D>(
      NEQ, *geo_dofs, *sol_dofs_fine, cell_reordering, P1, *local_ops);

  std::unique_ptr<ls::LocalTransferOps<Real>> prolongation_ops(new ls::LocalTransferOps<Real>());
  interpolation::CoarseScaleCorrectionOpBuilder::build_prolongation_ops<Cart2D>(
      NEQ, *geo_dofs, *sol_dofs_fine, cell_reordering, P1, *prolongation_ops);

  correction.build_block_matrix_sparsity_pattern(mesh_dual_graph_crs, block_sizes_coarse,
                                                 block_sizes_fine, cell_reordering,
                                                 std::move(local_ops), std::move(prolongation_ops));

  interpolation::VectorMeshFunction<Real> solution_coarse("", "solution_coarse");
  solution_coarse.resize(NEQ, nb_sol_dofs_coarse);

  interpolation::VectorMeshFunction<Real> solution_fine("", "solution_fine");
  solution_fine.resize(NEQ, nb_sol_dofs_fine);

  // --------------------------------------
  // Fill the solution field
  // --------------------------------------

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint ac = 0; ac < (*sol_dofs_fine).nb_active_cells(); ++ac)
  {
    const mesh::CellTopologyView<Cart2D> tcell_view = (*sol_dofs_fine).tcell(mesh::ActiveIdx(ac));
    const mesh::MeshEntity cell = (*sol_dofs_fine).active_cell(mesh::ActiveIdx(ac));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      const math::DenseConstVecView<Real> node_coord = cell_coords.row_transpose(n);
      const Real val =
          std::sin(2 * 3.14159 * node_coord[X0]) * std::sin(2 * 3.14159 * node_coord[X1]);

      interpolation::VectorMeshFunction<Real>::entry_type node_values =
          solution_fine.value(cell.vertex(n));

      for (Uint eq = 0; eq < NEQ; ++eq)
      {
        node_values[eq] = (eq + 1) * val;
      }
    }
  }

  // -----------------------------------------
  // Transform the solution to a Tpetra vector
  // Perform a restriction
  // Transform from vector to function
  // -----------------------------------------

  vec_function_to_trilinos_vector(solution_fine, tpetra_vec_fine);

  // Try to restrict the vector
  correction.restrict_vec(tpetra_vec_fine, tpetra_vec_coarse);

  // Prolongate the restriction back
  correction.prolongate_vec(tpetra_vec_coarse, tpetra_vec_fine);

  trilinos_vector_to_vec_function(tpetra_vec_coarse, NEQ, solution_coarse);
  trilinos_vector_to_vec_function(tpetra_vec_fine, NEQ, solution_fine);

  mesh::gmsh::GmshWriter gmsh_writer;
  const std::string outfilename_coarse = "output_coarse_scale_correction_2D_coarse.msh";
  const std::string outfilename_fine   = "output_coarse_scale_correction_2D_fine.msh";

  gmsh_writer.write_mesh_to_file(*input_mesh, "sol_dofs_fine", outfilename_fine);
  gmsh_writer.append_nodal_function_to_file(*input_mesh, outfilename_fine, solution_fine,
                                            "prolongated_solution");

  gmsh_writer.write_mesh_to_file(*input_mesh, "sol_dofs_coarse", outfilename_coarse);
  gmsh_writer.append_nodal_function_to_file(*input_mesh, outfilename_coarse, solution_coarse,
                                            "restricted_solution");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
