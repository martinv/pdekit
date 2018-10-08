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
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/scalar/RotationAdvection2D.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/DRDMethodImplicit.hpp"
#include "solver/rdm/PGRDMethodImplicit.hpp"
#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/bc/WeakDirichletBC.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/rdm/facetsplitters//FacetSplitters.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

// ----------------------------------------------------------------------------

class ComputeRotVector2D
{
  public:
  template <typename V1, typename V2>
  inline static void rotation_vector(const V1 &direction, V2 &rot_vector)
  {
    static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                  "'coord' must be a vector");
    static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                  "'coord' must be a vector");
    rot_vector[X0] = direction[X1];
    rot_vector[X1] = 1. - direction[X0];
  }
};

// ----------------------------------------------------------------------------

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::RotationAdvection2D<ComputeRotVector2D> phys_model;

int main(int argc, char *argv[])
{
  // ------------------------------------------------------
  // INITIALIZE ENVIRONMENT
  // ------------------------------------------------------

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // SET SCHEME TYPE
  // ------------------------------------------------------

  typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;
  // typedef rdm::DRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA,
  // rdm::FacetDG> scheme_type;
  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme = std::make_shared<scheme_type>();

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2d");

  const std::string infilename = "unit_square_tri_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  /*
  mesh2D->change_std_region_types([](const PointSetTag &tag_old) ->
  PointSetTag { return PointSetTag(tag_old.elem_shape(), P2,
  tag_old.ref_topology());
  });
  */

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  if (scheme->is_continuous())
  {
    MeshType::dof_storage_type::clone_continuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                                 PointSetID::Warpblend);
  }
  else
  {
    MeshType::dof_storage_type::clone_discontinuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                                    PointSetID::Warpblend);
  }

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P3;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  vector_function::ptr solution     = std::make_shared<vector_function>("", "solution");
  vector_function::ptr source_term  = std::make_shared<vector_function>("", "source_term");
  vector_function::ptr residual     = std::make_shared<vector_function>("", "residual");
  scalar_function::ptr update_coeff = std::make_shared<scalar_function>("", "update_coeff");
  // scalar_function::ptr blend_coeff_func =
  // std::make_shared<scalar_function>("", "blending_coeff");

  const Uint nb_nodes = (*sol_dofs).nb_nodes();

  solution->resize(phys_model::NEQ, nb_nodes);
  source_term->resize(phys_model::NEQ, nb_nodes);
  residual->resize(phys_model::NEQ, nb_nodes);
  update_coeff->resize(nb_nodes);
  // blend_coeff_func->resize(sol_topology.dof_storage().nb_cells());

  // Choose scheme, configure the function spaces (solution, residual,
  // update coefficients

  for (Uint c = 0; c < (*sol_dofs).nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = mesh2D->active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity active_cell             = (*sol_dofs).active_cell(mesh::ActiveIdx(c));
    const mesh::CellGeometry<MeshConfig::GDIM> cell_coords = tcell.coordinates();

    for (Uint n = 0; n < active_cell.nb_vert(); ++n)
    {
      const math::DenseVecView<const Real> node_coord = cell_coords.const_node_view(n);

      const Real dist = std::sqrt((node_coord[X0] - 0.5) * (node_coord[X0] - 0.5) +
                                  (node_coord[X1] - 0.3) * (node_coord[X1] - 0.3));

      interpolation::VectorMeshFunction<Real>::entry_type source_value =
          source_term->value(active_cell.vertex(n));

      if (dist <= 0.25)
      {
        source_value[0] = 10.0;
      }
      else
      {
        source_value[0] = 0.0;
      }
    }
  }

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  /*
  rdm::RDBlendingCoeff<MeshConfig> rd_blending_coeff_operator;
  rd_blending_coeff_operator.setup(geo_topology.dof_storage(),
  sol_topology.dof_storage(), "LF_Blend");
  */

  scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::sources, source_term);
  scheme->set_vec_function(SolverVecFn::residuals, residual);
  // scheme->set_blending_coeff(blend_coeff_func);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);
  solver::InitialCondition<MeshConfig> ic("zero_initial_condition", mesh2D);

  ic.set_domain(_2D, "interior");
  ic.set_expression(&Zero::value);
  ic.apply_function("sol_dofs", *solution);
  // ic.set_domain(_1D,"inlet");
  // ic.apply_function(hat,*solution);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  /*
  std::shared_ptr<scheme_type::bc_base_type> strong_bc_left =
      scheme->add_boundary_condition("StrongDirichlet", "strong_bc_left",
  "left"); strong_bc_left->set_expression(&Zero::value);
  */

  std::shared_ptr<scheme_type::bc_base_type> strong_inlet_bc =
      scheme->add_boundary_condition("StrongDirichlet", "strong_inlet", "bottom");
  strong_inlet_bc->set_expression(&InletJump2D::value);

  /*
  std::shared_ptr<scheme_type::bc_base_type> weak_inlet_bc =
      scheme->add_boundary_condition("WeakDirichlet", "weak_inlet", "bottom");
  weak_inlet_bc->set_expression(&InletJump2D::value);
  */

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;
  std::cout.precision(15);

  math::DenseDVec<Real> res_L1_norm(1);

  start = clock();
  for (Uint iter = 0; iter < 150; ++iter)
  {
    // Blending coefficient for nonlinear schemes
    /*
    rd_blending_coeff_operator.calculate(geo_topology.dof_storage(),
    sol_topology.dof_storage(), *solution, *blend_coeff_func); const Real dt
    = 3.0;
    */

    const Real CFL = iter < 10 ? 3.0 : std::min(5.0, 4.0 + std::pow(1.30, iter - 10));

    scheme->set_cfl(CFL);
    scheme->assemble_lhs_and_rhs(CFL);
    scheme->solve({SolverOption::RecomputePreconditioner});
    scheme->compute_residual_norm(res_L1_norm);

    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L1_norm);
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Write the output to file

  const std::string outfilename = "output_src_term_continuous_advection_p1_2D.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "solution");
  /*
  meshwriter.append_cell_function_to_file(*mesh2D, outfilename,
  *blend_coeff_func, "blending_coeff");
  */

  /*
  vtk::VtkWriter vtkwriter;
  vtkwriter.append_nodal_function_to_file(*sol_mesh, "output.vtu", *solution,
  "solution");
  */

  typedef mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_dof_iterator bdry_iter_type;

  const mesh::MeshBoundarySet<MeshConfig> &boundaries = mesh2D->all_boundaries();
  std::shared_ptr<mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> outlet_bdry =
      boundaries.domain("right");

  const bdry_iter_type outlet_begin = outlet_bdry->cbegin(*sol_dofs);
  const bdry_iter_type outlet_end   = outlet_bdry->cend(*sol_dofs);

  // Count nodes on boundary
  Uint nb_nodes_on_bdry = 0;
  for (bdry_iter_type bdry_it = outlet_begin; bdry_it != outlet_end; ++bdry_it)
  {
    const mesh::MeshEntity bdry_cell = bdry_it->mesh_entity();
    nb_nodes_on_bdry += bdry_cell.nb_vert();
  }

  std::vector<std::pair<Real, Real>> bdry_values;
  bdry_values.reserve(nb_nodes_on_bdry);
  bdry_values.resize(0);

  for (bdry_iter_type bdry_it = outlet_begin; bdry_it != outlet_end; ++bdry_it)
  {
    const mesh::PointSetTag geo_pt_set_id = bdry_it->geo_pt_set_id();
    const mesh::MeshEntity bdry_cell      = bdry_it->mesh_entity();
    const mesh::PointSetTag sol_pt_set_id = bdry_cell.pt_set_id();
    const auto facet_geo_coords           = bdry_it->cell_geometry();

    for (Uint v = 0; v < bdry_cell.nb_vert(); ++v)
    {
      const interpolation::VectorMeshFunction<Real>::const_entry_type sol_in_node =
          (*solution).const_value(bdry_cell.vertex(v));
      bdry_values.push_back(
          std::make_pair(facet_geo_coords.const_node_view(v)[X1], sol_in_node[0]));
    }
  }

  // Sort the values on boundary according to the first entry of each pair
  // (i.e. position on the outlet boundary)
  std::sort(bdry_values.begin(), bdry_values.end());

  // Save the solution on outlet boundary to a text file
  std::ofstream out_stream;
  out_stream.setf(std::ios::fixed);
  out_stream.precision(15);
  out_stream.open("rdm_src_adv_outlet_solution.txt");

  for (Uint i = 0; i < bdry_values.size(); ++i)
  {
    out_stream << std::setw(15) << bdry_values[i].first << " " << std::setw(15)
               << bdry_values[i].second << std::endl;
  }

  out_stream.close();

  mpi_env.finalize();

  return 0;
}
