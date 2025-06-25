#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "examples/testcases-rdm/PredefinedFunctions.hpp"
#include "interpolation/ErrorEstimator.hpp"
#include "interpolation/L2ProjectionGlobal.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "physics/scalar/RotationAdvection2D.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/SolverIO.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/RDSolver.hpp"
#include "solver/rdm/cellsplitters/CellSplitters.hpp"
#include "solver/time/ExplicitTimeStepper.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef interpolation::ScalarMeshFunction<Real> scalar_function;
typedef interpolation::VectorMeshFunction<Real> vector_function;
typedef physics::RotationAdvection2D<physics::AroundOrigin2D> phys_model;

// ----------------------------------------------------------------------------

const SFunc sol_basis = SFunc::Carnevali;

struct FwdTransformRule
{
  // Return: a tuple consisting of
  // 1) Shape function tag of source shape function (which we're projecting
  // FROM) 2) Shape function tag of target polynomial (i.e space we're
  // projecting TO)
  std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> sf_tags(const mesh::PointSetTag &geo_tag,
                                                       const mesh::PointSetTag &sol_tag_src,
                                                       const mesh::PointSetTag &sol_tag_tgt) const
  {
    const mesh::sf::SFTag src_poly_space_tag(sol_tag_src.elem_shape(), SFunc::Lagrange,
                                             sol_tag_src.poly_order(), ModalBasis::Modal);

    const mesh::sf::SFTag tgt_poly_space_tag(sol_tag_tgt.elem_shape(), sol_basis,
                                             sol_tag_tgt.poly_order(), ModalBasis::Modal);

    const std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> result =
        std::make_tuple(src_poly_space_tag, tgt_poly_space_tag);
    return result;
  }

  // Return: a tuple consisting of
  // 1) Quadrature order to be used in order to construct LHS matrix operator
  // 2) Quadrature order to be used in order to construct RHS
  std::tuple<Uint, Uint> quadrature_orders(const mesh::PointSetTag &geo_tag,
                                           const mesh::PointSetTag &sol_tag_src,
                                           const mesh::PointSetTag &sol_tag_tgt) const
  {
    const Uint src_quad_order =
        std::max(sol_tag_src.poly_order() + sol_tag_tgt.poly_order(), geo_tag.poly_order());
    const Uint tgt_quad_order = std::max(2 * sol_tag_tgt.poly_order(), geo_tag.poly_order());

    const std::tuple<Uint, Uint> result = std::make_tuple(src_quad_order, tgt_quad_order);
    return result;
  }
};

// ----------------------------------------------------------------------------

struct BwdTransformRule
{
  // Return: a tuple consisting of
  // 1) Shape function tag of source shape function (which we're projecting
  // FROM) 2) Shape function tag of target polynomial (i.e space we're
  // projecting TO)
  std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> sf_tags(const mesh::PointSetTag &geo_tag,
                                                       const mesh::PointSetTag &sol_tag_src,
                                                       const mesh::PointSetTag &sol_tag_tgt) const
  {
    const mesh::sf::SFTag src_poly_space_tag(sol_tag_src.elem_shape(), sol_basis,
                                             sol_tag_src.poly_order(), ModalBasis::Modal);

    const mesh::sf::SFTag tgt_poly_space_tag(sol_tag_tgt.elem_shape(), SFunc::Lagrange,
                                             sol_tag_tgt.poly_order(), ModalBasis::Modal);

    const std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> result =
        std::make_tuple(src_poly_space_tag, tgt_poly_space_tag);
    return result;
  }

  // Return: a tuple consisting of
  // 1) Quadrature order to be used in order to construct LHS matrix operator
  // 2) Quadrature order to be used in order to construct RHS
  std::tuple<Uint, Uint> quadrature_orders(const mesh::PointSetTag &geo_tag,
                                           const mesh::PointSetTag &sol_tag_src,
                                           const mesh::PointSetTag &sol_tag_tgt) const
  {
    const Uint src_quad_order =
        std::max(sol_tag_src.poly_order() + sol_tag_tgt.poly_order(), geo_tag.poly_order());
    const Uint tgt_quad_order = std::max(2 * sol_tag_tgt.poly_order(), geo_tag.poly_order());

    const std::tuple<Uint, Uint> result = std::make_tuple(src_quad_order, tgt_quad_order);
    return result;
  }
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void resample_solution_for_postprocessing(const mesh::Tria<MeshConfig> &mesh_in,
                                          const DofMap<MeshConfig> &dofs, const SFunc sfunc,
                                          interpolation::VectorMeshFunction<Real> &solution)
{
  const Uint nb_fields = solution.nb_fields();

  interpolation::VectorMeshFunction<Real> tmp_solution("", "tmp");
  tmp_solution.resize(nb_fields, solution.nb_entries());
  using Vandermonde_map_t = common::DataMap<mesh::PointSetTag, math::DenseDMat<Real>>;
  Vandermonde_map_t map_V_back_transf;

  std::vector<Real> tmp_values_in;
  std::vector<Real> tmp_values_out;

  for (typename DofMap<MeshConfig>::const_dof_iterator dof_it = dofs.cbegin();
       dof_it != dofs.cend(); ++dof_it)
  {
    const mesh::MeshEntity entity       = dof_it->mesh_entity();
    const mesh::PointSetTag std_reg_tag = entity.pt_set_id();

    common::PtrHandle<math::DenseDMat<Real>> V_handle =
        map_V_back_transf.std_region_data(std_reg_tag);

    if (V_handle.is_null())
    {
      V_handle = map_V_back_transf.create(std_reg_tag);

      mesh::sf::SFTag sf_tag(std_reg_tag.elem_shape(), SFunc::Carnevali, std_reg_tag.poly_order(),
                             ModalBasis::Modal);
      sf::ShapeFunction sfunc;
      sfunc.change_type(std_reg_tag, sf_tag);
      mesh::StdRegion std_reg(std_reg_tag);

      sfunc.get().compute_ref_values(std_reg.get().coordinates(), *V_handle);
    }

    const Uint nb_vert = entity.nb_vert();

    tmp_values_in.resize(nb_vert * nb_fields);
    tmp_values_out.resize(nb_vert * nb_fields);
    math::DenseMatView<Real> tmp_view_in(tmp_values_in.data(), nb_fields, nb_vert, nb_fields);
    math::DenseMatView<Real> tmp_view_out(tmp_values_out.data(), nb_fields, nb_vert, nb_fields);

    for (Uint v = 0; v < nb_vert; ++v)
    {
      const interpolation::VectorMeshFunction<Real>::const_entry_type dof_value =
          solution.const_value(entity.vertex(v));
      for (Uint f = 0; f < nb_fields; ++f)
      {
        tmp_view_in(v, f) = dof_value[f];
      }
    }

    tmp_view_out = (*V_handle) * tmp_view_in;

    for (Uint v = 0; v < nb_vert; ++v)
    {
      interpolation::VectorMeshFunction<Real>::entry_type dof_value =
          tmp_solution.value(entity.vertex(v));
      for (Uint f = 0; f < nb_fields; ++f)
      {
        dof_value[f] = tmp_view_out(v, f);
      }
    }
  } // Loop over all mesh entities

  solution = tmp_solution;
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh2D = std::make_shared<MeshType>("mesh2D");

  const std::string infilename = "advection-tri.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2D, "geo_dofs");

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2D->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh2D->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_continuous(*mesh2D, *geo_dofs, *sol_dofs, P1,
                                               PointSetID::Warpblend);

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

  typedef rdm::RDSolver<MeshConfig, phys_model> scheme_type;

  std::shared_ptr<scheme_type> scheme = std::make_shared<scheme_type>();
  // scheme_type scheme;

  /*
  std::shared_ptr<rdm::RDBlendingCoeff<MeshConfig>> blending_coeff =
      std::make_shared<rdm::RDBlendingCoeff<MeshConfig>>();
  blending_coeff->setup(topology.dof_storage(), sol_topology.dof_storage(),
  "LF_Blend");
  */

  scheme->configure_mesh_data(mesh2D, geo_dofs, sol_dofs);
  scheme->initialize_work_data(sol_basis, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::residuals, residual);

  time::ExplicitEuler time_stepper;
  time_stepper.setup(phys_model::NEQ, nb_nodes);

  // ------------------------------------------------------
  // Prepare transformation operators
  // ------------------------------------------------------

  common::PtrHandle<MeshType::dof_storage_type> sol_dofs_nodal =
      mesh2D->create_dof_storage("sol_dofs_nodal");

  MeshType::dof_storage_type::make_identical_copy(*sol_dofs, *sol_dofs_nodal);

  interpolation::L2ProjectionGlobal forward_trans;
  forward_trans.build_projection_operator<MeshConfig>(*mesh2D, *sol_dofs_nodal, *sol_dofs,
                                                      FwdTransformRule{});

  interpolation::L2ProjectionGlobal backward_trans;
  backward_trans.build_projection_operator<MeshConfig>(*mesh2D, *sol_dofs, *sol_dofs_nodal,
                                                       BwdTransformRule{});

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  vector_function::ptr tmp = std::make_shared<vector_function>("", "tmp");
  tmp->resize(phys_model::NEQ, nb_nodes);

  tmp->fill(0.0);
  solver::InitialCondition<MeshConfig> ic("hat_initial_condition", mesh2D);

  ic.set_domain(_2D, "fluid");
  ic.set_expression(&CosineHat2D::value);
  ic.apply_function("sol_dofs", *tmp);

  forward_trans.project<MeshConfig>(*mesh2D, *sol_dofs_nodal, *sol_dofs, FwdTransformRule{}, *tmp,
                                    *solution);

  // (*solution) = (*tmp);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  std::shared_ptr<scheme_type::bc_base_type> farfield =
      scheme->add_boundary_condition("WeakDirichlet", "weak_farfield", "farfield");
  farfield->set_expression(&Zero::value);

  std::shared_ptr<scheme_type::bc_base_type> strong_inlet =
      scheme->add_boundary_condition("WeakDirichlet", "weak_inlet", "inlet");
  strong_inlet->set_expression(&CosineHat2D::value);

  /*
  std::shared_ptr<scheme_type::bc_base_type> weak_inlet =
      scheme->add_boundary_condition("WeakDirichlet", "weak_inlet", "inlet");
  weak_inlet->set_expression(&CosineHat2D::value);
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
    const Real CFL = 0.5;

    /*
    residual->fill(0.0);
    time_update->reset_wave_speeds();

    // weak_inlet_bc.apply();
    // strong_inlet_bc.apply();
    // farfield.apply();

    scheme->iterate_by_std_region_type();
    scheme->apply_boundary_conditions();

    time_update->compute_local_time_step(0.5);
    const scalar_function &dt = time_update->time_update_function();
    (*solution) -= dt * (*residual);
    */

    // Blending coefficient for nonlinear schemes
    /*
    blending_coeff->calculate_blending_coeff(topology.dof_storage(),
    sol_topology.dof_storage(), *solution);
    */

    // time_stepper.advance_in_time(scheme, 0.5);

    scheme->assemble_lhs_and_rhs(CFL);
    scheme->solve({});

    norm_L1(*residual, res_L1_norm);

    solver::SolverIO::print_iter_and_res_norm(iter, CFL, res_L1_norm);
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  // --------------------------------------------

  // Transform from modal to nodal representation
  (*tmp) = (*solution);

  backward_trans.project<MeshConfig>(*mesh2D, *sol_dofs, *sol_dofs_nodal, BwdTransformRule{}, *tmp,
                                     *solution);

  /*
  resample_solution_for_postprocessing(*mesh2D, *sol_dofs, SFunc::Carnevali,
  *solution);
  */

  // --------------------------------------------

  /// Compute the error in the solution
  Rotation2DSolution exact_solution;
  interpolation::ErrorEstimator<MeshConfig> error_estim(mesh2D);
  error_estim.compute_pointwise_L2_error(exact_solution, *sol_dofs, *solution);
  error_estim.compute_L1_error(exact_solution, *mesh2D, *sol_dofs, *solution);
  error_estim.compute_L2_error(exact_solution, *mesh2D, *sol_dofs, *solution);
  error_estim.compute_infty_error(exact_solution, *mesh2D, *sol_dofs, *solution);

  /// Write the output to file
  const std::string outfilename = "output_advection_p1_2D_spectral.msh";
  gmsh::GmshWriter meshwriter;
  meshwriter.write_mesh_to_file(*mesh2D, "sol_dofs", outfilename);
  meshwriter.append_nodal_function_to_file(*mesh2D, outfilename, *solution, "solution");
  meshwriter.append_nodal_function_to_file(
      *mesh2D, outfilename, *scheme->scal_function(SolverScalFn::time_step), "time_step");

  vtk::VtkWriter vtkwriter;
  vtkwriter.append_nodal_function_to_file(*mesh2D, *sol_dofs, "output_advection_p1_2D_spectral.vtu",
                                          *solution, "solution");

  mpi_env.finalize();
  return 0;
}
