/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Euler2D_test
#include <boost/test/unit_test.hpp>

#include <ctime>
#include <iostream>
#include <memory>

#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "math/MathConstants.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/MeshManipulator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/Postprocessing.hpp"
#include "solver/rdm/PGRDMethodExplicit.hpp"
#include "solver/rdm/PGRDMethodImplicit.hpp"
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

void print_res_norm(const Uint iter, const Real CFL, const math::DenseDVec<Real> &norm)
{

  std::cout << "Iter = " << std::setw(5) << iter << "  CFL = " << std::setw(5)
            << std::setprecision(3) << CFL << " , " << std::setprecision(15) << " res =";
  for (Uint i = 0; i < norm.size(); ++i)
  {
    std::cout << " " << std::setw(20) << norm[i];
  }

  std::cout << " , log(res) =";
  for (Uint i = 0; i < norm.size(); ++i)
  {
    std::cout << " " << std::setw(20) << std::log(norm[i]);
  }
  std::cout << std::endl;
}

// ----------------------------------------------------------------------------

class InitCondMSEuler2D
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component)
  {
    const Real rho_init = 1.0;
    const Real u_init   = 77.2;
    const Real v_init   = 82.6;
    const Real p_init   = 2.94e5;

    Real result = 0.0;

    switch (component)
    {
      case 0:
        result = rho_init;
        break;
      case 1:
        result = rho_init * u_init;
        break;
      case 2:
        result = rho_init * v_init;
        break;
      case 3:
        const Real e_init =
            p_init / (1.4 - 1.) + 0.5 * rho_init * (u_init * u_init + v_init * v_init);
        result = e_init;
        break;
    };

    return result;
  }
};

// ----------------------------------------------------------------------------

class MSEuler2D
{
  public:
  // Constructor
  MSEuler2D();

  // Destructor
  ~MSEuler2D();

  // Given point coordinates, calculate the values of exact solution
  static void fill_exact_solution_values(const pdekit::result_of::dof_map_t<MeshConfig> &cell_dofs,
                                         vector_function &exact_solution);

  // Given point coordinates, calculate the values of exact solution
  static void fill_source_term_values(const pdekit::result_of::dof_map_t<MeshConfig> &cell_dofs,
                                      vector_function &source_term);

  static Real Dirichlet_bc_value(
      const math::DenseConstVecView<Real> &node,
      const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
      const Uint component);

  private:
  // Characteristic length
  static const Real L;
  // Ratio of specific heats
  static const Real gamma;

  // -------
  // Density
  // -------
  // Coefficients rho_x, rho_y
  static const Real rho0;
  static const Real coeff_rho[2];
  static const Real a_rho[2];

  // ------------
  // Velocity - x
  // ------------
  static const Real u0;
  static const Real coeff_u[2];
  static const Real a_u[2];

  // ------------
  // Velocity - y
  // ------------
  static const Real v0;
  static const Real coeff_v[2];
  static const Real a_v[2];

  // --------
  // Pressure
  // --------
  static const Real p0;
  static const Real coeff_p[2];
  static const Real a_p[2];
};

// ----------------------------------------------------------------------------
// Static variables
// ----------------------------------------------------------------------------

const Real MSEuler2D::L     = 1.0;
const Real MSEuler2D::gamma = 1.4;

const Real MSEuler2D::rho0         = 1.0;
const Real MSEuler2D::coeff_rho[2] = {0.15, -0.1};
const Real MSEuler2D::a_rho[2]     = {1.0, 0.5};

const Real MSEuler2D::u0         = 70.0;
const Real MSEuler2D::coeff_u[2] = {5.0, -7.0};
const Real MSEuler2D::a_u[2]     = {1.5, 0.6};

const Real MSEuler2D::v0         = 90.0;
const Real MSEuler2D::coeff_v[2] = {-15.0, 8.5};
const Real MSEuler2D::a_v[2]     = {0.5, 2. / 3.};

const Real MSEuler2D::p0         = 1.e5;
const Real MSEuler2D::coeff_p[2] = {0.2 * 1.e5, 0.5 * 1.e5};
const Real MSEuler2D::a_p[2]     = {2.0, 1.0};

// ----------------------------------------------------------------------------

void MSEuler2D::fill_exact_solution_values(
    const pdekit::result_of::dof_map_t<MeshConfig> &cell_dofs, vector_function &exact_solution)
{
  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell                         = cell_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      const math::DenseConstVecView<Real> node = cell_coords.row_transpose(n);

      vector_function::entry_type sol_at_pt = exact_solution.value(cell.vertex(n));

      const Real rho = coeff_rho[X0] * std::sin(a_rho[X0] * math::pi * node[X0] / L) +
                       coeff_rho[X1] * std::cos(a_rho[X1] * math::pi * node[X1] / L) + rho0;

      const Real u = coeff_u[X0] * std::sin(a_u[X0] * math::pi * node[X0] / L) +
                     coeff_u[X1] * std::cos(a_u[X1] * math::pi * node[X1] / L) + u0;

      const Real v = coeff_v[X0] * std::cos(a_v[X0] * math::pi * node[X0] / L) +
                     coeff_v[X1] * std::sin(a_v[X1] * math::pi * node[X1] / L) + v0;

      const Real p = coeff_p[X0] * std::cos(a_p[X0] * math::pi * node[X0] / L) +
                     coeff_p[X1] * std::sin(a_p[X1] * math::pi * node[X1] / L) + p0;

      const Real e = p / (gamma - 1.) + 0.5 * rho * (u * u + v * v);

      sol_at_pt[0] = rho;
      sol_at_pt[1] = rho * u;
      sol_at_pt[2] = rho * v;
      sol_at_pt[3] = e;
    }
  }
}

// ----------------------------------------------------------------------------

void MSEuler2D::fill_source_term_values(const pdekit::result_of::dof_map_t<MeshConfig> &cell_dofs,
                                        vector_function &source_term)
{
  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell                         = cell_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      const math::DenseConstVecView<Real> node = cell_coords.row_transpose(n);

      vector_function::entry_type src_at_pt = source_term.value(cell.vertex(n));

      // --------------------------------

      const Real rho = coeff_rho[X0] * std::sin(a_rho[X0] * math::pi * node[X0] / L) +
                       coeff_rho[X1] * std::cos(a_rho[X1] * math::pi * node[X1] / L) + rho0;

      const Real drho_dx = (a_rho[X0] * math::pi * coeff_rho[X0] / L) *
                           std::cos(a_rho[X0] * math::pi * node[X0] / L);

      const Real drho_dy = -(a_rho[X1] * math::pi * coeff_rho[X1] / L) *
                           std::sin(a_rho[X1] * math::pi * node[X1] / L);

      // --------------------------------

      const Real u = coeff_u[X0] * std::sin(a_u[X0] * math::pi * node[X0] / L) +
                     coeff_u[X1] * std::cos(a_u[X1] * math::pi * node[X1] / L) + u0;

      const Real du_dx =
          (a_u[X0] * math::pi * coeff_u[X0] / L) * std::cos(a_u[X0] * math::pi * node[X0] / L);
      const Real du_dy =
          -(a_u[X1] * math::pi * coeff_u[X1] / L) * std::sin(a_u[X1] * math::pi * node[X1] / L);

      // --------------------------------

      const Real v = coeff_v[X0] * std::cos(a_v[X0] * math::pi * node[X0] / L) +
                     coeff_v[X1] * std::sin(a_v[X1] * math::pi * node[X1] / L) + v0;

      const Real dv_dx =
          -(a_v[X0] * math::pi * coeff_v[X0] / L) * std::sin(a_v[X0] * math::pi * node[X0] / L);

      const Real dv_dy =
          (a_v[X1] * math::pi * coeff_v[X1] / L) * std::cos(a_v[X1] * math::pi * node[X1] / L);

      // --------------------------------

      const Real p = coeff_p[X0] * std::cos(a_p[X0] * math::pi * node[X0] / L) +
                     coeff_p[X1] * std::sin(a_p[X1] * math::pi * node[X1] / L) + p0;

      const Real dp_dx =
          -(a_p[X0] * math::pi * coeff_p[X0] / L) * std::sin(a_p[X0] * math::pi * node[X0] / L);

      const Real dp_dy =
          (a_p[X1] * math::pi * coeff_p[X1] / L) * std::cos(a_p[X1] * math::pi * node[X1] / L);

      // --------------------------------

      const Real e = p / (gamma - 1.) + 0.5 * rho * (u * u + v * v);

      const Real de_dx =
          dp_dx / (gamma - 1.) + 0.5 * drho_dx * (u * u + v * v) + rho * (u * du_dx + v * dv_dx);
      const Real de_dy =
          dp_dy / (gamma - 1.) + 0.5 * drho_dy * (u * u + v * v) + rho * (u * du_dy + v * dv_dy);

      // --------------------------------

      src_at_pt[0] = rho * (du_dx + dv_dy) + u * drho_dx + v * drho_dy;
      src_at_pt[1] = (drho_dx * u * u + 2.0 * rho * u * du_dx + dp_dx) +
                     (drho_dy * u * v + rho * du_dy * v + rho * u * dv_dy);
      src_at_pt[2] = (drho_dx * u * v + rho * du_dx * v + rho * u * dv_dx) +
                     (drho_dy * v * v + 2.0 * rho * v * dv_dy + dp_dy);
      src_at_pt[3] =
          ((de_dx + dp_dx) * u + (e + p) * du_dx) + ((de_dy + dp_dy) * v + (e + p) * dv_dy);
    }
  }
}

// ----------------------------------------------------------------------------

Real MSEuler2D::Dirichlet_bc_value(
    const math::DenseConstVecView<Real> &node,
    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution, const Uint component)
{

  const Real rho = coeff_rho[X0] * std::sin(a_rho[X0] * math::pi * node[X0] / L) +
                   coeff_rho[X1] * std::cos(a_rho[X1] * math::pi * node[X1] / L) + rho0;

  const Real u = coeff_u[X0] * std::sin(a_u[X0] * math::pi * node[X0] / L) +
                 coeff_u[X1] * std::cos(a_u[X1] * math::pi * node[X1] / L) + u0;

  const Real v = coeff_v[X0] * std::cos(a_v[X0] * math::pi * node[X0] / L) +
                 coeff_v[X1] * std::sin(a_v[X1] * math::pi * node[X1] / L) + v0;

  const Real p = coeff_p[X0] * std::cos(a_p[X0] * math::pi * node[X0] / L) +
                 coeff_p[X1] * std::sin(a_p[X1] * math::pi * node[X1] / L) + p0;

  Real result = 0.0;

  switch (component)
  {

    case 0:
      result = rho;
      break;

    case 1:
      result = rho * u;
      break;

    case 2:
      result = rho * v;
      break;

    case 3:
      result = p / (gamma - 1.) + 0.5 * rho * (u * u + v * v);
      break;
  };

  return result;
}

// ----------------------------------------------------------------------------

std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> build_solver(
    const std::shared_ptr<MeshType const> &mesh,
    common::PtrHandle<MeshType::dof_storage_type> geo_dofs,
    common::PtrHandle<MeshType::dof_storage_type> sol_dofs, const Uint quad_order,
    vector_function::ptr solution, vector_function::ptr sources, vector_function::ptr residuals,
    const std::string &scheme_type, const bool reorder_dofs)
{
  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  // typedef rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA>
  // scheme_type;

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> scheme;

  if (scheme_type == "LDA")
  {
    scheme = std::make_shared<rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGLDA>>();
    // scheme = std::make_shared<rdm::DRDMethodImplicit<MeshConfig,
    // phys_model, rdm::PGLDA, rdm::FacetDG>>();
  }
  else if (scheme_type == "N")
  {
    scheme = std::make_shared<rdm::PGRDMethodImplicit<MeshConfig, phys_model, rdm::PGN>>();
    // scheme = std::make_shared<rdm::DRDMethodImplicit<MeshConfig,
    // phys_model, rdm::PGN, FacetDG>>();
  }

  if (reorder_dofs)
  {
    // scheme->reorder_dofs(*mesh, *sol_dofs);

    std::vector<Int> reordering;
    scheme->compute_node_reordering(*mesh, *sol_dofs, reordering);
    (*sol_dofs).renumber_dofs(reordering);
  }

  scheme->configure_mesh_data(mesh, geo_dofs, sol_dofs);
  scheme->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quad_order, 1, 100);
  scheme->set_vec_function(SolverVecFn::solution, solution);
  scheme->set_vec_function(SolverVecFn::sources, sources);
  scheme->set_vec_function(SolverVecFn::residuals, residuals);

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  // Impose the BCs weakly:
  std::shared_ptr<rdm::RDMBCBase<MeshConfig, phys_model>> weak_Dirichlet_edge0 =
      scheme->add_boundary_condition("WeakDirichlet", "weak_dirichlet_edge0", "Edge0");
  weak_Dirichlet_edge0->set_expression(&MSEuler2D::Dirichlet_bc_value);

  std::shared_ptr<rdm::RDMBCBase<MeshConfig, phys_model>> weak_Dirichlet_edge1 =
      scheme->add_boundary_condition("WeakDirichlet", "weak_dirichlet_edge1", "Edge1");
  weak_Dirichlet_edge1->set_expression(&MSEuler2D::Dirichlet_bc_value);

  std::shared_ptr<rdm::RDMBCBase<MeshConfig, phys_model>> weak_Dirichlet_edge2 =
      scheme->add_boundary_condition("WeakDirichlet", "weak_dirichlet_edge2", "Edge2");
  weak_Dirichlet_edge2->set_expression(&MSEuler2D::Dirichlet_bc_value);

  std::shared_ptr<rdm::RDMBCBase<MeshConfig, phys_model>> weak_Dirichlet_edge3 =
      scheme->add_boundary_condition("WeakDirichlet", "weak_dirichlet_edge3", "Edge3");
  weak_Dirichlet_edge3->set_expression(&MSEuler2D::Dirichlet_bc_value);

  return scheme;
}

// ----------------------------------------------------------------------------

struct EulerUtestGlobalFixture
{
  /// common setup for each test case
  EulerUtestGlobalFixture()
  {
    /// arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~EulerUtestGlobalFixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(EulerTestSuite, EulerUtestGlobalFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(EulerTestSuite_ParenvInit)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(m_argc, m_argv);
  BOOST_CHECK_EQUAL(mpi_env.is_initialized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Euler2D_explicit_utest)
{
  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh = std::make_shared<MeshType>("mesh");

  MeshCreator::make_unit_quad(*mesh, "geo_dofs", 20, true);

  // The original mesh generated by MeshCreator is a square
  // <-1,1> x <-1,1>, but we need a square <0,1> x <0,1> - therefore
  // we will transform the original mesh
  math::DenseSVec<Real, _2D> transl_vector;
  transl_vector[X0] = 1.0;
  transl_vector[X1] = 1.0;
  MeshManipulator::translate_mesh(*mesh, "geo_dofs", transl_vector);
  MeshManipulator::scale_mesh(*mesh, "geo_dofs", 0.5);

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_continuous(*mesh, *geo_dofs, *sol_dofs, P1,
                                               PointSetID::Warpblend);

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P4;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  vector_function::ptr solution       = std::make_shared<vector_function>("", "solution");
  vector_function::ptr exact_solution = std::make_shared<vector_function>("", "exact_solution");
  vector_function::ptr source_term    = std::make_shared<vector_function>("", "source_term");
  vector_function::ptr residual       = std::make_shared<vector_function>("", "residual");
  vector_function::ptr sol_error      = std::make_shared<vector_function>("", "sol_error");

  const Uint nb_nodes = (*sol_dofs).nb_nodes();

  solution->resize(phys_model::NEQ, nb_nodes);
  exact_solution->resize(phys_model::NEQ, nb_nodes);
  source_term->resize(phys_model::NEQ, nb_nodes);
  residual->resize(phys_model::NEQ, nb_nodes);
  sol_error->resize(1, nb_nodes);

  MSEuler2D::fill_exact_solution_values(*sol_dofs, *exact_solution);
  MSEuler2D::fill_source_term_values(*sol_dofs, *source_term);

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  typedef rdm::PGRDMethodExplicit<MeshConfig, phys_model, rdm::PGLDA> scheme_type;
  scheme_type scheme;

  scheme.configure_mesh_data(mesh, mesh->dof_storage("geo_dofs"), mesh->dof_storage("sol_dofs"));
  scheme.initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1, 100);
  scheme.set_vec_function(SolverVecFn::solution, solution);
  scheme.set_vec_function(SolverVecFn::sources, source_term);
  scheme.set_vec_function(SolverVecFn::residuals, residual);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh);
  ic.set_domain(_2D, "Interior");
  ic.set_expression(&InitCondMSEuler2D::value);
  ic.apply_function("sol_dofs",
                    *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // BOUNDARY CONDITIONS
  // ------------------------------------------------------

  // Impose the BCs weakly:
  std::shared_ptr<scheme_type::bc_base_type> weak_Dirichlet_edge0 =
      scheme.add_boundary_condition("WeakDirichlet", "weak_dirichlet_edge0", "Edge0");
  weak_Dirichlet_edge0->set_expression(&MSEuler2D::Dirichlet_bc_value);

  std::shared_ptr<scheme_type::bc_base_type> weak_Dirichlet_edge1 =
      scheme.add_boundary_condition("WeakDirichlet", "weak_dirichlet_edge1", "Edge1");
  weak_Dirichlet_edge1->set_expression(&MSEuler2D::Dirichlet_bc_value);

  std::shared_ptr<scheme_type::bc_base_type> weak_Dirichlet_edge2 =
      scheme.add_boundary_condition("WeakDirichlet", "weak_dirichlet_edge2", "Edge2");
  weak_Dirichlet_edge2->set_expression(&MSEuler2D::Dirichlet_bc_value);

  std::shared_ptr<scheme_type::bc_base_type> weak_Dirichlet_edge3 =
      scheme.add_boundary_condition("WeakDirichlet", "weak_dirichlet_edge3", "Edge3");
  weak_Dirichlet_edge3->set_expression(&MSEuler2D::Dirichlet_bc_value);

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;

  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  const Real CFL = 0.5;

  start = clock();
  for (Uint iter = 0; iter < 2500; ++iter)
  {

    scheme.assemble_lhs_and_rhs(CFL);
    scheme.solve({});

    norm_L2(*residual, res_L2_norm);
    print_res_norm(iter, CFL, res_L2_norm);
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Compute the error
  ///
  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto quad_generator = [=](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order, PointSetID::Gauss);
  };

  interpolation::FunctionSpace<MeshConfig> geo_cell_space;
  geo_cell_space.set_reference_fe_values((*geo_dofs).as_range(), sf_generator, quad_generator);
  interpolation::FunctionSpace<MeshConfig> sol_cell_space;
  sol_cell_space.set_reference_fe_values((*sol_dofs).as_range(), sf_generator, quad_generator);

  for (Uint n = 0; n < nb_nodes; ++n)
  {
    interpolation::VectorMeshFunction<Real>::entry_type computed_value = solution->value(n);
    interpolation::VectorMeshFunction<Real>::entry_type exact_value    = exact_solution->value(n);
    interpolation::VectorMeshFunction<Real>::entry_type error_value    = sol_error->value(n);
    error_value[0] = computed_value[0] - exact_value[0];
  }

  const Real error_norm = solver::PostprocessingUtils<MeshConfig>::compute_function_norm(
      *geo_dofs, *sol_dofs, geo_cell_space, sol_cell_space, *sol_error);

  std::cout.precision(15);
  std::cout.setf(std::ios::fixed);
  std::cout << "L2 norm of error (explicit solver) = " << error_norm << std::endl;

  BOOST_CHECK_LE(error_norm, 5.e-3);

  /// Write the output to file

  gmsh::GmshWriter meshwriter;
  const std::string outfilename_exact = "solution_utest_Euler_continuous_rdm_reference.msh";

  meshwriter.write_mesh_to_file(*mesh, "sol_dofs", outfilename_exact);
  meshwriter.append_nodal_function_to_file(*mesh, outfilename_exact, *exact_solution, "u_exact");
  meshwriter.append_nodal_function_to_file(*mesh, outfilename_exact, *source_term, "source");

  const std::string outfilename_computed = "solution_utest_Euler_continuous_explicit_rdm.msh";
  meshwriter.write_mesh_to_file(*mesh, "sol_dofs", outfilename_computed);
  meshwriter.append_nodal_function_to_file(*mesh, outfilename_computed, *solution, "u");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Euler2D_implicit_utest)
{
  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh = std::make_shared<MeshType>("mesh");

  MeshCreator::make_unit_quad(*mesh, "geo_dofs", 20, true);

  // The original mesh generated by MeshCreator is a square
  // <-1,1> x <-1,1>, but we need a square <0,1> x <0,1> - therefore
  // we will transform the original mesh
  math::DenseSVec<Real, _2D> transl_vector;
  transl_vector[X0] = 1.0;
  transl_vector[X1] = 1.0;
  MeshManipulator::translate_mesh(*mesh, "geo_dofs", transl_vector);
  MeshManipulator::scale_mesh(*mesh, "geo_dofs", 0.5);

  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh->dof_storage("geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> sol_dofs = mesh->create_dof_storage("sol_dofs");

  MeshType::dof_storage_type::clone_continuous(*mesh, *geo_dofs, *sol_dofs, P1,
                                               PointSetID::Warpblend);

  // ------------------------------------------------------

  const PolyOrderID quadrature_order = P4;

  // ------------------------------------------------------
  // PREPARE MESH FUNCTIONS
  // ------------------------------------------------------

  vector_function::ptr solution       = std::make_shared<vector_function>("", "solution");
  vector_function::ptr exact_solution = std::make_shared<vector_function>("", "exact_solution");
  vector_function::ptr source_term    = std::make_shared<vector_function>("", "source_term");
  vector_function::ptr residual       = std::make_shared<vector_function>("", "residual");
  vector_function::ptr sol_error      = std::make_shared<vector_function>("", "sol_error");

  const Uint nb_nodes = (*sol_dofs).nb_nodes();

  solution->resize(phys_model::NEQ, nb_nodes);
  exact_solution->resize(phys_model::NEQ, nb_nodes);
  source_term->resize(phys_model::NEQ, nb_nodes);
  residual->resize(phys_model::NEQ, nb_nodes);
  sol_error->resize(1, nb_nodes);

  // ------------------------------------------------------
  // SOLVER SETUP
  // ------------------------------------------------------

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> n_scheme = build_solver(
      mesh, geo_dofs, sol_dofs, quadrature_order, solution, source_term, residual, "N", true);

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> lda_scheme = build_solver(
      mesh, geo_dofs, sol_dofs, quadrature_order, solution, source_term, residual, "LDA", false);

  // Source term values and exact solution have to be filled after we build
  // the solver, because the function 'build_solver' can reorder the degrees
  // of freedom in solution DofMap!

  MSEuler2D::fill_source_term_values(*sol_dofs, *source_term);
  MSEuler2D::fill_exact_solution_values(*sol_dofs, *exact_solution);

  // ------------------------------------------------------
  // INITIAL CONDITION
  // ------------------------------------------------------

  solution->fill(0.0);

  solver::InitialCondition<MeshConfig> ic("initial_cond", mesh);
  ic.set_domain(_2D, "Interior");
  ic.set_expression(&InitCondMSEuler2D::value);
  ic.apply_function("sol_dofs",
                    *solution); // The inlet function is also used as bc

  // ------------------------------------------------------
  // MAIN ITERATION LOOP
  // ------------------------------------------------------

  clock_t start, end;
  Real elapsed;

  math::DenseDVec<Real> res_L2_norm(phys_model::NEQ);

  std::shared_ptr<rdm::RDMethod<MeshConfig, phys_model>> active_scheme = n_scheme;

  start = clock();
  for (Uint iter = 0; iter < 100; ++iter)
  {
    const Real CFL = 10.0;

    active_scheme->set_cfl(CFL);

    if (iter % 5 == 0)
    {
      active_scheme->assemble_lhs_and_rhs(CFL);
      active_scheme->solve({SolverOption::RecomputePreconditioner});
    }
    else
    {
      active_scheme->assemble_rhs();
      active_scheme->solve({});
    }

    active_scheme->compute_residual_norm(res_L2_norm);
    print_res_norm(iter, CFL, res_L2_norm);

    if (iter == 20)
    {
      active_scheme = lda_scheme;
      active_scheme->assemble_lhs_and_rhs(CFL);
      active_scheme->solve({SolverOption::RecomputePreconditioner});
    }
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(3);
  std::cout << "CPU time (main iteration loop) = " << elapsed << " s" << std::endl;

  /// Compute the error
  interpolation::FunctionSpace<MeshConfig> geo_cell_space;

  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto quad_generator = [=](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order, PointSetID::Gauss);
  };

  geo_cell_space.set_reference_fe_values((*geo_dofs).as_range(), sf_generator, quad_generator);
  interpolation::FunctionSpace<MeshConfig> sol_cell_space;
  sol_cell_space.set_reference_fe_values((*sol_dofs).as_range(), sf_generator, quad_generator);

  for (Uint n = 0; n < nb_nodes; ++n)
  {
    interpolation::VectorMeshFunction<Real>::entry_type computed_value = solution->value(n);
    interpolation::VectorMeshFunction<Real>::entry_type exact_value    = exact_solution->value(n);
    interpolation::VectorMeshFunction<Real>::entry_type error_value    = sol_error->value(n);
    error_value[0] = computed_value[0] - exact_value[0];
  }

  const Real error_norm = solver::PostprocessingUtils<MeshConfig>::compute_function_norm(
      *geo_dofs, *sol_dofs, geo_cell_space, sol_cell_space, *sol_error);

  std::cout.precision(15);
  std::cout.setf(std::ios::fixed);
  std::cout << "L2 norm of error (implicit solver) = " << error_norm << std::endl;

  BOOST_CHECK_LE(error_norm, 5.e-3);

  /// Write the output to file

  gmsh::GmshWriter meshwriter;

  const std::string outfilename_computed = "solution_utest_Euler_continuous_implicit_rdm.msh";
  meshwriter.write_mesh_to_file(*mesh, "sol_dofs", outfilename_computed);
  meshwriter.append_nodal_function_to_file(*mesh, outfilename_computed, *solution, "u");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(EulerTestSuite_ParenvFinalize)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance();
  mpi_env.finalize();
  BOOST_CHECK_EQUAL(mpi_env.is_finalized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
