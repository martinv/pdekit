/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE flux_space_metric_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "interpolation/FluxSpaceMetric.hpp"
#include "interpolation/SolutionSpaceMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "mesh/adaptation/LocalInterpolator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

struct PhysicsMockup2D
{
  enum
  {
    NEQ = 4
  };

  enum
  {
    DIM = 2
  };

  struct Properties
  {
    math::DenseSVec<Real, _2D> coords;
    math::DenseSVec<Real, 4> vars;
    math::DenseSMat<Real, 4, _2D> grad_vars;
  };

  template <typename V1, typename V2, typename M1>
  static void compute_properties(const V1 &coord, const V2 &sol, const M1 &grad_vars, Properties &p)
  {
    p.coords    = coord;     // cache the coordinates locally
    p.vars      = sol;       // cache the variables locally
    p.grad_vars = grad_vars; // cache the gradient of variables locally
  }

  template <typename M1>
  static void flux(const Properties &p, M1 &flux)
  {
    flux(0, X) = p.vars[X] * p.vars[X]; // x^2
    flux(1, X) = p.vars[X] * p.vars[X];
    flux(2, X) = p.vars[X] * p.vars[X];
    flux(3, X) = p.vars[X] * p.vars[X];

    flux(0, Y) = p.vars[Y] * p.vars[Y]; // y^2
    flux(1, Y) = p.vars[Y] * p.vars[Y];
    flux(2, Y) = p.vars[Y] * p.vars[Y];
    flux(3, Y) = p.vars[Y] * p.vars[Y];
  }
};

// ----------------------------------------------------------------------------

struct PhysicsMockup3D
{
  enum
  {
    NEQ = 5
  };

  enum
  {
    DIM = 3
  };

  struct Properties
  {
    math::DenseSVec<Real, _3D> coords;
    math::DenseSVec<Real, 5> vars;
    math::DenseSMat<Real, 5, _3D> grad_vars;
  };

  template <typename V1, typename V2, typename M1>
  static void compute_properties(const V1 &coord, const V2 &sol, const M1 &grad_vars, Properties &p)
  {
    p.coords    = coord;     // cache the coordinates locally
    p.vars      = sol;       // cache the variables locally
    p.grad_vars = grad_vars; // cache the gradient of variables locally
  }

  template <typename M1>
  static void flux(const Properties &p, M1 &flux)
  {
    flux(0, X) = p.vars[X] * p.vars[X]; // x^2
    flux(1, X) = p.vars[X] * p.vars[X];
    flux(2, X) = p.vars[X] * p.vars[X];
    flux(3, X) = p.vars[X] * p.vars[X];
    flux(4, X) = p.vars[X] * p.vars[X];

    flux(0, Y) = p.vars[Y] * p.vars[Y]; // y^2
    flux(1, Y) = p.vars[Y] * p.vars[Y];
    flux(2, Y) = p.vars[Y] * p.vars[Y];
    flux(3, Y) = p.vars[Y] * p.vars[Y];
    flux(4, Y) = p.vars[Y] * p.vars[Y];

    flux(0, Z) = p.vars[Z] * p.vars[Z]; // z^2
    flux(1, Z) = p.vars[Z] * p.vars[Z];
    flux(2, Z) = p.vars[Z] * p.vars[Z];
    flux(3, Z) = p.vars[Z] * p.vars[Z];
    flux(4, Z) = p.vars[Z] * p.vars[Z];
  }
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void fill_solution(typename result_of::dof_map_t<MeshConfig> const &cell_dofs,
                   interpolation::VectorMeshFunction<Real> &solution, const SFunc shape_func,
                   const Uint quad_order, const Uint nb_fields)
{
  const Uint nb_nodes_in_mesh = cell_dofs.nb_nodes();
  solution.resize(nb_fields, nb_nodes_in_mesh);

  mesh::adapt::LocalInterpolator loc_interpolator;

  if (shape_func == SFunc::Lagrange)
  {
    for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
    {
      const mesh::MeshEntity cell = cell_dofs.active_cell(mesh::ActiveIdx(c));
      const mesh::CellTopologyView<MeshConfig> topo_view = cell_dofs.tcell(mesh::ActiveIdx(c));

      const math::DenseConstMatView<Real> interpolated_coord = loc_interpolator.transfer_coords(
          topo_view.pt_set_id(), cell.pt_set_id(), topo_view.coordinates());

      for (Uint n = 0; n < cell.nb_vert(); ++n)
      {
        interpolation::VectorMeshFunction<Real>::entry_type values_nodal =
            solution.value(cell.vertex(n));

        for (Uint f = 0; f < nb_fields; ++f)
        {
          // values_nodal[f] = node_coord[f % MeshConfig::GDIM];
          values_nodal[f] = interpolated_coord(n, f % MeshConfig::GDIM);
          // values_nodal[f] = node_coord[f % MeshConfig::GDIM] + 0.1
          // * (f % MeshConfig::GDIM);
        }
      }
    }
  }
  else if (sf_is_modal(shape_func))
  {
    common::DataMap<mesh::PointSetTag, math::DenseDMat<Real>> basis_transfer_mats;
    common::DataMap<mesh::PointSetTag, math::DenseDMat<Real>> V_modal_mats;
    common::DataMap<mesh::PointSetTag, math::DenseDMat<Real>> V_nodal_mats;

    std::vector<Real> wspace_nodal, wspace_modal;

    for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
    {
      const mesh::MeshEntity cell         = cell_dofs.active_cell(mesh::ActiveIdx(c));
      const mesh::PointSetTag std_reg_tag = cell.pt_set_id();

      const mesh::CellTopologyView<MeshConfig> topo_view = cell_dofs.tcell(mesh::ActiveIdx(c));

      common::PtrHandle<math::DenseDMat<Real>> transfer_mat =
          basis_transfer_mats.std_region_data(std_reg_tag);

      if (transfer_mat.is_null())
      {
        transfer_mat = basis_transfer_mats.create(std_reg_tag);

        mesh::StdRegion std_reg(std_reg_tag);

        const mesh::sf::SFTag sf_tag_nodal(std_reg_tag.elem_shape(), SFunc::Lagrange,
                                           std_reg_tag.poly_order(), ModalBasis::Modal);
        const mesh::sf::SFTag sf_tag_modal(std_reg_tag.elem_shape(), shape_func,
                                           std_reg_tag.poly_order(), ModalBasis::Modal);

        mesh::StdPointSet quadrature;
        quadrature.change_type(std_reg_tag.elem_shape(), quad_order, PointSetID::Gauss);

        interpolation::FEValues fe_values_nodal;
        fe_values_nodal.configure(std_reg_tag, sf_tag_nodal);
        fe_values_nodal.fill_Vandermonde(quadrature.get().coordinates(),
                                         quadrature.get().weights());

        interpolation::FEValues fe_values_modal;
        fe_values_modal.configure(std_reg_tag, sf_tag_modal);
        fe_values_modal.fill_Vandermonde(quadrature.get().coordinates(),
                                         quadrature.get().weights());

        const math::DenseDMat<Real> &V_nodal_tmp = fe_values_nodal.Vandermonde();
        const math::DenseDMat<Real> &V_modal_tmp = fe_values_modal.Vandermonde();

        common::PtrHandle<math::DenseDMat<Real>> V_modal = V_modal_mats.create(std_reg_tag);
        (*V_modal).resize(V_modal_tmp.rows(), V_modal_tmp.cols());
        (*V_modal) = V_modal_tmp;

        common::PtrHandle<math::DenseDMat<Real>> V_nodal = V_nodal_mats.create(std_reg_tag);
        (*V_nodal).resize(V_nodal_tmp.rows(), V_nodal_tmp.cols());
        (*V_nodal) = V_nodal_tmp;

        /*
        std::cout << "V(nodal) = \n" << *V_nodal << std::endl;
        std::cout << "V(modal) = \n" << *V_modal << std::endl;
        */

        const Uint nb_dof = (*V_nodal).cols();

        math::DenseDMat<Real> V_modal_inv;
        V_modal_inv.resize(nb_dof, nb_dof);
        (*V_modal).inv(V_modal_inv);

        (*transfer_mat).resize(nb_dof, nb_dof);
        (*transfer_mat) = V_modal_inv * (*V_nodal);

        /*
        std::cout << "Transfer matrix:\n" << std::endl;
        std::cout << *transfer_mat << std::endl;
        */
      }

      const math::DenseConstMatView<Real> interpolated_coord = loc_interpolator.transfer_coords(
          topo_view.pt_set_id(), std_reg_tag, topo_view.coordinates());

      wspace_nodal.resize(cell.nb_vert() * nb_fields);
      wspace_modal.resize(cell.nb_vert() * nb_fields);

      math::DenseMatView<Real> values_nodal(wspace_nodal.data(), nb_fields, cell.nb_vert(),
                                            nb_fields);
      math::DenseMatView<Real> values_modal(wspace_modal.data(), nb_fields, cell.nb_vert(),
                                            nb_fields);

      for (Uint n = 0; n < cell.nb_vert(); ++n)
      {
        for (Uint f = 0; f < nb_fields; ++f)
        {
          values_nodal(n, f) = interpolated_coord(n, f % MeshConfig::GDIM);
          // values_nodal(n, f) = node_coord[f % MeshConfig::GDIM] +
          // 0.1 * (f % MeshConfig::GDIM);
        }
      }

      // ---
      // Now we can express expansion coefficients in modal basis:
      values_modal = (*transfer_mat) * values_nodal;
      // ---

      common::PtrHandle<math::DenseDMat<Real>> V_nodal = V_nodal_mats.std_region_data(std_reg_tag);
      common::PtrHandle<math::DenseDMat<Real>> V_modal = V_modal_mats.std_region_data(std_reg_tag);

      const Uint nb_qd_pts = (*V_nodal).rows();

      // Express nodal values in quadrature points:
      std::vector<Real> buffer(nb_qd_pts * nb_fields);
      math::DenseMatView<Real> buffer_view(buffer.data(), nb_fields, nb_qd_pts, nb_fields);

      buffer_view = (*V_nodal) * values_nodal;
      // std::cout << "Nodal values in quadrature points:\n" <<
      // buffer_view << std::endl;

      // ---

      buffer_view = (*V_modal) * values_modal;
      // std::cout << "Modal values in quadrature points:\n" <<
      // buffer_view << std::endl;

      for (Uint n = 0; n < cell.nb_vert(); ++n)
      {
        solution.insert_value(cell.vertex(n), values_modal.row_transpose(n));
      }
    } // Loop over active cells
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
void test_metric(typename result_of::tria_t<MeshConfig> const &triangulation,
                 typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                 interpolation::VectorMeshFunction<Real> const &solution, const SFunc shape_func,
                 const PointSetID quadrature_type, const Uint quadrature_order, const Real tol)
{
  clock_t t1, t2, t3, t4;
  Real elapsed;

  typename interpolation::FunctionSpace<MeshConfig>::ptr fs_geometry =
      std::make_shared<interpolation::FunctionSpace<MeshConfig>>();

  const auto geo_sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  const auto geo_eval_pts_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::PointSetTag(shape, quadrature_order, quadrature_type);
  };

  const auto sol_sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, shape_func, order, ModalBasis::Modal);
  };

  const auto sol_eval_pts_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::PointSetTag(shape, quadrature_order, quadrature_type);
  };

  fs_geometry->set_reference_fe_values(
      common::make_iter_range(triangulation.cbegin_cells(), triangulation.cend_cells()),
      geo_sf_generator, geo_eval_pts_generator);

  typename interpolation::FunctionSpace<MeshConfig>::ptr fs_solution =
      std::make_shared<interpolation::FunctionSpace<MeshConfig>>();

  fs_solution->set_reference_fe_values(common::make_iter_range(sol_dofs.cbegin(), sol_dofs.cend()),
                                       sol_sf_generator, sol_eval_pts_generator);

  interpolation::GeometryCache<MeshConfig::GDIM> geometry_cache;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> geometry_metric;

  interpolation::SolutionCache solution_cache;
  interpolation::SolutionSpaceMetric<MeshConfig, DIM> solution_metric;
  interpolation::FluxSpaceMetric<MeshConfig, Physics, DIM> flux_metric;

  const mesh::MeshEntity first_geo_cell =
      (*sol_dofs.all_active_dof_groups()[0].begin()).mesh_entity();
  const mesh::PointSetTag first_geo_cell_tag = first_geo_cell.pt_set_id();
  const Uint geo_poly_order                  = first_geo_cell_tag.poly_order();

  geometry_cache.allocate(fs_geometry->discrete_elements().cbegin(),
                          fs_geometry->discrete_elements().cend(), sol_dofs.nb_active_cells());
  // We use the same function space for geometry and solution metric:
  geometry_metric.allocate_buffer(fs_geometry->discrete_elements().cbegin(),
                                  fs_geometry->discrete_elements().cend(),
                                  sol_dofs.nb_active_cells());

  solution_cache.allocate(fs_solution->reference_elements().cbegin(),
                          fs_solution->reference_elements().cend(), sol_dofs.nb_active_cells(),
                          Physics::NEQ);
  solution_metric.allocate_buffer(fs_solution->reference_elements().cbegin(),
                                  fs_solution->reference_elements().cend(),
                                  sol_dofs.nb_active_cells(), Physics::NEQ);
  flux_metric.allocate_buffer(SFunc::Lagrange, geo_poly_order,
                              fs_solution->reference_elements().cbegin(),
                              fs_solution->reference_elements().cend(), sol_dofs.nb_active_cells());

  BOOST_CHECK_EQUAL(geometry_metric.max_nb_blocks_in_buffer(), sol_dofs.nb_active_cells());

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity sol_cell                     = sol_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));

    const ElemShape cell_shape = tcell_view.pt_set_id().elem_shape();
    const Uint cell_order      = tcell_view.pt_set_id().poly_order();

    const mesh::PointSetTagExt geo_support_pt_set_ext =
        mesh::PointSetTagExt(tcell_view.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u);
    const mesh::sf::SFTag basis_tag         = geo_sf_generator(cell_shape, cell_order);
    const mesh::PointSetTag geo_eval_pt_set = geo_eval_pts_generator(cell_shape, cell_order);

    const mesh::PointSetTagExt geo_eval_pt_set_ext(geo_eval_pt_set, P0,
                                                   mesh::CellTransform::NO_TRANS, 0u);

    const mesh::DiscreteElemKey geo_key(geo_support_pt_set_ext, basis_tag, geo_eval_pt_set_ext);

    const mesh::CellGeometry<MeshConfig::GDIM> cell_coords = tcell_view.coordinates();

    geometry_cache.push_back_to_buffer(cell_coords, geo_key);

    const mesh::PointSetTagExt tmp(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u);

    solution_cache.push_back_to_buffer(
        sol_cell, solution,
        mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0u));
  }

  t1 = clock();

  geometry_metric.evaluate(geometry_cache, interpolation::RebuildMetricIndex{true});

  t2 = clock();

  solution_metric.evaluate(geometry_metric, solution_cache,
                           interpolation::ComputeMetricDerivs{true},
                           interpolation::RebuildMetricIndex{true});

  t3 = clock();

  flux_metric.evaluate(geometry_cache, geometry_metric, solution_cache, solution_metric,
                       interpolation::RebuildMetricIndex{true});

  t4 = clock();

  elapsed = ((Real)(t2 - t1)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (geometric metric terms) = " << elapsed << " s" << std::endl;
  elapsed = ((Real)(t3 - t2)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (solution metric terms) = " << elapsed << " s" << std::endl;
  elapsed = ((Real)(t4 - t3)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (flux metric terms) = " << elapsed << " s" << std::endl;

  typedef typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                                 DIM>::cellwise_metric geo_cell_metric_type;
  typedef typename interpolation::FluxSpaceMetric<MeshConfig, Physics, DIM>::cellwise_metric
      flux_cell_metric_type;

  // If the mesh is not a boundary mesh, check derivatives
  if (MeshConfig::GDIM == DIM)
  {
    math::DenseConstMatView<Real> flux_deriv[MeshConfig::GDIM];

    // Loop over all cell values and check computed flux derivatives against
    // expected values
    for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
    {
      flux_cell_metric_type flux_cell_met = flux_metric.cellwise_values(c);

      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        flux_deriv[d] = flux_cell_met.flux_derivatives(d);
      }

      geo_cell_metric_type geo_cell_met = geometry_metric.cellwise_values(c);

      const math::DenseConstMatView<Real> quad_pt_coords = geo_cell_met.interpolated_coords();

      for (Uint d = 0; d < Physics::DIM; ++d)
      {
        const math::DenseConstVecView<Real> coord_component = quad_pt_coords.col(d);
        const math::DenseConstVecView<Real> flux_component  = flux_deriv[d].col(d);

        for (Uint i = 0; i < coord_component.size(); ++i)
        {
          // Since the flux F0 = x^2, dF0/dx = 2 * x
          // Since the flux F1 = y^2, dF1/dy = 2 * y
          // Therefore the reference value is 2 * coord_component[i],
          // where coord_component is gradually the vector of all
          // x-values in quadrature points, then all y-values in
          // quadrature points etc.
          const Real reference = 2.0 * coord_component[i];
          BOOST_CHECK_CLOSE(flux_component[i], reference, tol);
        } // For loop to check flux components
      }
    }
  } // ( if MeshConfig::GDIM == DIM )

  /*
  else // Flux space metric computed on facets
  {
    for (Uint c = 0; c < cell_dofs.nb_cells(); ++c)
    {

      geo_cell_metric_type geo_cell_met = geometry_metric.cellwise_values(c);
      flux_cell_metric_type flux_cell_met = flux_metric.cellwise_values(c);


      const math::MatrixBlock<Real> quad_pt_coords =
  geo_cell_met.interpolated_coords();

      std::cout << "Cell = " << cell_dofs.cell(c) << std::endl;
      std::cout << "Quadrature points = " << quad_pt_coords << std::endl;

      for(Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        std::cout << "Fluxes[" << d << "]" << std::endl <<
  flux_cell_met.flux_values(d) << std::endl;
      }
      std::cout << "-------------------------------------------" << std::endl
  << std::endl;
    }
  } // else
  */

} // test_metric

typedef mesh::Cart2D MeshConfig2D;
typedef mesh::Cart3D MeshConfig3D;

typedef mesh::Tria<MeshConfig2D> MeshType2D;
typedef mesh::Tria<MeshConfig3D> MeshType3D;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_2D_Lagrange_utest)
{
  // std::cout << "DEBUG: Lagrange" << std::endl;
  const SFunc sf_type         = SFunc::Lagrange;
  const Uint quadrature_order = P4;

  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_tri_p2.msh", *mesh2d, "geo_dofs");

  common::PtrHandle<MeshType2D::dof_storage_type> geo_dofs = mesh2d->dof_storage("geo_dofs");
  common::PtrHandle<MeshType2D::dof_storage_type> sol_dofs = mesh2d->create_dof_storage("sol_dofs");

  MeshType2D::dof_storage_type::clone_discontinuous(*mesh2d, *geo_dofs, *sol_dofs, P2,
                                                    PointSetID::Equidist);

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> solution =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "solution");
  fill_solution(*sol_dofs, *solution, sf_type, quadrature_order, PhysicsMockup2D::NEQ);

  const Real tolerance = 2.e-9;
  test_metric<MeshConfig2D, PhysicsMockup2D, _2D>(*mesh2d, *sol_dofs, *solution, sf_type,
                                                  PointSetID::Gauss, quadrature_order, tolerance);
  std::cout << std::endl << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_2D_Carnevali_utest)
{
  std::cout << "DEBUG: Carnevali" << std::endl;
  const SFunc sf_type         = SFunc::Carnevali;
  const Uint quadrature_order = P4;

  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_tri_p2.msh", *mesh2d, "geo_dofs");

  common::PtrHandle<MeshType2D::dof_storage_type> geo_dofs = mesh2d->dof_storage("geo_dofs");
  common::PtrHandle<MeshType2D::dof_storage_type> sol_dofs = mesh2d->create_dof_storage("sol_dofs");

  MeshType2D::dof_storage_type::clone_discontinuous(*mesh2d, *geo_dofs, *sol_dofs, P2,
                                                    PointSetID::Equidist);

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> solution =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "solution");
  fill_solution(*sol_dofs, *solution, sf_type, quadrature_order, PhysicsMockup2D::NEQ);

  const Real tolerance = 2.e-9;
  test_metric<MeshConfig2D, PhysicsMockup2D, _2D>(*mesh2d, *sol_dofs, *solution, sf_type,
                                                  PointSetID::Gauss, quadrature_order, tolerance);
  std::cout << std::endl << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_3D_utest)
{
  const SFunc sf_type         = SFunc::Lagrange;
  const Uint quadrature_order = P4;

  MeshType3D::shared_ptr mesh3d = std::make_shared<MeshType3D>("mesh3D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_cube_tet_p3.msh", *mesh3d, "geo_dofs");

  common::PtrHandle<MeshType3D::dof_storage_type> geo_dofs = mesh3d->dof_storage("geo_dofs");
  common::PtrHandle<MeshType3D::dof_storage_type> sol_dofs = mesh3d->create_dof_storage("sol_dofs");

  MeshType3D::dof_storage_type::clone_discontinuous(*mesh3d, *geo_dofs, *sol_dofs, P3,
                                                    PointSetID::Equidist);

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> solution =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "solution");
  fill_solution(*sol_dofs, *solution, sf_type, quadrature_order, PhysicsMockup3D::NEQ);

  const Real tolerance = 2.e-9;
  test_metric<MeshConfig3D, PhysicsMockup3D, _3D>(*mesh3d, *sol_dofs, *solution, sf_type,
                                                  PointSetID::Gauss, quadrature_order, tolerance);

  std::cout << std::endl << std::endl;
}

// ----------------------------------------------------------------------------

/*
BOOST_AUTO_TEST_CASE(metric_1D_facets_utest)
{
  MeshManager<MeshConfig2D>::instance_type &meshes2d =
MeshManager<MeshConfig2D>::instance(); mesh::Tria<MeshConfig2D>::shared_ptr
mesh2d = meshes2d.mesh("square2D");

  interpolation::FunctionSpace<MeshConfig2D, _1D>::ptr fs_solution =
      std::make_shared<interpolation::FunctionSpace<MeshConfig2D, _1D>>(mesh2d);

  mesh::MeshBoundarySet<MeshConfig2D> const &boundary =
mesh2d->all_boundaries();

  mesh::TopologyAlgorithms::compute_connectivity(mesh2d->topology(), _1D, _0D,
true);

  test_metric<MeshConfig2D, PhysicsMockup2D, _1D>(mesh2d->topology().cells(_1D),
mesh2d->geometry(), fs_solution);

  std::cout << std::endl << std::endl;
}
*/
// ----------------------------------------------------------------------------

/*
BOOST_AUTO_TEST_CASE(metric_2D_boundary_utest)
{
  MeshManager<MeshConfig3D>::instance_type &meshes3d =
      MeshManager<MeshConfig3D>::instance();
  mesh::Tria<MeshConfig3D>::shared_ptr mesh3d =
      meshes3d.create_mesh("sphere2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_sphere_tet_p3.msh", *mesh3d);

  interpolation::FunctionSpace<MeshConfig3D, _2D>::ptr fs_solution =
      std::make_shared<interpolation::FunctionSpace<MeshConfig3D, _2D>>(mesh3d);

  mesh::MeshBoundarySet<MeshConfig3D> const &boundary =
      mesh3d->all_boundaries();

  const Real tolerance = 7.e-5;
  test_metric<MeshConfig3D, physics::Euler3DCons,
_2D>(*boundary.domain("boundary"),
                                 mesh3d->geometry(), fs_solution, 4. * PI,
                                 tolerance);
}

// ----------------------------------------------------------------------------

*/
