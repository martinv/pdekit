/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE solution_space_metric_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "interpolation/SolutionSpaceMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "mesh/adaptation/LocalInterpolator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"

using namespace pdekit;

template <typename MeshConfig>
void fill_solution(typename result_of::dof_map_t<MeshConfig> const &cell_dofs,
                   interpolation::VectorMeshFunction<Real> &solution, const SFunc shape_func,
                   const Uint quad_order)
{
  const Uint nb_fields        = MeshConfig::GDIM;
  const Uint nb_nodes_in_mesh = cell_dofs.nb_nodes();

  mesh::adapt::LocalInterpolator loc_interpolator;

  solution.resize(nb_fields, nb_nodes_in_mesh);

  if (shape_func == SFunc::Lagrange)
  {

    for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
      const mesh::MeshEntity cell = cell_dofs.active_cell(mesh::ActiveIdx(c));

      const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
          tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());
      for (Uint n = 0; n < cell_coords.rows(); ++n)
      {
        const math::DenseConstVecView<Real> node_coord = cell_coords.row_transpose(n);
        solution.insert_value(cell.vertex(n), node_coord);
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
      const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
      const mesh::MeshEntity cell = cell_dofs.active_cell(mesh::ActiveIdx(c));

      const mesh::PointSetTag std_reg_tag = cell.pt_set_id();

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

      const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
          tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

      wspace_nodal.resize(cell.nb_vert() * nb_fields);
      wspace_modal.resize(cell.nb_vert() * nb_fields);

      math::DenseMatView<Real> values_nodal(wspace_nodal.data(), nb_fields, cell.nb_vert(),
                                            nb_fields);
      math::DenseMatView<Real> values_modal(wspace_modal.data(), nb_fields, cell.nb_vert(),
                                            nb_fields);

      for (Uint n = 0; n < cell_coords.rows(); ++n)
      {
        const math::DenseConstVecView<Real> node_coord = cell_coords.row_transpose(n);
        for (Uint f = 0; f < nb_fields; ++f)
        {
          values_nodal(n, f) = node_coord[f];
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
// This function is supposed to be used with quadratures that place at
// least some of their points in element interior, such as standard
// Gauss quadrature
// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void test_metric(typename result_of::dof_map_t<MeshConfig> const &geo_dofs,
                 typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                 interpolation::VectorMeshFunction<Real> &solution, const SFunc shape_func,
                 const PointSetID quadrature_type, const Uint quadrature_order,
                 const Real reference_volume_value, const Real tol)
{
  std::cout << std::endl;
  std::cout << "**********************************************************" << std::endl;
  std::cout << " Test metric [" << MeshConfig::TDIM << "D, " << DIM << "D]" << std::endl;
  std::cout << " Shape function type: " << SFuncInfo::name(shape_func) << std::endl;
  std::cout << " Quadrature type:     " << PointSetInfo::name(quadrature_type) << std::endl;
  std::cout << "**********************************************************" << std::endl;

  clock_t t1, t2, t3;
  Real elapsed;

  typename interpolation::FunctionSpace<MeshConfig>::ptr fs_geometry =
      std::make_shared<interpolation::FunctionSpace<MeshConfig>>();

  const auto geo_sf_tag_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  const auto sol_sf_tag_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, shape_func, order, ModalBasis::Modal);
  };

  const auto quad_tag_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::PointSetTag(shape, quadrature_order, quadrature_type);
  };

  fs_geometry->set_reference_fe_values(common::make_iter_range(geo_dofs.cbegin(), geo_dofs.cend()),
                                       geo_sf_tag_generator, quad_tag_generator);

  typename interpolation::FunctionSpace<MeshConfig>::ptr fs_solution =
      std::make_shared<interpolation::FunctionSpace<MeshConfig>>();

  fs_solution->set_reference_fe_values(common::make_iter_range(sol_dofs.cbegin(), sol_dofs.cend()),
                                       sol_sf_tag_generator, quad_tag_generator);

  interpolation::GeometryCache<MeshConfig::GDIM> geometry_cache;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> geometry_metric;

  interpolation::SolutionCache solution_cache;
  interpolation::SolutionSpaceMetric<MeshConfig, DIM> solution_metric;

  geometry_cache.allocate(fs_geometry->discrete_elements().cbegin(),
                          fs_geometry->discrete_elements().cend(), geo_dofs.nb_active_cells());
  // For the moment, we use the same interpolation space for geometry and
  // solution
  geometry_metric.allocate_buffer(fs_geometry->discrete_elements().cbegin(),
                                  fs_geometry->discrete_elements().cend(),
                                  geo_dofs.nb_active_cells());

  solution_cache.allocate(fs_solution->reference_elements().cbegin(),
                          fs_solution->reference_elements().cend(), sol_dofs.nb_active_cells(),
                          MeshConfig::GDIM);

  solution_metric.allocate_buffer(fs_solution->reference_elements().cbegin(),
                                  fs_solution->reference_elements().cend(),
                                  sol_dofs.nb_active_cells(), MeshConfig::GDIM);

  BOOST_CHECK_EQUAL(geometry_metric.max_nb_blocks_in_buffer(), sol_dofs.nb_active_cells());

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity geo_cell                     = geo_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity sol_cell                     = sol_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), geo_cell.pt_set_id(), tcell_view.coordinates());

    const ElemShape elem_shape = geo_cell.pt_set_id().elem_shape();
    const Uint poly_order      = geo_cell.pt_set_id().poly_order();

    const mesh::PointSetTagExt geo_tag(geo_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0);
    const mesh::sf::SFTag sf_tag            = geo_sf_tag_generator(elem_shape, poly_order);
    const mesh::PointSetTag quad_pt_set_tag = quad_tag_generator(elem_shape, poly_order);
    const mesh::PointSetTagExt quad_pt_set_tag_ext(quad_pt_set_tag, P0,
                                                   mesh::CellTransform::NO_TRANS, 0);
    const mesh::DiscreteElemKey geo_key(geo_tag, sf_tag, quad_pt_set_tag_ext);

    geometry_cache.push_back_to_buffer(cell_coords, geo_key);
    solution_cache.push_back_to_buffer(
        sol_cell, solution,
        mesh::PointSetTagExt(sol_cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0));
  }

  t1 = clock();

  geometry_metric.evaluate(geometry_cache, interpolation::RebuildMetricIndex{true});

  t2 = clock();

  solution_metric.evaluate(geometry_metric, solution_cache,
                           interpolation::ComputeMetricDerivs{true},
                           interpolation::RebuildMetricIndex{true});

  t3 = clock();

  elapsed = ((Real)(t2 - t1)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (interpolation + computation of jacobians) = " << elapsed << " s"
            << std::endl;
  elapsed = ((Real)(t3 - t2)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (interpolation + reuse of jacobians) = " << elapsed << " s" << std::endl;

  Real domain_volume = 0.0;

  using geometry_metric_type =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             DIM>::cellwise_metric;
  using solution_metric_type =
      typename interpolation::SolutionSpaceMetric<MeshConfig, DIM>::cellwise_metric;

  math::DenseConstMatView<Real> sol_deriv[MeshConfig::GDIM];
  Real divergence;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    // Intentionally use the metric of the solution to check it correctly
    // reused the data from geometry metric cellwise_metric_type cell_met =
    // solution_metric.cellwise_values(c);
    geometry_metric_type geo_met = geometry_metric.cellwise_values(c);
    solution_metric_type sol_met = solution_metric.cellwise_values(c);

    math::DenseConstVecView<Real> jdet   = geo_met.jdet();
    math::DenseDVec<Real> const &weights = geo_met.pt_weights();

    if (MeshConfig::GDIM == DIM)
    {
      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        // std::cout << "Solution derivatives for field [" << d <<
        // "]:\n";
        sol_deriv[d] = sol_met.field_derivatives(d);
        // std::cout << sol_deriv[d] << std::endl;
      }
      // std::cout <<
      // "*****************************************************" <<
      // std::endl;
    }

    // As a test, we will integrate the divergence of the solution field,
    // which holds values (x,y,z) at every mesh node. This means that the
    // divergence is div(solution) = 1 + 1 = 2 (in 2D)
    //                                          and div(solution) = 1 + 1 +
    //                                          1 =
    // 3 (in 3D)
    // For the moment, this is done only for 'volume' metric, i.e. when
    // MeshConfig::GDIM == DIM
    // Otherwise we just integrate unity
    for (Uint q = 0; q < geo_met.nb_qd_pts(); ++q)
    {
      if (MeshConfig::GDIM == DIM)
      {
        divergence = 0.0;
        for (Uint d = 0; d < MeshConfig::GDIM; ++d)
        {
          divergence += sol_deriv[d](q, d);
        }
      }
      else
      {
        divergence = 1.0;
      }

      domain_volume += jdet[q] * weights[q] * divergence;
    }
  }

  std::cout.precision(10);
  std::cout.setf(std::ios::fixed);
  std::cout << "Domain measure = " << std::setw(10) << domain_volume << std::endl;

  BOOST_CHECK_CLOSE(domain_volume, reference_volume_value, tol);

  std::cout << "**********************************************************" << std::endl;
}

// ----------------------------------------------------------------------------
// This function is supposed to be used with quadratures that place points
// on element facets, such as FaceGauss quadrature
// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void test_metric_on_elem_bdry(
    typename result_of::dof_map_t<MeshConfig> const &cell_dofs,
    typename interpolation::FunctionSpace<MeshConfig, DIM>::ptr fs_solution,
    const PointSetID quadrature_type)
{
  std::cout << std::endl;
  std::cout << "**********************************************************" << std::endl;
  std::cout << " Test metric on element boundary [" << MeshConfig::TDIM << "D, " << DIM << "D]"
            << std::endl;
  std::cout << "**********************************************************" << std::endl;

  clock_t t1, t2;
  Real elapsed;

  const auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  const auto quad_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::PointSetTag(shape, P4, quadrature_type);
  };

  fs_solution->set_reference_fe_values(
      common::make_iter_range(cell_dofs.cbegin(), cell_dofs.cend()), sf_generator, quad_generator);

  interpolation::GeometryCache<MeshConfig::GDIM> geometry_cache;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> geometry_metric;

  interpolation::SolutionCache solution_cache;
  interpolation::SolutionSpaceMetric<MeshConfig, DIM> solution_metric;

  const Uint cache_size = 10u;

  geometry_cache.allocate(fs_solution->discrete_elements().cbegin(),
                          fs_solution->discrete_elements().cend(), cache_size);
  // For the moment, we use the same interpolation space for geometry and
  // solution
  geometry_metric.allocate_buffer(fs_solution->discrete_elements().cbegin(),
                                  fs_solution->discrete_elements().cend(), cache_size);

  solution_cache.allocate(fs_solution->reference_elements().cbegin(),
                          fs_solution->reference_elements().cend(), cache_size, MeshConfig::GDIM);
  solution_metric.allocate_buffer(fs_solution->reference_elements().cbegin(),
                                  fs_solution->reference_elements().cend(), cache_size,
                                  MeshConfig::GDIM);

  // geometry_cache.print_types();
  // solution_cache.print_types();

  BOOST_CHECK_EQUAL(geometry_metric.max_nb_blocks_in_buffer(), cache_size);

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> solution =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "solution");

  const Uint nb_fields        = MeshConfig::GDIM;
  const Uint nb_nodes_in_mesh = cell_dofs.nb_nodes();

  solution->resize(nb_fields, nb_nodes_in_mesh);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell                         = cell_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());
    for (Uint n = 0; n < cell_coords.rows(); ++n)
    {
      const math::DenseConstVecView<Real> node_coord = cell_coords.row_transpose(n);
      typename interpolation::VectorMeshFunction<Real>::entry_type nodal_value =
          solution->value(cell.vertex(n));

      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        //(*solution)(d, n) = node[d];
        nodal_value[d] = node_coord[d];
      }
    }
  }

  typedef typename interpolation::SolutionSpaceMetric<MeshConfig, DIM>::cellwise_metric
      solution_metric_type;

  t1 = clock();

  // common::StaticArray<math::ConstMatrixBlock<Real>, DIM> field_derivatives;

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell                         = cell_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    geometry_cache.flush();
    solution_cache.flush();

    geometry_metric.empty_buffer();
    solution_metric.empty_buffer();

    const ElemShape cell_shape = cell.pt_set_id().elem_shape();
    const Uint cell_order      = cell.pt_set_id().poly_order();

    for (Uint sub_cell = 0; sub_cell < cell.nb_sub_elements(MeshConfig::TDIM - 1); ++sub_cell)
    {
      const mesh::PointSetTagExt support_pts_ext(cell.pt_set_id(), P0,
                                                 mesh::CellTransform::NO_TRANS, sub_cell);
      const mesh::PointSetTag eval_pts = quad_generator(cell_shape, P4);
      const mesh::PointSetTagExt eval_pts_ext(eval_pts, P0, mesh::CellTransform::NO_TRANS,
                                              sub_cell);

      const mesh::sf::SFTag sf_tag = sf_generator(cell_shape, cell_order);

      const mesh::DiscreteElemKey geo_key(support_pts_ext, sf_tag, eval_pts_ext);

      geometry_cache.push_back_to_buffer(cell_coords, geo_key);
      solution_cache.push_back_to_buffer(
          cell, *solution,
          mesh::PointSetTagExt(cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, sub_cell));
    }

    geometry_metric.evaluate(geometry_cache, interpolation::RebuildMetricIndex{true});
    solution_metric.evaluate(geometry_metric, solution_cache,
                             interpolation::ComputeMetricDerivs{true},
                             interpolation::RebuildMetricIndex{true});

    // The solution field is u = [x,y]. This means that:
    // du_1/dx = 1, du_1/dy = 0, du_2/dx = 0, du_2/dy = 1
    // Now we're going to check that the d-th derivative of d-th component
    // of the solution field is one and all other derivatives of d-th
    // component of u are equal to zero
    for (Uint i = 0; i < solution_metric.nb_values_in_buffer(); ++i)
    {
      // geometry_metric_type geo_met =
      // geometry_metric.cellwise_values(i);
      solution_metric_type sol_met = solution_metric.cellwise_values(i);

      for (Uint d = 0; d < DIM; ++d)
      {
        const math::DenseConstMatView<Real> sol_derivatives = sol_met.field_derivatives(d);

        for (Uint q = 0; q < sol_met.nb_qd_pts(); ++q)
        {
          BOOST_CHECK_CLOSE(sol_derivatives(q, d), 1.0, 1.e-13);

          for (Uint n = 0; n < DIM; ++n)
          {
            if (n != d)
            {
              BOOST_CHECK_LE(std::abs(sol_derivatives(q, n)), 1.e-12);
            }
          }
        } // Loop over quadrature points

      } // Loop over dimensions

    } // Loop over values in buffer

  } // Loop over all cells in mesh

  t2 = clock();

  elapsed = ((Real)(t2 - t1)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (interpolation + computation of jacobians on cell "
               "facets) = "
            << elapsed << " s" << std::endl;

  std::cout << "**********************************************************" << std::endl;
}

// ----------------------------------------------------------------------------

typedef mesh::Cart2D MeshConfig2D;
typedef mesh::Cart3D MeshConfig3D;

typedef mesh::Tria<MeshConfig2D> MeshType2D;
typedef mesh::Tria<MeshConfig3D> MeshType3D;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(solution_cache_2D_utest)
{
  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_mixed_p1.msh", *mesh2d, "geo_dofs");

  const common::PtrHandle<mesh::DofMap<mesh::Cart2D> const> geo_dofs =
      mesh2d->dof_storage("geo_dofs");

  interpolation::VectorMeshFunction<Real> solution("", "solution");
  solution.resize(2, (*geo_dofs).nb_nodes());

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < (*geo_dofs).nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<mesh::Cart2D> tcell_view = (*geo_dofs).tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell = (*geo_dofs).active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < cell_coords.rows(); ++n)
    {
      const math::DenseConstVecView<Real> node_coord = cell_coords.row_transpose(n);

      typename interpolation::VectorMeshFunction<Real>::entry_type node_val =
          solution.value(cell.vertex(n));

      node_val[0] = node_coord[X0] + node_coord[X1];
      node_val[1] = node_coord[X0] * node_coord[X1];
    }
  }

  interpolation::FunctionSpace<mesh::Cart2D> sol_space;
  sol_space.set_reference_fe_values(
      common::make_iter_range((*geo_dofs).cbegin(), (*geo_dofs).cend()),
      [=](const ElemShape shape, const Uint order) {
        return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
      },
      [=](const ElemShape shape, const Uint order) {
        return mesh::PointSetTag(shape, P2, PointSetID::Gauss);
      });

  interpolation::SolutionCache sol_cache;
  sol_cache.allocate(sol_space.reference_elements().cbegin(), sol_space.reference_elements().cend(),
                     (*geo_dofs).nb_active_cells(), 2);

  const mesh::DofMap<mesh::Cart2D>::all_dof_type_ranges elem_ranges =
      (*geo_dofs).all_active_dof_groups();

  std::vector<common::Range1D<Uint>> num_ranges;

  Uint idx_start = 0;
  Uint idx_end   = 0;

  for (const mesh::DofMap<mesh::Cart2D>::const_dof_range_typed &elem_range : elem_ranges)
  {
    for (mesh::DofMap<mesh::Cart2D>::const_dof_iterator_typed eit = elem_range.begin();
         eit != elem_range.end(); ++eit)
    {
      const mesh::MeshEntity elem = eit->mesh_entity();
      const mesh::PointSetTagExt tag_ext(elem.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0);
      sol_cache.push_back_to_buffer(elem, solution, tag_ext);
      idx_end++;
    }
    num_ranges.push_back(common::Range1D<Uint>(idx_start, idx_end - 1));
    idx_start = idx_end;
  }

  std::cout << "Element ranges: " << std::endl;
  for (auto range : num_ranges)
  {
    std::cout << "[" << range.lbound() << "," << range.ubound() << "]" << std::endl;
  }

  const Uint local_node = 1;
  const Uint component  = 0;
  sol_cache.perturb_values(num_ranges[0], local_node, component, 1.e-1);
  const math::DenseConstVecView<Real> unperturbed_vals =
      sol_cache.unperturbed_values(num_ranges[0]);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_2D_Lagrange_utest)
{
  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_tri_p2.msh", *mesh2d, "geo_dofs");

  common::PtrHandle<MeshType2D::dof_storage_type> geo_dofs = mesh2d->dof_storage("geo_dofs");
  common::PtrHandle<MeshType2D::dof_storage_type> sol_dofs = mesh2d->create_dof_storage("sol_dofs");

  MeshType2D::dof_storage_type::clone_discontinuous(*mesh2d, *geo_dofs, *sol_dofs, P2,
                                                    PointSetID::Equidist);

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> solution =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "solution");
  fill_solution(*sol_dofs, *solution, SFunc::Lagrange, P4);

  const Real tolerance = 1.e-7;
  test_metric<MeshConfig2D, _2D>(*geo_dofs, *sol_dofs, *solution, SFunc::Lagrange,
                                 PointSetID::Gauss, P4, 2.0, tolerance);
}

// ----------------------------------------------------------------------------

#if 1
BOOST_AUTO_TEST_CASE(metric_2D_Carnevali_utest)
{
  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_tri_p2.msh", *mesh2d, "geo_dofs");

  common::PtrHandle<MeshType2D::dof_storage_type> geo_dofs = mesh2d->dof_storage("geo_dofs");
  common::PtrHandle<MeshType2D::dof_storage_type> sol_dofs = mesh2d->create_dof_storage("sol_dofs");

  MeshType2D::dof_storage_type::clone_discontinuous(*mesh2d, *geo_dofs, *sol_dofs, P2,
                                                    PointSetID::Equidist);

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> solution =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "solution");
  fill_solution(*sol_dofs, *solution, SFunc::Carnevali, P4);

  const Real tolerance = 1.e-7;
  test_metric<MeshConfig2D, _2D>(*geo_dofs, *sol_dofs, *solution, SFunc::Carnevali,
                                 PointSetID::Gauss, P4, 2.0, tolerance);
}
#endif

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_2D_face_quadrature_utest)
{
  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_tri_p2.msh", *mesh2d, "geo_dofs");

  interpolation::FunctionSpace<MeshConfig2D>::ptr fs_solution =
      std::make_shared<interpolation::FunctionSpace<MeshConfig2D>>();

  test_metric_on_elem_bdry<MeshConfig2D, _2D>(*(mesh2d->dof_storage("geo_dofs")), fs_solution,
                                              PointSetID::FaceGauss);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_3D_Lagrange_utest)
{
  MeshType3D::shared_ptr mesh3d = std::make_shared<MeshType3D>("mesh3D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_cube_tet_p3.msh", *mesh3d, "geo_dofs");

  common::PtrHandle<MeshType3D::dof_storage_type> geo_dofs = mesh3d->dof_storage("geo_dofs");
  common::PtrHandle<MeshType3D::dof_storage_type> sol_dofs = mesh3d->create_dof_storage("sol_dofs");

  MeshType3D::dof_storage_type::clone_discontinuous(*mesh3d, *geo_dofs, *sol_dofs, P3,
                                                    PointSetID::Equidist);

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> solution =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "solution");
  fill_solution(*sol_dofs, *solution, SFunc::Lagrange, P6);

  const Real tolerance = 1.e-6;
  test_metric<MeshConfig3D, _3D>(*geo_dofs, *sol_dofs, *solution, SFunc::Lagrange,
                                 PointSetID::Gauss, P6, 3.0, tolerance);
}

// ----------------------------------------------------------------------------

#if 0
BOOST_AUTO_TEST_CASE(metric_3D_Dubiner_utest)
{
  MeshType3D::shared_ptr mesh3d = MeshType3D::shared_ptr(new MeshType3D("mesh3D"));

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_cube_tet_p3.msh", *mesh3d, "geo_dofs");

  common::PtrHandle<MeshType3D::dof_storage_type> geo_dofs = mesh3d->dof_storage("geo_dofs");
  common::PtrHandle<MeshType3D::dof_storage_type> sol_dofs = mesh3d->create_dof_storage("sol_dofs");

  MeshType3D::dof_storage_type::clone_discontinuous(*mesh3d, *geo_dofs, *sol_dofs, P3,
                                                    PointSetID::Equidist);

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> solution =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "solution");
  fill_solution(*sol_dofs, *solution, SFunc::Modal, P6);

  const Real tolerance = 1.e-7;
  test_metric<MeshConfig3D, _3D>(*geo_dofs, *sol_dofs, *solution, SFunc::Modal, PointSetID::Gauss, P6,
                                 3.0, tolerance);
}
#endif

// ----------------------------------------------------------------------------

// This can only be enabled when face quadratures in 3D are available
/*
BOOST_AUTO_TEST_CASE(metric_3D_face_quadrature_utest)
{
  MeshManager<MeshConfig3D>::instance_type &meshes3d =
MeshManager<MeshConfig3D>::instance(); mesh::Tria<MeshConfig3D>::shared_ptr
mesh3d = meshes3d.mesh("mesh3D");

  interpolation::FunctionSpace<MeshConfig3D>::ptr fs_solution =
      std::make_shared<interpolation::FunctionSpace<MeshConfig3D>>(mesh3d);

  test_metric_on_elem_bdry<MeshConfig3D, _3D>(mesh3d->topology().cells(),
mesh3d->geometry(), fs_solution, FaceGauss);
}
*/

// ----------------------------------------------------------------------------

#if 0
BOOST_AUTO_TEST_CASE(metric_1D_boundary_utest)
{
  MeshType2D::shared_ptr mesh2d = MeshType2D::shared_ptr(new MeshType2D("circle2D"));

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_circle_tri_p3.msh", *mesh2d);

  interpolation::FunctionSpace<MeshConfig2D, _1D>::ptr fs_solution =
      std::make_shared<interpolation::FunctionSpace<MeshConfig2D, _1D>>(mesh2d);

  mesh::MeshBoundarySet<MeshConfig2D> const &boundary = mesh2d->topology().all_boundaries();

  const Real tolerance = 1.e-5;
  test_metric<MeshConfig2D, _1D>(*boundary.domain("boundary"), mesh2d->geometry(), fs_solution,
                                 Gauss, 2. * PI, tolerance);
}
#endif

// ----------------------------------------------------------------------------

#if 0
BOOST_AUTO_TEST_CASE(metric_2D_boundary_utest)
{
  MeshType3D::shared_ptr mesh3d = MeshType3D::shared_ptr(new MeshType3D("sphere2D"));

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_sphere_tet_p3.msh", *mesh3d);

  interpolation::FunctionSpace<MeshConfig3D, _2D>::ptr fs_solution =
      std::make_shared<interpolation::FunctionSpace<MeshConfig3D, _2D>>(mesh3d);

  mesh::MeshBoundarySet<MeshConfig3D> const &boundary = mesh3d->topology().all_boundaries();

  const Real tolerance = 7.e-5;
  test_metric<MeshConfig3D, _2D>(*boundary.domain("boundary"), mesh3d->geometry(), fs_solution,
                                 Gauss, 4. * PI, tolerance);
}
#endif

// ----------------------------------------------------------------------------
