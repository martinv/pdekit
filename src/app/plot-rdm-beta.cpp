#include <ctime>
#include <iostream>
#include <memory>

#include "common/PDEKit.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/DenseSVec.hpp"
#include "math/polynomials/JacobiPolynomial.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

typedef Tria<Cart2D> MeshType;

int main()
{
  // Create the mesh on which we plot the shape function
  MeshType::shared_ptr mesh = std::make_shared<MeshType>("mesh");

  const Uint nb_pts_on_edge = 50;
  MeshCreator::make_unit_triangle(*mesh, "geo_dofs", nb_pts_on_edge);

  const result_of::dof_map_t<Cart2D> &geo_dofs = (*mesh->dof_storage("geo_dofs"));

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> shapefunc =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "shape_function");

  std::shared_ptr<interpolation::VectorMeshFunction<Real>> testfunc =
      std::make_shared<interpolation::VectorMeshFunction<Real>>("", "beta_function");

  const Uint nb_nodes = geo_dofs.nb_nodes();

  // Create a shape function that we are going to evaluate

  sf::ShapeFunction sf;
  const mesh::PointSetTag std_reg_tag(ElemShape::Triag, P4, PointSetID::Equidist);
  const mesh::sf::SFTag sf_tag(ElemShape::Triag, SFunc::Lagrange, P4, ModalBasis::Modal);
  sf.change_type(std_reg_tag, sf_tag);

  math::DenseDMat<Real> mesh_coordinates;
  math::DenseDMat<Real> sf_values;
  std::vector<math::DenseDMat<Real>> sf_derivatives;
  mesh_coordinates.resize(nb_nodes, _2D);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (MeshType::dof_storage_type::const_dof_iterator it = geo_dofs.cbegin(); it != geo_dofs.cend();
       ++it)
  {
    const mesh::CellTopologyView<Cart2D> tcell_view = it->tcell();
    const mesh::MeshEntity active_cell              = it->mesh_entity();

    const math::DenseConstMatView<Real> active_cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), active_cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < active_cell.nb_vert(); ++n)
    {
      mesh_coordinates.insert_row(active_cell.vertex(n), active_cell_coords.row(n));
    }
  }

  sf.get().compute_ref_values(mesh_coordinates, sf_values);
  sf.get().compute_ref_derivatives(mesh_coordinates, sf_derivatives);

  shapefunc->resize(sf_values.cols(), nb_nodes);
  testfunc->resize(sf_values.cols(), nb_nodes);

  math::DenseDVec<Real> k_upwind;
  k_upwind.resize(sf_values.cols());

  math::DenseDVec<Real> k_upwind_plus;
  k_upwind_plus.resize(sf_values.cols());

  const std::vector<Real> advection_dir = {1.0, 1.0};

  for (Uint n = 0; n < nb_nodes; ++n)
  {
    interpolation::VectorMeshFunction<Real>::entry_type sf_entry_in_node    = shapefunc->value(n);
    interpolation::VectorMeshFunction<Real>::entry_type testf_entry_in_node = testfunc->value(n);

    k_upwind.fill(0.0);
    k_upwind_plus.fill(0.0);

    Real sum_k_upwind_plus = 0.0;

    // Loop over all polynomials computed by 'sf'
    for (Uint p = 0; p < sf_values.cols(); ++p)
    {
      sf_entry_in_node[p] = sf_values(n, p);

      for (Uint d = 0; d < _2D; ++d)
      {
        k_upwind[p] += sf_derivatives[d](n, p) * advection_dir[d];
      }
      k_upwind_plus[p] = std::max(0.0, k_upwind[p]);
      sum_k_upwind_plus += k_upwind_plus[p];
    }

    for (Uint p = 0; p < sf_values.cols(); ++p)
    {
      testf_entry_in_node[p] = k_upwind_plus[p] / sum_k_upwind_plus;
    }
  }

  // Save to gmsh and vtk file formats
  const std::string outfilename_base = "reference_test_function";

  gmsh::GmshWriter gmsh_writer;
  gmsh_writer.write_mesh_to_file(*mesh, "geo_dofs", outfilename_base + ".msh");
  gmsh_writer.append_nodal_function_to_file(*mesh, outfilename_base + ".msh", *shapefunc,
                                            "shape_function");
  gmsh_writer.append_nodal_function_to_file(*mesh, outfilename_base + ".msh", *testfunc,
                                            "beta_test_function");

  vtk::VtkWriter vtk_writer;

  vtk_writer.append_nodal_function_to_file(*mesh, geo_dofs, "reference_sf.vtu", *shapefunc,
                                           "shape_function");
  vtk_writer.append_nodal_function_to_file(*mesh, geo_dofs, "reference_beta.vtu", *testfunc,
                                           "beta_test_function");

  // ----------------------------------------------------------------------------

  return 0;
}
