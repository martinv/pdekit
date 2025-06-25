#include <ctime>
#include <iostream>
#include <memory>

#include "common/PDEKit.hpp"
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

  const Uint nb_nodes = geo_dofs.nb_nodes();

  // Create a shape function that we are going to evaluate

  sf::ShapeFunction sf;
  const mesh::PointSetTag std_reg_tag(ElemShape::Triag, P4, PointSetID::Equidist);
  const mesh::sf::SFTag sf_tag(ElemShape::Triag, SFunc::Lagrange, P4, ModalBasis::Modal);
  sf.change_type(std_reg_tag, sf_tag);

  math::DenseDMat<Real> mesh_coordinates;
  math::DenseDMat<Real> sf_values;
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

  shapefunc->resize(sf_values.cols(), nb_nodes);
  for (Uint n = 0; n < nb_nodes; ++n)
  {
    interpolation::VectorMeshFunction<Real>::entry_type sf_entry_in_node = shapefunc->value(n);

    // Loop over all polynomials computed by 'sf'
    for (Uint p = 0; p < sf_values.cols(); ++p)
    {
      sf_entry_in_node[p] = sf_values(n, p);
    }
  }

  // Save to gmsh and vtk file formats
  const std::string outfilename_base = "reference_interpolant";

  gmsh::GmshWriter gmsh_writer;
  gmsh_writer.write_mesh_to_file(*mesh, "geo_dofs", outfilename_base + ".msh");
  gmsh_writer.append_nodal_function_to_file(*mesh, outfilename_base + ".msh", *shapefunc,
                                            "shape_function");

  /*
  vtk::VtkWriter vtk_writer;
  vtk_writer.append_nodal_function_to_file(*mesh, outfilename_base + ".vtu",
  *shapefunc, "shape_function");
  */

  // ----------------------------------------------------------------------------

  math::JacobiPolynomial jp;

  const Uint nb_pts     = 201;
  const Real dx         = 2. / (nb_pts - 1);
  const Real alpha      = 0.;
  const Real beta       = 0.;
  const Uint poly_order = 4;

  math::DenseSVec<Real, poly_order + 1> values;
  math::DenseDMat<Real> ref_coords(nb_pts, 1);
  math::DenseDMat<Real> Lagrange_sf_values(nb_pts, poly_order + 1);

  std::ofstream outfile;
  outfile.open("Legendre_P4.dat");
  outfile.precision(14);

  for (Uint i = 0; i < nb_pts; ++i)
  {
    const Real x     = -1. + i * dx;
    ref_coords(i, 0) = x;

    for (Uint p = 0; p <= poly_order; ++p)
    {
      values[p] = jp(p, alpha, beta, x);
    }
    outfile << x << " " << values << std::endl;
  }
  outfile.close();

  mesh::sf::ShapeFunction Lagrange_sf;
  const mesh::PointSetTag std_reg_tag_line(ElemShape::Line, P4, PointSetID::Equidist);
  const mesh::sf::SFTag sf_tag_line(ElemShape::Line, SFunc::Lagrange, P4, ModalBasis::Modal);
  Lagrange_sf.change_type(std_reg_tag_line, sf_tag_line);

  Lagrange_sf.get().compute_ref_values(ref_coords, Lagrange_sf_values);

  outfile.open("Lagrange_P4.dat");

  for (Uint i = 0; i < nb_pts; ++i)
  {
    outfile << ref_coords(i, 0) << " " << Lagrange_sf_values.row(i) << std::endl;
  }

  outfile.close();

  return 0;
}
