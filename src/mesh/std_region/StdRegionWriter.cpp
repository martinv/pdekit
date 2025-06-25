#include <fstream>

#include "math/DenseConstVecView.hpp"
#include "mesh/io/gmsh/PDEKitToGmsh.hpp"
#include "mesh/std_region/StdRegion.hpp"
#include "mesh/std_region/StdRegionWriter.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

StdRegionWriter::StdRegionWriter()
{
}

// ----------------------------------------------------------------------------

StdRegionWriter::~StdRegionWriter()
{
}

// ----------------------------------------------------------------------------

void StdRegionWriter::write_to_vtk(PointSetTag const elem_type_tag, const std::string &filename)
{
  StdRegion std_region;
  std_region.change_type(elem_type_tag);

  const math::DenseDMat<Real> &node_coordinates = std_region.get().coordinates();

  const Uint topo_dim = node_coordinates.cols();

  std::string zero_coords = "";
  if (topo_dim == 1)
  {
    zero_coords = " 0.0 0.0";
  }
  else if (topo_dim == 2)
  {
    zero_coords = " 0.0";
  }

  Uint nb_pts_to_draw      = 0;
  Uint nb_segments_to_draw = 0;

  const Uint nb_edges = std_region.get().nb_entities(_1D);

  /*
  std::cout << "Edges" << std::endl;
  for (Uint e = 0; e < nb_edges; ++e)
  {
    std::cout << "E[" << e << "]" << (*std_region.get().elem_entity(_1D, e))
  << std::endl;
  }
  */

  for (Uint e = 0; e < nb_edges; ++e)
  {
    const std::shared_ptr<StdRegionEntity const> edge = std_region.get().elem_entity(_1D, e);
    nb_pts_to_draw += (*edge).nb_vert();
    nb_segments_to_draw += ((*edge).nb_vert() - 1);
  }

  std::ofstream outfile;
  outfile.precision(14);
  outfile.setf(std::ios::fixed);
  outfile.open(filename.c_str());

  outfile << "# vtk DataFile Version 2.0" << std::endl;
  outfile << "interpolation point set" << std::endl;
  outfile << "ASCII" << std::endl;
  outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
  outfile << "POINTS " << node_coordinates.rows() << " double" << std::endl;
  for (Uint n = 0; n < node_coordinates.rows(); ++n)
  {
    outfile << node_coordinates.const_row(n) << zero_coords << std::endl;
  }

  outfile << std::endl
          << "CELLS " << nb_segments_to_draw + node_coordinates.rows() << " "
          << 2 * node_coordinates.rows() + 3 * nb_segments_to_draw << std::endl;

  /// Print all nodes
  for (Uint c = 0; c < node_coordinates.rows(); ++c)
  {
    outfile << "1 " << c << std::endl;
  }

  /// Print all segments of boundary edges of this element
  /// This is just for pretty-printing
  // Helper vector to re-order the vertices for plotting
  std::vector<Uint> edge_vertices;
  for (Uint e = 0; e < nb_edges; ++e)
  {
    const std::shared_ptr<StdRegionEntity const> edge = std_region.get().elem_entity(_1D, e);
    edge_vertices.clear();
    // First edge vertex
    edge_vertices.push_back((*edge).vertex(0));
    // All interior vertices
    for (Uint n = 2; n < (*edge).nb_vert(); ++n)
    {
      edge_vertices.push_back((*edge).vertex(n));
    }
    // Last edge vertex (which is on position 1)
    edge_vertices.push_back((*edge).vertex(1));

    for (Uint n = 0; n < (edge_vertices.size() - 1); ++n)
    {
      outfile << "2 " << edge_vertices[n] << " " << edge_vertices[n + 1] << std::endl;
    }
  }

  outfile << std::endl
          << "CELL_TYPES " << node_coordinates.rows() + nb_segments_to_draw << std::endl;

  for (Uint c = 0; c < node_coordinates.rows(); ++c)
  {
    outfile << "1" << std::endl;
  }

  for (Uint s = 0; s < nb_segments_to_draw; ++s)
  {
    outfile << "3" << std::endl; // The type of all segments is a line
  }

  outfile << std::endl << "POINT_DATA " << node_coordinates.rows() << std::endl;
  outfile << "SCALARS node_xyz_coord float 1" << std::endl;
  outfile << "LOOKUP_TABLE default" << std::endl;

  for (Uint c = 0; c < node_coordinates.rows(); ++c)
  {
    double value = node_coordinates(c, 0);
    for (Uint d = 1; d < topo_dim; ++d)
    {
      value += node_coordinates(c, d);
    }
    outfile << value << std::endl;
  }

  outfile.close();
}

// ----------------------------------------------------------------------------

void StdRegionWriter::write_to_point3d(PointSetTag const elem_type_tag, const std::string &filename)
{
  /// TO PLOT THIS IN VISIT:
  /// Pseudocolor->point_xyz_value
  /// Label->operators->Connected components->Points to see node numbers

  StdRegion std_region;
  std_region.change_type(elem_type_tag);

  const math::DenseDMat<Real> &node_coordinates = std_region.get().coordinates();

  const Uint topo_dim = node_coordinates.cols();

  std::string zero_coords = " ";
  if (topo_dim == 1)
  {
    zero_coords = " 0.0 0.0 ";
  }
  else if (topo_dim == 2)
  {
    zero_coords = " 0.0 ";
  }

  std::ofstream outfile;
  outfile.precision(14);
  outfile.setf(std::ios::fixed);
  outfile.open(filename.c_str());

  outfile << "x y z point_xyz_value" << std::endl;

  for (Uint c = 0; c < node_coordinates.rows(); ++c)
  {
    double value = 0.0;
    for (Uint d = 0; d < topo_dim; ++d)
    {
      outfile << node_coordinates(c, d) << " ";
      value += node_coordinates(c, d);
    }
    outfile << zero_coords;
    outfile << value << std::endl;
  }

  outfile.close();
}

// ----------------------------------------------------------------------------

void StdRegionWriter::write_to_gmsh(PointSetTag const elem_type_tag, const std::string &filename)
{
  StdRegion std_region;
  std_region.change_type(elem_type_tag);

  const math::DenseDMat<Real> &node_coordinates = std_region.get().coordinates();

  const Uint topo_dim = node_coordinates.cols();

  std::string zero_coords = "";
  if (topo_dim == 1)
  {
    zero_coords = " 0.0 0.0";
  }
  else if (topo_dim == 2)
  {
    zero_coords = " 0.0";
  }

  std::ofstream outfile;
  outfile.precision(14);
  outfile.setf(std::ios::fixed);
  outfile.open(filename.c_str());

  outfile << "$MeshFormat" << std::endl;
  outfile << "2.2 0 8" << std::endl;
  outfile << "$EndMeshFormat" << std::endl;
  outfile << "$Nodes" << std::endl;
  outfile << node_coordinates.rows() << std::endl;

  for (Uint c = 0; c < node_coordinates.rows(); ++c)
  {
    outfile << c + 1 << " ";
    for (Uint d = 0; d < topo_dim; ++d)
    {
      outfile << node_coordinates(c, d) << " ";
    }
    outfile << zero_coords << std::endl;
  }

  outfile << "$EndNodes" << std::endl;

  gmsh::PDEKitToGmsh pdekit_to_gmsh;
  outfile << "$Elements" << std::endl;
  // In the future, we could perform complete decomposition of the element and
  // write
  // all of its edges and faces
  outfile << "1" << std::endl;
  const Uint gmsh_type =
      pdekit_to_gmsh.ref_topology_type_to_gmsh_type(std_region.get().pt_set_id());
  outfile << "1 " << gmsh_type << " 2 1 1";

  const std::shared_ptr<StdRegionEntity const> std_reg_connectivity =
      std_region.get().elem_entity(topo_dim, 0);

  for (Uint i = 0; i < (*std_reg_connectivity).nb_vert(); ++i)
  {
    outfile << " " << (*std_reg_connectivity).vertex(i) + 1u;
  }
  outfile << std::endl;

  outfile << "$EndElements" << std::endl;

  outfile.close();
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
