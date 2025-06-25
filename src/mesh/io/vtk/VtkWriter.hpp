#ifndef PDEKIT_VTK_Vtk_Writer_hpp
#define PDEKIT_VTK_Vtk_Writer_hpp

#include <ctime>
#include <fstream>

#include "common/StringUtils.hpp"
#include "mesh/MeshPredicates.hpp"
#include "mesh/Tria.hpp"
#include "mesh/io/vtk/PDEKitToVTK.hpp"

namespace pdekit
{

namespace interpolation
{
// Forward declaration
template <typename MFT>
class VectorMeshFunctionBase;
} // namespace interpolation

namespace mesh
{

namespace vtk
{

class VtkWriter
{
  public:
  /// Constructor
  VtkWriter();

  /// Destructor
  ~VtkWriter();

  /// Write a mesh
  /// @param mesh_root mesh which should be written
  /// @param filename to which file the mesh data should be written
  template <typename MeshConfig>
  void write_mesh_to_file(Tria<MeshConfig> const &mesh, const std::string &filename);

  /// Append field to file
  /// @param field pointer to the field which should be written
  /// @param name of the file to which the data is written
  template <typename MeshConfig, typename MFT>
  void append_nodal_function_to_file(
      const Tria<MeshConfig> &mesh, const typename result_of::dof_map_t<MeshConfig> &dof_handler,
      const std::string &filename, const interpolation::VectorMeshFunctionBase<MFT> &mesh_function,
      const std::string &field_name, const double time = 0.0, const Uint timestep = 0);

  private:
  /// METHODS

  template <typename DataMatrix>
  void append_nodal_function_to_vtu_file(const DataMatrix &data, const std::string &field_name,
                                         const double time, const int timestep = 0);

  template <typename MeshConfig>
  void write_coordinates(const Tria<MeshConfig> &mesh,
                         const typename result_of::dof_map_t<MeshConfig> &dof_handler);

  template <typename MeshConfig>
  void write_connectivity(const Tria<MeshConfig> &mesh,
                          const typename result_of::dof_map_t<MeshConfig> &dof_handler);

  /// DATA

  /// Input file stream
  std::ofstream m_outfile;
};

// ----------------------------------------------------------------------------==

template <typename MeshConfig, typename MFT>
void VtkWriter::append_nodal_function_to_file(
    const Tria<MeshConfig> &mesh, const typename result_of::dof_map_t<MeshConfig> &dof_handler,
    const std::string &filename, const interpolation::VectorMeshFunctionBase<MFT> &mesh_function,
    const std::string &field_name, const double time, const Uint timestep)
{
  clock_t start, end;
  Real elapsed;

  start = clock();

  m_outfile.open(filename.c_str());

  const Uint nb_cells = dof_handler.nb_active_cells();
  const Uint nb_nodes = dof_handler.nb_nodes();

  typename interpolation::VectorMeshFunctionBase<MFT>::MeshFunctionType const &data =
      mesh_function.wrapped_type();

  m_outfile.setf(std::ios::fixed);
  m_outfile.precision(17);

  m_outfile << "<?xml version=\"1.0\"?>" << std::endl;
  m_outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
               "byte_order=\"LittleEndian\">"
            << std::endl;
  m_outfile << "  <UnstructuredGrid>" << std::endl;

  m_outfile << "    <Piece NumberOfPoints=\"" << nb_nodes << "\" NumberOfCells=\"" << nb_cells
            << "\">" << std::endl;

  append_nodal_function_to_vtu_file(data, field_name, time, timestep);
  write_coordinates(mesh, dof_handler);
  write_connectivity(mesh, dof_handler);

  m_outfile << "    </Piece>" << std::endl;
  m_outfile << "  </UnstructuredGrid>" << std::endl;
  m_outfile << "</VTKFile>" << std::endl;

  // write_header( mesh );
  // write_coordinates( mesh );
  // write_connectivity( mesh );

  m_outfile.close();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "Writing " << filename << " took " << elapsed << " s" << std::endl;

  start = clock();
}

// ----------------------------------------------------------------------------==

template <typename DataMatrix>
void VtkWriter::append_nodal_function_to_vtu_file(const DataMatrix &data,
                                                  const std::string &field_name, const double time,
                                                  const int timestep)
{

  const std::string ascii_output_offset = "      ";
  std::vector<std::string> field_names;

  // Write the list of scalar data arrays
  m_outfile << ascii_output_offset << "<PointData Scalars=\"";
  for (Uint component = 0; component < data.nb_fields(); ++component)
  {
    field_names.push_back(field_name + "_" + common::StringUtils::to_string<Uint>(component));
    m_outfile << field_names[component];
    if (component < (data.nb_fields() - 1))
    {
      m_outfile << ",";
    }
  }
  m_outfile << "\">" << std::endl;

  // Write the actual data arrays
  for (Uint component = 0; component < data.nb_fields(); ++component)
  {
    m_outfile << ascii_output_offset << "  <DataArray Name=\"" << field_names[component]
              << R"(" type="Float32" format="ascii">)" << std::endl;
    for (Uint i = 0; i < data.nb_entries(); ++i)
    {
      typename DataMatrix::const_entry_type const nodal_value = data.const_value(i);
      m_outfile << ascii_output_offset << "    " << nodal_value[component] << std::endl;
    }
    m_outfile << ascii_output_offset << "  </DataArray>" << std::endl;
  }

  m_outfile << ascii_output_offset << "</PointData>" << std::endl;
  m_outfile << ascii_output_offset << "<CellData>" << std::endl;
  m_outfile << ascii_output_offset << "</CellData>" << std::endl;
}

// ----------------------------------------------------------------------------==

template <typename MeshConfig>
void VtkWriter::write_coordinates(const Tria<MeshConfig> &mesh,
                                  const typename result_of::dof_map_t<MeshConfig> &dof_handler)
{
  const std::string ascii_output_offset = "      ";
  m_outfile << ascii_output_offset << "<Points>" << std::endl;
  m_outfile << ascii_output_offset
            << "  <DataArray type=\"Float32\" "
               "Name=\"Nodes\" NumberOfComponents=\"3\" "
               "format=\"ascii\">"
            << std::endl;

  const Uint nb_nodes = dof_handler.nb_nodes();

  m_outfile.precision(17);
  //    m_outfile.setf(std::ios::fixed);

  std::vector<Real> raw_coords;
  raw_coords.resize(nb_nodes * MeshConfig::GDIM);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (typename result_of::dof_map_t<MeshConfig>::const_dof_iterator it = dof_handler.cbegin();
       it != dof_handler.cend(); ++it)
  {
    const mesh::MeshEntity active_cell             = it->mesh_entity();
    const mesh::CellTopologyView<MeshConfig> tcell = it->tcell();

    const math::DenseConstMatView<Real> active_cell_coords = loc_interpolator.transfer_coords(
        tcell.pt_set_id(), active_cell.pt_set_id(), tcell.coordinates());

    for (Uint v = 0; v < active_cell.nb_vert(); ++v)
    {
      const Uint offset = MeshConfig::GDIM * active_cell.vertex(v);
      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        raw_coords[offset + d] = active_cell_coords(v, d);
      }
    }
  }

  for (Uint inode = 0; inode < nb_nodes; ++inode)
  {
    m_outfile << ascii_output_offset << "    ";

    const Uint offset = MeshConfig::GDIM * inode;
    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      m_outfile << " " << raw_coords[offset + d];
    }
    if (MeshConfig::GDIM == _2D)
    {
      m_outfile << " 0.0";
    }

    m_outfile << std::endl;
  }

  m_outfile << ascii_output_offset << "  </DataArray>" << std::endl;
  m_outfile << ascii_output_offset << "</Points>" << std::endl;
}

// ----------------------------------------------------------------------------==

template <typename MeshConfig>
void VtkWriter::write_connectivity(const Tria<MeshConfig> &mesh,
                                   const typename result_of::dof_map_t<MeshConfig> &dof_handler)
{

  const std::string ascii_output_offset = "      ";
  m_outfile << ascii_output_offset << "<Cells>" << std::endl;
  m_outfile << ascii_output_offset
            << R"(  <DataArray type="Int64" Name="connectivity" format="ascii">)" << std::endl;

  // typename Tria<MeshConfig>::topology_type const &topology =
  // mesh.topology();

  // MeshBoundarySet const & boundaries  = mesh.all_boundaries();
  // typename MeshType::cell_connectivity const & conn = topology.cells();

  // const Uint total_nb_elements = boundaries.total_nb_active_cells() +
  // conn.nb_active_cells();

  /*
  Uint global_elem_id = 1;

  /// 1) WRITE ALL BOUNDARY CELLS
  for(const Cells::shared_ptr bdomain : boundaries.all_domains())
  {
    //const Uint dim = bdomain->dim();

    for(Uint c = 0; c < bdomain->nb_active_cells(); ++c)
    {
      //MeshEntity bcell = conn.cell(parent_cell_id);
      //bcell.transform_to_sub_entity(dim,local_id);
      const MeshEntity bcell = bdomain->cell(c);

      const Uint gmsh_type =
  Shared::ref_topology_type_to_gmsh_type(bcell.type_id());

      m_outfile << global_elem_id << " " << gmsh_type << " 2 " <<
  bdomain->tag(c)
                << " " << bdomain->tag(c);

      for(Uint node = 0; node < bcell.nb_vert(); ++node)
      {
        m_outfile << " " << bcell.vertex(node)+1;
      }
      m_outfile << std::endl;

      global_elem_id++;

    } // Loop over all cells of one boundary patch
  }
  */

  typedef typename result_of::dof_map_t<MeshConfig> cell_dofs_type;

  /// 2) WRITE ALL INTERNAL CELLS

  for (const typename cell_dofs_type::const_dof_range_typed &cells :
       dof_handler.all_active_dof_groups())
  {
    typename cell_dofs_type::const_dof_iterator_typed cell_iter = cells.begin();

    // const Uint gmsh_type =
    // Shared::ref_topology_type_to_gmsh_type(cell_iter->type_id());

    for (; cell_iter != cells.end(); ++cell_iter)
    {
      const MeshEntity one_cell = cell_iter->mesh_entity();
      m_outfile << ascii_output_offset << "  ";
      for (Uint node = 0; node < one_cell.nb_vert(); ++node)
      {
        m_outfile << "  " << one_cell.vertex(node);
      }
      m_outfile << std::endl;

    } // Loop over all cells of one group

  } // Loop over all element groups of the same dimension

  m_outfile << ascii_output_offset << "  </DataArray>" << std::endl;

  // --------------------------------------------

  // WRITE OFFSETS
  Uint connectivity_offset = 0;
  m_outfile << ascii_output_offset << R"(  <DataArray type="Int64" Name="offsets" format="ascii">)"
            << std::endl;

  for (const typename cell_dofs_type::const_dof_range_typed &cells :
       dof_handler.all_active_dof_groups())
  {
    typename cell_dofs_type::const_dof_iterator_typed cell_iter = cells.begin();

    for (; cell_iter != cells.end(); ++cell_iter)
    {
      const MeshEntity mesh_entity = cell_iter->mesh_entity();
      connectivity_offset += mesh_entity.nb_vert();
      m_outfile << ascii_output_offset << "    " << connectivity_offset << std::endl;

    } // Loop over all cells of one group

  } // Loop over all element groups of the same dimension
  m_outfile << ascii_output_offset << "  </DataArray>" << std::endl;

  // --------------------------------------------

  // WRITE ELEMENT TYPES
  m_outfile << ascii_output_offset << R"(  <DataArray type="UInt8" Name="types" format="ascii">)"
            << std::endl;

  for (const typename cell_dofs_type::const_dof_range_typed &cells :
       dof_handler.all_active_dof_groups())
  {
    typename cell_dofs_type::const_dof_iterator_typed cell_iter = cells.begin();

    for (; cell_iter != cells.end(); ++cell_iter)
    {
      const MeshEntity cell = cell_iter->mesh_entity();
      connectivity_offset += cell.nb_vert();
      m_outfile << ascii_output_offset << "     "
                << PDEKitToVTK::ref_topology_type_to_vtk_type(cell.pt_set_id()) << std::endl;

    } // Loop over all cells of one group

  } // Loop over all element groups of the same dimension
  m_outfile << ascii_output_offset << "  </DataArray>" << std::endl;

  // --------------------------------------------

  m_outfile << ascii_output_offset << "</Cells>" << std::endl;
}

// ----------------------------------------------------------------------------==

} // Namespace vtk

} // Namespace mesh

} // Namespace pdekit

#endif
