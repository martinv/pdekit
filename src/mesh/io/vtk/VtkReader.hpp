#ifndef PDEKIT_VTK_Vtk_Reader_hpp
#define PDEKIT_VTK_Vtk_Reader_hpp

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

class VtkReader
{
  public:
  /// Constructor
  VtkReader();

  /// Destructor
  ~VtkReader();

  /// Read a mesh from file. The type of file is detected based on filename
  /// extension
  template <typename MeshConfig>
  void read_mesh_from_file(const std::string &filename, Tria<MeshConfig> &mesh,
                           const std::string &geo_dofs_name);

  private:
  /// FUNCTIONS

  /// Re-initialize all member variables
  void reset();

  template <typename MeshConfig>
  void read_unstructured_grid(const std::string &filename, Tria<MeshConfig> &mesh,
                              const std::string &geo_dofs_name);

  template <typename MeshConfig>
  void read_polydata(const std::string &filename, Tria<MeshConfig> &mesh,
                     const std::string &geo_dofs_name);

  void read_coordinates(const Uint geo_dim);

  template <Uint GDIM, Uint TDIM>
  void read_connectivity_unstructured_grid(CellBuffer<GDIM, TDIM> &cell_buffer);

  template <Uint GDIM, Uint TDIM>
  void read_connectivity_polydata(CellBuffer<GDIM, TDIM> &cell_buffer);

  void read_celltypes_unstructured_grid(const Uint topo_dim);

  /// DATA

  /// Input file stream
  std::ifstream m_infile;

  /// Number of nodes
  Uint m_nb_nodes;

  /// Total number of elements
  Uint m_nb_elements;

  /// Total entries in the connectivity block
  Uint m_tot_nb_entries;

  /// Vector containing number of all nodes per dimension
  std::vector<Uint> m_nb_nodes_per_dim;

  /// Vector that remembers topological dimension for each element
  /// in the input file
  std::vector<Uint> m_elem_topo_dim;

  /// Vector of all possible element type exsiting in the mesh
  std::vector<std::unique_ptr<std::vector<StdRegion>>> m_element_type;

  /// Vector to hold coordinates
  std::vector<Real> m_coords;

  /// Array to hold a vector of connectivity for each topology dimension
  std::vector<std::unique_ptr<std::vector<Uint>>> m_connectivity;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void VtkReader::read_mesh_from_file(const std::string &filename, Tria<MeshConfig> &mesh,
                                    const std::string &geo_dofs_name)
{
  reset();

  clock_t start, end;
  Real elapsed;

  start = clock();

  std::string tempstr;
  std::vector<std::string> tempstrvector;
  common::StringUtils strtool;

  // 1st reading pass
  m_infile.open(filename.c_str());
  // Parse trhough the Mesh file

  tempstr = "";
  tempstrvector.resize(1);
  tempstrvector[0] = "";

  bool found_dataset_type = false;

  while (getline(m_infile, tempstr))
  {
    strtool.split_string(tempstr, ' ', tempstrvector);

    if ((!tempstrvector.empty()) && (tempstrvector[0] == "DATASET"))
    {
      std::cout << "Found the string DATASET" << std::endl;
      found_dataset_type = true;
      break;
    }
  }

  m_infile.close();

  if (!found_dataset_type || (tempstrvector.size() < 2))
  {
    std::cerr << "Did not find VTK data set type. Exiting." << std::endl;
    return;
  }

  if (tempstrvector[1] == "UNSTRUCTURED_GRID")
  {
    std::cout << "The VTK mesh is of UNSTRUCTURED GRID type" << std::endl;
    read_unstructured_grid(filename, mesh, geo_dofs_name);
  }

  if (tempstrvector[1] == "POLYDATA")
  {
    std::cout << "The VTK mesh is of POLYDATA type" << std::endl;
    read_polydata(filename, mesh, geo_dofs_name);
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "Reading " << filename << " took " << elapsed << " s" << std::endl;

} // End of read_mesh_from_file

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void VtkReader::read_unstructured_grid(const std::string &filename, Tria<MeshConfig> &mesh,
                                       const std::string &geo_dofs_name)
{
  std::string tempstr;
  std::vector<std::string> tempstrvector;
  common::StringUtils strtool;

  // 1st reading pass
  m_infile.open(filename.c_str());
  // Parse trhough the Mesh file

  tempstr = "";
  tempstrvector.resize(0);

  while (getline(m_infile, tempstr))
  {
    strtool.split_string(tempstr, ' ', tempstrvector);

    if ((!tempstrvector.empty()) && (tempstrvector[0] == "POINTS"))
    {
      std::cout << "Found the string POINTS" << std::endl;
      m_nb_nodes = strtool.from_string<Uint>(tempstrvector[1]);
      std::cout << "Number of points in the mesh: " << m_nb_nodes << std::endl;
      read_coordinates(MeshConfig::GDIM);
    }

    if ((!tempstrvector.empty()) && (tempstrvector[0] == "CELL_TYPES"))
    {
      std::cout << "Found the string CELL_TYPES" << std::endl;
      m_nb_elements = strtool.from_string<Uint>(tempstrvector[1]);
      read_celltypes_unstructured_grid(MeshConfig::TDIM);
    }
  }

  m_infile.close();

  // CellBuffer cell_buffer;

  // 2nd reading pass
  m_infile.open(filename.c_str());

  tempstr = "";
  tempstrvector.resize(0);

  CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> cell_buffer;

  while (getline(m_infile, tempstr))
  {
    strtool.split_string(tempstr, ' ', tempstrvector);

    if ((!tempstrvector.empty()) && (tempstrvector[0] == "CELLS"))
    {
      std::cout << "Found the string CELLS" << std::endl;
      m_tot_nb_entries = strtool.from_string<Uint>(tempstrvector[2]);
      std::cout << "Number of elements in the mesh: " << m_nb_elements << std::endl;
      read_connectivity_unstructured_grid(cell_buffer);
    }
  }

  m_infile.close();

  mesh.create_from_cells(cell_buffer);

  common::PtrHandle<typename Tria<MeshConfig>::dof_storage_type> geo_dofs =
      mesh.create_dof_storage(geo_dofs_name);

  (*geo_dofs).create_from_cells(cell_buffer);

  const std::vector<std::pair<Uint, std::string>> cell_tag_name = {{1, "interior"}};

  std::vector<Uint> material_ids;
  material_ids.resize(cell_buffer.nb_active_cells());
  material_ids.assign(cell_buffer.nb_active_cells(), 1);

  (*geo_dofs).tag_all_active_cells(cell_tag_name, material_ids);

} // End of read_unstructured_grid

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void VtkReader::read_polydata(const std::string &filename, Tria<MeshConfig> &mesh,
                              const std::string &geo_dofs_name)
{
  std::string tempstr;
  std::vector<std::string> tempstrvector;
  common::StringUtils strtool;

  // 1st reading pass
  m_infile.open(filename.c_str());
  // Parse trhough the Mesh file

  tempstr = "";
  tempstrvector.resize(0);

  CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> cell_buffer;

  while (getline(m_infile, tempstr))
  {
    strtool.split_string(tempstr, ' ', tempstrvector);

    if ((!tempstrvector.empty()) && (tempstrvector[0] == "POINTS"))
    {
      std::cout << "Found the string POINTS" << std::endl;
      m_nb_nodes = strtool.from_string<Uint>(tempstrvector[1]);
      std::cout << "Number of points in the mesh: " << m_nb_nodes << std::endl;
      read_coordinates(MeshConfig::GDIM);
    }

    if ((!tempstrvector.empty()) && (tempstrvector[0] == "POLYGONS"))
    {
      std::cout << "Found the string POLYGONS" << std::endl;
      m_nb_elements = strtool.from_string<Uint>(tempstrvector[1]);

      const Uint nb_conn_entries = strtool.from_string<Uint>(tempstrvector[2]);

      m_nb_nodes_per_dim.resize(MeshConfig::TDIM + 1);
      m_nb_nodes_per_dim.assign(MeshConfig::TDIM + 1, 0);
      m_nb_nodes_per_dim[_2D] = nb_conn_entries - m_nb_elements;
      read_connectivity_polydata(cell_buffer);
    }
  }
  m_infile.close();

  mesh.create_from_cells(cell_buffer);

  common::PtrHandle<typename Tria<MeshConfig>::dof_storage_type> geo_dofs =
      mesh.create_dof_storage(geo_dofs_name);

  (*geo_dofs).create_from_cells(cell_buffer);

  const std::vector<std::pair<Uint, std::string>> cell_tag_name = {{1, "interior"}};

  std::vector<Uint> material_ids;
  material_ids.resize(cell_buffer.nb_active_cells());
  material_ids.assign(cell_buffer.nb_active_cells(), 1);

  (*geo_dofs).tag_all_active_cells(cell_tag_name, material_ids);

} // End of read_polydata

// ----------------------------------------------------------------------------

template <Uint GDIM, Uint TDIM>
void VtkReader::read_connectivity_unstructured_grid(CellBuffer<GDIM, TDIM> &cell_buffer)
{
  m_connectivity.resize(TDIM + 1);

  for (Uint d = 0; d < (TDIM + 1); ++d)
  {
    m_connectivity[d] = std::unique_ptr<std::vector<Uint>>(new std::vector<Uint>());
  }

  for (Uint i_dim = 0; i_dim < m_nb_nodes_per_dim.size(); ++i_dim)
  {
    /*
    std::cout << "Dimension: " << i_dim << "D: total " <<
    m_nb_nodes_per_dim[i_dim] << " nodes"
              << std::endl;
    */
    m_connectivity[i_dim]->reserve(m_nb_nodes_per_dim[i_dim]);
    m_connectivity[i_dim]->resize(0);
  }

  Uint temp_var;

  std::vector<Uint> nb_counted_elem_per_dim(TDIM + 1);
  nb_counted_elem_per_dim.assign(TDIM + 1, 0);

  for (Uint i = 0; i < m_nb_elements; i++)
  {
    const Uint dim            = m_elem_topo_dim[i];
    const Uint elem_nr_in_dim = nb_counted_elem_per_dim[dim];

    const Uint nb_nodes = (*m_element_type[dim])[elem_nr_in_dim].get().nb_nodes();
    nb_counted_elem_per_dim[dim]++;

    m_infile >> temp_var;
    for (Uint i_entry = 0; i_entry < nb_nodes; i_entry++)
    {
      m_infile >> temp_var;
      m_connectivity[dim]->push_back(temp_var);
    }
  }

  /*
  for (Uint dim = 0; dim < (topo_dim + 1); ++dim)
  {
    std::cout << "Connectivity in dim " << dim << ": " << std::endl;
    for (Uint i = 0; i < m_connectivity[dim]->size(); ++i)
    {
      std::cout << (*m_connectivity[dim])[i] << " ";
    }
    std::cout << std::endl;
  }
  */

  std::unique_ptr<std::vector<Uint>> material_ids(new std::vector<Uint>());
  material_ids->resize(m_connectivity[TDIM]->size());
  material_ids->assign(m_connectivity[TDIM]->size(), 1);

  // Build a block array of cell coordinates
  std::unique_ptr<std::vector<Real>> coord_values(new std::vector<Real>());
  std::unique_ptr<std::vector<Uint>> coord_block_sizes(new std::vector<Uint>());

  // References for simpler notation ...
  std::vector<Uint> const &elem_connectivity = *m_connectivity[TDIM];
  std::vector<StdRegion> const &elem_types   = *m_element_type[TDIM];

  coord_values->resize(GDIM * elem_connectivity.size());
  coord_block_sizes->resize(elem_types.size());

  Uint next_coord_entry = 0;

  for (Uint i = 0; i < elem_connectivity.size(); ++i)
  {
    for (Uint d = 0; d < GDIM; ++d)
    {
      (*coord_values)[next_coord_entry + d] = m_coords[elem_connectivity[i] * GDIM + d];
    }
    next_coord_entry += GDIM;
  }

  for (Uint i = 0; i < elem_types.size(); ++i)
  {
    (*coord_block_sizes)[i] = GDIM * elem_types[i].get().nb_nodes();
  }

  std::unique_ptr<common::BlockArray<Real, Uint>> cell_coords(new common::BlockArray<Real, Uint>());
  cell_coords->build(std::move(coord_values), std::move(coord_block_sizes));

  cell_buffer.emplace_cell_data(std::move(m_connectivity[TDIM]), std::move(m_element_type[TDIM]),
                                std::move(material_ids), std::move(cell_coords));
}

// ----------------------------------------------------------------------------

template <Uint GDIM, Uint TDIM>
void VtkReader::read_connectivity_polydata(CellBuffer<GDIM, TDIM> &cell_buffer)
{
  std::cout << "[VTKReader] Reading Polygons ..." << std::endl;

  // m_elem_topo_dim.resize(m_nb_elements);

  m_element_type.resize(TDIM + 1);
  m_element_type[_2D] = std::unique_ptr<std::vector<StdRegion>>(new std::vector<StdRegion>());
  m_element_type[_2D]->resize(0);

  m_connectivity.resize(TDIM + 1);
  m_connectivity[_2D] = std::unique_ptr<std::vector<Uint>>(new std::vector<Uint>());
  m_connectivity[_2D]->reserve(m_nb_nodes_per_dim[_2D]);
  m_connectivity[_2D]->resize(0);

  Uint temp_var;
  StdRegion std_region;

  for (Uint i_elem = 0; i_elem < m_nb_elements; i_elem++)
  {
    m_infile >> temp_var;
    const Uint order = P1;
    if (temp_var == 3)
    {
      const ElemShape shape = ElemShape::Triag;
      const PointSetTag std_region_tag(shape, order, PointSetID::Equidist);
      std_region.change_type(std_region_tag);
    }
    else if (temp_var == 4)
    {
      const ElemShape shape = ElemShape::Quad;
      const PointSetTag std_region_tag(shape, order, PointSetID::Equidist);
      std_region.change_type(std_region_tag);
    }

    const Uint dim      = std_region.get().topo_dim();
    const Uint nb_nodes = std_region.get().nb_nodes();

    // m_elem_topo_dim[i_elem] = dim;

    m_element_type[dim]->push_back(std_region);
    m_nb_nodes_per_dim[dim] += nb_nodes;

    // m_infile >> temp_var;
    for (Uint i_entry = 0; i_entry < nb_nodes; i_entry++)
    {
      m_infile >> temp_var;
      m_connectivity[_2D]->push_back(temp_var);
    }
  }

  std::unique_ptr<std::vector<Uint>> material_ids(new std::vector<Uint>());
  material_ids->resize(m_connectivity[_2D]->size());
  material_ids->assign(m_connectivity[_2D]->size(), 1);

  // Build a block array of cell coordinates
  std::unique_ptr<std::vector<Real>> coord_values(new std::vector<Real>());
  std::unique_ptr<std::vector<Uint>> coord_block_sizes(new std::vector<Uint>());

  // References for simpler notation ...
  std::vector<Uint> const &elem_connectivity = *m_connectivity[TDIM];
  std::vector<StdRegion> const &elem_types   = *m_element_type[TDIM];

  coord_values->resize(GDIM * elem_connectivity.size());
  coord_block_sizes->resize(elem_types.size());

  Uint next_coord_entry = 0;

  for (Uint i = 0; i < elem_connectivity.size(); ++i)
  {
    for (Uint d = 0; d < GDIM; ++d)
    {
      (*coord_values)[next_coord_entry + d] = m_coords[elem_connectivity[i] * GDIM + d];
    }
    next_coord_entry += GDIM;
  }

  for (Uint i = 0; i < elem_types.size(); ++i)
  {
    (*coord_block_sizes)[i] = GDIM * elem_types[i].get().nb_nodes();
  }

  std::unique_ptr<common::BlockArray<Real, Uint>> cell_coords(new common::BlockArray<Real, Uint>());
  cell_coords->build(std::move(coord_values), std::move(coord_block_sizes));

  cell_buffer.emplace_cell_data(std::move(m_connectivity[_2D]), std::move(m_element_type[_2D]),
                                std::move(material_ids), std::move(cell_coords));

  std::cout << "[VTKReader] ... finished reading Polygons" << std::endl;
}

// ----------------------------------------------------------------------------

} // Namespace vtk

} // Namespace mesh

} // Namespace pdekit

#endif
