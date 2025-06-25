#ifndef PDEKIT_Mesh_Gmsh_Gmsh_Reader_hpp
#define PDEKIT_Mesh_Gmsh_Gmsh_Reader_hpp

#include <ctime>
#include <fstream>
#include <unordered_map>

#include "mesh/CellBuffer.hpp"
#include "mesh/Tria.hpp"
#include "mesh/containers/TopologyAlgorithms.hpp"
#include "mesh/io/gmsh/PDEKitToGmsh.hpp"

namespace pdekit
{

namespace mesh
{

namespace gmsh
{

class GmshReader
{
  public:
  /// Constructor
  GmshReader();

  /// Destructor
  ~GmshReader();

  /// Read mesh from file
  /// @param: filename from which the mesh should be read
  template <typename MeshConfig>
  void read_mesh_from_file(const std::string &filename, Tria<MeshConfig> &mesh,
                           const std::string &geo_dofs_name);

  template <typename MeshConfig>
  void read_mesh_from_file(const std::string &filename, Tria<MeshConfig> &mesh,
                           typename result_of::dof_map_t<MeshConfig> &cells_dofs);

  /// Read field from file
  template <typename MeshConfig, typename DataMatrix>
  void read_nodal_function_from_file(const Tria<MeshConfig> &mesh, const std::string &filename,
                                     DataMatrix &data, const std::vector<std::string> &field_names,
                                     const Real time = 0.0, const Uint timestep = 0);

  private:
  /// AUXILIARY TYPES AND TYPEDEFS:

  struct SubDomainHeaderData
  {
    std::string name;
    Uint idx;
    Uint topo_dim;
  };

  /// FUNCTIONS

  /// Re-initialize all member variables
  void reset();

  void get_record_positions();

  void read_coordinates(const Uint geo_dim);

  template <typename MeshConfig>
  void read_connectivity(Tria<MeshConfig> &mesh,
                         typename result_of::dof_map_t<MeshConfig> &interior_cells_dofs);

  /// DATA

  /// Input file stream
  std::ifstream m_infile;

  /// Buffer to read the file line by line
  std::string tempstr;

  /// Position of subdomain names block
  Uint m_subdomain_names_filepos;

  /// Position of node coordinates block
  Uint m_node_coord_filepos;

  /// Position of element connectivity block
  Uint m_elements_filepos;

  /// Maximum dimension of subdomain in the mesh
  Uint m_max_topo_dim;

  /// Number of nodes
  Uint m_nb_nodes;

  /// Total number of elements
  Uint m_nb_elements;

  /// Array to hold coordinates
  std::vector<Real> m_coords;

  /// Information about subdomains: dimension, physical tag and name for each
  /// subdomain
  std::vector<SubDomainHeaderData> m_inner_domains;

  /// Information about subdomains: dimension, physical tag and name for each
  /// subdomain
  std::vector<SubDomainHeaderData> m_boundary_domains;

  /// Information about subdomains: count how many elements are there for each
  /// subdomain We have four maps which associat a physical tag with number of
  /// elements in each dimension
  std::array<std::map<Uint, Uint>, 4> m_tag_to_elem_count;

  /// Number of elements of each type in each subdomain
  std::array<Uint, PDEKitToGmsh::NbElemTypes> m_nb_elem_of_type;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshReader::read_mesh_from_file(const std::string &filename, Tria<MeshConfig> &mesh,
                                     const std::string &geo_dofs_name)
{
  reset();

  clock_t start, end;
  Real elapsed;

  start = clock();

  m_infile.open(filename.c_str());

  get_record_positions();
  read_coordinates(MeshConfig::GDIM);

  mesh.create_dof_storage(geo_dofs_name);
  read_connectivity<MeshConfig>(mesh, *mesh.dof_storage(geo_dofs_name));

  m_infile.close();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "Reading " << filename << " took " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshReader::read_mesh_from_file(const std::string &filename, Tria<MeshConfig> &mesh,
                                     typename result_of::dof_map_t<MeshConfig> &cells_dofs)
{
  reset();

  clock_t start, end;
  Real elapsed;

  start = clock();

  m_infile.open(filename.c_str());

  get_record_positions();
  read_coordinates(MeshConfig::GDIM);
  read_connectivity<MeshConfig>(mesh, cells_dofs);

  m_infile.close();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "Reading " << filename << " took " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshReader::read_connectivity(Tria<MeshConfig> &mesh,
                                   typename result_of::dof_map_t<MeshConfig> &interior_cells_dofs)
{
  std::cout << "[GmshReader] Reading connectivity ..." << std::endl;

  const Uint GDIM = Tria<MeshConfig>::GDIM;
  const Uint TDIM = Tria<MeshConfig>::TDIM;

  typename result_of::mesh_boundary_set_t<MeshConfig> &boundaries = mesh.all_boundaries();

  interior_cells_dofs.init(TDIM);

  /// Total number of elements in 0D,1D,2D and 3D
  Uint total_nb_cells[TDIM + 1];

  /// Number of connectivity entries in 0D,1D,2D, and 3D
  /// Example: for 10 triangles, the number will be 10 x 3 = 30 connectivity
  /// entries in 2D for 20 hexahedra and 30 tetras, the number will be
  /// (20 x 8) + (30 x 4) = 160 + 120 = 280 connectivity entries in 3D
  Uint nb_connections[TDIM + 1];

  /// Initialize
  for (Uint dim = 0; dim < (TDIM + 1); ++dim)
  {
    total_nb_cells[dim] = 0;
    nb_connections[dim] = 0;
  }

  /// Loop over all element types present in the mesh. Resize corresponding
  /// connectivity tables

  PDEKitToGmsh pdekit_to_gmsh;

  for (Uint et = 0; et < pdekit_to_gmsh.nb_elem_types(); ++et)
  {
    if (m_nb_elem_of_type[et] > 0)
    {
      const Uint d = pdekit_to_gmsh.elem_dim(et);
      total_nb_cells[d] += m_nb_elem_of_type[et];
      nb_connections[d] += m_nb_elem_of_type[et] * pdekit_to_gmsh.nb_nodes_in_elem(et);
    }
  }

  /// Temporary data to store the connectivity of the cells
  /// (triangles, quads ... in 2D, tetra, hexas, pyramids ... in 3D)

  Uint current_vertex = 0;

  std::unique_ptr<std::vector<Uint>> cell_connections(new std::vector<Uint>(nb_connections[TDIM]));

  // std::unique_ptr<common::BlockArray<Real, Uint>> node_coords(new
  // common::BlockArray<Real, Uint>()); node_coords->resize(GDIM *
  // nb_connections[TDIM], total_nb_cells[TDIM]);

  std::unique_ptr<std::vector<StdRegion>> cell_types(
      new std::vector<StdRegion>(total_nb_cells[TDIM]));

  std::unique_ptr<std::vector<Uint>> cell_tags(new std::vector<Uint>(total_nb_cells[TDIM]));

  /// Temporary cell arrays to store the boundary cells (1D and 2D cells in 3D
  /// mesh, 1D cells in 2D mesh)

  // On i-th  position, we store one mesh entity that represents an interior
  // cell
  CellBuffer<GDIM, TDIM> interior_cells;
  interior_cells.reserve(nb_connections[TDIM], total_nb_cells[TDIM]);

  /// On i-th position, we store one mesh entity that represents a boundary
  /// cell
  std::vector<std::unique_ptr<CellBuffer<GDIM, TDIM>>> boundary_cells(TDIM);
  for (Uint d = 0; d < TDIM; ++d)
  {
    boundary_cells[d] = std::unique_ptr<CellBuffer<GDIM, TDIM>>(new CellBuffer<GDIM, TDIM>());
  }

  /// The vector of pairs 'correspondence_vector' has the following
  /// information: on i-th position, we have a pair [c,li], which stores a
  /// cell number [c] and a local index [li] which uniquely defines the
  /// sub-cell in the mesh which corresponds to the i-th boundary cell in the
  /// vector boundary_cells Example: if correspondence_vector[37] = < 274, 2
  /// >, it means that the 37th boundary cell in the gmsh file is in fact the
  /// second face of the volume cell number 274
  std::vector<IncidenceEntry> correspondence_vector;

  for (Uint d = 0; d < TDIM; ++d)
  {
    boundary_cells[d]->reserve(nb_connections[d], total_nb_cells[d]);
  }

  /// Vectors which store physical tags for boundary patches
  std::vector<Uint> domain_tag[TDIM + 1];

  /// Clear the information about how many cells have been read:
  /// we are going to read & count them again

  for (Uint d = 1; d < TDIM + 1; ++d)
  {
    domain_tag[d].resize(total_nb_cells[d]);
    total_nb_cells[d] = 0;
  }

  /// Temporary vector to load connectivity of one cell
  std::vector<Uint> vertices_in_one_cell;
  /// Temporary vector to load coordinates of one cell
  std::vector<Real> coordinates_in_one_cell;

  tempstr = "";

  m_infile.seekg(m_elements_filepos, std::ios::beg);
  getline(m_infile, tempstr);
  getline(m_infile, tempstr);

  Uint eidx, etype, nb_tags, phys_tag, other_tag;

  for (Uint ie = 0; ie < m_nb_elements; ++ie)
  {
    m_infile >> eidx;
    eidx--;
    m_infile >> etype;
    m_infile >> nb_tags;
    m_infile >> phys_tag;

    for (Uint i = 0; i < nb_tags - 1; ++i)
    {
      m_infile >> other_tag;
    }

    const Uint dim = pdekit_to_gmsh.elem_dim(etype);
    const PointSetTag cell_type_id =
        PointSetTag(pdekit_to_gmsh.element_shape(etype), pdekit_to_gmsh.elem_order(etype),
                    PointSetID::Equidist);

    const Uint nb_nodes_in_this_cell = pdekit_to_gmsh.nb_nodes_in_elem(etype);

    /// Case 1: we are reading an internal cell
    if (dim == TDIM)
    {

      for (Uint i = 0; i < nb_nodes_in_this_cell; ++i)
      {
        m_infile >> (*cell_connections)[current_vertex];
        (*cell_connections)[current_vertex]--; // The numbering in
                                               // PDEKit starts from 0,
                                               // not like in gmsh,
                                               // where it starts with 1
        current_vertex++;
      }

      (*cell_types)[total_nb_cells[dim]].change_type(cell_type_id);
    }

    /// Case 2: we are reading a boundary cell
    else
    {
      vertices_in_one_cell.resize(nb_nodes_in_this_cell);
      coordinates_in_one_cell.resize(nb_nodes_in_this_cell * GDIM);

      for (Uint i = 0; i < nb_nodes_in_this_cell; ++i)
      {
        m_infile >> vertices_in_one_cell[i];
        vertices_in_one_cell[i]--;
        for (Uint d = 0; d < GDIM; ++d)
        {
          coordinates_in_one_cell[i * GDIM + d] = m_coords[vertices_in_one_cell[i] * GDIM + d];
        }
      }

      boundary_cells[dim]->push_back_cell(total_nb_cells[dim], cell_type_id, vertices_in_one_cell,
                                          phys_tag, coordinates_in_one_cell);
    }

    if (dim > 0)
    {
      domain_tag[dim][total_nb_cells[dim]] = phys_tag;
    }

    total_nb_cells[dim]++;

    getline(m_infile, tempstr);
  } /// Loop over connectivity data in the gmsh file

  // Build a block array of cell coordinates
  std::unique_ptr<std::vector<Real>> coord_values(new std::vector<Real>());
  std::unique_ptr<std::vector<Uint>> coord_block_sizes(new std::vector<Uint>());

  coord_values->resize(MeshConfig::GDIM * cell_connections->size());
  coord_block_sizes->resize(cell_types->size());

  Uint next_coord_entry = 0;

  for (Uint i = 0; i < cell_connections->size(); ++i)
  {
    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      (*coord_values)[next_coord_entry + d] =
          m_coords[(*cell_connections)[i] * MeshConfig::GDIM + d];
    }
    next_coord_entry += MeshConfig::GDIM;
  }

  for (Uint i = 0; i < cell_types->size(); ++i)
  {
    (*coord_block_sizes)[i] = MeshConfig::GDIM * (*cell_types)[i].get().nb_nodes();
  }

  std::unique_ptr<common::BlockArray<Real, Uint>> cell_coords(new common::BlockArray<Real, Uint>());
  cell_coords->build(std::move(coord_values), std::move(coord_block_sizes));

  interior_cells.emplace_cell_data(std::move(cell_connections), std::move(cell_types),
                                   std::move(cell_tags), std::move(cell_coords));
  // interior_cells.print(true);

  /// Fill the main cell storage
  // topology.create_from_cells(interior_cells);
  mesh.create_from_cells(interior_cells);
  interior_cells_dofs.create_from_cells(interior_cells);

  std::vector<std::pair<Uint, std::string>> cell_tag_names;
  for (Uint cell_dom = 0; cell_dom < m_inner_domains.size(); ++cell_dom)
  {
    cell_tag_names.push_back(std::pair<Uint, std::string>(m_inner_domains[cell_dom].idx,
                                                          m_inner_domains[cell_dom].name));
  }

  interior_cells_dofs.tag_all_active_cells(cell_tag_names, domain_tag[TDIM]);

  // Debugging output:
  // topology.cells().list_types();

  /// Fill the boundary cells
  // boundary_cells[1].print(true);
  // boundary_cells[0].print(true);

  for (Uint idomain = 0; idomain < m_boundary_domains.size(); ++idomain)
  {
    boundaries.create(m_boundary_domains[idomain].topo_dim, m_boundary_domains[idomain].name);
  }

  // boundaries.list();

  std::unique_ptr<std::vector<IncidenceEntry>> parent_cell_id(new std::vector<IncidenceEntry>());

  /// Find all boundary edges and faces
  for (Uint d = 1; d < TDIM; ++d)
  {
    std::cout << "Looking for boundary cells of dimension " << d << ":" << std::endl;
#if 1
    TopologyAlgorithms::identify_subcells2(interior_cells_dofs, d, *boundary_cells[d],
                                           correspondence_vector);
#else
    using interior_dof_iter_t = typename result_of::dof_map_t<MeshConfig>::const_dof_iterator;
    using boundary_dof_iter_t = typename CellBuffer<GDIM, TDIM>::const_dof_iterator;

    common::IteratorRange<interior_dof_iter_t> range_interior(interior_cells_dofs.cbegin(),
                                                              interior_cells_dofs.cend());
    common::IteratorRange<boundary_dof_iter_t> range_bdry(boundary_cells[d]->cbegin(),
                                                          boundary_cells[d]->cend());

    TopologyAlgorithms::identify_subcells3(range_interior, range_bdry, d, correspondence_vector);
#endif

    const std::map<Uint, Uint> &tag_to_elem_count_map = m_tag_to_elem_count[d];

    for (Uint idomain = 0; idomain < m_boundary_domains.size(); ++idomain)
    {
      if (m_boundary_domains[idomain].topo_dim == d)
      {
        std::map<Uint, Uint>::const_iterator tag_to_elem_count_iter =
            tag_to_elem_count_map.find(m_boundary_domains[idomain].idx);

        const Uint nb_cells_in_bdry_patch = tag_to_elem_count_iter->second;

        parent_cell_id->reserve(nb_cells_in_bdry_patch);
        parent_cell_id->resize(0);

        for (Uint c = 0; c < correspondence_vector.size(); ++c)
        {
          /*
          if (domain_tag[d][c] == m_boundary_domains[idomain].idx)
          {
            parent_cell_id->push_back(IncidencePair(correspondence_vector[c].cell_idx,
                                                   correspondence_vector[c].local_id));
          }
          */
          if ((boundary_cells[d]->active_cell_tag(ActiveIdx(c)) ==
               m_boundary_domains[idomain].idx) &&
              (correspondence_vector[c].cell_idx != INVALID_CELL_ID))
          {
            parent_cell_id->push_back(correspondence_vector[c]);
          }
        }

        std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdry_dofs =
            boundaries.domain(m_boundary_domains[idomain].name, d);

        bdry_dofs->emplace_bdry_cell_ids(std::move(parent_cell_id),
                                         m_boundary_domains[idomain].idx);
      }
    } // Loop over boundary domains
  }

  std::cout << "Number of cells = " << mesh.nb_all_cells_in_all_levels() << std::endl;

  std::cout << "[GmshReader] ... finished reading connectivity" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename DataMatrix>
void GmshReader::read_nodal_function_from_file(const Tria<MeshConfig> &mesh,
                                               const std::string &filename, DataMatrix &data,
                                               const std::vector<std::string> &field_names,
                                               const Real time, const Uint timestep)
{

  const Uint nb_nodes = mesh.geometry().nb_nodes();
  data.resize(field_names.size(), nb_nodes);

  std::string current_field_name = "";
  Uint field_idx, nb_string_tags, nb_real_tags, nb_int_tags, int_tag, nb_node_entries, entry_idx;
  Real real_tag, data_value;

  m_infile.open(filename.c_str());

  tempstr = "";

  while (getline(m_infile, tempstr))
  {
    if (tempstr == "$NodeData")
    {
      m_infile >> nb_string_tags;
      m_infile >> current_field_name;

      bool this_field_should_be_read = false;
      field_idx                      = 0;

      for (Uint f = 0; f < field_names.size(); ++f)
      {
        if (current_field_name == field_names[f])
        {
          this_field_should_be_read = true;
          field_idx                 = f;
          break;
        }
      }

      if (this_field_should_be_read == true)
      {
        // Read real tags
        m_infile >> nb_real_tags;
        for (Uint rt = 0; rt < nb_real_tags; ++rt)
        {
          m_infile >> real_tag;
        }

        // Read integer tags
        m_infile >> nb_int_tags;
        for (Uint it = 0; it < nb_int_tags; ++it)
        {
          m_infile >> int_tag;
        }

        // The last integer tag is the number of node entries
        nb_node_entries = int_tag;

        for (Uint i = 0; i < nb_node_entries; ++i)
        {
          m_infile >> entry_idx;
          entry_idx--;

          typename DataMatrix::block_type node_data = data.value(entry_idx);

          m_infile >> data_value;
          node_data[field_idx] = data_value;
        }

      } // If this field should be read

    } // If this is a $NodeData block
  }

  m_infile.close();
}

// ----------------------------------------------------------------------------

} // Namespace gmsh

} // Namespace mesh

} // Namespace pdekit

#endif
