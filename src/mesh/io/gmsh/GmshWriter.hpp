#ifndef PDEKIT_GMSH_Gmsh_Writer_hpp
#define PDEKIT_GMSH_Gmsh_Writer_hpp

#include <ctime>
#include <fstream>

#include "common/StringUtils.hpp"
#include "mesh/MeshPredicates.hpp"
#include "mesh/Tria.hpp"
#include "mesh/io/gmsh/PDEKitToGmsh.hpp"
#include "mesh/iterators/CellDofIterator.hpp"

namespace pdekit
{

namespace interpolation
{
// Forward declaration
// Forward declaration
template <typename MFT>
class ScalarMeshFunctionBase;

template <typename MFT>
class VectorMeshFunctionBase;
} // namespace interpolation

namespace mesh
{

namespace gmsh
{

class GmshWriter
{
  public:
  /// Constructor
  GmshWriter();

  /// Destructor
  ~GmshWriter();

  /// Write aadaptive mesh to a file
  /// @param: mesh_root mesh which should be written
  /// @param: filename to which file the mesh data should be written
  template <typename MeshConfig>
  void write_mesh_to_file(const Tria<MeshConfig> &mesh, const std::string &filename);

  /// Write aadaptive mesh to a file
  /// @param: mesh_root mesh which should be written
  /// @param: filename to which file the mesh data should be written
  template <typename MeshConfig>
  void write_mesh_to_file(const Tria<MeshConfig> &mesh, const std::string &dof_handler_name,
                          const std::string &filename);

  /// Append field to file
  /// @param: field pointer to the field which should be written
  /// @param: name of the file to which the data is written
  template <typename MeshConfig, typename MFT>
  void append_nodal_function_to_file(
      const Tria<MeshConfig> &mesh, const std::string &filename,
      const interpolation::VectorMeshFunctionBase<MFT> &mesh_function,
      const std::string &field_name, const Real time = 0.0, const Uint timestep = 0);
  /// Append scalar function to file
  template <typename MeshConfig, typename MFT>
  void append_nodal_function_to_file(
      const Tria<MeshConfig> &mesh, const std::string &filename,
      const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function,
      const std::string &field_name, const Real time = 0.0, const Uint timestep = 0);

  /// Append field to file
  /// @param: field pointer to the field which should be written
  /// @param: name of the file to which the data is written
  template <typename MeshConfig, typename MFT>
  void append_cell_function_to_file(const Tria<MeshConfig> &mesh, const std::string &filename,
                                    const interpolation::VectorMeshFunctionBase<MFT> &mesh_function,
                                    const std::string &field_name, const Real time = 0.0,
                                    const Uint timestep = 0);

  /// Append field to file
  /// @param: field pointer to the field which should be written
  /// @param: name of the file to which the data is written
  template <typename MeshConfig, typename MFT>
  void append_cell_function_to_file(const Tria<MeshConfig> &mesh, const std::string &filename,
                                    const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function,
                                    const std::string &field_name, const Real time = 0.0,
                                    const Uint timestep = 0);

  /// Append field to file
  /// @param: field pointer to the field which should be written
  /// @param: name of the file to which the data is written
  template <typename MeshConfig, typename T>
  void append_cell_function_to_file(const Tria<MeshConfig> &mesh, const std::string &filename,
                                    const common::ArrayView<const T, _1D, Uint> &mesh_function,
                                    const std::string &field_name, const Real time = 0.0,
                                    const Uint timestep = 0);

  template <typename MeshConfig, typename MFT>
  void save_mesh_boundary(const Tria<MeshConfig> &mesh,
                          const typename result_of::dof_map_t<MeshConfig> &cell_dofs,
                          const std::string &filename,
                          const interpolation::VectorMeshFunctionBase<MFT> &mesh_function,
                          const std::string &field_name,
                          const std::vector<std::string> &domain_names, const Real time = 0.0,
                          const Uint timestep = 0);

  template <typename MeshConfig, typename MFT>
  void save_mesh_boundary(const Tria<MeshConfig> &mesh,
                          const typename result_of::dof_map_t<MeshConfig> &cell_dofs,
                          const std::string &filename,
                          const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function,
                          const std::string &field_name,
                          const std::vector<std::string> &domain_names, const Real time = 0.0,
                          const Uint timestep = 0);

  private:
  /// METHODS

  template <typename MeshConfig>
  void write_header(const Tria<MeshConfig> &mesh);

  template <typename MeshConfig>
  void write_header(typename result_of::dof_map_t<MeshConfig> const &mesh_interior_dofs,
                    MeshBoundarySet<MeshConfig> const &mesh_boundary);

  template <typename MeshConfig>
  void write_coordinates(const Tria<MeshConfig> &mesh);

  template <typename MeshConfig>
  void write_coordinates(typename result_of::dof_map_t<MeshConfig> const &cell_dofs);

  template <typename MeshConfig>
  void write_connectivity(const Tria<MeshConfig> &mesh);

  template <typename MeshConfig>
  void write_connectivity(typename result_of::dof_map_t<MeshConfig> const &cells,
                          MeshBoundarySet<MeshConfig> const &mesh_boundary);

  /// DATA

  /// Input file stream
  std::ofstream m_outfile;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshWriter::write_mesh_to_file(const Tria<MeshConfig> &mesh, const std::string &filename)
{
  clock_t start, end;
  Real elapsed;

  start = clock();

  m_outfile.open(filename.c_str());
  write_header(mesh);
  write_coordinates(mesh);
  write_connectivity(mesh);

  m_outfile.close();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "Writing " << filename << " took " << elapsed << " s" << std::endl;

  start = clock();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshWriter::write_mesh_to_file(const Tria<MeshConfig> &mesh,
                                    const std::string &dof_handler_name,
                                    const std::string &filename)
{
  clock_t start, end;
  Real elapsed;

  start = clock();

  m_outfile.open(filename.c_str());

  const common::PtrHandle<typename Tria<MeshConfig>::dof_storage_type const> dofs_ptr =
      mesh.dof_storage(dof_handler_name);

  write_header(*dofs_ptr, mesh.all_boundaries());
  write_coordinates<MeshConfig>(*dofs_ptr);
  write_connectivity(*dofs_ptr, mesh.all_boundaries());

  m_outfile.close();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "Writing " << filename << " took " << elapsed << " s" << std::endl;

  start = clock();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename MFT>
void GmshWriter::append_nodal_function_to_file(
    const Tria<MeshConfig> &mesh, const std::string &filename,
    const interpolation::VectorMeshFunctionBase<MFT> &mesh_function, const std::string &field_name,
    const Real time, const Uint timestep)
{
  typedef typename interpolation::VectorMeshFunctionBase<MFT>::MeshFunctionType MeshFunctionType;
  MeshFunctionType const &data = mesh_function.wrapped_type();

  m_outfile.open(filename.c_str(), std::ios::app);

  for (Uint component = 0; component < data.nb_fields(); ++component)
  {
    m_outfile << "$NodeData\n1" << std::endl;
    m_outfile << "\"" << field_name + "_" + common::StringUtils::to_string(component) << "\""
              << std::endl;

    m_outfile << "1" << std::endl; // Number of real tags
    m_outfile << time << std::endl;

    m_outfile << "3" << std::endl; // Number of integer tags
    m_outfile << timestep << std::endl;

    m_outfile << "1" << std::endl;
    m_outfile << data.nb_entries() << std::endl;

    m_outfile.precision(18);

    for (Uint n = 0; n < data.nb_entries(); ++n)
    {
      typename MeshFunctionType::const_entry_type const nodal_value = data.const_value(n);
      m_outfile << n + 1 << " " << nodal_value[component] << std::endl;
    }

    m_outfile << "$EndNodeData" << std::endl;
  }
  m_outfile.close();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename MFT>
void GmshWriter::append_nodal_function_to_file(
    const Tria<MeshConfig> &mesh, const std::string &filename,
    const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function, const std::string &field_name,
    const Real time, const Uint timestep)
{
  typedef typename interpolation::ScalarMeshFunctionBase<MFT>::MeshFunctionType MeshFunctionType;
  MeshFunctionType const &data = mesh_function.wrapped_type();

  m_outfile.open(filename.c_str(), std::ios::app);

  m_outfile << "$NodeData\n1" << std::endl;
  m_outfile << "\"" << field_name << "\"" << std::endl;

  m_outfile << "1" << std::endl; // Number of real tags
  m_outfile << time << std::endl;

  m_outfile << "3" << std::endl; // Number of integer tags
  m_outfile << timestep << std::endl;

  m_outfile << "1" << std::endl;
  m_outfile << data.nb_entries() << std::endl;

  m_outfile.precision(18);

  for (Uint n = 0; n < data.nb_entries(); ++n)
  {
    m_outfile << n + 1 << " " << data[n] << std::endl;
  }

  m_outfile << "$EndNodeData" << std::endl;

  m_outfile.close();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename MFT>
void GmshWriter::append_cell_function_to_file(
    const Tria<MeshConfig> &mesh, const std::string &filename,
    const interpolation::VectorMeshFunctionBase<MFT> &mesh_function, const std::string &field_name,
    const Real time, const Uint timestep)
{
  typedef typename interpolation::VectorMeshFunctionBase<MFT>::MeshFunctionType MeshFunctionType;
  MeshFunctionType const &data = mesh_function.wrapped_type();

  m_outfile.open(filename.c_str(), std::ios::app);

  for (Uint component = 0; component < data.nb_fields(); ++component)
  {
    m_outfile << "$ElementData\n1" << std::endl;
    m_outfile << "\"" << field_name + "_" + common::StringUtils::to_string(component) << "\""
              << std::endl;

    m_outfile << "1" << std::endl; // Number of real tags
    m_outfile << time << std::endl;

    m_outfile << "3" << std::endl; // Number of integer tags
    m_outfile << timestep << std::endl;

    m_outfile << "1" << std::endl;

    MeshBoundarySet<MeshConfig> const &mesh_boundary = mesh.all_boundaries();

    Uint nb_bdry_elements = 0;

    // Count all boundary cells
    for (const typename std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdomain :
         mesh_boundary.all_domains())
    {
      nb_bdry_elements += bdomain->nb_active_cells();
    }

    const Uint total_nb_elements = nb_bdry_elements + data.nb_entries();

    m_outfile << total_nb_elements << std::endl;
    m_outfile.precision(18);

    // Write dummy data for boundary elements
    for (Uint n = 0; n < nb_bdry_elements; ++n)
    {
      m_outfile << n + 1 << " " << 0.0 << std::endl;
    }

    /*
    for (Uint n = 0; n < data.nb_entries(); ++n)
    {
      typename DataMatrix::const_entry_type const cell_value =
    data.const_value(n); m_outfile << n + 1 << " " << cell_value[component]
    << std::endl;
    }
    */

    for (Uint n = 0; n < data.nb_entries(); ++n)
    {
      typename MeshFunctionType::const_entry_type const cell_value = data.const_value(n);
      m_outfile << nb_bdry_elements + n + 1 << " " << cell_value[component] << std::endl;
    }

    m_outfile << "$EndElementData" << std::endl;
  }
  m_outfile.close();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename MFT>
void GmshWriter::append_cell_function_to_file(
    const Tria<MeshConfig> &mesh, const std::string &filename,
    const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function, const std::string &field_name,
    const Real time, const Uint timestep)
{
  typename interpolation::ScalarMeshFunctionBase<MFT>::MeshFunctionType const &data =
      mesh_function.wrapped_type();

  m_outfile.open(filename.c_str(), std::ios::app);

  for (Uint component = 0; component < data.nb_fields(); ++component)
  {
    m_outfile << "$ElementData\n1" << std::endl;
    m_outfile << "\"" << field_name + "_" + common::StringUtils::to_string(component) << "\""
              << std::endl;

    m_outfile << "1" << std::endl; // Number of real tags
    m_outfile << time << std::endl;

    m_outfile << "3" << std::endl; // Number of integer tags
    m_outfile << timestep << std::endl;

    m_outfile << "1" << std::endl;

    MeshBoundarySet<MeshConfig> const &mesh_boundary = mesh.all_boundaries();

    Uint nb_bdry_elements = 0;

    // Count all boundary cells
    for (const typename std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdomain :
         mesh_boundary.all_domains())
    {
      nb_bdry_elements += bdomain->nb_active_cells();
    }

    const Uint total_nb_elements = nb_bdry_elements + data.nb_entries();

    m_outfile << total_nb_elements << std::endl;
    m_outfile.precision(18);

    // Write dummy data for boundary elements
    for (Uint n = 0; n < nb_bdry_elements; ++n)
    {
      m_outfile << n + 1 << " " << 0.0 << std::endl;
    }

    /*
    for (Uint n = 0; n < data.nb_entries(); ++n)
    {
      typename DataMatrix::const_entry_type const cell_value =
    data.const_value(n); m_outfile << n + 1 << " " << cell_value[component]
    << std::endl;
    }
    */

    for (Uint n = 0; n < data.nb_entries(); ++n)
    {
      m_outfile << nb_bdry_elements + n + 1 << " " << data[n] << std::endl;
    }

    m_outfile << "$EndElementData" << std::endl;
  }
  m_outfile.close();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename T>
void GmshWriter::append_cell_function_to_file(
    const Tria<MeshConfig> &mesh, const std::string &filename,
    const common::ArrayView<const T, _1D, Uint> &mesh_function, const std::string &field_name,
    const Real time, const Uint timestep)
{
  m_outfile.open(filename.c_str(), std::ios::app);

  m_outfile << "$ElementData\n1" << std::endl;
  m_outfile << "\"" << field_name + "\"" << std::endl;

  m_outfile << "1" << std::endl; // Number of real tags
  m_outfile << time << std::endl;

  m_outfile << "3" << std::endl; // Number of integer tags
  m_outfile << timestep << std::endl;

  m_outfile << "1" << std::endl;

  MeshBoundarySet<MeshConfig> const &mesh_boundary = mesh.all_boundaries();

  Uint nb_bdry_elements = 0;

  // Count all boundary cells
  for (const typename std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdomain :
       mesh_boundary.all_domains())
  {
    nb_bdry_elements += bdomain->nb_active_cells();
  }

  const Uint total_nb_elements = nb_bdry_elements + mesh_function.size();

  m_outfile << total_nb_elements << std::endl;
  m_outfile.precision(18);

  // Write dummy data for boundary elements
  for (Uint n = 0; n < nb_bdry_elements; ++n)
  {
    m_outfile << n + 1 << " " << 0.0 << std::endl;
  }

  /*
  for (Uint n = 0; n < data.nb_entries(); ++n)
  {
    typename DataMatrix::const_entry_type const cell_value =
  data.const_value(n); m_outfile << n + 1 << " " << cell_value[component] <<
  std::endl;
  }
  */

  for (Uint n = 0; n < mesh_function.size(); ++n)
  {
    m_outfile << nb_bdry_elements + n + 1 << " " << mesh_function[n] << std::endl;
  }

  m_outfile << "$EndElementData" << std::endl;

  m_outfile.close();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename MFT>
void GmshWriter::save_mesh_boundary(const Tria<MeshConfig> &mesh,
                                    const typename result_of::dof_map_t<MeshConfig> &cell_dofs,
                                    const std::string &filename,
                                    const interpolation::VectorMeshFunctionBase<MFT> &mesh_function,
                                    const std::string &field_name,
                                    const std::vector<std::string> &domain_names, const Real time,
                                    const Uint timestep)
{
  typename interpolation::VectorMeshFunctionBase<MFT>::MeshFunctionType const &data =
      mesh_function.wrapped_type();

  clock_t start, end;
  Real elapsed;

  start = clock();

  std::map<Uint, std::string> bdry_domain_tags;
  Uint total_nb_elems = 0;

  const MeshBoundarySet<MeshConfig> &mesh_boundary = mesh.all_boundaries();

  Uint tot_nb_vert_on_boundary = 0;

  for (const auto &name : domain_names)
  {
    typename result_of::mesh_boundary_set_t<MeshConfig>::bdry_facets_shared_ptr bdomain =
        mesh_boundary.domain(name);

    total_nb_elems += bdomain->nb_active_cells();
    bdry_domain_tags.insert(std::pair<Uint, std::string>(bdomain->material_id(), bdomain->name()));

    for (typename BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_dof_iterator it =
             bdomain->cbegin(cell_dofs);
         it != bdomain->cend(cell_dofs); ++it)
    {
      const MeshEntity cell = it->mesh_entity();
      tot_nb_vert_on_boundary += cell.nb_vert();
    }
  } // Loop over boundary domains

  m_outfile.open(filename.c_str());

  m_outfile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames" << std::endl;

  m_outfile << domain_names.size() << std::endl;

  for (const auto &name : domain_names)
  {
    typename result_of::mesh_boundary_set_t<MeshConfig>::bdry_facets_shared_ptr bdomain =
        mesh_boundary.domain(name);

    m_outfile << bdomain->dim() << " " << bdomain->material_id() << " \"" << bdomain->name() << "\""
              << std::endl;
  }

  m_outfile << "$EndPhysicalNames" << std::endl;
  m_outfile << "$Nodes" << std::endl;
  m_outfile << tot_nb_vert_on_boundary << std::endl;
  m_outfile.precision(14);

  Uint node_id = 0;

  adapt::LocalInterpolator loc_interpolator;

  for (const auto &name : domain_names)
  {
    typename result_of::mesh_boundary_set_t<MeshConfig>::bdry_facets_shared_ptr bdomain =
        mesh_boundary.domain(name);
    {
      for (typename BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_dof_iterator it =
               bdomain->cbegin(cell_dofs);
           it != bdomain->cend(cell_dofs); ++it)
      {
        const auto facet_geom        = it->cell_geometry();
        const PointSetTag geo_pt_set = it->geo_pt_set_id();
        const PointSetTag dof_pt_set = it->pt_set_id();

        const math::DenseConstMatView<Real> active_cell_coords =
            loc_interpolator.transfer_coords(geo_pt_set, dof_pt_set, facet_geom);

        for (Uint n = 0; n < active_cell_coords.rows(); ++n)
        {
          const math::DenseConstVecView<Real> node = active_cell_coords.row_transpose(n);
          m_outfile << node_id + 1;
          node_id++;

          for (Uint d = 0; d < MeshConfig::GDIM; ++d)
          {
            m_outfile << " " << node[d];
          }
          for (Uint d = MeshConfig::GDIM; d < _3D; ++d)
          {
            m_outfile << " 0.0";
          }
          m_outfile << std::endl;
        }
      }
    } // Loop over boundary domains
  }

  m_outfile << "$EndNodes" << std::endl;

  PDEKitToGmsh pdekit_to_gmsh;

  m_outfile << "$Elements" << std::endl;
  m_outfile << total_nb_elems << std::endl;

  total_nb_elems = 0;
  node_id        = 0;

  for (const auto &name : domain_names)
  {
    typename result_of::mesh_boundary_set_t<MeshConfig>::bdry_facets_shared_ptr bdomain =
        mesh_boundary.domain(name);
    {
      for (typename BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_dof_iterator it =
               bdomain->cbegin(cell_dofs);
           it != bdomain->cend(cell_dofs); ++it)
      {
        const MeshEntity active_cell = it->mesh_entity();
        const Uint gmsh_elem_type =
            pdekit_to_gmsh.ref_topology_type_to_gmsh_type(active_cell.pt_set_id());

        m_outfile << total_nb_elems + 1 << " " << gmsh_elem_type << " 2 " << bdomain->material_id()
                  << " " << bdomain->material_id();
        total_nb_elems++;

        for (Uint v = 0; v < active_cell.nb_vert(); ++v)
        {
          m_outfile << " " << node_id + 1;
          node_id++;
        }
        m_outfile << std::endl;
      }
    }
  }

  m_outfile << "$EndElements" << std::endl;

  typedef typename interpolation::VectorMeshFunctionBase<MFT>::MeshFunctionType function_type;
  typedef typename function_type::const_entry_type const_entry_type;

  for (Uint component = 0; component < data.nb_fields(); ++component)
  {
    m_outfile << "$NodeData\n1" << std::endl;
    m_outfile << "\"" << field_name + "_" + common::StringUtils::to_string(component) << "\""
              << std::endl;

    m_outfile << "1" << std::endl; // Number of real tags
    m_outfile << time << std::endl;

    m_outfile << "3" << std::endl; // Number of integer tags
    m_outfile << timestep << std::endl;

    m_outfile << "1" << std::endl;
    m_outfile << tot_nb_vert_on_boundary << std::endl;

    m_outfile.precision(18);

    node_id = 0;

    for (Uint idom = 0; idom < domain_names.size(); ++idom)
    {
      typename result_of::mesh_boundary_set_t<MeshConfig>::bdry_facets_shared_ptr bdomain =
          mesh_boundary.domain(domain_names[idom]);
      {
        for (typename BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_dof_iterator it =
                 bdomain->cbegin(cell_dofs);
             it != bdomain->cend(cell_dofs); ++it)
        {
          const MeshEntity active_cell = it->mesh_entity();

          for (Uint n = 0; n < active_cell.nb_vert(); ++n)
          {
            const_entry_type node_data = data.const_value(active_cell.vertex(n));
            m_outfile << node_id + 1 << " " << node_data[component] << std::endl;
            node_id++;
          } // Loop over vertices of one active cell

        } // Loop over one boundary domain
      }
    } // Loop over domain names

    m_outfile << "$EndNodeData" << std::endl;
  } // Loop over all components

  m_outfile.close();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "Writing " << filename << " took " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename MFT>
void GmshWriter::save_mesh_boundary(const Tria<MeshConfig> &mesh,
                                    const typename result_of::dof_map_t<MeshConfig> &cell_dofs,
                                    const std::string &filename,
                                    const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function,
                                    const std::string &field_name,
                                    const std::vector<std::string> &domain_names, const Real time,
                                    const Uint timestep)
{
  typename interpolation::ScalarMeshFunctionBase<MFT>::MeshFunctionType const &data =
      mesh_function.wrapped_type();

  clock_t start, end;
  Real elapsed;

  start = clock();

  std::map<Uint, std::string> bdry_domain_tags;
  Uint total_nb_elems = 0;

  const MeshBoundarySet<MeshConfig> &mesh_boundary = mesh.all_boundaries();

  Uint tot_nb_vert_on_boundary = 0;

  for (const auto &name : domain_names)
  {
    typename result_of::mesh_boundary_set_t<MeshConfig>::bdry_facets_shared_ptr bdomain =
        mesh_boundary.domain(name);

    total_nb_elems += bdomain->nb_active_cells();
    bdry_domain_tags.insert(std::pair<Uint, std::string>(bdomain->material_id(), bdomain->name()));

    for (typename BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_dof_iterator it =
             bdomain->cbegin(cell_dofs);
         it != bdomain->cend(cell_dofs); ++it)
    {
      const MeshEntity cell = it->mesh_entity();
      tot_nb_vert_on_boundary += cell.nb_vert();
    }
  } // Loop over boundary domains

  m_outfile.open(filename.c_str());

  m_outfile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames" << std::endl;

  m_outfile << domain_names.size() << std::endl;

  for (const auto &name : domain_names)
  {
    typename result_of::mesh_boundary_set_t<MeshConfig>::bdry_facets_shared_ptr bdomain =
        mesh_boundary.domain(name);

    m_outfile << bdomain->dim() << " " << bdomain->material_id() << " \"" << bdomain->name() << "\""
              << std::endl;
  }

  m_outfile << "$EndPhysicalNames" << std::endl;
  m_outfile << "$Nodes" << std::endl;
  m_outfile << tot_nb_vert_on_boundary << std::endl;
  m_outfile.precision(14);

  Uint node_id = 0;

  adapt::LocalInterpolator loc_interpolator;

  for (const auto &name : domain_names)
  {
    typename result_of::mesh_boundary_set_t<MeshConfig>::bdry_facets_shared_ptr bdomain =
        mesh_boundary.domain(name);
    {
      for (typename BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_dof_iterator it =
               bdomain->cbegin(cell_dofs);
           it != bdomain->cend(cell_dofs); ++it)
      {
        const auto facet_geom        = it->cell_geometry();
        const PointSetTag geo_pt_set = it->geo_pt_set_id();
        const PointSetTag dof_pt_set = it->pt_set_id();

        const math::DenseConstMatView<Real> active_cell_coords =
            loc_interpolator.transfer_coords(geo_pt_set, dof_pt_set, facet_geom);

        for (Uint n = 0; n < active_cell_coords.rows(); ++n)
        {
          const math::DenseConstVecView<Real> node = active_cell_coords.row_transpose(n);
          m_outfile << node_id + 1;
          node_id++;

          for (Uint d = 0; d < MeshConfig::GDIM; ++d)
          {
            m_outfile << " " << node[d];
          }
          for (Uint d = MeshConfig::GDIM; d < _3D; ++d)
          {
            m_outfile << " 0.0";
          }
          m_outfile << std::endl;
        }
      }
    } // Loop over boundary domains
  }

  m_outfile << "$EndNodes" << std::endl;

  PDEKitToGmsh pdekit_to_gmsh;

  m_outfile << "$Elements" << std::endl;
  m_outfile << total_nb_elems << std::endl;

  total_nb_elems = 0;
  node_id        = 0;

  for (const auto &name : domain_names)
  {
    typename result_of::mesh_boundary_set_t<MeshConfig>::bdry_facets_shared_ptr bdomain =
        mesh_boundary.domain(name);
    {
      for (typename BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_dof_iterator it =
               bdomain->cbegin(cell_dofs);
           it != bdomain->cend(cell_dofs); ++it)
      {
        const MeshEntity active_cell = it->mesh_entity();
        const Uint gmsh_elem_type =
            pdekit_to_gmsh.ref_topology_type_to_gmsh_type(active_cell.pt_set_id());

        m_outfile << total_nb_elems + 1 << " " << gmsh_elem_type << " 2 " << bdomain->material_id()
                  << " " << bdomain->material_id();
        total_nb_elems++;

        for (Uint v = 0; v < active_cell.nb_vert(); ++v)
        {
          m_outfile << " " << node_id + 1;
          node_id++;
        }
        m_outfile << std::endl;
      }
    }
  }

  m_outfile << "$EndElements" << std::endl;

  // typedef typename
  // interpolation::VectorMeshFunctionBase<MFT>::MeshFunctionType
  // function_type; typedef typename function_type::const_entry_type
  // const_entry_type;

  m_outfile << "$NodeData\n1" << std::endl;
  m_outfile << "\"" << field_name << "\"" << std::endl;

  m_outfile << "1" << std::endl; // Number of real tags
  m_outfile << time << std::endl;

  m_outfile << "3" << std::endl; // Number of integer tags
  m_outfile << timestep << std::endl;

  m_outfile << "1" << std::endl;
  m_outfile << tot_nb_vert_on_boundary << std::endl;

  m_outfile.precision(18);

  node_id = 0;

  for (Uint idom = 0; idom < domain_names.size(); ++idom)
  {
    typename result_of::mesh_boundary_set_t<MeshConfig>::bdry_facets_shared_ptr bdomain =
        mesh_boundary.domain(domain_names[idom]);
    {
      for (typename BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_dof_iterator it =
               bdomain->cbegin(cell_dofs);
           it != bdomain->cend(cell_dofs); ++it)
      {
        const MeshEntity active_cell = it->mesh_entity();

        for (Uint n = 0; n < active_cell.nb_vert(); ++n)
        {
          m_outfile << node_id + 1 << " " << data[active_cell.vertex(n)] << std::endl;
          node_id++;
        } // Loop over vertices of one active cell

      } // Loop over one boundary domain
    }
  } // Loop over domain names

  m_outfile << "$EndNodeData" << std::endl;

  m_outfile.close();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "Writing " << filename << " took " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshWriter::write_header(const Tria<MeshConfig> &mesh)
{
  // typedef typename result_of::dof_map_t<MeshConfig> cells_type;

  m_outfile << "$MeshFormat\n2.1 0 8\n$EndMeshFormat\n$PhysicalNames" << std::endl;

  const MeshBoundarySet<MeshConfig> &mesh_boundary = mesh.all_boundaries();

  m_outfile << 1 + mesh_boundary.all_domains().size() << std::endl;
  Uint bdry_idx = 1;

  for (std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdry_ptr :
       mesh_boundary.all_domains())
  {
    m_outfile << bdry_ptr->dim() << " " << bdry_idx << " \"" << bdry_idx << "\"" << std::endl;
    bdry_idx++;
  }

  m_outfile << mesh.topo_dim() << " " << bdry_idx << " \"" << bdry_idx << "\"" << std::endl;
  m_outfile << "$EndPhysicalNames" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshWriter::write_header(typename result_of::dof_map_t<MeshConfig> const &mesh_interior_dofs,
                              MeshBoundarySet<MeshConfig> const &mesh_boundary)
{
  // typedef typename result_of::dof_map_t<MeshConfig> cells_type;

  m_outfile << "$MeshFormat\n2.1 0 8\n$EndMeshFormat\n$PhysicalNames" << std::endl;

  m_outfile << mesh_interior_dofs.nb_cell_tag_types() + mesh_boundary.all_domains().size()
            << std::endl;

  for (std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdry_ptr :
       mesh_boundary.all_domains())
  {
    m_outfile << bdry_ptr->dim() << " " << bdry_ptr->material_id() << " \"" << bdry_ptr->name()
              << "\"" << std::endl;
  }

  std::map<Uint, std::string> cell_tag_names;
  mesh_interior_dofs.all_tag_names(cell_tag_names);

  for (std::map<Uint, std::string>::const_iterator cell_tag_iter = cell_tag_names.begin();
       cell_tag_iter != cell_tag_names.end(); ++cell_tag_iter)
  {
    m_outfile << mesh_interior_dofs.dim() << " " << cell_tag_iter->first << " \""
              << cell_tag_iter->second << "\"" << std::endl;
  }

  m_outfile << "$EndPhysicalNames" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshWriter::write_coordinates(const Tria<MeshConfig> &mesh)
{
  using node_view_t = typename CellGeometry<MeshConfig::GDIM - 1>::node_view_t;

  std::vector<Uint> node_counts;
  const MeshBoundarySet<MeshConfig> &mesh_boundary = mesh.all_boundaries();

  for (std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdry_ptr :
       mesh_boundary.all_domains())
  {
    Uint bdry_node_count = 0;

    for (Uint facet_id = 0; facet_id < bdry_ptr->nb_active_cells(); ++facet_id)
    {
      const FlatIdx cell_id                    = bdry_ptr->parent_cell_id(ActiveIdx(facet_id));
      const Uint local_id                      = bdry_ptr->local_id(ActiveIdx(facet_id));
      const CellTopologyView<MeshConfig> tcell = mesh.cell(cell_id);
      const std::shared_ptr<const StdRegionEntity> std_reg_entity =
          tcell.sub_entity(mesh.topo_dim() - 1, local_id);
      bdry_node_count += std_reg_entity->nb_vert();
    }
    node_counts.push_back(bdry_node_count);
  }

  // Count interior nodes
  Uint nb_interior_nodes = 0;

  for (Uint ac = 0; ac < mesh.nb_active_cells(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh.active_cell(ActiveIdx(ac));
    nb_interior_nodes += tcell.nb_vert();
  }

  node_counts.push_back(nb_interior_nodes);
  const Uint tot_nb_nodes = std::accumulate(node_counts.cbegin(), node_counts.cend(), 0);

  // Write boundary node coordinates to file
  Uint node_id = 1;

  m_outfile.precision(17);
  m_outfile << "$Nodes" << std::endl;
  m_outfile << tot_nb_nodes << std::endl;

  for (std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdry_ptr :
       mesh_boundary.all_domains())
  {
    for (Uint facet_id = 0; facet_id < bdry_ptr->nb_active_cells(); ++facet_id)
    {
      const FlatIdx cell_id                    = bdry_ptr->parent_cell_id(ActiveIdx(facet_id));
      const Uint local_id                      = bdry_ptr->local_id(ActiveIdx(facet_id));
      const CellTopologyView<MeshConfig> tcell = mesh.cell(cell_id);
      const CellGeometry<MeshConfig::GDIM> facet_coords =
          tcell.coordinates(mesh.topo_dim() - 1, local_id);
      for (Uint i = 0; i < facet_coords.size(); ++i)
      {
        m_outfile << node_id++;
        const node_view_t node_coord = facet_coords.const_node_view(i);
        for (Uint d = 0; d < node_coord.size(); ++d)
        {
          m_outfile << " " << node_coord[d];
        }
        for (Uint d = node_coord.size(); d < _3D; ++d)
        {
          m_outfile << " 0.0";
        }
        m_outfile << std::endl;
      }
    }
  } // Loop over all mesh boundaries

  // Write interior node coordinates to file
  for (Uint ac = 0; ac < mesh.nb_active_cells(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell         = mesh.active_cell(ActiveIdx(ac));
    const CellGeometry<MeshConfig::GDIM> cell_coords = tcell.coordinates();

    for (Uint i = 0; i < cell_coords.size(); ++i)
    {
      m_outfile << node_id++;
      const node_view_t node_coord = cell_coords.const_node_view(i);
      for (Uint d = 0; d < node_coord.size(); ++d)
      {
        m_outfile << " " << node_coord[d];
      }
      for (Uint d = node_coord.size(); d < _3D; ++d)
      {
        m_outfile << " 0.0";
      }
      m_outfile << std::endl;
    }
  } // Loop over all active mesh cells
  m_outfile << "$EndNodes" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshWriter::write_coordinates(typename result_of::dof_map_t<MeshConfig> const &cell_dofs)
{
  m_outfile << "$Nodes" << std::endl;

  const Uint nb_nodes = cell_dofs.nb_nodes();

  m_outfile << nb_nodes << std::endl;

  m_outfile.precision(17);
  //    m_outfile.setf(std::ios::fixed);

  std::vector<Real> raw_coords;
  raw_coords.resize(nb_nodes * MeshConfig::GDIM);

  adapt::LocalInterpolator loc_interpolator;

  for (typename result_of::dof_map_t<MeshConfig>::const_dof_iterator it = cell_dofs.cbegin();
       it != cell_dofs.cend(); ++it)
  {
    const CellTopologyView<MeshConfig> tcell_view = it->tcell();
    const MeshEntity active_cell                  = it->mesh_entity();

    const math::DenseConstMatView<Real> active_cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), active_cell.pt_set_id(), tcell_view.coordinates());

    for (Uint v = 0; v < active_cell_coords.rows(); ++v)
    {
      const math::DenseConstVecView<Real> node_coord = active_cell_coords.row_transpose(v);

      const Uint offset = MeshConfig::GDIM * active_cell.vertex(v);
      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        raw_coords[offset + d] = node_coord[d];
      }
    }
  }

  for (Uint node = 0; node < nb_nodes; ++node)
  {
    m_outfile << node + 1;
    const Uint offset = MeshConfig::GDIM * node;
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

  m_outfile << "$EndNodes" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshWriter::write_connectivity(const Tria<MeshConfig> &mesh)
{
  // Write boundary node coordinates to file
  Uint elem_id = 1;
  Uint node_id = 1;

  Uint nb_active_elems = 0;

  const MeshBoundarySet<MeshConfig> &mesh_boundary = mesh.all_boundaries();

  for (std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdry_ptr :
       mesh_boundary.all_domains())
  {
    nb_active_elems += bdry_ptr->nb_active_cells();
  }

  nb_active_elems += mesh.nb_active_cells();

  m_outfile << "$Elements" << std::endl;
  m_outfile << nb_active_elems << std::endl;

  PDEKitToGmsh pdekit_to_gmsh;

  Uint bdry_id = 1;

  for (std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdry_ptr :
       mesh_boundary.all_domains())
  {
    for (Uint facet_id = 0; facet_id < bdry_ptr->nb_active_cells(); ++facet_id)
    {
      const FlatIdx cell_id                    = bdry_ptr->parent_cell_id(ActiveIdx(facet_id));
      const Uint local_id                      = bdry_ptr->local_id(ActiveIdx(facet_id));
      const CellTopologyView<MeshConfig> tcell = mesh.cell(cell_id);

      const std::shared_ptr<const StdRegionEntity> std_reg_entity =
          tcell.sub_entity(mesh.topo_dim() - 1, local_id);

      const Uint nb_facet_vert = std_reg_entity->nb_vert();

      const Uint gmsh_type =
          pdekit_to_gmsh.ref_topology_type_to_gmsh_type(std_reg_entity->pt_set_id());
      m_outfile << elem_id++ << " " << gmsh_type << " 2 " << bdry_id << " " << bdry_id;

      for (Uint v = 0; v < nb_facet_vert; ++v)
      {
        m_outfile << " " << node_id++;
      }
      m_outfile << std::endl;
    }

    bdry_id++;
  }

  // Write interior node coordinates to file
  const Uint interior_id = mesh_boundary.nb_domains() + 1;

  for (Uint ac = 0; ac < mesh.nb_active_cells(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh.active_cell(ActiveIdx(ac));
    /// Return the cell type of this cell
    const StdRegion active_cell_type = tcell.std_region();

    const Uint nb_cell_vert = tcell.nb_vert();
    const Uint gmsh_type =
        pdekit_to_gmsh.ref_topology_type_to_gmsh_type(active_cell_type.get().pt_set_id());

    m_outfile << elem_id++ << " " << gmsh_type << " 2 " << interior_id << " " << interior_id;

    for (Uint v = 0; v < nb_cell_vert; ++v)
    {
      m_outfile << " " << node_id++;
    }
    m_outfile << std::endl;
  } // Loop over all active mesh cells
  m_outfile << "$EndElements" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void GmshWriter::write_connectivity(typename result_of::dof_map_t<MeshConfig> const &cells,
                                    MeshBoundarySet<MeshConfig> const &mesh_boundary)
{

  PDEKitToGmsh pdekit_to_gmsh;

  m_outfile << "$Elements" << std::endl;

  const Uint total_nb_elements = mesh_boundary.total_nb_cells() + cells.nb_active_cells();

  m_outfile << total_nb_elements << std::endl;

  Uint global_elem_id = 1;

  // 1) WRITE ALL BOUNDARY CELLS
  for (const std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>> bdomain :
       mesh_boundary.all_domains())
  {
    // const Uint dim = bdomain->dim();

    for (Uint c = 0; c < bdomain->nb_active_cells(); ++c)
    {
      /*
      const IncidencePair bcell_id = bdomain->bdry_cell_id(c);
      MeshEntity bcell = cells.active_cell(bcell_id.cell_idx);
      bcell.local_transform(bdomain->dim(), bcell_id.local_id);
      */

      const MeshEntity bcell = bdomain->active_cell(cells, ActiveIdx(c));

      const Uint gmsh_type = pdekit_to_gmsh.ref_topology_type_to_gmsh_type(bcell.pt_set_id());

      m_outfile << global_elem_id << " " << gmsh_type << " 2 " << bdomain->material_id() << " "
                << bdomain->material_id();

      for (Uint node = 0; node < bcell.nb_vert(); ++node)
      {
        m_outfile << " " << bcell.vertex(node) + 1;
      }
      m_outfile << std::endl;

      global_elem_id++;

    } // Loop over all cells of one boundary patch
  }
#if 0
  // 2) WRITE ALL INTERNAL CELLS
  typedef typename result_of::dof_map_t<MeshConfig> cells_dof_storage;

  for (const typename cells_dof_storage::const_dof_range_typed &cell_group :
       cells.all_active_dof_groups())
  {
    typename cells_dof_storage::const_dof_iterator_typed cell_iter = cell_group.begin();

    const Uint gmsh_type =
        pdekit_to_gmsh.ref_topology_type_to_gmsh_type(cell_iter->std_region_id());

    for (; cell_iter != cell_group.end(); ++cell_iter)
    {
      m_outfile << global_elem_id << " " << gmsh_type << " 2 "
                << cells.active_cell_tag(cell_iter->idx()) << " "
                << cells.active_cell_tag(cell_iter->idx());

      const MeshEntity &one_cell = *cell_iter;
      for (Uint node = 0; node < one_cell.nb_vert(); ++node)
      {
        m_outfile << " " << one_cell.vertex(node) + 1;
      }
      m_outfile << std::endl;

      global_elem_id++;

    } // Loop over all cells of one group

  } // Loop over all element groups of the same dimension
#endif

  // 2) WRITE ALL INTERNAL CELLS
  for (Uint ac = 0; ac < cells.nb_active_cells(); ++ac)
  {
    const MeshEntity one_cell = cells.active_cell(ActiveIdx(ac));

    const Uint gmsh_type = pdekit_to_gmsh.ref_topology_type_to_gmsh_type(one_cell.pt_set_id());

    m_outfile << global_elem_id << " " << gmsh_type << " 2 "
              << cells.active_cell_tag(ActiveIdx(one_cell.idx())) << " "
              << cells.active_cell_tag(ActiveIdx(one_cell.idx()));

    for (Uint node = 0; node < one_cell.nb_vert(); ++node)
    {
      m_outfile << " " << one_cell.vertex(node) + 1;
    }
    m_outfile << std::endl;

    global_elem_id++;

  } // Loop over all element groups of the same dimension

  m_outfile << "$EndElements" << std::endl;
}

// ----------------------------------------------------------------------------

} // Namespace gmsh

} // Namespace mesh

} // Namespace pdekit

#endif
