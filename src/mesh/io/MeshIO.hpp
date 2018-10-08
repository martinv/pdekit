#ifndef PDEKIT_MESH_IO_hpp
#define PDEKIT_MESH_IO_hpp

#include "common/StringUtils.hpp"

// Include file reader headers
#include "mesh/io/gmsh/GmshReader.hpp"

// Include file writer headers
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"

namespace pdekit
{

namespace mesh
{

class MeshIO
{
  public:
  /// Read a mesh from file. The type of file is detected based on filename
  /// extension
  template <typename MeshType>
  static void read_mesh_from_file(const std::string &filename, MeshType &mesh);

  /// Write a mesh to file. The type of file is detected based on filename
  /// extension
  template <typename MeshType>
  static void write_mesh_to_file(const MeshType &mesh, const std::string &filename);

  /// Append field to file
  /// @param name of the file to which the data is written
  template <typename MeshType, typename DataMatrix>
  static void append_nodal_function_to_file(const MeshType &mesh, const std::string &filename,
                                            const DataMatrix &data, const std::string &field_name,
                                            const double time = 0.0, const int timestep = 0);
};

// ----------------------------------------------------------------------------

template <typename MeshType>
void MeshIO::read_mesh_from_file(const std::string &filename, MeshType &mesh)
{
  std::vector<std::string> words;
  common::StringUtils::split_string(filename, '.', words);

  if (words.back() == "msh")
  {
    gmsh::GmshReader reader;
    reader.read_mesh_from_file(filename, mesh);
  }
  else
  {
    std::cerr << "Error: unknown file type for mesh reading (unrecognized "
                 "extension \""
              << words.back() << "\")" << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshType>
void MeshIO::write_mesh_to_file(const MeshType &mesh, const std::string &filename)
{
  std::vector<std::string> words;
  common::StringUtils::split_string(filename, '.', words);

  if (words.back() == "msh")
  {
    gmsh::GmshWriter writer;
    writer.write_mesh_to_file(mesh, filename);
  }
  else if (words.back() == "vtk")
  {
    vtk::VtkWriter writer;
    writer.write_mesh_to_file(mesh, filename);
  }
  else
  {
    std::cerr << "Error: unknown file type for mesh writing (unrecognized "
                 "extension \""
              << words.back() << "\")" << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshType, typename DataMatrix>
void MeshIO::append_nodal_function_to_file(const MeshType &mesh, const std::string &filename,
                                           const DataMatrix &data, const std::string &field_name,
                                           const double time, const int timestep)
{
  std::vector<std::string> words;
  common::StringUtils::split_string(filename, '.', words);

  if (words.back() == "msh")
  {
    gmsh::GmshWriter writer;
    writer.append_nodal_function_to_file(mesh, filename, data, field_name, time, timestep);
  }
  else if (words.back() == "vtk")
  {
    vtk::VtkWriter writer;
    writer.append_nodal_function_to_file(mesh, filename, data, field_name, time, timestep);
  }
  else
  {
    std::cerr << "Error: unknown file type to append nodal data (unrecognized "
                 "extension \""
              << words.back() << "\")" << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
