#ifndef PDEKIT_Mesh_Std_Region_Writer_hpp
#define PDEKIT_Mesh_Std_Region_Writer_hpp

#include "mesh/std_region/StdRegionBuilder.hpp"

namespace pdekit
{

namespace mesh
{

class StdRegionWriter
{
  public:
  /// Default constructor
  StdRegionWriter();

  /// Destructor
  ~StdRegionWriter();

  /// Write the interpolation point set to a vtk file
  static void write_to_vtk(PointSetTag const elem_type_tag, const std::string &filename);

  /// Write the interpolation point set to a point 3D file
  static void write_to_point3d(PointSetTag const elem_type_tag, const std::string &filename);

  /// Write the interpolation point set to a gmsh file
  static void write_to_gmsh(PointSetTag const elem_type_tag, const std::string &filename);

  private:
};

} // namespace mesh

} // namespace pdekit

#endif
