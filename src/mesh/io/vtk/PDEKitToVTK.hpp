#ifndef PDEKIT_Mesh_IO_VTK_PDEKitToVTK_hpp
#define PDEKIT_Mesh_IO_VTK_PDEKitToVTK_hpp

#include <map>
#include <string>

#include "common/PDEKit.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

namespace vtk
{

class PDEKitToVTK
{

  public:
  /// Number of element types that vtk currently supports.
  enum
  {
    NbElemTypes = 26
  };

  PDEKitToVTK();

  ~PDEKitToVTK();

  static Uint nb_elem_types();

  static ElemShape element_shape(const Uint vtk_elem_id);

  static Uint nb_nodes_in_elem(const Uint vtk_elem_id);

  static Uint elem_dim(const Uint vtk_elem_id);

  static Uint elem_order(const Uint vtk_elem_id);

  static const std::string elem_geo_name(const Uint vtk_elem_id);

  /// Convert an element type given as string to a number
  static Uint ref_topology_type_to_vtk_type(const PointSetTag rt_type_id);

  private:
  /// Store the shape of each element type
  static const ElemShape ElementShape[NbElemTypes];

  /// Store the number of nodes for each element type
  static const Uint NodesInVtkElem[NbElemTypes];

  /// 'Default' dimension of each element type
  static const Uint VtkElemDim[NbElemTypes];

  /// Polynomial order of Lagrange shape functions of each element type
  static const Uint VtkElemOrder[NbElemTypes];

  /// Strings holding the name of the geometrical shape of each element type
  static const std::string VtkElemGeoName[NbElemTypes];

  /// A map which assigns to a give string holding element name a number of
  /// corresponding element type number in vtk
  static const std::map<PointSetTag, Uint> PDEKitToVTKTypeMap;
};

} // Namespace vtk

} // Namespace mesh

} // Namespace pdekit

#endif
