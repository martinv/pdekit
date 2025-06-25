#ifndef PDEKIT_Mesh_IO_PDEKitToGmsh_hpp
#define PDEKIT_Mesh_IO_PDEKitToGmsh_hpp

#include <map>
#include <string>
#include <vector>

#include "common/PDEKit.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

namespace gmsh
{

class PDEKitToGmsh
{

  public:
  /// Number of element types that gmsh currently supports. Gmsh supports many
  /// element types, we create an array of length 40 so that the types can be
  /// indexed from 1 to 39 as they are stored
  /// in the msh file format. Index 0 does not refer to any type.
  /// See http://www.geuz.org/gmsh/doc/texinfo/gmsh.html#MSH-ASCII-file-format
  enum
  {
    NbElemTypes = 150
  };

  PDEKitToGmsh();

  ~PDEKitToGmsh();

  Uint nb_elem_types() const;

  ElemShape element_shape(const Uint gmsh_elem_id) const;

  Uint nb_nodes_in_elem(const Uint gmsh_elem_id) const;

  Uint elem_dim(const Uint gmsh_elem_id) const;

  Uint elem_order(const Uint gmsh_elem_id) const;

  const std::string elem_geo_name(const Uint gmsh_elem_id) const;

  /// Convert an element type given as string to a number
  Uint ref_topology_type_to_gmsh_type(const PointSetTag rt_type_id) const;

  private:
  /// METHODS
  static bool static_data_initialized;
  static void initialize_static_data();

  /// DATA
  /// Store the shape of each element type
  static std::vector<ElemShape> ElementShape;

  /// Store the number of nodes for each element type
  static std::vector<Uint> NodesInGmshElem;

  /// 'Default' dimension of each element type
  static std::vector<Uint> GmshElemDim;

  /// Polynomial order of Lagrange shape functions of each element type
  static std::vector<Uint> GmshElemOrder;

  /// Strings holding the name of the geometrical shape of each element type
  static std::vector<std::string> GmshElemGeoName;

  /// A map which assigns to a give string holding element name a number of
  /// corresponding element type number in gmsh
  static std::map<PointSetTag, Uint> PDEKitToGmshId;
};

} // Namespace gmsh

} // Namespace mesh

} // Namespace pdekit

#endif
