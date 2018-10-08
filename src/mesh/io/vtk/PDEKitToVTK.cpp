#include "mesh/io/vtk/PDEKitToVTK.hpp"
#include "boost/assign/list_of.hpp"
#include "mesh/std_region/StdRegionEntity.hpp"

namespace pdekit
{

namespace mesh
{

namespace vtk
{

const ElemShape PDEKitToVTK::ElementShape[NbElemTypes] = {
    ElemShape::Undefined, ElemShape::Point,     ElemShape::Undefined,
    ElemShape::Line,                                                  //  0 - 3
    ElemShape::Undefined, ElemShape::Triag,     ElemShape::Undefined, //  4 - 6
    ElemShape::Undefined, ElemShape::Undefined, ElemShape::Quad,      //  7 - 9
    ElemShape::Tetra,     ElemShape::Undefined, ElemShape::Hexa,      // 10 - 12
    ElemShape::Undefined, ElemShape::Pyramid,   ElemShape::Undefined, // 13 - 15
    ElemShape::Undefined, ElemShape::Undefined, ElemShape::Undefined, // 16 - 18
    ElemShape::Undefined, ElemShape::Undefined, ElemShape::Line,      // 19 - 21
    ElemShape::Triag,     ElemShape::Quad,      ElemShape::Tetra,     // 22 - 24
    ElemShape::Hexa                                                   // 25
};

const Uint PDEKitToVTK::NodesInVtkElem[NbElemTypes] = {0, 1, 0,  2, 0, 3, 0, 0, 0, 4, 4, //  0 - 10
                                                       0, 8, 0,  5, 0, 0, 0, 0, 0, 0, 3, // 11 - 21
                                                       6, 8, 10, 20};                    // 22 - 25

const Uint PDEKitToVTK::VtkElemDim[NbElemTypes] = {
    _0D, _0D, _0D, _1D, _0D, _2D, _0D, _0D, _0D, _2D, _3D, //  0 - 10
    _0D, _3D, _0D, _3D, _0D, _0D, _0D, _0D, _0D, _0D, _1D, // 11 - 21
    _2D, _2D, _3D, _3D};                                   // 22 - 25

const Uint PDEKitToVTK::VtkElemOrder[NbElemTypes] = {
    P0, P0, P0, P1, P0, P1, P0, P0, P0, P1, P1, //  0 - 10
    P0, P1, P0, P1, P0, P0, P0, P0, P0, P0, P2, // 11 - 21
    P2, P2, P2, P2};                            // 22 - 25

const std::string PDEKitToVTK::VtkElemGeoName[NbElemTypes] = {
    "Undefined", "Point",     "Undefined", "Line", //  0 - 3
    "Undefined", "Triag",     "Undefined",         //  4 - 6
    "Undefined", "Undefined", "Quad",              //  7 - 9
    "Tetra",     "Undefined", "Hexa",              // 10 - 12
    "Undefined", "Pyramid",   "Undefined",         // 13 - 15
    "Undefined", "Undefined", "Undefined",         // 16 - 18
    "Undefined", "Undefined", "Line",              // 19 - 21
    "Triag",     "Quad",      "Tetra",             // 22 - 24
    "Hexa"};                                       // 25

const std::map<PointSetTag, Uint> PDEKitToVTK::PDEKitToVTKTypeMap =
    boost::assign::map_list_of(PointSetTag(ElemShape::Line, P1, PointSetID::Equidist),
                               3)(PointSetTag(ElemShape::Line, P2, PointSetID::Equidist),
                                  21)(PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist), 5)(
        PointSetTag(ElemShape::Triag, P2, PointSetID::Equidist),
        22)(PointSetTag(ElemShape::Quad, P1, PointSetID::Equidist),
            9)(PointSetTag(ElemShape::Quad, P2, PointSetID::Equidist),
               23)(PointSetTag(ElemShape::Tetra, P1, PointSetID::Equidist),
                   10)(PointSetTag(ElemShape::Tetra, P2, PointSetID::Equidist),
                       24)(PointSetTag(ElemShape::Hexa, P1, PointSetID::Equidist),
                           12)(PointSetTag(ElemShape::Hexa, P2, PointSetID::Equidist),
                               25)(PointSetTag(ElemShape::Prism, P1, PointSetID::Equidist),
                                   13)(PointSetTag(ElemShape::Line, P1, PointSetID::Warpblend), 3)(
        PointSetTag(ElemShape::Line, P2, PointSetID::Warpblend),
        21)(PointSetTag(ElemShape::Triag, P1, PointSetID::Warpblend),
            5)(PointSetTag(ElemShape::Triag, P2, PointSetID::Warpblend),
               22)(PointSetTag(ElemShape::Quad, P1, PointSetID::Warpblend),
                   9)(PointSetTag(ElemShape::Quad, P2, PointSetID::Warpblend),
                      23)(PointSetTag(ElemShape::Tetra, P1, PointSetID::Warpblend),
                          10)(PointSetTag(ElemShape::Tetra, P2, PointSetID::Warpblend),
                              24)(PointSetTag(ElemShape::Hexa, P1, PointSetID::Warpblend),
                                  12)(PointSetTag(ElemShape::Hexa, P2, PointSetID::Warpblend), 25)(
        PointSetTag(ElemShape::Prism, P1, PointSetID::Warpblend), 13);

// --------------------------------------------------------------------

PDEKitToVTK::PDEKitToVTK()
{
}

// --------------------------------------------------------------------

PDEKitToVTK::~PDEKitToVTK()
{
}

// --------------------------------------------------------------------

Uint PDEKitToVTK::nb_elem_types()
{
  return NbElemTypes;
}

// --------------------------------------------------------------------

ElemShape PDEKitToVTK::element_shape(const Uint vtk_elem_id)
{
  return ElementShape[vtk_elem_id];
}

// --------------------------------------------------------------------

Uint PDEKitToVTK::nb_nodes_in_elem(const Uint vtk_elem_id)
{
  return NodesInVtkElem[vtk_elem_id];
}

// --------------------------------------------------------------------

Uint PDEKitToVTK::elem_dim(const Uint vtk_elem_id)
{
  return VtkElemDim[vtk_elem_id];
}

// --------------------------------------------------------------------

Uint PDEKitToVTK::elem_order(const Uint vtk_elem_id)
{
  return VtkElemOrder[vtk_elem_id];
}

// --------------------------------------------------------------------

const std::string PDEKitToVTK::elem_geo_name(const Uint vtk_elem_id)
{
  return VtkElemGeoName[vtk_elem_id];
}

// --------------------------------------------------------------------

Uint PDEKitToVTK::ref_topology_type_to_vtk_type(const PointSetTag rt_type_id)
{
  std::map<PointSetTag, Uint>::const_iterator vtk_type_iterator;

  vtk_type_iterator = PDEKitToVTKTypeMap.find(rt_type_id);

  if (vtk_type_iterator != PDEKitToVTKTypeMap.end())
  {
    return vtk_type_iterator->second;
  }

  //  for(vtk_type_iterator = SCFDToGmsh.begin(); vtk_type_iterator !=
  // SCFDToGmsh.end(); ++vtk_type_iterator)
  //  {
  //    if ( vtk_type_iterator->first == name )
  //    {
  //      return vtk_type_iterator->second;
  //    }
  //  }

  // Default return value
  return 0;
}

// --------------------------------------------------------------------

} // Namespace vtk

} // Namespace mesh

} // Namespace pdekit
