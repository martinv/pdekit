#include "mesh/io/gmsh/PDEKitToGmsh.hpp"
#include "boost/assign/list_of.hpp"
#include "mesh/std_region/StdRegionEntity.hpp"

namespace pdekit
{

namespace mesh
{

namespace gmsh
{

// ----------------------------------------------------------------------------

// Static variables:

bool PDEKitToGmsh::static_data_initialized = false;

std::vector<ElemShape> PDEKitToGmsh::ElementShape = std::vector<ElemShape>();

std::vector<Uint> PDEKitToGmsh::NodesInGmshElem = std::vector<Uint>();

std::vector<Uint> PDEKitToGmsh::GmshElemDim = std::vector<Uint>();

std::vector<Uint> PDEKitToGmsh::GmshElemOrder = std::vector<Uint>();

std::vector<std::string> PDEKitToGmsh::GmshElemGeoName = std::vector<std::string>();

std::map<PointSetTag, Uint> PDEKitToGmsh::PDEKitToGmshId = std::map<PointSetTag, Uint>();

// ----------------------------------------------------------------------------

PDEKitToGmsh::PDEKitToGmsh()
{
  if (!static_data_initialized)
  {
    initialize_static_data();
    static_data_initialized = true;
  }
}

// ----------------------------------------------------------------------------

PDEKitToGmsh::~PDEKitToGmsh()
{
}

// ----------------------------------------------------------------------------

Uint PDEKitToGmsh::nb_elem_types() const
{
  return NbElemTypes;
}

// ----------------------------------------------------------------------------

ElemShape PDEKitToGmsh::element_shape(const Uint gmsh_elem_id) const
{
  return ElementShape[gmsh_elem_id];
}

// ----------------------------------------------------------------------------

Uint PDEKitToGmsh::nb_nodes_in_elem(const Uint gmsh_elem_id) const
{
  return NodesInGmshElem[gmsh_elem_id];
}

// ----------------------------------------------------------------------------

Uint PDEKitToGmsh::elem_dim(const Uint gmsh_elem_id) const
{
  return GmshElemDim[gmsh_elem_id];
}

// ----------------------------------------------------------------------------

Uint PDEKitToGmsh::elem_order(const Uint gmsh_elem_id) const
{
  return GmshElemOrder[gmsh_elem_id];
}

// ----------------------------------------------------------------------------

const std::string PDEKitToGmsh::elem_geo_name(const Uint gmsh_elem_id) const
{
  return GmshElemGeoName[gmsh_elem_id];
}

// ----------------------------------------------------------------------------

Uint PDEKitToGmsh::ref_topology_type_to_gmsh_type(const PointSetTag rt_type_id) const
{
  std::map<PointSetTag, Uint>::const_iterator gmsh_type_iterator;

  gmsh_type_iterator = PDEKitToGmshId.find(rt_type_id);

  if (gmsh_type_iterator != PDEKitToGmshId.end())
  {
    return gmsh_type_iterator->second;
  }

  //  for(gmsh_type_iterator = SCFDToGmsh.begin(); gmsh_type_iterator !=
  // SCFDToGmsh.end(); ++gmsh_type_iterator)
  //  {
  //    if ( gmsh_type_iterator->first == name )
  //    {
  //      return gmsh_type_iterator->second;
  //    }
  //  }

  // Default return value
  return 0;
}

// ----------------------------------------------------------------------------
// For element type numbers defined in GMSH, look in
// Common/GmshDefines.h in Gmsh source code
// ----------------------------------------------------------------------------

void PDEKitToGmsh::initialize_static_data()
{
  ElementShape.resize(NbElemTypes, ElemShape::Undefined);
  NodesInGmshElem.resize(NbElemTypes, 0u);
  GmshElemDim.resize(NbElemTypes, _0D);
  GmshElemOrder.resize(NbElemTypes, P0);
  GmshElemGeoName.resize(NbElemTypes, "");

  // --------------------------------------------------------------------------

  // Line elements - orders P1 - P6
  // The number of nodes in line of order P is
  // N = (P+1)
  // GMSH element type 1 is P1 line
  // GMSH element type 8 is P2 line
  // GMSH element type 26 is P3 line etc.
  std::vector<Uint> line_elem_idx = {1, 8, 26, 27, 28, 62, 63, 64, 65, 66};
  for (Uint i = 0; i < line_elem_idx.size(); ++i)
  {
    ElementShape[line_elem_idx[i]] = ElemShape::Line;
    // P1 element corresponds to i = 0, i.e. (i+2) == (P+1)
    NodesInGmshElem[line_elem_idx[i]] = 2 + i;
    GmshElemDim[line_elem_idx[i]]     = _1D;
    GmshElemOrder[line_elem_idx[i]]   = 1 + i;
    GmshElemGeoName[line_elem_idx[i]] = "Line";

    PDEKitToGmshId[PointSetTag(ElemShape::Line, 1 + i, PointSetID::Equidist)]  = line_elem_idx[i];
    PDEKitToGmshId[PointSetTag(ElemShape::Line, 1 + i, PointSetID::Warpblend)] = line_elem_idx[i];
  }

  // --------------------------------------------------------------------------

  // Triangle elements - orders P1 - P6
  // The number of nodes in triag of order P is
  // N = (P+1)(P+2)/2
  // GMSH element type 2 is P1 triag
  // GMSH element type 9 is P2 triag
  // GMSH element type 21 is P3 triag etc.
  std::vector<Uint> triag_elem_idx = {2, 9, 21, 23, 25, 42, 43, 44, 45, 46};
  for (Uint i = 0; i < triag_elem_idx.size(); ++i)
  {
    ElementShape[triag_elem_idx[i]] = ElemShape::Triag;
    // P1 element corresponds to i = 0, i.e. (i+2) == (P+1)
    NodesInGmshElem[triag_elem_idx[i]] = (i + 2) * (i + 3) / 2;
    GmshElemDim[triag_elem_idx[i]]     = _2D;
    GmshElemOrder[triag_elem_idx[i]]   = 1 + i;
    GmshElemGeoName[triag_elem_idx[i]] = "Triag";

    PDEKitToGmshId[PointSetTag(ElemShape::Triag, 1 + i, PointSetID::Equidist)]  = triag_elem_idx[i];
    PDEKitToGmshId[PointSetTag(ElemShape::Triag, 1 + i, PointSetID::Warpblend)] = triag_elem_idx[i];
  }

  // --------------------------------------------------------------------------

  // Quadrilateral elements - orders P1 - P4
  // The number of nodes in quadrilateral of order P is
  // N = (P+1)(P+1)
  // GMSH element type 3 is P1 quadrilateral
  // GMSH element type 10 is P2 quadrilateral
  // GMSH element type 36 is P3 quadrilateral
  // GMSH element type 37 is P4 quadrilateral
  // GMSH element type 38 is P5 quadrilateral
  // GMSH element type 47 is P6 quadrilateral
  // GMSH element type 48 is P7 quadrilateral
  // GMSH element type 49 is P8 quadrilateral
  // GMSH element type 50 is P9 quadrilateral
  // GMSH element type 51 is P10 quadrilateral
  std::vector<Uint> quad_elem_idx = {3, 10, 36, 37, 38, 47, 48, 49, 50, 51};
  for (Uint i = 0; i < quad_elem_idx.size(); ++i)
  {
    ElementShape[quad_elem_idx[i]] = ElemShape::Quad;
    // P1 element corresponds to i = 0, i.e. (i+2) == (P+1)
    NodesInGmshElem[quad_elem_idx[i]] = (i + 2) * (i + 2);
    GmshElemDim[quad_elem_idx[i]]     = _2D;
    GmshElemOrder[quad_elem_idx[i]]   = 1 + i;
    GmshElemGeoName[quad_elem_idx[i]] = "Quad";

    PDEKitToGmshId[PointSetTag(ElemShape::Quad, 1 + i, PointSetID::Equidist)]  = quad_elem_idx[i];
    PDEKitToGmshId[PointSetTag(ElemShape::Quad, 1 + i, PointSetID::Warpblend)] = quad_elem_idx[i];
  }

  // --------------------------------------------------------------------------

  // Tetrahedral elements - orders P1 - P10
  // The number of nodes in tetra of order P is
  // N = (P+1)(P+2)(P+3)/6
  // GMSH element type 3 is P1 quadrilateral
  // GMSH element type 10 is P2 quadrilateral
  // GMSH element type 36 is P3 quadrilateral etc.
  std::vector<Uint> tetra_elem_idx = {4, 11, 29, 30, 31, 71, 72, 73, 74, 75};
  for (Uint i = 0; i < tetra_elem_idx.size(); ++i)
  {
    ElementShape[tetra_elem_idx[i]] = ElemShape::Tetra;
    // P1 element corresponds to i = 0, i.e. (i+2) == (P+1)
    NodesInGmshElem[tetra_elem_idx[i]] = (i + 2) * (i + 3) * (i + 4) / 6;
    GmshElemDim[tetra_elem_idx[i]]     = _3D;
    GmshElemOrder[tetra_elem_idx[i]]   = 1 + i;
    GmshElemGeoName[tetra_elem_idx[i]] = "Tetra";

    PDEKitToGmshId[PointSetTag(ElemShape::Tetra, 1 + i, PointSetID::Equidist)]  = tetra_elem_idx[i];
    PDEKitToGmshId[PointSetTag(ElemShape::Tetra, 1 + i, PointSetID::Warpblend)] = tetra_elem_idx[i];
  }

  // --------------------------------------------------------------------------

  // Hexahedral elements - orders P1 - P2
  // The number of nodes in hexahedron of order P is
  // N = (P+1)(P+1)(P+1)
  // GMSH element type  5 is P1 hexahedron
  // GMSH element type 12 is P2 hexahedron
  // GMSH element type 98 is P9 hexahedron
  std::vector<Uint> hexa_elem_idx = {5, 12, 92, 93, 94, 95, 96, 97, 98};
  for (Uint i = 0; i < hexa_elem_idx.size(); ++i)
  {
    ElementShape[hexa_elem_idx[i]] = ElemShape::Hexa;
    // P1 element corresponds to i = 0, i.e. (i+2) == (P+1)
    NodesInGmshElem[hexa_elem_idx[i]] = (i + 2) * (i + 2) * (i + 2);
    GmshElemDim[hexa_elem_idx[i]]     = _3D;
    GmshElemOrder[hexa_elem_idx[i]]   = 1 + i;
    GmshElemGeoName[hexa_elem_idx[i]] = "Hexa";

    PDEKitToGmshId[PointSetTag(ElemShape::Hexa, 1 + i, PointSetID::Equidist)]  = hexa_elem_idx[i];
    PDEKitToGmshId[PointSetTag(ElemShape::Hexa, 1 + i, PointSetID::Warpblend)] = hexa_elem_idx[i];
  }

  // --------------------------------------------------------------------------

  // Pyramidal elements - orders P1 - P6
  // The number of nodes in hexahedron of order P is
  // N = (P+1)(P+2)(2*P+3)/6 { = 1 + 2^2 + 3^2 + 4^2 ... (p+1)^2 }
  // GMSH element type 7 is P1 pyramid
  // GMSH element type 14 is P2 pyramid
  std::vector<Uint> pyramid_elem_idx = {7, 14, 118, 119, 120, 121};
  for (Uint i = 0; i < pyramid_elem_idx.size(); ++i)
  {
    ElementShape[pyramid_elem_idx[i]] = ElemShape::Pyramid;
    // P1 element corresponds to i = 0, i.e. (i+2) == (P+1)
    NodesInGmshElem[pyramid_elem_idx[i]] = (i + 2) * (i + 3) * (2 * i + 5) / 6;
    GmshElemDim[pyramid_elem_idx[i]]     = _3D;
    GmshElemOrder[pyramid_elem_idx[i]]   = 1 + i;
    GmshElemGeoName[pyramid_elem_idx[i]] = "Pyramid";

    PDEKitToGmshId[PointSetTag(ElemShape::Pyramid, 1 + i, PointSetID::Equidist)] =
        pyramid_elem_idx[i];
  }

  // --------------------------------------------------------------------------

  // Prismatic elements - orders P1 - P2
  // The number of nodes in prism of order P is
  // N = (P+1)(P+2)(P+1)/2
  // GMSH element type 6 is P1 prism
  // GMSH element type 13 is P2 prism
  std::vector<Uint> prism_elem_idx = {6, 13};
  for (Uint i = 0; i < prism_elem_idx.size(); ++i)
  {
    ElementShape[prism_elem_idx[i]] = ElemShape::Prism;
    // P1 element corresponds to i = 0, i.e. (i+2) == (P+1)
    NodesInGmshElem[prism_elem_idx[i]] = (i + 2) * (i + 3) * (i + 2) / 2;
    GmshElemDim[prism_elem_idx[i]]     = _3D;
    GmshElemOrder[prism_elem_idx[i]]   = 1 + i;
    GmshElemGeoName[prism_elem_idx[i]] = "Prism";

    PDEKitToGmshId[PointSetTag(ElemShape::Prism, 1 + i, PointSetID::Equidist)]  = prism_elem_idx[i];
    PDEKitToGmshId[PointSetTag(ElemShape::Prism, 1 + i, PointSetID::Warpblend)] = prism_elem_idx[i];
  }
}

// ----------------------------------------------------------------------------

} // Namespace gmsh

} // Namespace mesh

} // Namespace pdekit
