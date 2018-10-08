#ifndef PDEKIT_Mesh_Element_Topology_hpp
#define PDEKIT_Mesh_Element_Topology_hpp

#include "common/Constants.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
//                         GENERIC INFORMATION ABOUT ELEMENT TYPES
// ----------------------------------------------------------------------------

template <ElemShape ElementShape>
struct ElementTopology;

// ----------------------------------------------------------------------------

template <>
struct ElementTopology<ElemShape::Line>
{
  static const Uint TopoDim   = _1D;
  static const Uint NbVerts   = 2;
  static const Uint NbEdges   = 1;
  static const Uint NbFaces   = 0;
  static const Uint NbVolumes = 0;
};

// ----------------------------------------------------------------------------

template <>
struct ElementTopology<ElemShape::Triag>
{
  static const Uint TopoDim   = _2D;
  static const Uint NbVerts   = 3;
  static const Uint NbEdges   = 3;
  static const Uint NbFaces   = 1;
  static const Uint NbVolumes = 0;
};

// ----------------------------------------------------------------------------

template <>
struct ElementTopology<ElemShape::Quad>
{
  static const Uint TopoDim   = _2D;
  static const Uint NbVerts   = 4;
  static const Uint NbEdges   = 4;
  static const Uint NbFaces   = 1;
  static const Uint NbVolumes = 0;
};

// ----------------------------------------------------------------------------

template <>
struct ElementTopology<ElemShape::Tetra>
{
  static const Uint TopoDim   = _3D;
  static const Uint NbVerts   = 4;
  static const Uint NbEdges   = 6;
  static const Uint NbFaces   = 4;
  static const Uint NbVolumes = 1;
};

// ----------------------------------------------------------------------------

template <>
struct ElementTopology<ElemShape::Hexa>
{
  static const Uint TopoDim   = _3D;
  static const Uint NbVerts   = 8;
  static const Uint NbEdges   = 12;
  static const Uint NbFaces   = 6;
  static const Uint NbVolumes = 1;
};

// ----------------------------------------------------------------------------

template <>
struct ElementTopology<ElemShape::Pyramid>
{
  static const Uint TopoDim   = _3D;
  static const Uint NbVerts   = 6;
  static const Uint NbEdges   = 9;
  static const Uint NbFaces   = 5;
  static const Uint NbVolumes = 1;
};

// ----------------------------------------------------------------------------

} // Namespace mesh

} // Namespace pdekit

#endif
