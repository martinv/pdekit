#ifndef PDEKIT_Common_Constants_hpp
#define PDEKIT_Common_Constants_hpp

// This file defines generic 'constants' - aliases and enums
// for element shape, polynomial order, shape function type and dimension

#include <string>

#include "common/PDEKit.hpp"

namespace pdekit
{

// ----------------------------------------------------------------------------
//                            ORIENTATION
// ----------------------------------------------------------------------------

enum LeftRightOrientation : SUint
{
  LEFT  = 0,
  RIGHT = 1
};

// ----------------------------------------------------------------------------
//                            ELEMENT SHAPE
// ----------------------------------------------------------------------------

// Label for element shape
enum class ElemShape : Uint
{
  Undefined = 0,
  Point     = 1,
  Line      = 2,
  Triag     = 3,
  Quad      = 4,
  Tetra     = 10,
  Hexa      = 11,
  Prism     = 12,
  Pyramid   = 13
};

class ElemShapeInfo
{
  public:
  static constexpr Uint nb_instances()
  {
    return NbInstances;
  }

  static const std::string name(const ElemShape shape);

  static const ElemShape value(const Uint idx);

  private:
  enum
  {
    NbInstances = 16
  };

  static const std::string Names[NbInstances];
  static const ElemShape Value[NbInstances];
};

// ----------------------------------------------------------------------------
//                           POLYNOMIAL ORDER
// ----------------------------------------------------------------------------

// Label for polynomial order

enum PolyOrderID
{
  P0  = 0,
  P1  = 1,
  P2  = 2,
  P3  = 3,
  P4  = 4,
  P5  = 5,
  P6  = 6,
  P7  = 7,
  P8  = 8,
  P9  = 9,
  P10 = 10,
  P11 = 11,
  P12 = 12,
  P13 = 13,
  P14 = 14,
  P15 = 15,
  P16 = 16,
  P17 = 17,
  P18 = 18,
  P19 = 19,
  P20 = 20,
  P21 = 21,
  P22 = 22,
  P23 = 23,
  P24 = 24,
  P25 = 25,
  P26 = 26,
  P27 = 27,
  P28 = 28,
  P29 = 29,
  P30 = 30
};

class PolyOrder
{
  public:
  static constexpr Uint nb_instances()
  {
    return NbInstances;
  }

  static const std::string name(const Uint order_id);

  static const Uint value(const Uint idx);

  private:
  enum
  {
    NbInstances = 31
  };
  static const std::string Names[NbInstances];
  static const Uint Value[NbInstances];
};

// ----------------------------------------------------------------------------
//                         SHAPE FUNCTION TYPE
// ----------------------------------------------------------------------------

enum class SFunc : Uint
{
  Undefined       = 0,
  Lagrange        = 1,
  BezierBernstein = 2,
  Modal           = 3,
  Carnevali       = 4,
  C0Modal         = 5,
  ThomasRaviart   = 6,
  BDFM            = 7
};

// Label for shape functions
class SFuncInfo
{
  public:
  static constexpr Uint nb_instances()
  {
    return NbInstances;
  }

  static const std::string name(const SFunc sf);

  static const SFunc value(const Uint idx);

  private:
  enum
  {
    NbInstances = 8
  };

  static const std::string Names[NbInstances];
  static const SFunc Value[NbInstances];
};

// ----------------------------------------------------------------------------
//                           MODAL BASIS TYPE
// ----------------------------------------------------------------------------

enum class ModalBasis : Uint
{
  Undefined = 0,
  Modal     = 1,
  Carnevali = 2,
  C0Modal   = 3
};

// Label for shape functions
class ModalBasisInfo
{
  public:
  static constexpr Uint nb_instances()
  {
    return NbInstances;
  }

  static const std::string name(const ModalBasis mb);

  static const ModalBasis value(const Uint idx);

  private:
  enum
  {
    NbInstances = 4
  };
  static const std::string Names[NbInstances];
  static const ModalBasis Value[NbInstances];
};

// ----------------------------------------------------------------------------
//                             DIMENSION
// ----------------------------------------------------------------------------

// Label for dimension
enum DimID
{
  _0D = 0,
  _1D = 1,
  _2D = 2,
  _3D = 3
};

struct Dim
{
  public:
  static constexpr Uint nb_instances()
  {
    return NbInstances;
  }

  static const std::string name(const DimID dim_id);

  static const Uint value(const Uint idx);

  private:
  enum
  {
    NbInstances = 4
  };
  static const std::string Names[NbInstances];
  static const Uint Value[NbInstances];
};

// ----------------------------------------------------------------------------
//                   Interpolation point set
// ----------------------------------------------------------------------------

enum class PointSetID : Uint
{
  Undefined       = 0,
  Equidist        = 1,
  Warpblend       = 2,
  Fekete          = 3,
  Gauss           = 4,
  GaussLobatto    = 5,
  NewtonCotes     = 6, // Equidistantly spaced quadrature points
  FaceGauss       = 7,
  FaceLobatto     = 8,
  FaceNewtonCotes = 9
};

// Label for interpolation point set
class PointSetInfo
{
  public:
  static constexpr Uint nb_instances()
  {
    return NbInstances;
  }

  static const std::string name(const PointSetID id);

  static const PointSetID value(const Uint idx);

  private:
  enum
  {
    NbInstances = 10
  };
  static const std::string Names[NbInstances];
  static const PointSetID Value[NbInstances];
};

// ----------------------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------------------

bool sf_is_nodal(const SFunc sf);

bool sf_is_modal(const SFunc sf);

ModalBasis cast_modal_sf_to_modal_basis(const SFunc sf);

SFunc cast_modal_basis_to_modal_sf(const ModalBasis basis);

} // Namespace pdekit

#endif
