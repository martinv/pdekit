#include "common/Constants.hpp"

namespace pdekit
{

// ----------------------------------------------------------------------------

const std::string ElemShapeInfo::name(const ElemShape shape)
{
  return Names[static_cast<std::underlying_type<ElemShape>::type>(shape)];
}

const ElemShape ElemShapeInfo::value(const Uint idx)
{
  return Value[idx];
}

const std::string ElemShapeInfo::Names[ElemShapeInfo::NbInstances] = {
    "Undefined", "Point",     "Line",  "Triag", "Quad",  "Undefined", "Undefined", "Undefined",
    "Undefined", "Undefined", "Tetra", "Hexa",  "Prism", "Pyramid",   "Undefined", "Undefined"};

const ElemShape ElemShapeInfo::Value[ElemShapeInfo::NbInstances] = {
    ElemShape::Undefined, ElemShape::Point,     ElemShape::Line,      ElemShape::Triag,
    ElemShape::Quad,      ElemShape::Undefined, ElemShape::Undefined, ElemShape::Undefined,
    ElemShape::Undefined, ElemShape::Undefined, ElemShape::Tetra,     ElemShape::Hexa,
    ElemShape::Prism,     ElemShape::Pyramid,   ElemShape::Undefined, ElemShape::Undefined};

// ----------------------------------------------------------------------------

const std::string PolyOrder::name(const Uint order_id)
{
  return Names[order_id];
}

const Uint PolyOrder::value(const Uint idx)
{
  return Value[idx];
}

const std::string PolyOrder::Names[PolyOrder::NbInstances] = {
    "P0",  "P1",  "P2",  "P3",  "P4",  "P5",  "P6",  "P7",  "P8",  "P9",  "P10",
    "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20", "P21",
    "P22", "P23", "P24", "P25", "P26", "P27", "P28", "P29", "P30"};
const Uint PolyOrder::Value[PolyOrder::NbInstances] = {
    P0,  P1,  P2,  P3,  P4,  P5,  P6,  P7,  P8,  P9,  P10, P11, P12, P13, P14, P15,
    P16, P17, P18, P19, P20, P21, P22, P23, P24, P25, P26, P27, P28, P29, P30};

// ----------------------------------------------------------------------------

const std::string SFuncInfo::name(const SFunc sf)
{
  return Names[static_cast<std::underlying_type<SFunc>::type>(sf)];
}

const SFunc SFuncInfo::value(const Uint idx)
{
  return Value[idx];
}

const std::string SFuncInfo::Names[SFuncInfo::NbInstances] = {
    "Undefined", "Lagrange", "BezierBernstein", "Modal",
    "Carnevali", "C0Modal",  "ThomasRaviart",   "BDFM"};

const SFunc SFuncInfo::Value[SFuncInfo::NbInstances] = {
    SFunc::Undefined, SFunc::Lagrange, SFunc::BezierBernstein, SFunc::Modal,
    SFunc::Carnevali, SFunc::C0Modal,  SFunc::ThomasRaviart,   SFunc::BDFM};

// ----------------------------------------------------------------------------

const std::string ModalBasisInfo::name(const ModalBasis mb)
{
  return Names[static_cast<std::underlying_type<ModalBasis>::type>(mb)];
}

const ModalBasis ModalBasisInfo::value(const Uint idx)
{
  return Value[idx];
}

const std::string ModalBasisInfo::Names[ModalBasisInfo::NbInstances] = {"Undefined", "Modal",
                                                                        "Carnevali", "C0Modal"};

const ModalBasis ModalBasisInfo::Value[ModalBasisInfo::NbInstances] = {
    ModalBasis::Undefined, ModalBasis::Modal, ModalBasis::Carnevali, ModalBasis::C0Modal};

// ----------------------------------------------------------------------------

const std::string Dim::name(const DimID dim_id)
{
  return Names[static_cast<std::underlying_type<DimID>::type>(dim_id)];
}

const Uint Dim::value(const Uint idx)
{
  return Value[idx];
}

const std::string Dim::Names[Dim::NbInstances] = {"0D", "1D", "2D", "3D"};

const Uint Dim::Value[Dim::NbInstances] = {_0D, _1D, _2D, _3D};

// ----------------------------------------------------------------------------

const std::string PointSetInfo::name(const PointSetID id)
{
  return Names[static_cast<std::underlying_type<PointSetID>::type>(id)];
}

const PointSetID PointSetInfo::value(const Uint idx)
{
  return Value[idx];
}

const std::string PointSetInfo::Names[PointSetInfo::NbInstances] = {
    "Undefined",   "Equidist",    "Warpblend", "Fekete",      "Gauss",
    "GaussLobato", "NewtonCotes", "FaceGauss", "FaceLobatto", "FaceNewtonCotes"};

const PointSetID PointSetInfo::Value[PointSetInfo::NbInstances] = {
    PointSetID::Undefined,      PointSetID::Equidist,  PointSetID::Warpblend,
    PointSetID::Fekete,         PointSetID::Gauss,     PointSetID::GaussLobatto,
    PointSetID::NewtonCotes,    PointSetID::FaceGauss, PointSetID::FaceLobatto,
    PointSetID::FaceNewtonCotes};

// ----------------------------------------------------------------------------

bool sf_is_nodal(const SFunc sf)
{
  if ((sf == SFunc::Lagrange) || (sf == SFunc::BezierBernstein))
  {
    return true;
  }
  return false;
}

// ----------------------------------------------------------------------------

bool sf_is_modal(const SFunc sf)
{
  if ((sf == SFunc::Modal) || (sf == SFunc::C0Modal) || (sf == SFunc::Carnevali))
  {
    return true;
  }
  return false;
}

// ----------------------------------------------------------------------------

ModalBasis cast_modal_sf_to_modal_basis(const SFunc sf)
{
  if (sf == SFunc::Modal)
  {
    return ModalBasis::Modal;
  }
  else if (sf == SFunc::Carnevali)
  {
    return ModalBasis::Carnevali;
  }
  else if (sf == SFunc::C0Modal)
  {
    return ModalBasis::C0Modal;
  }
  return ModalBasis::Undefined;
}

// ----------------------------------------------------------------------------

SFunc cast_modal_basis_to_modal_sf(const ModalBasis basis)
{
  if (basis == ModalBasis::Modal)
  {
    return SFunc::Modal;
  }
  else if (basis == ModalBasis::Carnevali)
  {
    return SFunc::Carnevali;
  }
  else if (basis == ModalBasis::C0Modal)
  {
    return SFunc::C0Modal;
  }
  return SFunc::Undefined;
}

// ----------------------------------------------------------------------------

} // Namespace pdekit
