#include "mesh/CellTransform.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

const std::string CellTransformName::value[CellTransformName::NbInstances] = {
    "NO_TRANS",
    "RESTRICT_TO_CODIM_1",
    "RESTRICT_TO_CODIM_2",
    "COARSEN",
    "UNIFORM_REFINE",
    "ANISO_REFINE_ORTHO_FACE_0",
    "ANISO_REFINE_ORTHO_FACE_1",
    "ANISO_REFINE_ORTHO_FACE_2",
    "ANISO_REFINE_ORTHO_FACE_3",
    "ANISO_REFINE_ORTHO_FACE_4",
    "ANISO_REFINE_ORTHO_FACE_5"};

// ----------------------------------------------------------------------------

const CellTransform CellTransformValue::value[CellTransformValue::NbInstances] = {
    CellTransform::NO_TRANS,
    CellTransform::RESTRICT_TO_CODIM_1,
    CellTransform::RESTRICT_TO_CODIM_2,
    CellTransform::COARSEN,
    CellTransform::UNIFORM_REFINE,
    CellTransform::ANISO_REFINE_ORTHO_FACE_0,
    CellTransform::ANISO_REFINE_ORTHO_FACE_1,
    CellTransform::ANISO_REFINE_ORTHO_FACE_2,
    CellTransform::ANISO_REFINE_ORTHO_FACE_3,
    CellTransform::ANISO_REFINE_ORTHO_FACE_4,
    CellTransform::ANISO_REFINE_ORTHO_FACE_5};

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const CellTransform adapt_op_id)
{
  os << CellTransformName::value[static_cast<std::underlying_type<CellTransform>::type>(
      adapt_op_id)];
  return os;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
