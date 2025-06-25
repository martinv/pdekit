#ifndef PDEKIT_Mesh_Cell_Adapt_Case_hpp
#define PDEKIT_Mesh_Cell_Adapt_Case_hpp

#include "common/Constants.hpp"

namespace pdekit
{

namespace mesh
{

enum class CellTransform : Uint
{
  NO_TRANS                  = 0,
  RESTRICT_TO_CODIM_1       = 1,
  RESTRICT_TO_CODIM_2       = 2,
  COARSEN                   = 3,
  UNIFORM_REFINE            = 4,
  ANISO_REFINE_ORTHO_FACE_0 = 5,  // Aniso refinement in direction orthogonal to face 0
  ANISO_REFINE_ORTHO_FACE_1 = 6,  // Aniso refinement in direction orthogonal to face 1
  ANISO_REFINE_ORTHO_FACE_2 = 7,  // Aniso refinement in direction orthogonal to face 2
  ANISO_REFINE_ORTHO_FACE_3 = 8,  // Aniso refinement in direction orthogonal to face 3
  ANISO_REFINE_ORTHO_FACE_4 = 9,  // Aniso refinement in direction orthogonal to face 4
  ANISO_REFINE_ORTHO_FACE_5 = 10, // Aniso refinement in direction orthogonal to face 5
};

struct CellTransformName
{
  static const Uint NbInstances = 11;
  static const std::string value[NbInstances];
};

struct CellTransformValue
{
  static const Uint NbInstances = 11;
  static const CellTransform value[NbInstances];
  static const Uint FirstRefine      = 4;
  static const Uint LastRefine       = 10;
  static const Uint FirstAnisoRefine = 5;
  static const Uint LastAnisoRefine  = 10;
};

std::ostream &operator<<(std::ostream &os, const CellTransform adapt_op_id);

class CellTransformTraits
{
  public:
  inline static bool is_refinement(const CellTransform transform)
  {
    const std::underlying_type<CellTransform>::type transform_idx =
        static_cast<std::underlying_type<CellTransform>::type>(transform);
    return (CellTransformValue::FirstRefine <= transform_idx) &&
           (transform_idx <= CellTransformValue::LastRefine);
  }

  inline static bool is_aniso_refinement(const CellTransform transform)
  {
    /*
    return (transform == CellTransform::ANISO_REFINE_ORTHO_FACE_0) ||
           (transform == CellTransform::ANISO_REFINE_ORTHO_FACE_1) ||
           (transform == CellTransform::ANISO_REFINE_ORTHO_FACE_2) ||
           (transform == CellTransform::ANISO_REFINE_ORTHO_FACE_3) ||
           (transform == CellTransform::ANISO_REFINE_ORTHO_FACE_4) ||
           (transform == CellTransform::ANISO_REFINE_ORTHO_FACE_5);
    */

    const std::underlying_type<CellTransform>::type transform_idx =
        static_cast<std::underlying_type<CellTransform>::type>(transform);
    return (CellTransformValue::FirstAnisoRefine <= transform_idx) &&
           (transform_idx <= CellTransformValue::LastAnisoRefine);
  }

  inline static CellTransform aniso_refinement_ortho_face(const Uint face_id)
  {
    return CellTransformValue::value[CellTransformValue::FirstAnisoRefine + face_id];
  }

  // Return the index of face to which this refinement is orthogonal
  inline static Uint aniso_refinement_face_id(const CellTransform transform)
  {
    const std::underlying_type<CellTransform>::type transform_idx =
        static_cast<std::underlying_type<CellTransform>::type>(transform);
    return transform_idx - CellTransformValue::FirstAnisoRefine;
  }
};

} // namespace mesh

} // namespace pdekit

#endif
