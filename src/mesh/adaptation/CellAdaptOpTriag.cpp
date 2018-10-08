#include "mesh/adaptation/CellAdaptOpTriag.hpp"
#include "mesh/CellTransform.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

namespace internal
{

// ----------------------------------------------------------------------------
//                   CellAdaptOpTriagUniformRefine
// ----------------------------------------------------------------------------

namespace
{

void child_coords_tri_uniform_refine_child0(math::DenseDMat<Real> const &ref_coords,
                                            math::DenseConstVecView<Real> const &coord_in,
                                            math::DenseVecView<Real> &coord_out)
{
  /*
  coord_out[X0] = 0.5 * coord_in[X0] - 0.5;
  coord_out[X1] = 0.5 * coord_in[X1] - 0.5;
  */

  coord_out[X0] = ref_coords(0, X0) + 0.5 * (coord_in[X0] - ref_coords(0, X0));
  coord_out[X1] = ref_coords(0, X1) + 0.5 * (coord_in[X1] - ref_coords(0, X1));
}

void child_coords_tri_uniform_refine_child1(math::DenseDMat<Real> const &ref_coords,
                                            math::DenseConstVecView<Real> const &coord_in,
                                            math::DenseVecView<Real> &coord_out)
{
  /*
  coord_out[X0] = 0.5 * coord_in[X0] + 0.5;
  coord_out[X1] = 0.5 * coord_in[X1] - 0.5;
  */

  coord_out[X0] =
      0.5 * (ref_coords(0, X0) + ref_coords(1, X0)) + 0.5 * (coord_in[X0] - ref_coords(0, X0));
  coord_out[X1] = ref_coords(0, X1) + 0.5 * (coord_in[X1] - ref_coords(0, X1));
}

void child_coords_tri_uniform_refine_child2(math::DenseDMat<Real> const &ref_coords,
                                            math::DenseConstVecView<Real> const &coord_in,
                                            math::DenseVecView<Real> &coord_out)
{
  /*
  coord_out[X0] = 0.5 * coord_in[X0] - 0.5;
  coord_out[X1] = 0.5 * coord_in[X1] + 0.5;
  */

  coord_out[X0] = ref_coords(0, X0) + 0.5 * (coord_in[X0] - ref_coords(0, X0));
  coord_out[X1] =
      0.5 * (ref_coords(0, X1) + ref_coords(2, X1)) + 0.5 * (coord_in[X1] - ref_coords(0, X1));
}

void child_coords_tri_uniform_refine_child3(math::DenseDMat<Real> const &ref_coords,
                                            math::DenseConstVecView<Real> const &coord_in,
                                            math::DenseVecView<Real> &coord_out)
{
  /*
  coord_out[X0] = -0.5 * coord_in[X0] - 0.5;
  coord_out[X1] = -0.5 * coord_in[X1] - 0.5;
  */

  coord_out[X0] =
      -0.5 * (coord_in[X0] - ref_coords(0, X0)) + 0.5 * (ref_coords(1, X0) + ref_coords(2, X0));
  coord_out[X1] =
      -0.5 * (coord_in[X1] - ref_coords(0, X1)) + 0.5 * (ref_coords(1, X1) + ref_coords(2, X1));
}

} // End of anonymous namespace

CellAdaptOpTriagUniformRefine::CellAdaptOpTriagUniformRefine() : CellAdaptOpBase()
{
}

CellAdaptOpTriagUniformRefine::~CellAdaptOpTriagUniformRefine()
{
}

void CellAdaptOpTriagUniformRefine::set_facet_adapt_op_ids(
    CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Triag, CellTransform::UNIFORM_REFINE);

  facet_adapt_op_ids.resize(3);
  facet_adapt_op_ids[0] = CellTransform::UNIFORM_REFINE;
  facet_adapt_op_ids[1] = CellTransform::UNIFORM_REFINE;
  facet_adapt_op_ids[2] = CellTransform::UNIFORM_REFINE;
}

Uint CellAdaptOpTriagUniformRefine::nb_parent_facets() const
{
  return 3u;
}

Uint CellAdaptOpTriagUniformRefine::nb_child_elems() const
{
  return 4u;
}

/// Fill the matrix of reference coordinates (of principal vertices)
void CellAdaptOpTriagUniformRefine::fill_parent_ref_coords(math::DenseDMat<Real> &ref_coords) const
{
  ref_coords.resize(3, _2D);

  ref_coords(0, X0) = -1.0;
  ref_coords(0, X1) = -1.0;
  ref_coords(1, X0) = 1.0;
  ref_coords(1, X1) = -1.0;
  ref_coords(2, X0) = -1.0;
  ref_coords(2, X1) = 1.0;
}

void CellAdaptOpTriagUniformRefine::set_child_elem_shapes(
    std::vector<ElemShape> &child_elem_shapes) const
{
  child_elem_shapes.resize(4);
  child_elem_shapes[0] = ElemShape::Triag;
  child_elem_shapes[1] = ElemShape::Triag;
  child_elem_shapes[2] = ElemShape::Triag;
  child_elem_shapes[3] = ElemShape::Triag;
}

void CellAdaptOpTriagUniformRefine::fill_local_child_incidences(
    std::vector<IncidenceEntry> &incidences, std::vector<EntityDofRealign> &permutations,
    std::vector<Uint> &local_offsets) const
{
  incidences.resize(6);

  incidences[0].cell_idx = 0;
  incidences[0].local_id = 1;
  incidences[1].cell_idx = 3;
  incidences[1].local_id = 1;

  incidences[2].cell_idx = 1;
  incidences[2].local_id = 2;
  incidences[3].cell_idx = 3;
  incidences[3].local_id = 2;

  incidences[4].cell_idx = 2;
  incidences[4].local_id = 0;
  incidences[5].cell_idx = 3;
  incidences[5].local_id = 0;

  permutations.resize(6);

  PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

  EntityRealignCode identity_p =
      EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0);
  EntityRealignCode flip_p =
      EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0);

  permutations[0].change_type(elem_tag, identity_p);
  permutations[1].change_type(elem_tag, flip_p);

  permutations[2].change_type(elem_tag, identity_p);
  permutations[3].change_type(elem_tag, flip_p);

  permutations[4].change_type(elem_tag, identity_p);
  permutations[5].change_type(elem_tag, flip_p);

  local_offsets.resize(4);
  local_offsets[0] = 0;
  local_offsets[1] = 2;
  local_offsets[2] = 4;
  local_offsets[3] = 6;
}

// ----------------------------------------------------------------------------

void CellAdaptOpTriagUniformRefine::fill_parent_child_incidences_on_facet(
    const Uint facet, std::vector<IncidenceEntry> &incidences,
    std::vector<EntityDofRealign> &permutations) const
{
  incidences.resize(2);

  switch (facet)
  {
    case 0u:
      incidences[0].cell_idx = 0;
      incidences[0].local_id = 0;
      incidences[1].cell_idx = 1;
      incidences[1].local_id = 0;
      break;

    case 1u:
      incidences[0].cell_idx = 1;
      incidences[0].local_id = 1;
      incidences[1].cell_idx = 2;
      incidences[1].local_id = 1;
      break;

    case 2u:
      incidences[0].cell_idx = 2;
      incidences[0].local_id = 2;
      incidences[1].cell_idx = 0;
      incidences[1].local_id = 2;
      break;
  };

  permutations.resize(2);

  PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

  const EntityRealignCode pcode_1 =
      EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0, ElemShape::Line, 0, 0);
  const EntityRealignCode pcode_2 =
      EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1, ElemShape::Line, 0, 0);

  permutations[0].change_type(elem_tag, pcode_1);
  permutations[1].change_type(elem_tag, pcode_2);
}

// ----------------------------------------------------------------------------

void CellAdaptOpTriagUniformRefine::fill_child_coord_transformers(
    std::vector<std::function<void(math::DenseDMat<Real> const &ref_coords,
                                   math::DenseConstVecView<Real> const &coord_in,
                                   math::DenseVecView<Real> &coord_out)>> &transformers) const
{
  transformers.resize(4);
  transformers[0] = child_coords_tri_uniform_refine_child0;
  transformers[1] = child_coords_tri_uniform_refine_child1;
  transformers[2] = child_coords_tri_uniform_refine_child2;
  transformers[3] = child_coords_tri_uniform_refine_child3;
}

// ----------------------------------------------------------------------------
//                   CellAdaptOpTriagAnisoRefineOrthoFace0
// ----------------------------------------------------------------------------

namespace
{

void child_coords_tri_aniso_refine_face0_child0(math::DenseDMat<Real> const &ref_coords,
                                                math::DenseConstVecView<Real> const &coord_in,
                                                math::DenseVecView<Real> &coord_out)
{
  /*
  coord_out[X0] = 0.5 * coord_in[X0] - 0.5;
  coord_out[X1] = coord_in[X1];
  */

  coord_out[X0] = ref_coords(0, X0) + 0.5 * (coord_in[X0] - ref_coords(0, X0));
  coord_out[X1] = coord_in[X1];
}

void child_coords_tri_aniso_refine_face0_child1(math::DenseDMat<Real> const &ref_coords,
                                                math::DenseConstVecView<Real> const &coord_in,
                                                math::DenseVecView<Real> &coord_out)
{
  /*
  coord_out[X0] = 0.5 * coord_in[X0] + 0.5 - 0.5 * (1.0 + coord_in[X1]);
  coord_out[X1] = coord_in[X1];
  */

  coord_out[X0] = 0.5 * (ref_coords(0, X0) + ref_coords(1, X0)) +
                  0.5 * (coord_in[X0] - ref_coords(0, X0)) -
                  0.5 * (ref_coords(1, X0) - ref_coords(0, X0)) /
                      (ref_coords(2, X1) - ref_coords(0, X1)) * (coord_in[X1] - ref_coords(0, X1));
  coord_out[X1] = coord_in[X1];
}

} // End of anonymous namespace

CellAdaptOpTriagAnisoRefineOrthoFace0::CellAdaptOpTriagAnisoRefineOrthoFace0() : CellAdaptOpBase()
{
}

CellAdaptOpTriagAnisoRefineOrthoFace0::~CellAdaptOpTriagAnisoRefineOrthoFace0()
{
}

void CellAdaptOpTriagAnisoRefineOrthoFace0::set_facet_adapt_op_ids(
    CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Triag, CellTransform::ANISO_REFINE_ORTHO_FACE_0);

  facet_adapt_op_ids.resize(3);
  facet_adapt_op_ids[0] = CellTransform::UNIFORM_REFINE;
  facet_adapt_op_ids[1] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[2] = CellTransform::NO_TRANS;
}

Uint CellAdaptOpTriagAnisoRefineOrthoFace0::nb_parent_facets() const
{
  return 3u;
}

Uint CellAdaptOpTriagAnisoRefineOrthoFace0::nb_child_elems() const
{
  return 2u;
}

void CellAdaptOpTriagAnisoRefineOrthoFace0::fill_parent_ref_coords(
    math::DenseDMat<Real> &ref_coords) const
{
  ref_coords.resize(3, _2D);

  ref_coords(0, X0) = -1.0;
  ref_coords(0, X1) = -1.0;
  ref_coords(1, X0) = 1.0;
  ref_coords(1, X1) = -1.0;
  ref_coords(2, X0) = -1.0;
  ref_coords(2, X1) = 1.0;
}

void CellAdaptOpTriagAnisoRefineOrthoFace0::set_child_elem_shapes(
    std::vector<ElemShape> &child_elem_shapes) const
{
  child_elem_shapes.resize(2);
  child_elem_shapes[0] = ElemShape::Triag;
  child_elem_shapes[1] = ElemShape::Triag;
}

void CellAdaptOpTriagAnisoRefineOrthoFace0::fill_local_child_incidences(
    std::vector<IncidenceEntry> &incidences, std::vector<EntityDofRealign> &permutations,
    std::vector<Uint> &local_offsets) const
{
  incidences.resize(2);

  incidences[0].cell_idx = 0;
  incidences[0].local_id = 1;
  incidences[1].cell_idx = 1;
  incidences[1].local_id = 2;

  permutations.resize(2);

  PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

  EntityRealignCode identity_p =
      EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0);
  EntityRealignCode flip_p =
      EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0);

  permutations[0].change_type(elem_tag, identity_p);
  permutations[1].change_type(elem_tag, flip_p);

  local_offsets.resize(2);
  local_offsets[0] = 0;
  local_offsets[1] = 2;
}

// ----------------------------------------------------------------------------

void CellAdaptOpTriagAnisoRefineOrthoFace0::fill_parent_child_incidences_on_facet(
    const Uint facet, std::vector<IncidenceEntry> &incidences,
    std::vector<EntityDofRealign> &permutations) const
{

  switch (facet)
  {
    case 0u:
    {
      incidences.resize(2);
      incidences[0].cell_idx = 0;
      incidences[0].local_id = 0;
      incidences[1].cell_idx = 1;
      incidences[1].local_id = 0;

      permutations.resize(2);

      const PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

      const EntityRealignCode pcode_1 = EntityRealignCode(
          ElemShape::Line, CellTransform::UNIFORM_REFINE, 0, ElemShape::Line, 0, 0);
      const EntityRealignCode pcode_2 = EntityRealignCode(
          ElemShape::Line, CellTransform::UNIFORM_REFINE, 1, ElemShape::Line, 0, 0);

      permutations[0].change_type(elem_tag, pcode_1);
      permutations[1].change_type(elem_tag, pcode_2);
    }

    break;

    case 1u:
    {
      incidences.resize(1);
      incidences[0].cell_idx = 1;
      incidences[0].local_id = 1;

      permutations.resize(1);

      const PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

      const EntityRealignCode pcode =
          EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0);

      permutations[0].change_type(elem_tag, pcode);
    }

    break;

    case 2u:
    {
      incidences.resize(1);
      incidences[0].cell_idx = 0;
      incidences[0].local_id = 2;

      permutations.resize(1);

      const PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

      const EntityRealignCode pcode =
          EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0);

      permutations[0].change_type(elem_tag, pcode);
    }

    break;
  };
}

// ----------------------------------------------------------------------------

void CellAdaptOpTriagAnisoRefineOrthoFace0::fill_child_coord_transformers(
    std::vector<std::function<void(math::DenseDMat<Real> const &ref_coords,
                                   math::DenseConstVecView<Real> const &coord_in,
                                   math::DenseVecView<Real> &coord_out)>> &transformers) const
{
  transformers.resize(2);
  transformers[0] = child_coords_tri_aniso_refine_face0_child0;
  transformers[1] = child_coords_tri_aniso_refine_face0_child1;
}

// ----------------------------------------------------------------------------
//                   CellAdaptOpTriagAnisoRefineOrthoFace1
// ----------------------------------------------------------------------------

namespace
{

void child_coords_tri_aniso_refine_face1_child0(math::DenseDMat<Real> const &ref_coords,
                                                math::DenseConstVecView<Real> const &coord_in,
                                                math::DenseVecView<Real> &coord_out)
{
  /*
  coord_out[X0] = coord_in[X0] + 0.5 * (1.0 + coord_in[X1]);
  coord_out[X1] = 0.5 * coord_in[X1] - 0.5;
  */

  coord_out[X0] = coord_in[X0] + 0.5 * (ref_coords(1, X0) - ref_coords(0, X0)) /
                                     (ref_coords(2, X1) - ref_coords(0, X1)) *
                                     (coord_in[X1] - ref_coords(0, X1));
  coord_out[X1] = ref_coords(0, X1) + 0.5 * (coord_in[X1] - ref_coords(0, X1));
}

void child_coords_tri_aniso_refine_face1_child1(math::DenseDMat<Real> const &ref_coords,
                                                math::DenseConstVecView<Real> const &coord_in,
                                                math::DenseVecView<Real> &coord_out)
{
  /*
  coord_out[X0] = 0.5 * coord_in[X0] - 0.5;
  coord_out[X1] = coord_in[X1] + 0.5 * (1.0 + coord_in[X0]);
  */

  coord_out[X0] = ref_coords(0, X0) + 0.5 * (coord_in[X0] - ref_coords(0, X0));
  coord_out[X1] = coord_in[X1] +
                  (0.5 * (ref_coords(1, X1) + ref_coords(2, X1)) - ref_coords(0, X1)) /
                      (ref_coords(1, X0) - ref_coords(0, X0)) * (coord_in[X0] - ref_coords(0, X0));
}

} // End of anonymous namespace

CellAdaptOpTriagAnisoRefineOrthoFace1::CellAdaptOpTriagAnisoRefineOrthoFace1() : CellAdaptOpBase()
{
}

CellAdaptOpTriagAnisoRefineOrthoFace1::~CellAdaptOpTriagAnisoRefineOrthoFace1()
{
}

void CellAdaptOpTriagAnisoRefineOrthoFace1::set_facet_adapt_op_ids(
    CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Triag, CellTransform::ANISO_REFINE_ORTHO_FACE_1);

  facet_adapt_op_ids.resize(3);
  facet_adapt_op_ids[0] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[1] = CellTransform::UNIFORM_REFINE;
  facet_adapt_op_ids[2] = CellTransform::NO_TRANS;
}

Uint CellAdaptOpTriagAnisoRefineOrthoFace1::nb_parent_facets() const
{
  return 3u;
}

Uint CellAdaptOpTriagAnisoRefineOrthoFace1::nb_child_elems() const
{
  return 2u;
}

void CellAdaptOpTriagAnisoRefineOrthoFace1::fill_parent_ref_coords(
    math::DenseDMat<Real> &ref_coords) const
{
  ref_coords.resize(3, _2D);

  ref_coords(0, X0) = -1.0;
  ref_coords(0, X1) = -1.0;
  ref_coords(1, X0) = 1.0;
  ref_coords(1, X1) = -1.0;
  ref_coords(2, X0) = -1.0;
  ref_coords(2, X1) = 1.0;
}

void CellAdaptOpTriagAnisoRefineOrthoFace1::set_child_elem_shapes(
    std::vector<ElemShape> &child_elem_shapes) const
{
  child_elem_shapes.resize(2);
  child_elem_shapes[0] = ElemShape::Triag;
  child_elem_shapes[1] = ElemShape::Triag;
}

void CellAdaptOpTriagAnisoRefineOrthoFace1::fill_local_child_incidences(
    std::vector<IncidenceEntry> &incidences, std::vector<EntityDofRealign> &permutations,
    std::vector<Uint> &local_offsets) const
{
  incidences.resize(2);

  incidences[0].cell_idx = 0;
  incidences[0].local_id = 2;
  incidences[1].cell_idx = 1;
  incidences[1].local_id = 0;

  permutations.resize(2);

  PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

  EntityRealignCode identity_p =
      EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0);
  EntityRealignCode flip_p =
      EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0);

  permutations[0].change_type(elem_tag, identity_p);
  permutations[1].change_type(elem_tag, flip_p);

  local_offsets.resize(2);
  local_offsets[0] = 0;
  local_offsets[1] = 2;
}

// ----------------------------------------------------------------------------

void CellAdaptOpTriagAnisoRefineOrthoFace1::fill_parent_child_incidences_on_facet(
    const Uint facet, std::vector<IncidenceEntry> &incidences,
    std::vector<EntityDofRealign> &permutations) const
{

  switch (facet)
  {
    case 0u:
    {
      incidences.resize(1);
      incidences[0].cell_idx = 0;
      incidences[0].local_id = 0;

      permutations.resize(1);

      const PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

      const EntityRealignCode pcode =
          EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0);

      permutations[0].change_type(elem_tag, pcode);
    }

    break;

    case 1u:
    {
      incidences.resize(2);
      incidences[0].cell_idx = 0;
      incidences[0].local_id = 1;
      incidences[1].cell_idx = 1;
      incidences[1].local_id = 1;

      permutations.resize(2);

      const PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

      const EntityRealignCode pcode_0 = EntityRealignCode(
          ElemShape::Line, CellTransform::UNIFORM_REFINE, 0, ElemShape::Line, 0, 0);

      const EntityRealignCode pcode_1 = EntityRealignCode(
          ElemShape::Line, CellTransform::UNIFORM_REFINE, 1, ElemShape::Line, 0, 0);

      permutations[0].change_type(elem_tag, pcode_0);
      permutations[1].change_type(elem_tag, pcode_1);
    }

    break;

    case 2u:
    {
      incidences.resize(1);
      incidences[0].cell_idx = 1;
      incidences[0].local_id = 2;

      permutations.resize(1);

      const PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

      const EntityRealignCode pcode =
          EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0);

      permutations[0].change_type(elem_tag, pcode);
    }

    break;
  };
}

// ----------------------------------------------------------------------------

void CellAdaptOpTriagAnisoRefineOrthoFace1::fill_child_coord_transformers(
    std::vector<std::function<void(math::DenseDMat<Real> const &ref_coords,
                                   math::DenseConstVecView<Real> const &coord_in,
                                   math::DenseVecView<Real> &coord_out)>> &transformers) const
{
  transformers.resize(2);
  transformers[0] = child_coords_tri_aniso_refine_face1_child0;
  transformers[1] = child_coords_tri_aniso_refine_face1_child1;
}

// ----------------------------------------------------------------------------
//                   CellAdaptOpTriagAnisoRefineOrthoFace2
// ----------------------------------------------------------------------------

namespace
{

void child_coords_tri_aniso_refine_face2_child0(math::DenseDMat<Real> const &ref_coords,
                                                math::DenseConstVecView<Real> const &coord_in,
                                                math::DenseVecView<Real> &coord_out)
{
  /*
  coord_out[X0] = coord_in[X0];
  coord_out[X1] = 0.5 * coord_in[X1] - 0.5;
  */

  coord_out[X0] = coord_in[X0];
  coord_out[X1] = ref_coords(0, X1) + 0.5 * (coord_in[X1] - ref_coords(0, X1));
}

void child_coords_tri_aniso_refine_face2_child1(math::DenseDMat<Real> const &ref_coords,
                                                math::DenseConstVecView<Real> const &coord_in,
                                                math::DenseVecView<Real> &coord_out)
{
  /*
  coord_out[X0] = coord_in[X0];
  coord_out[X1] = 0.5 * coord_in[X1] - 0.5 * coord_in[X0];
  */

  coord_out[X0] = coord_in[X0];
  coord_out[X1] = ref_coords(0, X1) + 0.5 * (coord_in[X1] - ref_coords(0, X1)) +
                  0.5 * (ref_coords(2, X1) - ref_coords(0, X1)) -
                  0.5 * (ref_coords(2, X1) - ref_coords(0, X1)) /
                      (ref_coords(1, X0) - ref_coords(0, X0)) * (coord_in[X0] - ref_coords(0, X0));
}

} // End of anonymous namespace

CellAdaptOpTriagAnisoRefineOrthoFace2::CellAdaptOpTriagAnisoRefineOrthoFace2() : CellAdaptOpBase()
{
}

CellAdaptOpTriagAnisoRefineOrthoFace2::~CellAdaptOpTriagAnisoRefineOrthoFace2()
{
}

void CellAdaptOpTriagAnisoRefineOrthoFace2::set_facet_adapt_op_ids(
    CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Triag, CellTransform::ANISO_REFINE_ORTHO_FACE_2);

  facet_adapt_op_ids.resize(3);
  facet_adapt_op_ids[0] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[1] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[2] = CellTransform::UNIFORM_REFINE;
}

Uint CellAdaptOpTriagAnisoRefineOrthoFace2::nb_parent_facets() const
{
  return 3u;
}

Uint CellAdaptOpTriagAnisoRefineOrthoFace2::nb_child_elems() const
{
  return 2u;
}

void CellAdaptOpTriagAnisoRefineOrthoFace2::fill_parent_ref_coords(
    math::DenseDMat<Real> &ref_coords) const
{
  ref_coords.resize(3, _2D);

  ref_coords(0, X0) = -1.0;
  ref_coords(0, X1) = -1.0;
  ref_coords(1, X0) = 1.0;
  ref_coords(1, X1) = -1.0;
  ref_coords(2, X0) = -1.0;
  ref_coords(2, X1) = 1.0;
}

void CellAdaptOpTriagAnisoRefineOrthoFace2::set_child_elem_shapes(
    std::vector<ElemShape> &child_elem_shapes) const
{
  child_elem_shapes.resize(2);
  child_elem_shapes[0] = ElemShape::Triag;
  child_elem_shapes[1] = ElemShape::Triag;
}

void CellAdaptOpTriagAnisoRefineOrthoFace2::fill_local_child_incidences(
    std::vector<IncidenceEntry> &incidences, std::vector<EntityDofRealign> &permutations,
    std::vector<Uint> &local_offsets) const
{
  incidences.resize(2);

  incidences[0].cell_idx = 0;
  incidences[0].local_id = 1;
  incidences[1].cell_idx = 1;
  incidences[1].local_id = 0;

  permutations.resize(2);

  PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

  EntityRealignCode identity_p =
      EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0);
  EntityRealignCode flip_p =
      EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0);

  permutations[0].change_type(elem_tag, identity_p);
  permutations[1].change_type(elem_tag, flip_p);

  local_offsets.resize(2);
  local_offsets[0] = 0;
  local_offsets[1] = 2;
}

// ----------------------------------------------------------------------------

void CellAdaptOpTriagAnisoRefineOrthoFace2::fill_parent_child_incidences_on_facet(
    const Uint facet, std::vector<IncidenceEntry> &incidences,
    std::vector<EntityDofRealign> &permutations) const
{

  switch (facet)
  {
    case 0u:
    {
      incidences.resize(1);
      incidences[0].cell_idx = 0;
      incidences[0].local_id = 0;

      permutations.resize(1);

      const PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

      const EntityRealignCode pcode =
          EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0);

      permutations[0].change_type(elem_tag, pcode);
    }

    break;

    case 1u:
    {
      incidences.resize(1);
      incidences[0].cell_idx = 1;
      incidences[0].local_id = 1;

      permutations.resize(1);

      const PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

      const EntityRealignCode pcode =
          EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0);

      permutations[0].change_type(elem_tag, pcode);
    }

    break;

    case 2u:
    {
      incidences.resize(2);
      incidences[0].cell_idx = 1;
      incidences[0].local_id = 2;
      incidences[1].cell_idx = 0;
      incidences[1].local_id = 2;

      permutations.resize(2);

      const PointSetTag elem_tag(ElemShape::Line, P1, PointSetID::Equidist);

      const EntityRealignCode pcode_0 = EntityRealignCode(
          ElemShape::Line, CellTransform::UNIFORM_REFINE, 0, ElemShape::Line, 0, 0);
      const EntityRealignCode pcode_1 = EntityRealignCode(
          ElemShape::Line, CellTransform::UNIFORM_REFINE, 1, ElemShape::Line, 0, 0);

      permutations[0].change_type(elem_tag, pcode_0);
      permutations[1].change_type(elem_tag, pcode_1);
    }

    break;
  };
}

// ----------------------------------------------------------------------------

void CellAdaptOpTriagAnisoRefineOrthoFace2::fill_child_coord_transformers(
    std::vector<std::function<void(math::DenseDMat<Real> const &ref_coords,
                                   math::DenseConstVecView<Real> const &coord_in,
                                   math::DenseVecView<Real> &coord_out)>> &transformers) const
{
  transformers.resize(2);
  transformers[0] = child_coords_tri_aniso_refine_face2_child0;
  transformers[1] = child_coords_tri_aniso_refine_face2_child1;
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace adapt

} // namespace mesh

} // namespace pdekit
