#include "mesh/adaptation/CellAdaptOpQuad.hpp"
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
//                   CellAdaptOpQuadUniformRefine
// ----------------------------------------------------------------------------

namespace
{

void child_coords_quad_uniform_refine_child0(math::DenseDMat<Real> const &ref_coords,
                                             math::DenseConstVecView<Real> const &coord_in,
                                             math::DenseVecView<Real> &coord_out)
{
  coord_out[X0] = 0.5 * coord_in[X0] - 0.5;
  coord_out[X1] = 0.5 * coord_in[X1] - 0.5;
}

void child_coords_quad_uniform_refine_child1(math::DenseDMat<Real> const &ref_coords,
                                             math::DenseConstVecView<Real> const &coord_in,
                                             math::DenseVecView<Real> &coord_out)
{
  coord_out[X0] = 0.5 * coord_in[X0] + 0.5;
  coord_out[X1] = 0.5 * coord_in[X1] - 0.5;
}

void child_coords_quad_uniform_refine_child2(math::DenseDMat<Real> const &ref_coords,
                                             math::DenseConstVecView<Real> const &coord_in,
                                             math::DenseVecView<Real> &coord_out)
{
  coord_out[X0] = 0.5 * coord_in[X0] + 0.5;
  coord_out[X1] = 0.5 * coord_in[X1] + 0.5;
}

void child_coords_quad_uniform_refine_child3(math::DenseDMat<Real> const &ref_coords,
                                             math::DenseConstVecView<Real> const &coord_in,
                                             math::DenseVecView<Real> &coord_out)
{
  coord_out[X0] = 0.5 * coord_in[X0] - 0.5;
  coord_out[X1] = 0.5 * coord_in[X1] + 0.5;
}

} // End of anonymous namespace

CellAdaptOpQuadUniformRefine::CellAdaptOpQuadUniformRefine() : CellAdaptOpBase()
{
}

CellAdaptOpQuadUniformRefine::~CellAdaptOpQuadUniformRefine()
{
}

void CellAdaptOpQuadUniformRefine::set_facet_adapt_op_ids(
    CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Quad, CellTransform::UNIFORM_REFINE);

  facet_adapt_op_ids.resize(4);
  facet_adapt_op_ids[0] = CellTransform::UNIFORM_REFINE;
  facet_adapt_op_ids[1] = CellTransform::UNIFORM_REFINE;
  facet_adapt_op_ids[2] = CellTransform::UNIFORM_REFINE;
  facet_adapt_op_ids[3] = CellTransform::UNIFORM_REFINE;
}

Uint CellAdaptOpQuadUniformRefine::nb_parent_facets() const
{
  return 4u;
}

Uint CellAdaptOpQuadUniformRefine::nb_child_elems() const
{
  return 4u;
}

void CellAdaptOpQuadUniformRefine::fill_parent_ref_coords(math::DenseDMat<Real> &ref_coords) const
{
  ref_coords.resize(4, _2D);

  ref_coords(0, X0) = -1.0;
  ref_coords(0, X1) = -1.0;
  ref_coords(1, X0) = 1.0;
  ref_coords(1, X1) = -1.0;
  ref_coords(2, X0) = 1.0;
  ref_coords(2, X1) = 1.0;
  ref_coords(3, X0) = -1.0;
  ref_coords(3, X1) = 1.0;
}

void CellAdaptOpQuadUniformRefine::set_child_elem_shapes(
    std::vector<ElemShape> &child_elem_shapes) const
{
  child_elem_shapes.resize(4);
  child_elem_shapes[0] = ElemShape::Quad;
  child_elem_shapes[1] = ElemShape::Quad;
  child_elem_shapes[2] = ElemShape::Quad;
  child_elem_shapes[3] = ElemShape::Quad;
}

void CellAdaptOpQuadUniformRefine::fill_local_child_incidences(
    std::vector<IncidenceEntry> &incidences, std::vector<EntityDofRealign> &permutations,
    std::vector<Uint> &local_offsets) const
{
  incidences.resize(8);

  incidences[0].cell_idx = 0;
  incidences[0].local_id = 1;
  incidences[1].cell_idx = 1;
  incidences[1].local_id = 3;

  incidences[2].cell_idx = 1;
  incidences[2].local_id = 2;
  incidences[3].cell_idx = 2;
  incidences[3].local_id = 0;

  incidences[4].cell_idx = 2;
  incidences[4].local_id = 3;
  incidences[5].cell_idx = 3;
  incidences[5].local_id = 1;

  incidences[6].cell_idx = 0;
  incidences[6].local_id = 2;
  incidences[7].cell_idx = 3;
  incidences[7].local_id = 0;

  permutations.resize(8);

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

  permutations[6].change_type(elem_tag, identity_p);
  permutations[7].change_type(elem_tag, flip_p);

  local_offsets.resize(5);
  local_offsets[0] = 0;
  local_offsets[1] = 2;
  local_offsets[2] = 4;
  local_offsets[3] = 6;
  local_offsets[4] = 8;
}

// ----------------------------------------------------------------------------

void CellAdaptOpQuadUniformRefine::fill_parent_child_incidences_on_facet(
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
      incidences[1].cell_idx = 3;
      incidences[1].local_id = 2;
      break;

    case 3u:
      incidences[0].cell_idx = 3;
      incidences[0].local_id = 3;
      incidences[1].cell_idx = 0;
      incidences[1].local_id = 3;
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

void CellAdaptOpQuadUniformRefine::fill_child_coord_transformers(
    std::vector<std::function<void(math::DenseDMat<Real> const &ref_coords,
                                   math::DenseConstVecView<Real> const &coord_in,
                                   math::DenseVecView<Real> &coord_out)>> &transformers) const
{
  transformers.resize(4);
  transformers[0] = child_coords_quad_uniform_refine_child0;
  transformers[1] = child_coords_quad_uniform_refine_child1;
  transformers[2] = child_coords_quad_uniform_refine_child2;
  transformers[3] = child_coords_quad_uniform_refine_child3;
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace adapt

} // namespace mesh

} // namespace pdekit
