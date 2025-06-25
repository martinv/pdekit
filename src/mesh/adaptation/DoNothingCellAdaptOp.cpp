#include "mesh/adaptation/DoNothingCellAdaptOp.hpp"
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
//                   DoNothingCellSplitStrategy
// ----------------------------------------------------------------------------

DoNothingCellAdaptOp::DoNothingCellAdaptOp() : CellAdaptOpBase()
{
}

DoNothingCellAdaptOp::~DoNothingCellAdaptOp()
{
}

void DoNothingCellAdaptOp::set_facet_adapt_op_ids(CellAdaptOpTag &tag,
                                                  std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Undefined, CellTransform::NO_TRANS);

  facet_adapt_op_ids.resize(0);
}

Uint DoNothingCellAdaptOp::nb_parent_facets() const
{
  return 0u;
}

Uint DoNothingCellAdaptOp::nb_child_elems() const
{
  return 0u;
}

void DoNothingCellAdaptOp::fill_parent_ref_coords(math::DenseDMat<Real> &ref_coords) const
{
  ref_coords.resize(0, 0);
}

void DoNothingCellAdaptOp::set_child_elem_shapes(std::vector<ElemShape> &child_elem_shapes) const
{
  child_elem_shapes.resize(0);
}

void DoNothingCellAdaptOp::fill_local_child_incidences(std::vector<IncidenceEntry> &incidences,
                                                       std::vector<EntityDofRealign> &permutations,
                                                       std::vector<Uint> &local_offsets) const
{
  incidences.resize(0);
  permutations.resize(0);
  local_offsets.resize(0);
}

// ----------------------------------------------------------------------------

void DoNothingCellAdaptOp::fill_parent_child_incidences_on_facet(
    const Uint facet, std::vector<IncidenceEntry> &incidences,
    std::vector<EntityDofRealign> &permutations) const
{
  incidences.resize(0);
  permutations.resize(0);
}

// ----------------------------------------------------------------------------

void DoNothingCellAdaptOp::fill_child_coord_transformers(
    std::vector<std::function<void(math::DenseDMat<Real> const &ref_coords,
                                   math::DenseConstVecView<Real> const &coord_in,
                                   math::DenseVecView<Real> &coord_out)>> &transformers) const
{
}

// ----------------------------------------------------------------------------
//                          ELEMENT SHAPE: LINE
// ----------------------------------------------------------------------------

void DoNothingCellAdaptOpLine::set_tag(CellAdaptOpTag &tag,
                                       std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Line, CellTransform::NO_TRANS);
  facet_adapt_op_ids.resize(0);
}

// ----------------------------------------------------------------------------
//                        ELEMENT SHAPE: TRIANGLE
// ----------------------------------------------------------------------------

void DoNothingCellAdaptOpTriag::set_tag(CellAdaptOpTag &tag,
                                        std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Triag, CellTransform::NO_TRANS);

  facet_adapt_op_ids.resize(3);
  facet_adapt_op_ids[0] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[1] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[2] = CellTransform::NO_TRANS;
}

// ----------------------------------------------------------------------------
//                      ELEMENT SHAPE: QUADRILATERAL
// ----------------------------------------------------------------------------

void DoNothingCellAdaptOpQuad::set_tag(CellAdaptOpTag &tag,
                                       std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Quad, CellTransform::NO_TRANS);

  facet_adapt_op_ids.resize(4);
  facet_adapt_op_ids[0] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[1] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[2] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[3] = CellTransform::NO_TRANS;
}

// ----------------------------------------------------------------------------
//                       ELEMENT SHAPE: TETRAHEDRON
// ----------------------------------------------------------------------------

void DoNothingCellAdaptOpTetra::set_tag(CellAdaptOpTag &tag,
                                        std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Tetra, CellTransform::NO_TRANS);

  facet_adapt_op_ids.resize(4);
  facet_adapt_op_ids[0] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[1] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[2] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[3] = CellTransform::NO_TRANS;
}

// ----------------------------------------------------------------------------
//                       ELEMENT SHAPE: HEXAHEDRON
// ----------------------------------------------------------------------------

void DoNothingCellAdaptOpHexa::set_tag(CellAdaptOpTag &tag,
                                       std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Hexa, CellTransform::NO_TRANS);

  facet_adapt_op_ids.resize(6);
  facet_adapt_op_ids[0] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[1] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[2] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[3] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[4] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[5] = CellTransform::NO_TRANS;
}

// ----------------------------------------------------------------------------
//                        ELEMENT SHAPE: PRISM
// ----------------------------------------------------------------------------

void DoNothingCellAdaptOpPrism::set_tag(CellAdaptOpTag &tag,
                                        std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Prism, CellTransform::NO_TRANS);

  facet_adapt_op_ids.resize(5);
  facet_adapt_op_ids[0] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[1] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[2] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[3] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[4] = CellTransform::NO_TRANS;
}

// ----------------------------------------------------------------------------
//                        ELEMENT SHAPE: PYRAMID
// ----------------------------------------------------------------------------

void DoNothingCellAdaptOpPyramid::set_tag(CellAdaptOpTag &tag,
                                          std::vector<CellTransform> &facet_adapt_op_ids)
{
  tag = CellAdaptOpTag(ElemShape::Pyramid, CellTransform::NO_TRANS);

  facet_adapt_op_ids.resize(5);
  facet_adapt_op_ids[0] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[1] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[2] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[3] = CellTransform::NO_TRANS;
  facet_adapt_op_ids[4] = CellTransform::NO_TRANS;
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace adapt

} // namespace mesh

} // namespace pdekit
