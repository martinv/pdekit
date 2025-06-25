#include <cassert>

#include "mesh/CellTransform.hpp"
#include "mesh/local_topology/CellSubdomainTag.hpp"
#include "mesh/local_topology/TraceTupleLine.hpp"

namespace pdekit
{

namespace mesh
{

namespace internal
{

// ----------------------------------------------------------------------------
// I) Two incident line segments, neither is a subdomain (i.e. obtained by
//    local refinement) of another segment
// ----------------------------------------------------------------------------

LineL_L0_LineR_L0::LineL_L0_LineR_L0() : TraceTupleBase()
{
}

LineL_L0_LineR_L0::~LineL_L0_LineR_L0()
{
}

void LineL_L0_LineR_L0::set_subdomain_tags(CellSubdomainTag &sub_tag_L, CellSubdomainTag &sub_tag_R)
{
  sub_tag_L = CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line);
  sub_tag_R = CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line);
}

void LineL_L0_LineR_L0::set_incidence_comp_functions(
    std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
{
  functions.push_back(
      std::make_tuple(CellTransform::NO_TRANS, CellTransform::NO_TRANS, do_nothing));
  functions.push_back(std::make_tuple(CellTransform::UNIFORM_REFINE, CellTransform::NO_TRANS,
                                      split_left_side_two_segments));
  functions.push_back(std::make_tuple(CellTransform::NO_TRANS, CellTransform::UNIFORM_REFINE,
                                      split_right_side_two_segments));
  functions.push_back(std::make_tuple(CellTransform::UNIFORM_REFINE, CellTransform::UNIFORM_REFINE,
                                      split_both_sides_two_segments));
}

// ----------------------------------------------------------------------------

void LineL_L0_LineR_L0::fill_incident_p1_nodes(const EntityRealignCode input_facet_perm_L,
                                               const EntityRealignCode input_facet_perm_R,
                                               std::vector<std::pair<Uint, Uint>> &p1_verts) const
{
  const Uint nb_transforms_left =
      (input_facet_perm_L.nb_flips() + input_facet_perm_L.nb_rotations()) % 2;
  const Uint nb_transforms_right =
      (input_facet_perm_R.nb_flips() + input_facet_perm_R.nb_rotations()) % 2;

  p1_verts.resize(2);

  // Case where the coordinate systems of two adjacent line faces are aligned
  // without transforming either of them or performing one rotation/flip on
  // BOTH of them. In this case, local node 0 of left line is incident to
  // local node 0 of right line
  //               local node 1 of left line is incident to local node 1 of
  //               right line
  // Note: This should actually neven occur in 2D mesh!
  if (nb_transforms_left == nb_transforms_right)
  {
    p1_verts[0].first  = 0;
    p1_verts[0].second = 0;
    p1_verts[1].first  = 1;
    p1_verts[1].second = 1;
  }
  else // The only other option is that one line is transformed (rotated or
       // flipped once) and the other is not
  {
    p1_verts[0].first  = 0;
    p1_verts[0].second = 1;
    p1_verts[1].first  = 1;
    p1_verts[1].second = 0;
  }
}

// ----------------------------------------------------------------------------

void LineL_L0_LineR_L0::do_nothing(const EntityRealignCode input_facet_perm_L,
                                   const EntityRealignCode input_facet_perm_R,
                                   const std::vector<IncidenceEntry> &incidences_subcells_L,
                                   const std::vector<IncidenceEntry> &incidences_subcells_R,
                                   std::vector<std::vector<IncidenceEntry>> &incidences,
                                   std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(1);
  incidences[0].resize(2);
  incidences[0][0] = incidences_subcells_L[0];
  incidences[0][1] = incidences_subcells_R[0];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  facet_permutations.resize(1);
  facet_permutations[0].resize(2);
  facet_permutations[0][0].change_type(line_std_reg_tag, input_facet_perm_L);
  facet_permutations[0][1].change_type(line_std_reg_tag, input_facet_perm_R);
}

// ----------------------------------------------------------------------------

void LineL_L0_LineR_L0::split_left_side_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  // Compute the number of flips needed to make the RIGHT (outer) facet
  // permutation an identity const Uint nb_flips_r = (
  // input_facet_perm_R.nb_flips() % 2 == 0 ) ? 0 : 1; const Uint nb_rot_r = (
  // input_facet_perm_R.nb_rotations() % 2 == 0 ) ? 0 : 1;

  const EntityRealignCode pcode_id   = EntityRealignCode::identity(ElemShape::Line);
  const EntityRealignCode pcode_flip = EntityRealignCode::single_flip(ElemShape::Line);

  incidences.resize(2);
  incidences[0].resize(2);
  incidences[1].resize(2);

  facet_permutations.resize(2);
  facet_permutations[0].resize(2);
  facet_permutations[1].resize(2);

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped
  // The new permutations have to be generated so that identity is ON THE
  // RIGHT
  if ((input_facet_perm_L == pcode_id) && (input_facet_perm_R == pcode_flip))
  {
    // New child subcell on position 0 on the left is incident to the
    // 'parent level' cell on the right
    incidences[0][LEFT]  = incidences_subcells_L[0];
    incidences[0][RIGHT] = incidences_subcells_R[0];

    // New child subcell on position 1 on the left is incident to the
    // 'parent level' cell on the right
    incidences[1][LEFT]  = incidences_subcells_L[1];
    incidences[1][RIGHT] = incidences_subcells_R[0];

    // The following means that the the first left entity (one sub-triangle
    // face after the splitting of the left master element) is incident to
    // the subdomain 0 of right facet (which was not split) and the second
    // left entity (one facet of another child triangle on the left) is
    // incident to the subdomain 1 of right facet

    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][0].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
    facet_permutations[0][1].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 1, 0));

    // Second row: second pair of incidences
    facet_permutations[1][0].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
    facet_permutations[1][1].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 1, 0));
  }
  // Second case: entity permutation on the left is one flip
  // Right facet permutation is identity
  // The new permutations have to be generated so that identity is ON THE
  // RIGHT
  else if ((input_facet_perm_L == pcode_flip) && (input_facet_perm_R == pcode_id))
  {
    // New child subcell on position 0 on the left is incident to the
    // 'parent level' cell on the right
    incidences[0][LEFT]  = incidences_subcells_L[1];
    incidences[0][RIGHT] = incidences_subcells_R[0];

    // New child subcell on position 1 on the left is incident to the
    // 'parent level' cell on the right
    incidences[1][LEFT]  = incidences_subcells_L[0];
    incidences[1][RIGHT] = incidences_subcells_R[0];

    // The following means that the the first left entity (one sub-triangle
    // face after the splitting of the left master element) is incident to
    // the subdomain 1 of right facet (which was not split) and the second
    // left entity (one facet of another child triangle on the left) is
    // incident to the subdomain 0 of right facet

    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 0, 0));

    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "LocalIndicencePatternLine::line_no_split_line_no_split_left_"
                 "side_two_segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations of incident\n"
              << " facets" << std::endl;
    return;
  }
}

// ----------------------------------------------------------------------------

void LineL_L0_LineR_L0::split_right_side_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  // Compute the number of flips needed to make the RIGHT (outer) facet
  // permutation an identity const Uint nb_flips_r = (
  // input_facet_perm_R.nb_flips() % 2 == 0 ) ? 0 : 1; const Uint nb_rot_r = (
  // input_facet_perm_R.nb_rotations() % 2 == 0 ) ? 0 : 1;

  const EntityRealignCode pcode_id   = EntityRealignCode::identity(ElemShape::Line);
  const EntityRealignCode pcode_flip = EntityRealignCode::single_flip(ElemShape::Line);

  /*
  std::cout << "Identity:    " << pcode_id.as_string() << std::endl;
  std::cout << "Flip:        " << pcode_flip.as_string() << std::endl;
  std::cout << "Input left:  " << input_facet_perm_L.as_string() << std::endl;
  std::cout << "Input right: " << input_facet_perm_R.as_string() << std::endl;
  */

  incidences.resize(2);
  incidences[0].resize(2);
  incidences[1].resize(2);

  facet_permutations.resize(2);
  facet_permutations[0].resize(2);
  facet_permutations[1].resize(2);

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped
  // Expected result for incidences:
  // (Left subdomain 0)  - (Right child cell 1)
  // (Left subdomain 1)  - (Right child cell 0)
  if ((input_facet_perm_L == pcode_id) && (input_facet_perm_R == pcode_flip))
  {
    // Cell on the left is incident to the 'child level' cell with index 0
    // on the right
    incidences[0][LEFT]  = incidences_subcells_L[0];
    incidences[0][RIGHT] = incidences_subcells_R[1];

    // New child subcell on position 1 on the left is incident to the
    // 'parent level' cell on the right
    incidences[1][LEFT]  = incidences_subcells_L[0];
    incidences[1][RIGHT] = incidences_subcells_R[0];

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 0, 0));
    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    // The following means that the the first right entity (one sub-triangle
    // face after the splitting of the right master element) is incident to
    // the subdomain 0 of left facet (which was not split) and the second
    // right entity (one facet of another child triangle on the right) is
    // incident to the subdomain 1 of left facet

    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 0, 0));
    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
  }
  // Second case: entity permutation on the left is identity
  // Right facet is flipped
  // Expected result for incidences:
  // (Left subdomain 0)  - (Right child cell 1)
  // (Left subdomain 1)  - (Right child cell 0)
  else if ((input_facet_perm_L == pcode_flip) && (input_facet_perm_R == pcode_id))
  {
    // New child subcell on position 0 on the right is incident to the
    // 'parent level' cell on the left
    incidences[0][LEFT]  = incidences_subcells_L[0];
    incidences[0][RIGHT] = incidences_subcells_R[0];

    // New child subcell on position 1 on the right is incident to the
    // 'parent level' cell on the left
    incidences[1][LEFT]  = incidences_subcells_L[0];
    incidences[1][RIGHT] = incidences_subcells_R[1];

    // The following means that the the first right entity (one sub-triangle
    // face after the splitting of the right master element) is incident to
    // the subdomain 0 of left facet (which was not split) and the second
    // right entity (one facet of another child triangle on the right) is
    // incident to the subdomain 1 of left facet

    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 1, 0));
    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 1, 0));
    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "CellSubdomainTupleLine::line_no_split_line_no_split_right_"
                 "side_two_segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations of incident\n"
              << " facets" << std::endl;
    return;
  }
}

// ----------------------------------------------------------------------------

void LineL_L0_LineR_L0::split_both_sides_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  // Compute the number of flips needed to make the RIGHT (outer) facet
  // permutation an identity const Uint nb_flips_r = (
  // input_facet_perm_R.nb_flips() % 2 == 0 ) ? 0 : 1; const Uint nb_rot_r = (
  // input_facet_perm_R.nb_rotations() % 2 == 0 ) ? 0 : 1;

  const EntityRealignCode pcode_id   = EntityRealignCode::identity(ElemShape::Line);
  const EntityRealignCode pcode_flip = EntityRealignCode::single_flip(ElemShape::Line);

  incidences.resize(2);
  incidences[0].resize(2);
  incidences[1].resize(2);

  facet_permutations.resize(2);
  facet_permutations[0].resize(2);
  facet_permutations[1].resize(2);

  // Child cell 0 on the left is incident to the child cell 1 on the right
  incidences[0][0] = incidences_subcells_L[0];
  incidences[0][1] = incidences_subcells_R[1];

  // Child cell 1 on the left is incident to the child cell 0 on the right
  incidences[1][0] = incidences_subcells_L[1];
  incidences[1][1] = incidences_subcells_R[0];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped
  // Expected result for incidences:
  // (Left child cell 0)  - (Right child cell 1)
  // (Left child cell 1)  - (Right child cell 0)
  if ((input_facet_perm_L == pcode_id) && (input_facet_perm_R == pcode_flip))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent
    ///

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
  }
  // Second case: entity permutation on the left is identity
  // Right facet is flipped
  // Expected result for incidences:
  // (Left subdomain 0)  - (Right child cell 1)
  // (Left subdomain 1)  - (Right child cell 0)
  else if ((input_facet_perm_L == pcode_flip) && (input_facet_perm_R == pcode_id))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "LocalIndicencePatternLine::line_no_split_line_no_split_both_"
                 "sides_two_segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations of incident\n"
              << " facets" << std::endl;
    return;
  }
}

// ----------------------------------------------------------------------------
// IIa) Two incident line segments, right one is a subdomain of another segment
//      that had already been refined
// ----------------------------------------------------------------------------

LineL_L1_LineR_L0PartType1Pos0::LineL_L1_LineR_L0PartType1Pos0() : TraceTupleBase()
{
}

LineL_L1_LineR_L0PartType1Pos0::~LineL_L1_LineR_L0PartType1Pos0()
{
}

void LineL_L1_LineR_L0PartType1Pos0::set_subdomain_tags(CellSubdomainTag &sub_tag_L,
                                                        CellSubdomainTag &sub_tag_R)
{
  sub_tag_L = CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line);
  sub_tag_R = CellSubdomainTag(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0, ElemShape::Line);
}

void LineL_L1_LineR_L0PartType1Pos0::set_incidence_comp_functions(
    std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
{
  functions.push_back(
      std::make_tuple(CellTransform::NO_TRANS, CellTransform::NO_TRANS, do_nothing));
  functions.push_back(std::make_tuple(CellTransform::NO_TRANS, CellTransform::UNIFORM_REFINE,
                                      split_right_side_two_segments));
  functions.push_back(std::make_tuple(CellTransform::UNIFORM_REFINE, CellTransform::UNIFORM_REFINE,
                                      split_both_sides_two_segments));
}

void LineL_L1_LineR_L0PartType1Pos0::fill_incident_p1_nodes(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    std::vector<std::pair<Uint, Uint>> &p1_verts) const
{
  const Uint nb_transforms_left =
      (input_facet_perm_L.nb_flips() + input_facet_perm_L.nb_rotations()) % 2;
  const Uint nb_transforms_right =
      (input_facet_perm_R.nb_flips() + input_facet_perm_R.nb_rotations()) % 2;

  p1_verts.resize(1);

  // Left face interfaces a sub-segment of the right face. Since this is
  // 'Partition type 1' and 'Position 0', it means that only the local vertex
  // 0 of the right (larger) face can be in contact with the left face. Hence
  // p1_verts[0].second MUST BE 0 in all cases. Now depending on whether both
  // systems are aligned or not, this vertex touches the local vertex 0 or
  // local vertex 1 of the left face.

  p1_verts[0].second = 0;

  if (nb_transforms_left == nb_transforms_right)
  {
    p1_verts[0].first = 0;
  }
  else
  {
    p1_verts[0].first = 1;
  }
}

void LineL_L1_LineR_L0PartType1Pos0::do_nothing(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(1);
  incidences[0].resize(2);
  incidences[0][LEFT]  = incidences_subcells_L[0];
  incidences[0][RIGHT] = incidences_subcells_R[0];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  facet_permutations.resize(1);
  facet_permutations[0].resize(2);
  facet_permutations[0][LEFT].change_type(line_std_reg_tag, input_facet_perm_L);
  facet_permutations[0][RIGHT].change_type(line_std_reg_tag, input_facet_perm_R);
}

void LineL_L1_LineR_L0PartType1Pos0::split_right_side_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(1);
  incidences[0].resize(2);

  facet_permutations.resize(1);
  facet_permutations[0].resize(2);

  assert(incidences_subcells_L.size() == 1);
  assert(incidences_subcells_R.size() == 2);

  // Child facet 0 on the left is incident to the child facet 0 on the right
  // Note that this is BY DEFINITION: line on left is on level 0 and not split
  // - only one entity to match The line on the right is larger/coarser: one
  // level above, but we're splitting it here This means that after it's
  // split, one of its children will match left facet and the other child will
  // not be considered here at all We also know that the left line is incident
  // to subregion 1 (SEE THE NAME OF THE FUNCTION!) of the right line
  incidences[0][LEFT]  = incidences_subcells_L[0];
  incidences[0][RIGHT] = incidences_subcells_R[0];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped

  if ((input_facet_perm_L.nb_flips() == 0) && (input_facet_perm_R.nb_flips() == 1))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
  }

  // Second case: left facet is flipped
  // Right facet is identity

  else if ((input_facet_perm_L.nb_flips() == 1) && (input_facet_perm_R.nb_flips() == 0))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "LineL_Level0_LineR_Level1Pos0::split_right_side_two_"
                 "segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations "
                 "of incident\n"
              << " facets" << std::endl;
    return;
  }
}

void LineL_L1_LineR_L0PartType1Pos0::split_both_sides_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(2);
  incidences[0].resize(2);
  incidences[1].resize(2);

  facet_permutations.resize(2);
  facet_permutations[0].resize(2);
  facet_permutations[1].resize(2);

  assert(incidences_subcells_L.size() == 2);
  assert(incidences_subcells_R.size() == 2);

  // Note that the following two identities are BY DEFINITION
  // line on left is on level 0 and split into two sub-entities on level 1
  // The line on the right is already on level 1, but will be split as well
  // Only its child corresponding to to subregion 0 (SEE THE NAME OF THE
  // FUNCTION!) will appear in the incidence identities!

  // Child facet 0 on the left is incident to the child facet 0 on the right
  incidences[0][LEFT]  = incidences_subcells_L[0];
  incidences[0][RIGHT] = incidences_subcells_R[0];

  // Child facet 1 on the left is incident to the child facet 1 on the right
  incidences[1][LEFT]  = incidences_subcells_L[1];
  incidences[1][RIGHT] = incidences_subcells_R[0];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped

  if ((input_facet_perm_L.nb_flips() == 0) && (input_facet_perm_R.nb_flips() == 1))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 1, 0));

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 1, 0));
  }

  // Second case: left facet is flipped
  // Right facet is identity

  else if ((input_facet_perm_L.nb_flips() == 1) && (input_facet_perm_R.nb_flips() == 0))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 0, 0));

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "LineL_Level0_LineR_Level1Pos0::split_both_sides_two_"
                 "segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations "
                 "of incident\n"
              << " facets" << std::endl;
    return;
  }
}

// ----------------------------------------------------------------------------
// IIb) Two incident line segments, right one is a subdomain of another segment
//      that had already been refined
// ----------------------------------------------------------------------------

LineL_L1_LineR_L0PartType1Pos1::LineL_L1_LineR_L0PartType1Pos1() : TraceTupleBase()
{
}

LineL_L1_LineR_L0PartType1Pos1::~LineL_L1_LineR_L0PartType1Pos1()
{
}

void LineL_L1_LineR_L0PartType1Pos1::set_subdomain_tags(CellSubdomainTag &sub_tag_L,
                                                        CellSubdomainTag &sub_tag_R)
{
  sub_tag_L = CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line);
  sub_tag_R = CellSubdomainTag(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1, ElemShape::Line);
}

void LineL_L1_LineR_L0PartType1Pos1::set_incidence_comp_functions(
    std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
{
  functions.push_back(
      std::make_tuple(CellTransform::NO_TRANS, CellTransform::NO_TRANS, do_nothing));
  functions.push_back(std::make_tuple(CellTransform::NO_TRANS, CellTransform::UNIFORM_REFINE,
                                      split_right_side_two_segments));
  functions.push_back(std::make_tuple(CellTransform::UNIFORM_REFINE, CellTransform::UNIFORM_REFINE,
                                      split_both_sides_two_segments));
}

void LineL_L1_LineR_L0PartType1Pos1::fill_incident_p1_nodes(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    std::vector<std::pair<Uint, Uint>> &p1_verts) const
{
  const Uint nb_transforms_left =
      (input_facet_perm_L.nb_flips() + input_facet_perm_L.nb_rotations()) % 2;
  const Uint nb_transforms_right =
      (input_facet_perm_R.nb_flips() + input_facet_perm_R.nb_rotations()) % 2;

  p1_verts.resize(1);

  // Left face interfaces a sub-segment of the right face. Since this is
  // 'Partition type 1' and 'Position 1', it means that only the local vertex
  // 1 of the right (larger) face can be in contact with the left face. Hence
  // p1_verts[0].second MUST BE 1 in all cases. Now depending on whether both
  // systems are aligned or not, this vertex touches the local vertex 0 or
  // local vertex 1 of the left face.

  p1_verts[0].second = 1;

  if (nb_transforms_left == nb_transforms_right)
  {
    p1_verts[0].first = 1;
  }
  else
  {
    p1_verts[0].first = 0;
  }
}

void LineL_L1_LineR_L0PartType1Pos1::do_nothing(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(1);
  incidences[0].resize(2);
  incidences[0][LEFT]  = incidences_subcells_L[0];
  incidences[0][RIGHT] = incidences_subcells_R[0];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  facet_permutations.resize(1);
  facet_permutations[0].resize(2);
  facet_permutations[0][LEFT].change_type(line_std_reg_tag, input_facet_perm_L);
  facet_permutations[0][RIGHT].change_type(line_std_reg_tag, input_facet_perm_R);
}

void LineL_L1_LineR_L0PartType1Pos1::split_right_side_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(1);
  incidences[0].resize(2);

  facet_permutations.resize(1);
  facet_permutations[0].resize(2);

  assert(incidences_subcells_L.size() == 1);
  assert(incidences_subcells_R.size() == 2);

  // Child facet 0 on the left is incident to the child facet 0 on the right
  // Note that this is BY DEFINITION: line on left is on level 0 and not split
  // - only one entity to match The line on the right is larger/coarser: one
  // level above, but we're splitting it here This means that after it's
  // split, one of its children will match left facet and the other child will
  // not be considered here at all We also know that the left line is incident
  // to subregion 1 (SEE THE NAME OF THE FUNCTION!) of the right line
  incidences[0][LEFT]  = incidences_subcells_L[0];
  incidences[0][RIGHT] = incidences_subcells_R[1];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped

  if ((input_facet_perm_L.nb_flips() == 0) && (input_facet_perm_R.nb_flips() == 1))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
  }

  // Second case: left facet is flipped
  // Right facet is identity

  else if ((input_facet_perm_L.nb_flips() == 1) && (input_facet_perm_R.nb_flips() == 0))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "LineL_Level0_LineR_Level1Pos1::split_both_sides_two_"
                 "segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations "
                 "of incident\n"
              << " facets" << std::endl;
    return;
  }
}

void LineL_L1_LineR_L0PartType1Pos1::split_both_sides_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(2);
  incidences[0].resize(2);
  incidences[1].resize(2);

  facet_permutations.resize(2);
  facet_permutations[0].resize(2);
  facet_permutations[1].resize(2);

  assert(incidences_subcells_L.size() == 2);
  assert(incidences_subcells_R.size() == 2);

  // Note that the following two identities are BY DEFINITION
  // line on left is on level 0 and split into two sub-entities on level 1
  // The line on the right is already on level 1, but will be split as well
  // Only its child corresponding to to subregion 0 (SEE THE NAME OF THE
  // FUNCTION!) will appear in the incidence identities!

  // Child facet 0 on the left is incident to the child facet 0 on the right
  incidences[0][LEFT]  = incidences_subcells_L[0];
  incidences[0][RIGHT] = incidences_subcells_R[1];

  // Child facet 1 on the left is incident to the child facet 1 on the right
  incidences[1][LEFT]  = incidences_subcells_L[1];
  incidences[1][RIGHT] = incidences_subcells_R[1];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped

  if ((input_facet_perm_L.nb_flips() == 0) && (input_facet_perm_R.nb_flips() == 1))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 1, 0));

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 1, 0));
  }

  // Second case: left facet is flipped
  // Right facet is identity

  else if ((input_facet_perm_L.nb_flips() == 1) && (input_facet_perm_R.nb_flips() == 0))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 0, 0));

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "LineL_Level0_LineR_Level1Pos1::split_both_sides_two_"
                 "segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations "
                 "of incident\n"
              << " facets" << std::endl;
    return;
  }
}

// ----------------------------------------------------------------------------
// IIIa) Two incident line segments, left one is a subdomain of another segment
//       and incident to right facet that had already been refined
// ----------------------------------------------------------------------------

LineL_L0PartType1Pos0_LineR_L1::LineL_L0PartType1Pos0_LineR_L1() : TraceTupleBase()
{
}

LineL_L0PartType1Pos0_LineR_L1::~LineL_L0PartType1Pos0_LineR_L1()
{
}

void LineL_L0PartType1Pos0_LineR_L1::set_subdomain_tags(CellSubdomainTag &sub_tag_L,
                                                        CellSubdomainTag &sub_tag_R)
{
  sub_tag_L = CellSubdomainTag(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0, ElemShape::Line);
  sub_tag_R = CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line);
}

void LineL_L0PartType1Pos0_LineR_L1::set_incidence_comp_functions(
    std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
{
  functions.push_back(
      std::make_tuple(CellTransform::NO_TRANS, CellTransform::NO_TRANS, do_nothing));
  functions.push_back(std::make_tuple(CellTransform::UNIFORM_REFINE, CellTransform::NO_TRANS,
                                      split_left_side_two_segments));
  functions.push_back(std::make_tuple(CellTransform::UNIFORM_REFINE, CellTransform::UNIFORM_REFINE,
                                      split_both_sides_two_segments));
}

void LineL_L0PartType1Pos0_LineR_L1::fill_incident_p1_nodes(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    std::vector<std::pair<Uint, Uint>> &p1_verts) const
{
  const Uint nb_transforms_left =
      (input_facet_perm_L.nb_flips() + input_facet_perm_L.nb_rotations()) % 2;
  const Uint nb_transforms_right =
      (input_facet_perm_R.nb_flips() + input_facet_perm_R.nb_rotations()) % 2;

  p1_verts.resize(1);

  // Right face interfaces a sub-segment of the left face. Since this is
  // 'Partition type 1' and 'Position 0', it means that only the local vertex
  // 0 of the left (larger) face can be in contact with the right face. Hence
  // p1_verts[0].first MUST BE 0 in all cases. Now depending on whether both
  // systems are aligned or not, this vertex touches the local vertex 0 or
  // local vertex 1 of the right face.

  p1_verts[0].first = 0;

  if (nb_transforms_left == nb_transforms_right)
  {
    p1_verts[0].second = 0;
  }
  else
  {
    p1_verts[0].second = 1;
  }
}

void LineL_L0PartType1Pos0_LineR_L1::do_nothing(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(1);
  incidences[0].resize(2);
  incidences[0][LEFT]  = incidences_subcells_L[0];
  incidences[0][RIGHT] = incidences_subcells_R[0];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  facet_permutations.resize(1);
  facet_permutations[0].resize(2);
  facet_permutations[0][LEFT].change_type(line_std_reg_tag, input_facet_perm_L);
  facet_permutations[0][RIGHT].change_type(line_std_reg_tag, input_facet_perm_R);
}

void LineL_L0PartType1Pos0_LineR_L1::split_left_side_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(1);
  incidences[0].resize(2);

  facet_permutations.resize(1);
  facet_permutations[0].resize(2);

  assert(incidences_subcells_L.size() == 1);
  assert(incidences_subcells_R.size() == 2);

  // Child facet 0 on the left is incident to the child facet 0 on the right
  // Note that this is BY DEFINITION: line on the right is on level 0 and not
  // split - only one entity to match The line on the left is larger/coarser:
  // one level above, but we're using only a part of it - the right line is
  // incident to subregion 1 (SEE THE NAME OF THE FUNCTION!) of the left line
  incidences[0][LEFT]  = incidences_subcells_L[0];
  incidences[0][RIGHT] = incidences_subcells_R[0];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped

  if ((input_facet_perm_L.nb_flips() == 0) && (input_facet_perm_R.nb_flips() == 1))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
  }

  // Second case: left facet is flipped
  // Right facet is identity

  else if ((input_facet_perm_L.nb_flips() == 1) && (input_facet_perm_R.nb_flips() == 0))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "LineL_Level1Pos0_LineR_Level0::split_left_side_two_segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations "
                 "of incident\n"
              << " facets" << std::endl;
    return;
  }
}

void LineL_L0PartType1Pos0_LineR_L1::split_both_sides_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(2);
  incidences[0].resize(2);
  incidences[1].resize(2);

  facet_permutations.resize(2);
  facet_permutations[0].resize(2);
  facet_permutations[1].resize(2);

  assert(incidences_subcells_L.size() == 2);
  assert(incidences_subcells_R.size() == 2);

  // Note that the following two identities are BY DEFINITION
  // line on left is on level 0 and split into two sub-entities on level 1
  // The line on the right is already on level 1, but will be split as well
  // Only its child corresponding to to subregion 0 (SEE THE NAME OF THE
  // FUNCTION!) will appear in the incidence identities!

  // Child facet 0 on the left is incident to the child facet 0 on the right
  incidences[0][LEFT]  = incidences_subcells_L[0];
  incidences[0][RIGHT] = incidences_subcells_R[0];

  // Child facet 1 on the left is incident to the child facet 1 on the right
  incidences[1][LEFT]  = incidences_subcells_L[0];
  incidences[1][RIGHT] = incidences_subcells_R[1];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped

  if ((input_facet_perm_L.nb_flips() == 0) && (input_facet_perm_R.nb_flips() == 1))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 0, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 0, 0));

    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
  }

  // Second case: left facet is flipped
  // Right facet is identity

  else if ((input_facet_perm_L.nb_flips() == 1) && (input_facet_perm_R.nb_flips() == 0))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 1, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 1, 0));

    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "LineL_Level1Pos0_LineR_Level0::split_both_sides_two_"
                 "segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations "
                 "of incident\n"
              << " facets" << std::endl;
    return;
  }
}

// ----------------------------------------------------------------------------
// IIIb) Two incident line segments, left one is a subdomain of another segment
//       and incident to right facet that had already been refined
// ----------------------------------------------------------------------------

LineL_L0PartType1Pos1_LineR_L1::LineL_L0PartType1Pos1_LineR_L1() : TraceTupleBase()
{
}

LineL_L0PartType1Pos1_LineR_L1::~LineL_L0PartType1Pos1_LineR_L1()
{
}

void LineL_L0PartType1Pos1_LineR_L1::set_subdomain_tags(CellSubdomainTag &sub_tag_L,
                                                        CellSubdomainTag &sub_tag_R)
{
  sub_tag_L = CellSubdomainTag(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1, ElemShape::Line);
  sub_tag_R = CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line);
}

void LineL_L0PartType1Pos1_LineR_L1::set_incidence_comp_functions(
    std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
{
  functions.push_back(
      std::make_tuple(CellTransform::NO_TRANS, CellTransform::NO_TRANS, do_nothing));
  functions.push_back(std::make_tuple(CellTransform::UNIFORM_REFINE, CellTransform::NO_TRANS,
                                      split_left_side_two_segments));
  functions.push_back(std::make_tuple(CellTransform::UNIFORM_REFINE, CellTransform::UNIFORM_REFINE,
                                      split_both_sides_two_segments));
}

void LineL_L0PartType1Pos1_LineR_L1::fill_incident_p1_nodes(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    std::vector<std::pair<Uint, Uint>> &p1_verts) const
{
  const Uint nb_transforms_left =
      (input_facet_perm_L.nb_flips() + input_facet_perm_L.nb_rotations()) % 2;
  const Uint nb_transforms_right =
      (input_facet_perm_R.nb_flips() + input_facet_perm_R.nb_rotations()) % 2;

  p1_verts.resize(1);

  // Right face interfaces a sub-segment of the left face. Since this is
  // 'Partition type 1' and 'Position 1', it means that only the local vertex
  // 1 of the left (larger) face can be in contact with the right face. Hence
  // p1_verts[0].first MUST BE 1 in all cases. Now depending on whether both
  // systems are aligned or not, this vertex touches the local vertex 0 or
  // local vertex 1 of the right face.

  p1_verts[0].first = 1;

  if (nb_transforms_left == nb_transforms_right)
  {
    p1_verts[0].second = 1;
  }
  else
  {
    p1_verts[0].second = 0;
  }
}

void LineL_L0PartType1Pos1_LineR_L1::do_nothing(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(1);
  incidences[0].resize(2);
  incidences[0][LEFT]  = incidences_subcells_L[0];
  incidences[0][RIGHT] = incidences_subcells_R[0];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  facet_permutations.resize(1);
  facet_permutations[0].resize(2);
  facet_permutations[0][LEFT].change_type(line_std_reg_tag, input_facet_perm_L);
  facet_permutations[0][RIGHT].change_type(line_std_reg_tag, input_facet_perm_R);
}

void LineL_L0PartType1Pos1_LineR_L1::split_left_side_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(1);
  incidences[0].resize(2);

  facet_permutations.resize(1);
  facet_permutations[0].resize(2);

  assert(incidences_subcells_L.size() == 1);
  assert(incidences_subcells_R.size() == 2);

  // Child facet 0 on the left is incident to the child facet 0 on the right
  // Note that this is BY DEFINITION: line on the right is on level 0 and not
  // split - only one entity to match The line on the left is larger/coarser:
  // one level above, but we're using only a part of it - the right line is
  // incident to subregion 1 (SEE THE NAME OF THE FUNCTION!) of the left line
  incidences[0][LEFT]  = incidences_subcells_L[1];
  incidences[0][RIGHT] = incidences_subcells_R[0];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped

  if ((input_facet_perm_L.nb_flips() == 0) && (input_facet_perm_R.nb_flips() == 1))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
  }

  // Second case: left facet is flipped
  // Right facet is identity

  else if ((input_facet_perm_L.nb_flips() == 1) && (input_facet_perm_R.nb_flips() == 0))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "LineL_Level1Pos1_LineR_Level0::split_left_side_two_segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations "
                 "of incident\n"
              << " facets" << std::endl;
    return;
  }
}

void LineL_L0PartType1Pos1_LineR_L1::split_both_sides_two_segments(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    const std::vector<IncidenceEntry> &incidences_subcells_L,
    const std::vector<IncidenceEntry> &incidences_subcells_R,
    std::vector<std::vector<IncidenceEntry>> &incidences,
    std::vector<std::vector<EntityDofRealign>> &facet_permutations)
{
  incidences.resize(2);
  incidences[0].resize(2);
  incidences[1].resize(2);

  facet_permutations.resize(2);
  facet_permutations[0].resize(2);
  facet_permutations[1].resize(2);

  assert(incidences_subcells_L.size() == 2);
  assert(incidences_subcells_R.size() == 2);

  // Note that the following two identities are BY DEFINITION
  // line on left is on level 1 and split into two sub-entities
  // The line on the right is on level 0, but will be split as well
  // Only the child of left facet child corresponding to subregion 1
  // (SEE THE NAME OF THE FUNCTION!) will appear in the incidence identities!

  // Child facet 0 on the left is incident to the child facet 0 on the right
  incidences[0][LEFT]  = incidences_subcells_L[1];
  incidences[0][RIGHT] = incidences_subcells_R[0];

  // Child facet 1 on the left is incident to the child facet 1 on the right
  incidences[1][LEFT]  = incidences_subcells_L[1];
  incidences[1][RIGHT] = incidences_subcells_R[1];

  const PointSetTag line_std_reg_tag = PointSetTag(ElemShape::Line, P1, PointSetID::Equidist);

  // First case: entity permutation on the left is identity
  // Right facet is flipped

  if ((input_facet_perm_L.nb_flips() == 0) && (input_facet_perm_R.nb_flips() == 1))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 0, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 0, 0));

    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 1, 0));
  }

  // Second case: left facet is flipped
  // Right facet is identity

  else if ((input_facet_perm_L.nb_flips() == 1) && (input_facet_perm_R.nb_flips() == 0))
  {
    /// @remark MAKE SURE THAT THE NUMBER OF FLIPS AND ROTATIONS IN THE
    /// ENTITY
    ///         PERMUTATION CODE IS __PRECISELY__ THE SAME AS IN THE
    ///         PARENT!!! Child entities 'inherit' the coordinate system of
    ///         parent - to align with incident entities on the 'other side'
    ///         of the trace, they have to be permuted in the same manner as
    ///         parent

    // First row: first pair of incidences
    facet_permutations[0][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1,
                                            ElemShape::Line, 1, 0));

    facet_permutations[0][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));

    // Second row: second pair of incidences
    facet_permutations[1][LEFT].change_type(
        line_std_reg_tag, EntityRealignCode(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0,
                                            ElemShape::Line, 1, 0));

    facet_permutations[1][RIGHT].change_type(
        line_std_reg_tag,
        EntityRealignCode(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line, 0, 0));
  }
  else
  {
    std::cerr << "LineL_Level1Pos1_LineR_Level0::split_both_sides_two_"
                 "segments::\n"
              << " don't know how to match entities due to incompatible "
                 "permutations "
                 "of incident\n"
              << " facets" << std::endl;
    return;
  }
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace mesh

} // namespace pdekit
