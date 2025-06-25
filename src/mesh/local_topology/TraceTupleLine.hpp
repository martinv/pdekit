#ifndef PDEKIT_Mesh_Trace_Tuple_Line_hpp
#define PDEKIT_Mesh_Trace_Tuple_Line_hpp

#include <vector>

#include "mesh/local_topology/TraceTupleBase.hpp"

namespace pdekit
{

namespace mesh
{

namespace internal
{

// ----------------------------------------------------------------------------

// Basic version of line facet incidences:
// 1 line segment on the left
// 1 line segment on the right
// The left facet is not partitioned in order to match the right facet
// The right facet is not partitioned in order to match the left facet

class LineL_L0_LineR_L0 : public TraceTupleBase
{
  public:
  /// Default constructor
  LineL_L0_LineR_L0();

  /// Default destructor
  ~LineL_L0_LineR_L0() override;

  /// Set the facet incidence tag
  void set_subdomain_tags(CellSubdomainTag &sub_tag_L, CellSubdomainTag &sub_tag_R) override;

  /// Set incidence functions
  void set_incidence_comp_functions(
      std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
      override;

  /// Set function that is able to find incident P1 nodes on nonconforming
  /// facets
  void fill_incident_p1_nodes(const EntityRealignCode input_facet_perm_L,
                              const EntityRealignCode input_facet_perm_R,
                              std::vector<std::pair<Uint, Uint>> &p1_verts) const override;

  private:
  /// Functions to compute incidences
  static void do_nothing(const EntityRealignCode input_facet_perm_L,
                         const EntityRealignCode input_facet_perm_R,
                         const std::vector<IncidenceEntry> &incidences_subcells_L,
                         const std::vector<IncidenceEntry> &incidences_subcells_R,
                         std::vector<std::vector<IncidenceEntry>> &incidences,
                         std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_left_side_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_right_side_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_both_sides_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);
};

// ----------------------------------------------------------------------------

// IIa) Incidences when left line is on lower refinement level than
//      the right line. This means that the left line is smaller
//      than the right line and interfaces only a part of the right line.
//      In this case, we consider partition type 1 (the right line is
//      'virtually' halved into two sub-segments) and it is the 'first' half
//      (measured from origin of the local coordinate system of the right line),
//      which interfaces the left line. This is represented by 'Pos0' (as
//      'position 0) in name of the class.
//
//      In case of refinement, either only the right line is refined, or both
//      left and right lines are refined. Any other refinement operation
//      would mean that the refinement levels of newly created entities
//      on the left and right sides would differ by more than 1, which is
//      currently not supported.

class LineL_L1_LineR_L0PartType1Pos0 : public TraceTupleBase
{
  public:
  /// Default constructor
  LineL_L1_LineR_L0PartType1Pos0();

  /// Default destructor
  ~LineL_L1_LineR_L0PartType1Pos0() override;

  /// Set the facet incidence tag
  void set_subdomain_tags(CellSubdomainTag &sub_tag_L, CellSubdomainTag &sub_tag_R) override;

  /// Set incidence functions
  void set_incidence_comp_functions(
      std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
      override;

  /// Find incident nodes of two adjacent nonconforming facets
  void fill_incident_p1_nodes(const EntityRealignCode input_facet_perm_L,
                              const EntityRealignCode input_facet_perm_R,
                              std::vector<std::pair<Uint, Uint>> &p1_verts) const override;

  private:
  /// Functions to compute incidences
  static void do_nothing(const EntityRealignCode input_facet_perm_L,
                         const EntityRealignCode input_facet_perm_R,
                         const std::vector<IncidenceEntry> &incidences_subcells_L,
                         const std::vector<IncidenceEntry> &incidences_subcells_R,
                         std::vector<std::vector<IncidenceEntry>> &incidences,
                         std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_right_side_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_both_sides_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);
};

// ----------------------------------------------------------------------------

// IIb) Incidences when left line is on lower refinement level than
//      the right line. This means that the left line is smaller
//      than the right line and interfaces only a part of the right line.
//      In this case, we consider partition type 1 (the right line is
//      'virtually' halved into two sub-segments) and it is the 'second' half
//      (measured from origin of the local coordinate system of the right line),
//      which interfaces the left line. This is represented by 'Pos1' (as
//      'position 1) in name of the class.
//
//      In case of refinement, either only the right line is refined, or both
//      left and right lines are refined. Any other refinement operation
//      would mean that the refinement levels of newly created entities
//      on the left and right sides would differ by more than 1, which is
//      currently not supported.

class LineL_L1_LineR_L0PartType1Pos1 : public TraceTupleBase
{
  public:
  /// Default constructor
  LineL_L1_LineR_L0PartType1Pos1();

  /// Default destructor
  ~LineL_L1_LineR_L0PartType1Pos1() override;

  /// Set the facet incidence tag
  void set_subdomain_tags(CellSubdomainTag &sub_tag_L, CellSubdomainTag &sub_tag_R) override;

  /// Set incidence functions
  void set_incidence_comp_functions(
      std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
      override;

  /// Set function that is able to find incident P1 nodes on nonconforming
  /// facets
  void fill_incident_p1_nodes(const EntityRealignCode input_facet_perm_L,
                              const EntityRealignCode input_facet_perm_R,
                              std::vector<std::pair<Uint, Uint>> &p1_verts) const override;

  private:
  /// Functions to compute incidences
  static void do_nothing(const EntityRealignCode input_facet_perm_L,
                         const EntityRealignCode input_facet_perm_R,
                         const std::vector<IncidenceEntry> &incidences_subcells_L,
                         const std::vector<IncidenceEntry> &incidences_subcells_R,
                         std::vector<std::vector<IncidenceEntry>> &incidences,
                         std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_right_side_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_both_sides_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);
};

// ----------------------------------------------------------------------------

// IIIa) Incidences when left right is on lower refinement level than
//       the left line. This means that the right line is smaller
//       than the left line and interfaces only a part of the left line.
//       In this case, we consider partition type 1 (the left line is
//       'virtually' halved into two sub-segments) and it is the 'first' half
//       (measured from origin of the local coordinate system of the left line),
//       which interfaces the right line. This is represented by 'Pos0' (as
//       'position 0) in name of the class.
//
//       In case of refinement, either only the left line is refined, or both
//       left and right lines are refined. Any other refinement operation
//       would mean that the refinement levels of newly created entities
//       on the left and right sides would differ by more than 1, which is
//       currently not supported.

class LineL_L0PartType1Pos0_LineR_L1 : public TraceTupleBase
{
  public:
  /// Default constructor
  LineL_L0PartType1Pos0_LineR_L1();

  /// Default destructor
  ~LineL_L0PartType1Pos0_LineR_L1() override;

  /// Set the facet incidence tag
  void set_subdomain_tags(CellSubdomainTag &sub_tag_L, CellSubdomainTag &sub_tag_R) override;

  /// Set incidence functions
  void set_incidence_comp_functions(
      std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
      override;

  /// Set function that is able to find incident P1 nodes on nonconforming
  /// facets
  void fill_incident_p1_nodes(const EntityRealignCode input_facet_perm_L,
                              const EntityRealignCode input_facet_perm_R,
                              std::vector<std::pair<Uint, Uint>> &p1_verts) const override;

  private:
  /// Functions to compute incidences
  static void do_nothing(const EntityRealignCode input_facet_perm_L,
                         const EntityRealignCode input_facet_perm_R,
                         const std::vector<IncidenceEntry> &incidences_subcells_L,
                         const std::vector<IncidenceEntry> &incidences_subcells_R,
                         std::vector<std::vector<IncidenceEntry>> &incidences,
                         std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_left_side_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_both_sides_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);
};

// ----------------------------------------------------------------------------

// IIIb) Incidences when left right is on lower refinement level than
//       the left line. This means that the right line is smaller
//       than the left line and interfaces only a part of the left line.
//       In this case, we consider partition type 1 (the left line is
//       'virtually' halved into two sub-segments) and it is the 'second' half
//       (measured from origin of the local coordinate system of the left line),
//       which interfaces the right line. This is represented by 'Pos1' (as
//       'position 1) in name of the class.
//
//       In case of refinement, either only the left line is refined, or both
//       left and right lines are refined. Any other refinement operation
//       would mean that the refinement levels of newly created entities
//       on the left and right sides would differ by more than 1, which is
//       currently not supported.

class LineL_L0PartType1Pos1_LineR_L1 : public TraceTupleBase
{
  public:
  /// Default constructor
  LineL_L0PartType1Pos1_LineR_L1();

  /// Default destructor
  ~LineL_L0PartType1Pos1_LineR_L1() override;

  /// Set the facet incidence tag
  void set_subdomain_tags(CellSubdomainTag &sub_tag_L, CellSubdomainTag &sub_tag_R) override;

  /// Set incidence functions
  void set_incidence_comp_functions(
      std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
      override;

  /// Set function that is able to find incident P1 nodes on nonconforming
  /// facets
  void fill_incident_p1_nodes(const EntityRealignCode input_facet_perm_L,
                              const EntityRealignCode input_facet_perm_R,
                              std::vector<std::pair<Uint, Uint>> &p1_verts) const override;

  private:
  /// Functions to compute incidences
  static void do_nothing(const EntityRealignCode input_facet_perm_L,
                         const EntityRealignCode input_facet_perm_R,
                         const std::vector<IncidenceEntry> &incidences_subcells_L,
                         const std::vector<IncidenceEntry> &incidences_subcells_R,
                         std::vector<std::vector<IncidenceEntry>> &incidences,
                         std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_left_side_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);

  static void split_both_sides_two_segments(
      const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
      const std::vector<IncidenceEntry> &incidences_subcells_L,
      const std::vector<IncidenceEntry> &incidences_subcells_R,
      std::vector<std::vector<IncidenceEntry>> &incidences,
      std::vector<std::vector<EntityDofRealign>> &facet_permutations);
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace mesh

} // namespace pdekit

#endif
