#ifndef PDEKIT_Mesh_Adaptation_Do_Nothing_Cell_Adapt_Op_hpp
#define PDEKIT_Mesh_Adaptation_Do_Nothing_Cell_Adapt_Op_hpp

#include <vector>

#include "mesh/adaptation/CellAdaptOpBase.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

namespace internal
{

// ----------------------------------------------------------------------------

// This is the basic 'isotropic' version of triangle split:
// the triangle is split into four triangles of equal shape and size

class DoNothingCellAdaptOp : public CellAdaptOpBase
{
  public:
  /// Default constructor
  DoNothingCellAdaptOp();

  /// Default destructor
  ~DoNothingCellAdaptOp() override;

  /// Set the facet incidence tag
  void set_facet_adapt_op_ids(CellAdaptOpTag &tag,
                              std::vector<CellTransform> &facet_adapt_op_ids) override;

  /// Number of facets of the parent
  Uint nb_parent_facets() const override;

  /// Set number of child elements after adaptation
  Uint nb_child_elems() const override;

  /// Fill the matrix of reference coordinates (of principal vertices)
  void fill_parent_ref_coords(math::DenseDMat<Real> &ref_coords) const override;

  /// Fill vector of child element shapes (i.e. after the parent element is
  /// split)
  void set_child_elem_shapes(std::vector<ElemShape> &child_elem_shapes) const override;

  /// Fill vector of child element shapes (i.e. after the parent element is
  /// split)
  void fill_local_child_incidences(std::vector<IncidenceEntry> &incidences,
                                   std::vector<EntityDofRealign> &permutations,
                                   std::vector<Uint> &local_offsets) const override;

  /// Fill a table of incidences of child facets that cover given facet of
  /// parent element Fill a table of permutations of child facets that cover
  /// given facet of parent element
  void fill_parent_child_incidences_on_facet(
      const Uint facet, std::vector<IncidenceEntry> &incidences,
      std::vector<EntityDofRealign> &permutations) const override;

  /// Fill a vector of funtions which know how to transform coordinates of
  /// parent element into the coordinates of sub-element Each function is of
  /// the type std::function<void(math::ConstVectorBlock<Real> const&,
  ///                                                 math::VectorBlock<Real>&)
  void fill_child_coord_transformers(
      std::vector<std::function<void(
          math::DenseDMat<Real> const &ref_coords, math::DenseConstVecView<Real> const &coord_in,
          math::DenseVecView<Real> &coord_out)>> &transformers) const override;
};

// ----------------------------------------------------------------------------
//                          ELEMENT SHAPE: LINE
// ----------------------------------------------------------------------------

class DoNothingCellAdaptOpLine : public DoNothingCellAdaptOp
{
  public:
  /// Set the facet incidence tag
  void set_tag(CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids);
};

// ----------------------------------------------------------------------------
//                        ELEMENT SHAPE: TRIANGLE
// ----------------------------------------------------------------------------

class DoNothingCellAdaptOpTriag : public DoNothingCellAdaptOp
{
  public:
  /// Set the facet incidence tag
  void set_tag(CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids);
};

// ----------------------------------------------------------------------------
//                      ELEMENT SHAPE: QUADRILATERAL
// ----------------------------------------------------------------------------

class DoNothingCellAdaptOpQuad : public DoNothingCellAdaptOp
{
  public:
  /// Set the facet incidence tag
  void set_tag(CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids);
};

// ----------------------------------------------------------------------------
//                       ELEMENT SHAPE: TETRAHEDRON
// ----------------------------------------------------------------------------

class DoNothingCellAdaptOpTetra : public DoNothingCellAdaptOp
{
  public:
  /// Set the facet incidence tag
  void set_tag(CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids);
};

// ----------------------------------------------------------------------------
//                       ELEMENT SHAPE: HEXAHEDRON
// ----------------------------------------------------------------------------

class DoNothingCellAdaptOpHexa : public DoNothingCellAdaptOp
{
  public:
  /// Set the facet incidence tag
  void set_tag(CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids);
};

// ----------------------------------------------------------------------------
//                        ELEMENT SHAPE: PRISM
// ----------------------------------------------------------------------------

class DoNothingCellAdaptOpPrism : public DoNothingCellAdaptOp
{
  public:
  /// Set the facet incidence tag
  void set_tag(CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids);
};

// ----------------------------------------------------------------------------
//                        ELEMENT SHAPE: PYRAMID
// ----------------------------------------------------------------------------

class DoNothingCellAdaptOpPyramid : public DoNothingCellAdaptOp
{
  public:
  /// Set the facet incidence tag
  void set_tag(CellAdaptOpTag &tag, std::vector<CellTransform> &facet_adapt_op_ids);
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace adapt

} // namespace mesh

} // namespace pdekit

#endif
