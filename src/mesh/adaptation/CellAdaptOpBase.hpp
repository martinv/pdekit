#ifndef PDEKIT_Mesh_Adapt_Cell_Adapt_Op_Base_hpp
#define PDEKIT_Mesh_Adapt_Cell_Adapt_Op_Base_hpp

#include <functional>
#include <vector>

#include "common/Constants.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseVecView.hpp"
#include "mesh/adaptation/CellAdaptOpTag.hpp"
#include "mesh/local_topology/CellSubdomainTag.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"
#include "mesh/std_region/StdRegion.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

namespace internal
{

// ----------------------------------------------------------------------------

class CellAdaptOpBase
{
  public:
  /// Default constructor
  CellAdaptOpBase();

  /// Default destructor
  virtual ~CellAdaptOpBase();

  /// type name of this base class
  static std::string type_name()
  {
    return "CellSplitStrategyBase";
  }

  /// Set the facet incidence tag
  virtual void set_facet_adapt_op_ids(CellAdaptOpTag &tag,
                                      std::vector<CellTransform> &facet_adapt_op_ids) = 0;

  /// Number of facets of the parent
  virtual Uint nb_parent_facets() const = 0;

  /// Set number of child elements after adaptation
  virtual Uint nb_child_elems() const = 0;

  /// Fill the matrix of reference coordinates (of principal vertices)
  virtual void fill_parent_ref_coords(math::DenseDMat<Real> &ref_coords) const = 0;

  /// Fill vector of child element shapes (i.e. after the parent element is
  /// split)
  virtual void set_child_elem_shapes(std::vector<ElemShape> &child_elem_shapes) const = 0;

  /// Fill vector of child element shapes (i.e. after the parent element is
  /// split)
  virtual void fill_local_child_incidences(std::vector<IncidenceEntry> &incidences,
                                           std::vector<EntityDofRealign> &permutations,
                                           std::vector<Uint> &local_offsets) const = 0;

  /// Fill a table of incidences of child facets that cover given facet of
  /// parent element Fill a table of permutations of child facets that cover
  /// given facet of parent element
  virtual void fill_parent_child_incidences_on_facet(
      const Uint facet, std::vector<IncidenceEntry> &incidences,
      std::vector<EntityDofRealign> &permutations) const = 0;

  /// Fill a vector of funtions which know how to transform coordinates of
  /// parent element into the coordinates of sub-element Each function is of
  /// the type std::function<void(math::ConstVectorBlock<Real> const&,
  ///                                                 math::VectorBlock<Real>&)
  virtual void fill_child_coord_transformers(
      std::vector<std::function<void(
          math::DenseDMat<Real> const &ref_coords, math::DenseConstVecView<Real> const &coord_in,
          math::DenseVecView<Real> &coord_out)>> &transformers) const = 0;
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace adapt

} // namespace mesh

} // namespace pdekit

#endif
