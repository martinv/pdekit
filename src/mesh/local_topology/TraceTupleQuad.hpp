#ifndef PDEKIT_Mesh_Trace_Tuple_Quad_hpp
#define PDEKIT_Mesh_Trace_Tuple_Quad_hpp

#include <vector>

#include "mesh/local_topology/TraceTupleBase.hpp"

namespace pdekit
{

namespace mesh
{

namespace internal
{

// ----------------------------------------------------------------------------

// This is the basic one-sided version of triangle facet incidences:
// Left: a quadrilateral facet, which is not split
// Right: a quadrilateral facet which is split into four subtriangles

class TraceTupleQuadVariant0 : public TraceTupleBase
{
  public:
  /// Default constructor
  TraceTupleQuadVariant0();

  /// Default destructor
  ~TraceTupleQuadVariant0() override;

  /// Set the facet incidence tag
  void set_subdomain_tags(CellSubdomainTag &sub_tag_left, CellSubdomainTag &sub_tag_right) override;

  /// Set incidence functions
  void set_incidence_comp_functions(
      std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
      override;

  /// Find incident nodes of two adjacent nonconforming facets
  void fill_incident_p1_nodes(const EntityRealignCode input_facet_perm_L,
                              const EntityRealignCode input_facet_perm_R,
                              std::vector<std::pair<Uint, Uint>> &p1_verts) const override;
};

// This is the basic one-sided version of triangle facet incidences:
// Left: a quadrilateral facet which is split into four subtriangles
// Right: a quadrilateral facet, which is not split

class TraceTupleQuadVariant1 : public TraceTupleBase
{
  public:
  /// Default constructor
  TraceTupleQuadVariant1();

  /// Default destructor
  ~TraceTupleQuadVariant1() override;

  /// Set the facet incidence tag
  void set_subdomain_tags(CellSubdomainTag &sub_tag_left, CellSubdomainTag &sub_tag_right) override;

  /// Set incidence functions
  void set_incidence_comp_functions(
      std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
      override;

  /// Find incident nodes of two adjacent nonconforming facets
  void fill_incident_p1_nodes(const EntityRealignCode input_facet_perm_L,
                              const EntityRealignCode input_facet_perm_R,
                              std::vector<std::pair<Uint, Uint>> &p1_verts) const override;
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace mesh

} // namespace pdekit

#endif
