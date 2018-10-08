#include "mesh/local_topology/TraceTupleQuad.hpp"
#include "mesh/local_topology/CellSubdomainTag.hpp"

namespace pdekit
{

namespace mesh
{

namespace internal
{

// ----------------------------------------------------------------------------
//                   FacetIncidenceQuadVariant0
// ----------------------------------------------------------------------------

TraceTupleQuadVariant0::TraceTupleQuadVariant0() : TraceTupleBase()
{
}

TraceTupleQuadVariant0::~TraceTupleQuadVariant0()
{
}

void TraceTupleQuadVariant0::set_subdomain_tags(CellSubdomainTag &sub_tag_left,
                                                CellSubdomainTag &sub_tag_right)
{
  // adapt_op_left = NO_SPLIT;
  // adapt_op_right = QUAD_TO_FOUR_SUBQUADS;
}

void TraceTupleQuadVariant0::set_incidence_comp_functions(
    std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
{
}

void TraceTupleQuadVariant0::fill_incident_p1_nodes(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    std::vector<std::pair<Uint, Uint>> &p1_verts) const
{
}

// ----------------------------------------------------------------------------
//                   FacetIncidenceQuadVariant1
// ----------------------------------------------------------------------------

TraceTupleQuadVariant1::TraceTupleQuadVariant1() : TraceTupleBase()
{
}

TraceTupleQuadVariant1::~TraceTupleQuadVariant1()
{
}

void TraceTupleQuadVariant1::set_subdomain_tags(CellSubdomainTag &sub_tag_left,
                                                CellSubdomainTag &sub_tag_right)
{
  // adapt_op_left = QUAD_TO_FOUR_SUBQUADS;
  // adapt_op_right = NO_SPLIT;
}

void TraceTupleQuadVariant1::set_incidence_comp_functions(
    std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
{
}

void TraceTupleQuadVariant1::fill_incident_p1_nodes(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    std::vector<std::pair<Uint, Uint>> &p1_verts) const
{
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace mesh

} // namespace pdekit
