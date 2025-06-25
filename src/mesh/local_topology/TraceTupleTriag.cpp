#include "mesh/local_topology/TraceTupleTriag.hpp"
#include "mesh/CellTransform.hpp"
#include "mesh/local_topology/CellSubdomainTag.hpp"

namespace pdekit
{

namespace mesh
{

namespace internal
{

// ----------------------------------------------------------------------------
//                   TraceTupleTriagVariant0
// ----------------------------------------------------------------------------

TraceTupleTriagVariant0::TraceTupleTriagVariant0() : TraceTupleBase()
{
}

TraceTupleTriagVariant0::~TraceTupleTriagVariant0()
{
}

void TraceTupleTriagVariant0::set_subdomain_tags(CellSubdomainTag &sub_tag_left,
                                                 CellSubdomainTag &sub_tag_right)
{
  // adapt_op_left = NO_SPLIT;
  // adapt_op_right = UNIFORM_REFINE;
}

void TraceTupleTriagVariant0::set_incidence_comp_functions(
    std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
{
}

void TraceTupleTriagVariant0::fill_incident_p1_nodes(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    std::vector<std::pair<Uint, Uint>> &p1_verts) const
{
}

// ----------------------------------------------------------------------------
//                   TraceTupleTriagVariant1
// ----------------------------------------------------------------------------

TraceTupleTriagVariant1::TraceTupleTriagVariant1() : TraceTupleBase()
{
}

TraceTupleTriagVariant1::~TraceTupleTriagVariant1()
{
}

void TraceTupleTriagVariant1::set_subdomain_tags(CellSubdomainTag &sub_tag_left,
                                                 CellSubdomainTag &sub_tag_right)
{
  // adapt_op_left = UNIFORM_REFINE;
  // adapt_op_right = NO_SPLIT;
}

void TraceTupleTriagVariant1::set_incidence_comp_functions(
    std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions)
{
}

void TraceTupleTriagVariant1::fill_incident_p1_nodes(
    const EntityRealignCode input_facet_perm_L, const EntityRealignCode input_facet_perm_R,
    std::vector<std::pair<Uint, Uint>> &p1_verts) const
{
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace mesh

} // namespace pdekit
