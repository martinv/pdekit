#include <iostream>

#include "mesh/local_topology/CellSubdomainTag.hpp"
#include "mesh/local_topology/TraceEntityTuple.hpp"
#include "mesh/local_topology/TraceTupleFactory.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

const std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint>
    TraceEntityTupleInstance::undefined =
        std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint>(
            EntityRealignCode(ElemShape::Undefined, CellTransform::NO_TRANS, 0,
                              ElemShape::Undefined, 0, 0),
            0,
            EntityRealignCode(ElemShape::Undefined, CellTransform::NO_TRANS, 0,
                              ElemShape::Undefined, 0, 0),
            0);

// ----------------------------------------------------------------------------

void TraceEntityTupleInstance::construct(
    const std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint> &key,
    TraceEntityTupleInstance &inc_pattern)
{
  TraceTupleFactory::instance_type &inc_factory = TraceTupleFactory::instance();

  const CellSubdomainTag sub_tag_L(std::get<0>(key).elem_shape(), std::get<0>(key).adapt_op_id(),
                                   std::get<0>(key).local_pos_in_parent(),
                                   std::get<0>(key).parent_shape());
  const CellSubdomainTag sub_tag_R(std::get<2>(key).elem_shape(), std::get<2>(key).adapt_op_id(),
                                   std::get<2>(key).local_pos_in_parent(),
                                   std::get<2>(key).parent_shape());

  const std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint> cell_subdomain_tag_key =
      std::make_tuple(sub_tag_L, std::get<1>(key), sub_tag_R, std::get<3>(key));

  const TraceTupleFactory::instance_type::const_product_base_ptr facet_inc_base =
      inc_factory.create(cell_subdomain_tag_key);

  /*
  std::cout << "CellSubdomainTupleInstance:: Subdomain tag key L: " <<
  key.first.as_string() << std::endl; std::cout <<
  "CellSubdomainTupleInstance:: Subdomain tag key R: " <<
  key.second.as_string() << std::endl;
  */
  /*
  facet_inc_base->set_subdomain_tags(inc_pattern.m_subdomain_tags.first,
                                     inc_pattern.m_subdomain_tags.second);
  */

  inc_pattern.m_facet_perm_tags = key;

  /*
  std::cout << "CellSubdomainTupleInstance:: Setting subdomain tags " <<
  inc_pattern.m_subdomain_tags.first.as_string() << " "
            << " and " << inc_pattern.m_subdomain_tags.second.as_string() <<
  std::endl;
  */

  facet_inc_base->set_incidence_comp_functions(inc_pattern.m_incidences_update_functions);

  facet_inc_base->fill_incident_p1_nodes(std::get<0>(key), std::get<2>(key),
                                         inc_pattern.m_incident_p1_nodes);
}

// ----------------------------------------------------------------------------

TraceEntityTupleInstance::TraceEntityTupleInstance()
{
  m_facet_perm_tags = undefined;
}

// ----------------------------------------------------------------------------

TraceEntityTupleInstance::TraceEntityTupleInstance(
    std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint> const &key)
{
  construct(key, *this);
}

// ----------------------------------------------------------------------------

TraceEntityTupleInstance::TraceEntityTupleInstance(
    const TraceEntityTupleInstance &other_inc_pattern)
{
  m_facet_perm_tags = other_inc_pattern.m_facet_perm_tags;
}

// ----------------------------------------------------------------------------

TraceEntityTupleInstance &TraceEntityTupleInstance::operator=(
    const TraceEntityTupleInstance &other_inc_pattern)
{
  m_facet_perm_tags = other_inc_pattern.m_facet_perm_tags;

  return *this;
}

// ----------------------------------------------------------------------------

TraceEntityTupleInstance::~TraceEntityTupleInstance()
{
}

// ----------------------------------------------------------------------------

void TraceEntityTupleInstance::combine_incidences_on_adjacent_facets(
    const CellTransform input_facet_transform_L, const CellTransform input_facet_transform_R,
    const std::vector<IncidenceEntry> &child_incidences_L,
    const std::vector<IncidenceEntry> &child_incidences_R,
    std::vector<std::vector<IncidenceEntry>> &combined_child_incidences,
    std::vector<std::vector<EntityDofRealign>> &combined_child_permutations) const
{
  for (Uint i = 0; i < m_incidences_update_functions.size(); ++i)
  {
    if ((std::get<0>((m_incidences_update_functions[i])) == input_facet_transform_L) &&
        (std::get<1>((m_incidences_update_functions[i])) == input_facet_transform_R))
    {
      std::get<2>(m_incidences_update_functions[i])(
          std::get<0>(m_facet_perm_tags), std::get<2>(m_facet_perm_tags), child_incidences_L,
          child_incidences_R, combined_child_incidences, combined_child_permutations);

      // For debugging:
      // Force the generated permutations to have the same number of
      // rotations and flips as the parent entities
      const Uint nb_flips_L     = std::get<0>(m_facet_perm_tags).nb_flips();
      const Uint nb_rotations_L = std::get<0>(m_facet_perm_tags).nb_rotations();

      const Uint nb_flips_R     = std::get<2>(m_facet_perm_tags).nb_flips();
      const Uint nb_rotations_R = std::get<2>(m_facet_perm_tags).nb_rotations();

      for (Uint p = 0; p < combined_child_permutations.size(); ++p)
      {
        EntityRealignCode pcode_L = combined_child_permutations[p][LEFT].get().code();
        EntityRealignCode pcode_R = combined_child_permutations[p][RIGHT].get().code();

        pcode_L.set_nb_flips(nb_flips_L);
        pcode_L.set_nb_rotations(nb_rotations_L);
        pcode_R.set_nb_flips(nb_flips_R);
        pcode_R.set_nb_rotations(nb_rotations_R);

        const PointSetTag tag_L = combined_child_permutations[p][LEFT].get().type_id();
        const PointSetTag tag_R = combined_child_permutations[p][RIGHT].get().type_id();

        combined_child_permutations[p][LEFT].change_type(tag_L, pcode_L);
        combined_child_permutations[p][RIGHT].change_type(tag_R, pcode_R);
      }

      return;
    }
  }
  std::cerr << "TraceEntityTuple::update_incidences_on_adjacent_facets: didn't "
               "find "
               "suitable\n"
            << "operation to perform for pair \n  " << std::get<0>(m_facet_perm_tags).as_string()
            << "," << std::get<2>(m_facet_perm_tags).as_string() << "] with ops {"
            << CellTransformName::value[static_cast<Uint>(input_facet_transform_L)] << " - "
            << CellTransformName::value[static_cast<Uint>(input_facet_transform_R)] << "}"
            << std::endl;

  std::cerr << "  TraceEntityTuple::update_incidences_on_adjacent_facets: have "
               "the following ops:\n";
  for (Uint i = 0; i < m_incidences_update_functions.size(); ++i)
  {
    std::cout << "    " << i << ") "
              << CellTransformName::value[static_cast<Uint>(
                     std::get<0>((m_incidences_update_functions[i])))]
              << " - "
              << CellTransformName::value[static_cast<Uint>(
                     std::get<1>((m_incidences_update_functions[i])))]
              << std::endl;
  }
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const TraceEntityTupleInstance &inc_pattern)
{
  /*
  os << "Local incidence pattern: {" <<
  ElemShape::Names[inc_pattern.m_shape_L]
  << "}{"
     << ElemShape::Names[inc_pattern.m_shape_L] << "}" << std::endl;
  */
  os << "(" << std::get<0>(inc_pattern.m_facet_perm_tags).as_string() << ","
     << std::get<2>(inc_pattern.m_facet_perm_tags).as_string() << ")"; // << std::endl;
  return os;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os,
                         const std::tuple<EntityRealignCode, EntityRealignCode> &tuple_key)
{
  os << std::get<0>(tuple_key).as_string() << "-" << std::get<1>(tuple_key).as_string();
  return os;
}

// ----------------------------------------------------------------------------

std::tuple<Uint, Uint> get_relative_refinement_levels(const Uint level_L, const Uint level_R)
{
  const Uint min_level = std::min(level_L, level_R);
  return std::tuple<Uint, Uint>(level_L - min_level, level_R - min_level);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
