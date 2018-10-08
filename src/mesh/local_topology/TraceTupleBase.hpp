#ifndef PDEKIT_Mesh_Trace_Tuple_Base_hpp
#define PDEKIT_Mesh_Trace_Tuple_Base_hpp

#include <functional>
#include <vector>

#include "common/Constants.hpp"
#include "mesh/EntityRealignCode.hpp"
#include "mesh/local_topology/CellSubdomainTag.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"

namespace pdekit
{

namespace mesh
{

namespace internal
{

// ----------------------------------------------------------------------------

class TraceTupleBase
{
  public:
  /// TYPEDEFS

  typedef std::function<void(const EntityRealignCode input_facet_perm_L,
                             const EntityRealignCode input_facet_perm_R,
                             const std::vector<IncidenceEntry> &incidences_subcells_left,
                             const std::vector<IncidenceEntry> &incidences_subcells_right,
                             std::vector<std::vector<IncidenceEntry>> &incidences,
                             std::vector<std::vector<EntityDofRealign>> &facet_permutations)>
      update_incidences_fct;

  /// Default constructor
  TraceTupleBase();

  /// Default destructor
  virtual ~TraceTupleBase();

  /// type name of this base class
  static std::string type_name()
  {
    return "TraceTupleBase";
  }

  /// Set the facet incidence tag
  virtual void set_subdomain_tags(CellSubdomainTag &sub_tag_L, CellSubdomainTag &sub_tag_R) = 0;

  /// Set incidence functions
  virtual void set_incidence_comp_functions(
      std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>> &functions) = 0;

  /// Set function that is able to find incident P1 nodes on nonconforming
  /// facets
  virtual void fill_incident_p1_nodes(const EntityRealignCode input_facet_perm_L,
                                      const EntityRealignCode input_facet_perm_R,
                                      std::vector<std::pair<Uint, Uint>> &p1_verts) const = 0;
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace mesh

} // namespace pdekit

#endif
