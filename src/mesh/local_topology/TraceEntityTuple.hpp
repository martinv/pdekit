#ifndef PDEKIT_Mesh_Local_Topology_Trace_Entity_Tuple_hpp
#define PDEKIT_Mesh_Local_Topology_Trace_Entity_Tuple_hpp

#include <array>
#include <functional>
#include <tuple>
#include <vector>

#include "common/Flyweight.hpp"
#include "mesh/CellTransform.hpp"
#include "mesh/EntityRealignCode.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"

namespace pdekit
{

namespace mesh
{

class TraceEntityTupleInstance
{
  public:
  /// STATIC VARIABLES

  static const std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint> undefined;

  /// Static method to fill a permutation
  static void construct(const std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint> &key,
                        TraceEntityTupleInstance &inc_pattern);

  /// Default constructor
  TraceEntityTupleInstance();

  /// Constructor permutation codes and refinement levels
  TraceEntityTupleInstance(std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint> const &key);

  /// Copy constructor
  TraceEntityTupleInstance(const TraceEntityTupleInstance &other_inc_pattern);

  /// Assignment operator
  TraceEntityTupleInstance &operator=(const TraceEntityTupleInstance &other_inc_pattern);

  /// Destructor
  ~TraceEntityTupleInstance();

  /// Return the interpolation point set type that
  /// this reference element represents
  const std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint> subdomain_tags() const;

  /// Return the pairs holding indices of p1 nodes that coincide
  const std::vector<std::pair<Uint, Uint>> &incident_p1_nodes() const;

  void combine_incidences_on_adjacent_facets(
      const CellTransform input_facet_transform_L, const CellTransform input_facet_transform_R,
      const std::vector<IncidenceEntry> &child_incidences_L,
      const std::vector<IncidenceEntry> &child_incidences_R,
      std::vector<std::vector<IncidenceEntry>> &combined_child_incidences,
      std::vector<std::vector<EntityDofRealign>> &combined_child_permutations) const;

  /// Print the internal data
  friend std::ostream &operator<<(std::ostream &os, const TraceEntityTupleInstance &inc_pattern);

  private:
  /// Type of reference topology of the L and R incident
  /// facet. Each tag knows
  /// - the element shape
  /// - adaptation operation that was applied to parent facet
  /// - position in parent facet
  /// - shape of parent facet
  std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint> m_facet_perm_tags;

  typedef std::function<void(const EntityRealignCode input_facet_perm_L,
                             const EntityRealignCode input_facet_perm_R,
                             const std::vector<IncidenceEntry> &incidences_subcells_L,
                             const std::vector<IncidenceEntry> &incidences_subcells_R,
                             std::vector<std::vector<IncidenceEntry>> &incidences,
                             std::vector<std::vector<EntityDofRealign>> &facet_permutations)>
      update_incidences_fct;

  typedef std::function<void(const EntityRealignCode input_facet_perm_L,
                             const EntityRealignCode input_facet_perm_R,
                             std::vector<std::pair<Uint, Uint>> &incident_p1_nodes)>
      detect_incident_p1_nodes_fct;

  std::vector<std::tuple<CellTransform, CellTransform, update_incidences_fct>>
      m_incidences_update_functions;

  std::vector<std::pair<Uint, Uint>> m_incident_p1_nodes;
};

// ----------------------------------------------------------------------------

inline const std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint> TraceEntityTupleInstance::
    subdomain_tags() const
{
  return m_facet_perm_tags;
}

// ----------------------------------------------------------------------------

inline const std::vector<std::pair<Uint, Uint>> &TraceEntityTupleInstance::incident_p1_nodes() const
{
  return m_incident_p1_nodes;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(
    std::ostream &os,
    const std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint> &tuple_key);

// ----------------------------------------------------------------------------

struct CellSubdomainTupleFlyweightPolicy
{
  typedef std::tuple<EntityRealignCode, Uint, EntityRealignCode, Uint> key_type;
};

// ----------------------------------------------------------------------------

typedef common::Flyweight<TraceEntityTupleInstance, CellSubdomainTupleFlyweightPolicy>
    TraceEntityTuple;

// ----------------------------------------------------------------------------

std::tuple<Uint, Uint> get_relative_refinement_levels(const Uint level_L, const Uint level_R);

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
