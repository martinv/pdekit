#ifndef PDEKIT_Mesh_Std_Region_hpp
#define PDEKIT_Mesh_Std_Region_hpp

#include <unordered_map>

#include "common/BlockArray.hpp"
#include "common/Flyweight.hpp"
#include "common/PtrHandle.hpp"
#include "math/DenseConstVecView.hpp"
#include "mesh/std_region/StdRegionEntity.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

class StdRegionInstance
{
  public:
  /// TYPEDEFS

  using incidence_table_type = common::BlockArray<SUint, SUint>;
  /// Type of list that stores all entities incident to given entity
  using incidence_list_type = common::ArrayView<const SUint, _1D, SUint>;
  using coordinates_type    = StdRegionEntity::coordinates_type;

  /// STATIC VARIABLES

  static const PointSetTag undefined;

  /// Static method to fill a reference element
  static void construct(const PointSetTag pt_set_id, StdRegionInstance &std_reg_instance);

  /// Constructor
  StdRegionInstance();

  /// Constructor which directly sets the
  /// reference element type
  StdRegionInstance(const PointSetTag pt_set_id);

  /// Copy constructor
  StdRegionInstance(const StdRegionInstance &std_reg_instance);

  /// Assignement operator
  StdRegionInstance &operator=(const StdRegionInstance &std_reg_instance);

  /// Destructor
  ~StdRegionInstance();

  /// Return the interpolation point set type that
  /// this reference element represents
  PointSetTag pt_set_id() const;

  /// Return the topological dimension of this reference element
  Uint topo_dim() const;

  /// Return the number of entities of given dimension
  /// Example: nb_entities(2) returns the number of faces
  /// (=2D entities) of this reference element
  Uint nb_entities(const Uint dim) const;

  /// Return the number of nodes in this reference element
  Uint nb_nodes() const;

  /// Return the number of p1 nodes in this standard region
  Uint nb_p1_nodes() const;

  /// Return the reference coordinates of this reference element
  coordinates_type const &coordinates() const;

  /// Return the sub-entity of dimension dim with index id
  /// Example: entity(1,3) would return the third
  /// edge (= 1D entity) of this reference element
  std::shared_ptr<StdRegionEntity const> elem_entity(const Uint dim, const Uint id) const;

  /// Return ids of all entities of dimension dim which are incident to an
  /// entity
  /// @param my_dim    -  dimension of entity whose incidences we are looking
  /// for
  /// @param my_id     -  local index of entity whose incidences we are
  /// looking for
  /// @param other_dim -  dimension of incident entities we are searching for
  /// @example            incident_entities(2,3,1) would mean that we want a
  /// list
  ///                     of all 1D entities (edges) that are incident to the
  ///                     third 2D entity (face)
  /// This method just forwards the call to the underlying reference element
  /// instance (m_instance_ptr)
  incidence_list_type incident_entities(const Uint my_dim, const Uint my_id,
                                        const Uint other_dim) const;

  /// Get a facet normal (in reference space)
  const math::DenseConstVecView<Real> facet_normal(const Uint facet_id) const;

  /// Return the name of this reference element
  const std::string type_name() const;

  /// Print all the hierarchy of reference entities
  void print_complete_topology() const;

  private:
  /// DATA

  /// Interpolation point set id
  PointSetTag m_type_id;

  /// Topological dimension
  Uint m_topo_dim;

  /// Number of P1 nodes
  Uint m_nb_p1_nodes;

  /// Indidences (d0,d1) are stored on position (d0*3+d1)
  std::array<incidence_table_type, (_3D + 1) * (_3D + 1)> m_incidences;
  common::BlockArray<std::shared_ptr<StdRegionEntity>, Uint> m_entity_storage;

  /// Storage for facet normals
  math::DenseDMat<Real> m_facet_normals;
};

// ----------------------------------------------------------------------------

inline PointSetTag StdRegionInstance::pt_set_id() const
{
  return m_type_id;
}

// ----------------------------------------------------------------------------

inline Uint StdRegionInstance::topo_dim() const
{
  return m_topo_dim;
}

// ----------------------------------------------------------------------------

inline Uint StdRegionInstance::nb_nodes() const
{
  return m_entity_storage.const_block(m_topo_dim)[0]->nb_vert();
}

// ----------------------------------------------------------------------------

inline Uint StdRegionInstance::nb_p1_nodes() const
{
  return m_nb_p1_nodes;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

struct StdRegionFlyweightPolicy
{
  typedef PointSetTag key_type;
};

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

typedef common::Flyweight<StdRegionInstance, StdRegionFlyweightPolicy> StdRegion;

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
