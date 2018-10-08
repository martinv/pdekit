#ifndef PDEKIT_Mesh_Std_Region_Builder_hpp
#define PDEKIT_Mesh_Std_Region_Builder_hpp

#include <array>
#include <memory>

#include "common/BlockArray.hpp"
#include "common/Constants.hpp"
#include "math/DenseDMat.hpp"
#include "mesh/EntityRealignCode.hpp"
#include "mesh/std_region/StdRegionEntity.hpp"

namespace pdekit
{

namespace mesh
{

class StdRegionBuilder
{
  public:
  /// TYPEDEFS
  typedef std::shared_ptr<StdRegionBuilder> ptr;
  typedef std::shared_ptr<StdRegionBuilder const> const_ptr;

  /// Default constructor
  StdRegionBuilder();

  /// Destructor
  virtual ~StdRegionBuilder() = 0;

  static std::string type_name()
  {
    return "StdRegionBuilder";
  }

  virtual std::string name() const = 0;

  /// Return the interpolation point set id
  virtual PointSetTag ref_topology_id() const = 0;

  /// Return the number of nodes (degrees of freedom) for this ips
  virtual Uint nb_dof() const = 0;

  /// Return the number of p1 nodes (degrees of freedom)
  virtual Uint nb_p1_dof() const = 0;

  /// Return the topological dimension of this point set
  virtual Uint topo_dim() const = 0;

  /// Fill vectors with local connectivities
  virtual void local_entities(const Uint dim, std::vector<StdRegionEntity> &entity_list) = 0;

  /// Fill local connectivity in reference space
  virtual void fill_reference_topology(
      std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &ref_incidences,
      common::BlockArray<std::shared_ptr<StdRegionEntity>, Uint> &ref_entities) = 0;

  /// Fill a matrix with local coordinates
  virtual void coordinates(math::DenseDMat<Real> &coordinates) const = 0;

  /// Fill a matrix with face normals
  /// Each row of the matrix corresponds to one facet normal
  virtual void facet_normals(math::DenseDMat<Real> &facet_normals) const = 0;

  /// Fill a vector which represents a permutation of the nodes of a point set
  virtual void permutation(const EntityRealignCode &permutation_code,
                           std::vector<Uint> &permutation_vec) = 0;

  protected:
};

} // namespace mesh

} // namespace pdekit

#endif
