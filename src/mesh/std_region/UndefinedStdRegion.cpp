#include "mesh/std_region/UndefinedStdRegion.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

// Undefined StdRegion has 0 edges, 0 faces and 0 volumes, but we allocate 1 for
// each since the
// static array has to hold default entities of type 'invalid'
const std::array<Uint, 4> UndefinedStdRegion::TopologyStorage = {1, 1, 1, 1};

// ----------------------------------------------------------------------------

Uint UndefinedStdRegion::nb_dof(const Uint poly_order)
{
  return 0;
}

// ----------------------------------------------------------------------------

void UndefinedStdRegion::fill_topology_relations(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &incidences)
{
}

// ----------------------------------------------------------------------------

void UndefinedStdRegion::fill_edge_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  StdRegionEntity &current_entity = *entity_vec[0];

  common::ArrayView<const Uint, _1D, Uint> vert_ids;
  common::ArrayView<const Uint, _1D, Uint> p1_vert_flags;
  current_entity.construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined), vert_ids,
                           p1_vert_flags, 0);
}

// ----------------------------------------------------------------------------

void UndefinedStdRegion::fill_face_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  StdRegionEntity &current_entity = *entity_vec[0];

  common::ArrayView<const Uint, _1D, Uint> vert_ids;
  common::ArrayView<const Uint, _1D, Uint> p1_vert_flags;
  current_entity.construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined), vert_ids,
                           p1_vert_flags, 0);
}

// ----------------------------------------------------------------------------

void UndefinedStdRegion::fill_volume_to_node_connectivity(
    const Uint poly_order,
    common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> &entity_vec)
{
  StdRegionEntity &current_entity = *entity_vec[0];

  common::ArrayView<const Uint, _1D, Uint> vert_ids;
  common::ArrayView<const Uint, _1D, Uint> p1_vert_flags;
  current_entity.construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined), vert_ids,
                           p1_vert_flags, 0);
}

// ----------------------------------------------------------------------------

void UndefinedStdRegion::fill_coordinates(const Uint poly_order, math::DenseDMat<Real> &coord)
{
  coord.resize(0, 0);
}

// ----------------------------------------------------------------------------

void UndefinedStdRegion::fill_facet_normals(math::DenseDMat<Real> &normals)
{
  normals.resize(0, 0);
}

// ----------------------------------------------------------------------------

void UndefinedStdRegion::fill_permutation(const Uint poly_order,
                                          const EntityRealignCode &permutation_code,
                                          std::vector<Uint> &permutation_vec)
{
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------

void UndefinedStdRegion::fill_rotation_permutation(const Uint poly_order,
                                                   std::vector<Uint> &permutation_vec)
{
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------

void UndefinedStdRegion::fill_flip_permutation(const Uint poly_order,
                                               std::vector<Uint> &permutation_vec)
{
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
