#ifndef PDEKIT_Mesh_Std_Region_Builder_T_hpp
#define PDEKIT_Mesh_Std_Region_Builder_T_hpp

#include "mesh/std_region/PointSetTag.hpp"
#include "mesh/std_region/StdRegionBuilder.hpp"

namespace pdekit
{

namespace mesh
{

template <typename ConcreteStdR>
class StdRegionBuilderT : public StdRegionBuilder
{
  public:
  /// Constructor
  StdRegionBuilderT(const std::tuple<Uint> &params);

  /// Destructor
  ~StdRegionBuilderT() override;

  /// Return the name of this point set
  std::string name() const override;

  /// Return the id of this interpolation piont set
  PointSetTag ref_topology_id() const override;

  /// Return the number of nodes (degrees of freedom) for this ps
  Uint nb_dof() const override;

  /// Return the number of p1 nodes (degrees of freedom)
  Uint nb_p1_dof() const override;

  /// Return the topological dimension of this point set
  Uint topo_dim() const override;

  /// Fill vectors with local connectivities
  void local_entities(const Uint dim, std::vector<StdRegionEntity> &entity_list) override;

  /// Fill local connectivity in reference space
  void fill_reference_topology(
      std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &ref_topology,
      common::BlockArray<std::shared_ptr<StdRegionEntity>, Uint> &ref_entities) override;

  /// Fill a matrix with local coordinates
  void coordinates(math::DenseDMat<Real> &coordinates) const override;

  /// Fill a matrix with face normals
  /// Each row of the matrix corresponds to one facet normal
  void facet_normals(math::DenseDMat<Real> &facet_normals) const override;

  /// Fill a vector representing permutation of nodes of this point set
  void permutation(const EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
  const Uint m_poly_order;
};

// ----------------------------------------------------------------------------
// Methods of the template class StdRegionT
// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
StdRegionBuilderT<ConcreteStdR>::StdRegionBuilderT(const std::tuple<Uint> &params)
    : m_poly_order(std::get<0>(params))
{
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
StdRegionBuilderT<ConcreteStdR>::~StdRegionBuilderT()
{
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
std::string StdRegionBuilderT<ConcreteStdR>::name() const
{
  return PointSetTag::fields_to_string(ConcreteStdR::GeomShape, m_poly_order, ConcreteStdR::Type);
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
PointSetTag StdRegionBuilderT<ConcreteStdR>::ref_topology_id() const
{
  return PointSetTag(ConcreteStdR::GeomShape, m_poly_order, ConcreteStdR::Type);
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
Uint StdRegionBuilderT<ConcreteStdR>::nb_dof() const
{
  return ConcreteStdR::nb_dof(m_poly_order);
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
Uint StdRegionBuilderT<ConcreteStdR>::nb_p1_dof() const
{
  return ConcreteStdR::NbVert;
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
Uint StdRegionBuilderT<ConcreteStdR>::topo_dim() const
{
  return ConcreteStdR::TopoDim;
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
void StdRegionBuilderT<ConcreteStdR>::local_entities(const Uint dim,
                                                     std::vector<StdRegionEntity> &entity_list)
{
  /*
  switch(dim)
  {
    case _1D:
    {
      ConcreteStdR::fill_edge_to_node_connectivity(entity_list);
      break;
    }
    case _2D:
    {
      ConcreteStdR::fill_face_to_node_connectivity(entity_list);
      break;
    }
    case _3D:
    {
      ConcreteStdR::fill_volume_to_node_connectivity(entity_list);
      break;
    }

  }
  */
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
void StdRegionBuilderT<ConcreteStdR>::fill_reference_topology(
    std::array<common::BlockArray<SUint, SUint>, (_3D + 1) * (_3D + 1)> &ref_incidences,
    common::BlockArray<std::shared_ptr<StdRegionEntity>, Uint> &ref_entities)
{
  // Fill the relations between topological entities
  ConcreteStdR::fill_topology_relations(ref_incidences);

  // Fill all entities in reference element, where each entity contains
  // reference nodes and
  // knows about its type and topological dimension
  std::unique_ptr<std::vector<std::shared_ptr<StdRegionEntity>>> entities(
      new std::vector<std::shared_ptr<StdRegionEntity>>());
  std::unique_ptr<std::vector<Uint>> offsets(new std::vector<Uint>());

  Uint tot_num_entities = 0;
  offsets->resize(ConcreteStdR::TopologyStorage.size() + 1);
  (*offsets)[0] = 0;

  for (Uint d = 0; d < ConcreteStdR::TopologyStorage.size(); ++d)
  {
    (*offsets)[d + 1] = ConcreteStdR::TopologyStorage[d];
    tot_num_entities += ConcreteStdR::TopologyStorage[d];
  }

  entities->resize(tot_num_entities);
  for (Uint i = 0; i < tot_num_entities; ++i)
  {
    (*entities)[i] = std::make_shared<StdRegionEntity>();
  }

  for (Uint i = 1; i < offsets->size(); ++i)
  {
    (*offsets)[i] += (*offsets)[i - 1];
  }

  ref_entities.build_from_offsets(std::move(entities), std::move(offsets));

  // Nodes - will have just one default value
  // StdRegionEntity &_0D_entity = ref_entities.row(1)[_0D];
  common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> _0D_entities =
      ref_entities.block(_0D);

  const common::ArrayView<const Uint, _1D, Uint> node_ids_view;
  const common::ArrayView<const Uint, _1D, Uint> p1_vert_flags_view;

  _0D_entities[0]->construct(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined),
                             node_ids_view, p1_vert_flags_view, 0);

  // Edges
  common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> edges = ref_entities.block(_1D);
  ConcreteStdR::fill_edge_to_node_connectivity(m_poly_order, edges);

  // Faces
  common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> faces = ref_entities.block(_2D);
  ConcreteStdR::fill_face_to_node_connectivity(m_poly_order, faces);

  // Volumes
  common::ArrayView<std::shared_ptr<StdRegionEntity>, _1D, Uint> volumes = ref_entities.block(_3D);
  ConcreteStdR::fill_volume_to_node_connectivity(m_poly_order, volumes);
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
void StdRegionBuilderT<ConcreteStdR>::coordinates(math::DenseDMat<Real> &coordinates) const
{
  ConcreteStdR::fill_coordinates(m_poly_order, coordinates);
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
void StdRegionBuilderT<ConcreteStdR>::facet_normals(math::DenseDMat<Real> &facet_normals) const
{
  ConcreteStdR::fill_facet_normals(facet_normals);
}

// ----------------------------------------------------------------------------

template <typename ConcreteStdR>
void StdRegionBuilderT<ConcreteStdR>::permutation(const EntityRealignCode &permutation_code,
                                                  std::vector<Uint> &permutation_vec)
{
  ConcreteStdR::fill_permutation(m_poly_order, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
