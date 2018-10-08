#ifndef PDEKIT_Interpolation_Function_Space_Traces_hpp
#define PDEKIT_Interpolation_Function_Space_Traces_hpp

#include <functional>
#include <memory>
#include <set>

#include "common/DataMap.hpp"
#include "interpolation/FEValues.hpp"
#include "mesh/DiscreteElemKey.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/Tria.hpp"
#include "mesh/point_set/QuadratureAdaptTransformAlgoFactory.hpp"
#include "mesh/point_set/StdPointSet.hpp"

namespace pdekit
{

namespace interpolation
{

template <typename MeshConfig, Uint DIM = MeshConfig::TDIM>
class FunctionSpaceTraces
{
  public:
  /// TYPEDEFS:
  using ref_elem_map = common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>;
  // One discrete element has
  // 1) A tag of point set (reference element) where it is defined
  // 2) A tag of polynomial (Finite Element basis) used in this element
  // 3) A tag of point set where the basis is evaluated (quadrature points)
  using discrete_elem_store = std::vector<mesh::DiscreteElemKey>;
  using sf_generator_fcn = std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)>;
  using quad_generator_fcn =
      std::function<mesh::PointSetTag(const ElemShape, const Uint elem_order)>;

  /// Constructor
  FunctionSpaceTraces();

  /// Deleted copy constructor
  FunctionSpaceTraces(const FunctionSpaceTraces &other) = delete;

  /// Destructor
  ~FunctionSpaceTraces();

  /// Deleted assignment operator
  FunctionSpaceTraces &operator=(const FunctionSpaceTraces &other) = delete;

  const sf_generator_fcn &sf_generator() const;

  const quad_generator_fcn &quad_generator() const;

  /// Configure the reference elements for all element types on this part of
  /// mesh boundary
  void set_reference_fe_values(
      mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry_facets,
      typename mesh::Tria<MeshConfig> const &tria,
      const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
      const std::function<mesh::PointSetTag(const ElemShape shape, const Uint elem_order)>
          &quad_generator,
      const bool use_filter = false);

  /// Configure the reference elements for all element types on this part of
  /// mesh boundary
  void set_reference_fe_values(
      mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry_facets,
      typename result_of::dof_map_t<MeshConfig> const &cell_dofs,
      const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
      const std::function<mesh::PointSetTag(const ElemShape shape, const Uint elem_order)>
          &quad_generator,
      const bool use_filter = false);

  /// Configure the reference elements for all element types on the skeleton
  /// of this mesh
  void set_reference_fe_values_on_skeleton(
      typename mesh::Tria<MeshConfig> const &cell_topology,
      const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
      const std::function<mesh::PointSetTag(const ElemShape shape, const Uint elem_order)>
          &quad_generator,
      const bool use_filter = false);

  /// Configure the reference elements for all element types on the skeleton
  /// of this mesh
  void set_reference_fe_values_on_skeleton(
      typename mesh::Tria<MeshConfig> const &cell_topology,
      typename result_of::dof_map_t<MeshConfig> const &cell_dofs,
      const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
      const std::function<mesh::PointSetTag(const ElemShape shape, const Uint elem_order)>
          &quad_generator,
      const bool use_filter = false);

  private:
  /// METHODS
  inline static bool discrete_elem_found(const discrete_elem_store &elems,
                                         const mesh::DiscreteElemKey candidate)
  {
    for (const auto &elem : elems)
    {
      if (candidate == elem)
      {
        return true;
      }
    }
    return false;
  }

  protected:
  /// DATA
  /// Function that generates suitable shape functions
  sf_generator_fcn m_sf_generator;

  /// Function that generates appropriate quadratures
  quad_generator_fcn m_quad_generator;

  /// Shape functions for all reference elements in this function space
  ref_elem_map m_elements;

  discrete_elem_store m_discrete_elems;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
FunctionSpaceTraces<MeshConfig, DIM>::FunctionSpaceTraces()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
FunctionSpaceTraces<MeshConfig, DIM>::~FunctionSpaceTraces()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
const typename FunctionSpaceTraces<MeshConfig, DIM>::sf_generator_fcn &FunctionSpaceTraces<
    MeshConfig, DIM>::sf_generator() const
{
  return m_sf_generator;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
const typename FunctionSpaceTraces<MeshConfig, DIM>::quad_generator_fcn &FunctionSpaceTraces<
    MeshConfig, DIM>::quad_generator() const
{
  return m_quad_generator;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void FunctionSpaceTraces<MeshConfig, DIM>::set_reference_fe_values(
    mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry_facets,
    typename mesh::Tria<MeshConfig> const &tria,
    const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
    const std::function<mesh::PointSetTag(const ElemShape shape, const Uint elem_order)>
        &quad_generator,
    const bool use_filter)
{
  m_elements.clear();
  m_discrete_elems.clear();

  m_sf_generator   = sf_generator;
  m_quad_generator = quad_generator;

  // Loop over the boundary and find all element types
  std::set<mesh::PointSetTag> boundary_cell_types;
  for (Uint c = 0; c < bdry_facets.nb_active_cells(); ++c)
  {
    const mesh::FlatIdx parent_cell_id = bdry_facets.parent_cell_id(mesh::ActiveIdx(c));
    const Uint local_id                = bdry_facets.local_id(mesh::ActiveIdx(c));

    const mesh::CellTopologyView<MeshConfig> tcell          = tria.cell(parent_cell_id);
    std::shared_ptr<const mesh::StdRegionEntity> sub_entity = tcell.sub_entity(DIM, local_id);

    boundary_cell_types.insert(sub_entity->pt_set_id());
  }

  mesh::StdPointSet quadrature;

  for (std::set<mesh::PointSetTag>::const_iterator iter = boundary_cell_types.cbegin();
       iter != boundary_cell_types.cend(); ++iter)
  {
    const mesh::PointSetTag std_reg_tag = *iter;

    const ElemShape shape = std_reg_tag.elem_shape();
    const Uint order      = std_reg_tag.poly_order();

    mesh::sf::SFTag const sf_type_tag = sf_generator(shape, order);

    const mesh::PointSetTag quad_tag = quad_generator(shape, order);

    quadrature.change_type(quad_tag);

    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      common::PtrHandle<FEValues> elem_ptr = m_elements.create(
          mesh::PointSetTagExt(std_reg_tag, P0, mesh::CellTransform::NO_TRANS, ie));
      (*elem_ptr).configure(std_reg_tag, sf_type_tag);

      (*elem_ptr).fill_Vandermonde(quadrature.get().coordinates(ie), quadrature.get().weights(ie),
                                   use_filter);
    }

    // Loop that generates new discrete elements
    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      const mesh::PointSetTagExt ext_std_reg_tag(std_reg_tag, P0, mesh::CellTransform::NO_TRANS,
                                                 ie);
      const mesh::PointSetTagExt ext_quad_tag(quad_tag, P0, mesh::CellTransform::NO_TRANS, 0);
      const mesh::DiscreteElemKey candidate(ext_std_reg_tag, sf_type_tag, ext_quad_tag);
      if (!discrete_elem_found(m_discrete_elems, candidate))
      {
        m_discrete_elems.push_back(candidate);
      }
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void FunctionSpaceTraces<MeshConfig, DIM>::set_reference_fe_values(
    mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry_facets,
    typename result_of::dof_map_t<MeshConfig> const &cell_dofs,
    const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
    const std::function<mesh::PointSetTag(const ElemShape shape, const Uint elem_order)>
        &quad_generator,
    const bool use_filter)
{
  m_elements.clear();
  m_discrete_elems.clear();

  m_sf_generator   = sf_generator;
  m_quad_generator = quad_generator;

  // Loop over the boundary and find all element types
  std::set<mesh::PointSetTag> boundary_cell_types;
  for (Uint c = 0; c < bdry_facets.nb_active_cells(); ++c)
  {
    /*
    const mesh::IncidencePair bcell_id = bdry_facets.bdry_cell_id(c);
    mesh::MeshEntity bcell = cell_dofs.active_cell(bcell_id.cell_idx);
    bcell.local_transform(bdry_facets.dim(), bcell_id.local_id);
    */

    const mesh::MeshEntity bcell = bdry_facets.active_cell(cell_dofs, mesh::ActiveIdx(c));
    boundary_cell_types.insert(bcell.pt_set_id());
  }

  mesh::StdPointSet quadrature;

  for (std::set<mesh::PointSetTag>::const_iterator iter = boundary_cell_types.cbegin();
       iter != boundary_cell_types.cend(); ++iter)
  {
    const mesh::PointSetTag std_reg_tag = *iter;

    const ElemShape shape = std_reg_tag.elem_shape();
    const Uint order      = std_reg_tag.poly_order();

    mesh::sf::SFTag const sf_type_tag = sf_generator(shape, order);

    const mesh::PointSetTag quad_tag = quad_generator(shape, order);

    quadrature.change_type(quad_tag);

    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      common::PtrHandle<FEValues> elem_ptr = m_elements.create(
          mesh::PointSetTagExt(std_reg_tag, P0, mesh::CellTransform::NO_TRANS, ie));
      (*elem_ptr).configure(std_reg_tag, sf_type_tag);

      (*elem_ptr).fill_Vandermonde(quadrature.get().coordinates(ie), quadrature.get().weights(ie),
                                   use_filter);
    }

    // Loop that generates new discrete elements
    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      const mesh::PointSetTagExt ext_std_reg_tag(std_reg_tag, P0, mesh::CellTransform::NO_TRANS,
                                                 ie);
      const mesh::PointSetTagExt ext_quad_tag(quad_tag, P0, mesh::CellTransform::NO_TRANS, 0);
      const mesh::DiscreteElemKey candidate(ext_std_reg_tag, sf_type_tag, ext_quad_tag);
      if (!discrete_elem_found(m_discrete_elems, candidate))
      {
        m_discrete_elems.push_back(candidate);
      }
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void FunctionSpaceTraces<MeshConfig, DIM>::set_reference_fe_values_on_skeleton(
    typename mesh::Tria<MeshConfig> const &cell_topology,
    const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
    const std::function<mesh::PointSetTag(const ElemShape shape, const Uint elem_order)>
        &quad_generator,
    const bool use_filter)
{
  mesh::StdPointSet quadrature;

  m_elements.clear();
  m_discrete_elems.clear();

  m_sf_generator   = sf_generator;
  m_quad_generator = quad_generator;

  const Uint facet_dim = MeshConfig::TDIM - 1;

  for (Uint f = 0; f < cell_topology.active_skeleton_size(); ++f)
  {
    const mesh::TraceIncidences facet_block =
        cell_topology.active_skeleton_entry(facet_dim, mesh::ActiveIdx(f));

    // In case this is an interior facet, the block size is 2 and there is
    // one facet in block on the 'left' side and one facet on the 'right'
    // side If the size of the block is 1, then this is a boundary facet

    const Uint local_idx_L = 0;
    const Uint local_idx_R = (facet_block.size() == 2) ? 1 : 0;

    // In order to get the 'active' cell indexes, we need to retrieve
    // topology cells and ask for active cell indexes through topology cells
    const mesh::CellTopologyView<MeshConfig> tcell_L =
        cell_topology.cell(mesh::FlatIdx(facet_block.cell_id(local_idx_L)));
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        cell_topology.cell(mesh::FlatIdx(facet_block.cell_id(local_idx_R)));

    const mesh::PointSetTag pt_set_id_L = tcell_L.pt_set_id(DIM, facet_block.local_id(local_idx_L));
    const mesh::PointSetTag pt_set_id_R = tcell_R.pt_set_id(DIM, facet_block.local_id(local_idx_R));

    const mesh::EntityRealignCode local_sub_tag_L =
        facet_block.permutation(local_idx_L).get().code();
    mesh::PointSetTagExt key_L(pt_set_id_L, P0, local_sub_tag_L.adapt_op_id(),
                               local_sub_tag_L.local_pos_in_parent());

    common::PtrHandle<FEValues> fe_values_L = m_elements.std_region_data(key_L);

    if (fe_values_L.is_null())
    {
      std::cout << "Adding " << key_L << " to collection" << std::endl;

      const ElemShape shape = pt_set_id_L.elem_shape();
      const Uint order      = pt_set_id_L.poly_order();

      mesh::sf::SFTag const facet_sf_type = sf_generator(shape, order);
      const mesh::PointSetTag quad_tag    = quad_generator(shape, order);

      fe_values_L = m_elements.create(key_L);
      (*fe_values_L).configure(pt_set_id_L, facet_sf_type);

      quadrature.change_type(quad_tag);

      const mesh::CellTransform adapt_op_id_L = key_L.cell_transform_id();

      if (adapt_op_id_L == mesh::CellTransform::NO_TRANS)
      {
        (*fe_values_L)
            .fill_Vandermonde(quadrature.get().coordinates(), quadrature.get().weights(),
                              use_filter);
      }
      else
      {
        QuadratureAdaptTransformAlgoFactory::instance_type &quad_adapt_algo_fact =
            QuadratureAdaptTransformAlgoFactory::instance();
        const QuadratureAdaptTransformAlgoFactory::instance_type::const_product_base_ptr
            algo_transform = quad_adapt_algo_fact.create(adapt_op_id_L);

        math::DenseDMat<Real> transf_quad_coords;
        math::DenseDVec<Real> transf_quad_weights;

        algo_transform->compute_transformed_coords(quadrature.get().coordinates(), key_L.local_id(),
                                                   transf_quad_coords);
        algo_transform->compute_transformed_weights(quadrature.get().weights(), key_L.local_id(),
                                                    transf_quad_weights);

        (*fe_values_L).fill_Vandermonde(transf_quad_coords, transf_quad_weights, use_filter);
      }
    }

    // NEW
    // Compute discrete element candidate
    const ElemShape shape_L = pt_set_id_L.elem_shape();
    const Uint order_L      = pt_set_id_L.poly_order();

    mesh::sf::SFTag const facet_sf_type_L = sf_generator(shape_L, order_L);
    const mesh::PointSetTag quad_tag_L    = quad_generator(shape_L, order_L);
    const mesh::PointSetTagExt ext_quad_tag_L(quad_tag_L, P0, local_sub_tag_L.adapt_op_id(),
                                              local_sub_tag_L.local_pos_in_parent());

    const mesh::DiscreteElemKey candidate_L(key_L, facet_sf_type_L, ext_quad_tag_L);
    if (!discrete_elem_found(m_discrete_elems, candidate_L))
    {
      m_discrete_elems.push_back(candidate_L);
    }

    if (facet_block.size() == 2)
    {
      const mesh::EntityRealignCode local_sub_tag_R =
          facet_block.permutation(local_idx_R).get().code();
      mesh::PointSetTagExt key_R(pt_set_id_R, P0, local_sub_tag_R.adapt_op_id(),
                                 local_sub_tag_R.local_pos_in_parent());

      common::PtrHandle<FEValues> fe_values_R = m_elements.std_region_data(key_R);

      if (fe_values_R.is_null())
      {
        std::cout << "Adding " << key_R << " to collection" << std::endl;
        const ElemShape shape = pt_set_id_R.elem_shape();
        const Uint order      = pt_set_id_R.poly_order();

        mesh::sf::SFTag const facet_sf_type = sf_generator(shape, order);
        const mesh::PointSetTag quad_tag    = quad_generator(shape, order);

        fe_values_R = m_elements.create(key_R);
        (*fe_values_R).configure(pt_set_id_R, facet_sf_type);

        quadrature.change_type(quad_tag);

        const mesh::CellTransform adapt_op_id_R = key_R.cell_transform_id();

        if (adapt_op_id_R == mesh::CellTransform::NO_TRANS)
        {
          (*fe_values_R)
              .fill_Vandermonde(quadrature.get().coordinates(), quadrature.get().weights(),
                                use_filter);
        }
        else
        {
          QuadratureAdaptTransformAlgoFactory::instance_type &quad_adapt_algo_fact =
              QuadratureAdaptTransformAlgoFactory::instance();
          const QuadratureAdaptTransformAlgoFactory::instance_type::const_product_base_ptr
              algo_transform = quad_adapt_algo_fact.create(adapt_op_id_R);

          math::DenseDMat<Real> transf_quad_coords;
          math::DenseDVec<Real> transf_quad_weights;

          algo_transform->compute_transformed_coords(quadrature.get().coordinates(),
                                                     key_R.local_id(), transf_quad_coords);
          algo_transform->compute_transformed_weights(quadrature.get().weights(), key_R.local_id(),
                                                      transf_quad_weights);

          (*fe_values_R).fill_Vandermonde(transf_quad_coords, transf_quad_weights, use_filter);
        }
      }

      // NEW
      // Compute discrete element candidate
      const ElemShape shape_R = pt_set_id_R.elem_shape();
      const Uint order_R      = pt_set_id_R.poly_order();

      mesh::sf::SFTag const facet_sf_type_R = sf_generator(shape_R, order_R);
      const mesh::PointSetTag quad_tag_R    = quad_generator(shape_R, order_R);
      const mesh::PointSetTagExt ext_quad_tag_R(quad_tag_R, P0, local_sub_tag_R.adapt_op_id(),
                                                local_sub_tag_R.local_pos_in_parent());

      const mesh::DiscreteElemKey candidate_R(key_R, facet_sf_type_R, ext_quad_tag_R);
      if (!discrete_elem_found(m_discrete_elems, candidate_R))
      {
        m_discrete_elems.push_back(candidate_R);
      }

    } // If size of this block is 2

  } // Loop over the mesh skeleton
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void FunctionSpaceTraces<MeshConfig, DIM>::set_reference_fe_values_on_skeleton(
    typename mesh::Tria<MeshConfig> const &cell_topology,
    typename result_of::dof_map_t<MeshConfig> const &cell_dofs,
    const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
    const std::function<mesh::PointSetTag(const ElemShape shape, const Uint elem_order)>
        &quad_generator,
    const bool use_filter)
{
  mesh::StdPointSet quadrature;

  m_elements.clear();
  m_discrete_elems.clear();

  m_sf_generator   = sf_generator;
  m_quad_generator = quad_generator;

  const Uint facet_dim = MeshConfig::TDIM - 1;

  for (Uint f = 0; f < cell_topology.active_skeleton_size(); ++f)
  {
    const mesh::TraceIncidences facet_block =
        cell_topology.active_skeleton_entry(facet_dim, mesh::ActiveIdx(f));

    // In case this is an interior facet, the block size is 2 and there is
    // one facet in block on the 'left' side and one facet on the 'right'
    // side If the size of the block is 1, then this is a boundary facet

    const Uint local_idx_L = 0;
    const Uint local_idx_R = (facet_block.size() == 2) ? 1 : 0;

    // In order to get the 'active' cell indexes, we need to retrieve
    // topology cells and ask for active cell indexes through topology cells
    const mesh::CellTopologyView<MeshConfig> tcell_L =
        cell_topology.cell(mesh::FlatIdx(facet_block.cell_id(local_idx_L)));
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        cell_topology.cell(mesh::FlatIdx(facet_block.cell_id(local_idx_R)));

    // Get the actual 'active' cell indices
    const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
    const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

    mesh::MeshEntity facet_L = cell_dofs.active_cell(active_cell_id_L);
    facet_L.local_transform(facet_dim, facet_block.local_id(local_idx_L));

    const mesh::EntityRealignCode local_sub_tag_L =
        facet_block.permutation(local_idx_L).get().code();
    mesh::PointSetTagExt key_L(facet_L.pt_set_id(), P0, local_sub_tag_L.adapt_op_id(),
                               local_sub_tag_L.local_pos_in_parent());

    common::PtrHandle<FEValues> fe_values_L = m_elements.std_region_data(key_L);

    if (fe_values_L.is_null())
    {
      const mesh::PointSetTag ref_topo_id = facet_L.pt_set_id();

      std::cout << "Adding " << key_L << " to collection" << std::endl;

      const ElemShape shape = ref_topo_id.elem_shape();
      const Uint order      = ref_topo_id.poly_order();

      mesh::sf::SFTag const facet_sf_type = sf_generator(shape, order);
      const mesh::PointSetTag quad_tag    = quad_generator(shape, order);

      fe_values_L = m_elements.create(key_L);
      (*fe_values_L).configure(ref_topo_id, facet_sf_type);

      quadrature.change_type(quad_tag);

      const mesh::CellTransform adapt_op_id_L = key_L.cell_transform_id();

      if (adapt_op_id_L == mesh::CellTransform::NO_TRANS)
      {
        (*fe_values_L)
            .fill_Vandermonde(quadrature.get().coordinates(), quadrature.get().weights(),
                              use_filter);
      }
      else
      {
        QuadratureAdaptTransformAlgoFactory::instance_type &quad_adapt_algo_fact =
            QuadratureAdaptTransformAlgoFactory::instance();
        const QuadratureAdaptTransformAlgoFactory::instance_type::const_product_base_ptr
            algo_transform = quad_adapt_algo_fact.create(adapt_op_id_L);

        math::DenseDMat<Real> transf_quad_coords;
        math::DenseDVec<Real> transf_quad_weights;

        algo_transform->compute_transformed_coords(quadrature.get().coordinates(), key_L.local_id(),
                                                   transf_quad_coords);
        algo_transform->compute_transformed_weights(quadrature.get().weights(), key_L.local_id(),
                                                    transf_quad_weights);

        (*fe_values_L).fill_Vandermonde(transf_quad_coords, transf_quad_weights, use_filter);
      }
    }

    // NEW
    // Compute discrete element candidate
    const mesh::PointSetTag ref_topo_id_L = facet_L.pt_set_id();
    const ElemShape shape_L               = ref_topo_id_L.elem_shape();
    const Uint order_L                    = ref_topo_id_L.poly_order();

    mesh::sf::SFTag const facet_sf_type_L = sf_generator(shape_L, order_L);
    const mesh::PointSetTag quad_tag_L    = quad_generator(shape_L, order_L);
    const mesh::PointSetTagExt ext_quad_tag_L(quad_tag_L, P0, local_sub_tag_L.adapt_op_id(),
                                              local_sub_tag_L.local_pos_in_parent());

    const mesh::DiscreteElemKey candidate_L(key_L, facet_sf_type_L, ext_quad_tag_L);
    if (!discrete_elem_found(m_discrete_elems, candidate_L))
    {
      m_discrete_elems.push_back(candidate_L);
    }

    if (facet_block.size() == 2)
    {
      mesh::MeshEntity facet_R = cell_dofs.active_cell(active_cell_id_R);
      facet_R.local_transform(facet_dim, facet_block.local_id(local_idx_R));

      const mesh::EntityRealignCode local_sub_tag_R =
          facet_block.permutation(local_idx_R).get().code();
      mesh::PointSetTagExt key_R(facet_R.pt_set_id(), P0, local_sub_tag_R.adapt_op_id(),
                                 local_sub_tag_R.local_pos_in_parent());

      common::PtrHandle<FEValues> fe_values_R = m_elements.std_region_data(key_R);

      if (fe_values_R.is_null())
      {
        const mesh::PointSetTag ref_topo_id = facet_R.pt_set_id();

        std::cout << "Adding " << key_R << " to collection" << std::endl;
        const ElemShape shape = ref_topo_id.elem_shape();
        const Uint order      = ref_topo_id.poly_order();

        mesh::sf::SFTag const facet_sf_type = sf_generator(shape, order);
        const mesh::PointSetTag quad_tag    = quad_generator(shape, order);

        fe_values_R = m_elements.create(key_R);
        (*fe_values_R).configure(ref_topo_id, facet_sf_type);

        quadrature.change_type(quad_tag);

        const mesh::CellTransform adapt_op_id_R = key_R.cell_transform_id();

        if (adapt_op_id_R == mesh::CellTransform::NO_TRANS)
        {
          (*fe_values_R)
              .fill_Vandermonde(quadrature.get().coordinates(), quadrature.get().weights(),
                                use_filter);
        }
        else
        {
          QuadratureAdaptTransformAlgoFactory::instance_type &quad_adapt_algo_fact =
              QuadratureAdaptTransformAlgoFactory::instance();
          const QuadratureAdaptTransformAlgoFactory::instance_type::const_product_base_ptr
              algo_transform = quad_adapt_algo_fact.create(adapt_op_id_R);

          math::DenseDMat<Real> transf_quad_coords;
          math::DenseDVec<Real> transf_quad_weights;

          algo_transform->compute_transformed_coords(quadrature.get().coordinates(),
                                                     key_R.local_id(), transf_quad_coords);
          algo_transform->compute_transformed_weights(quadrature.get().weights(), key_R.local_id(),
                                                      transf_quad_weights);

          (*fe_values_R).fill_Vandermonde(transf_quad_coords, transf_quad_weights, use_filter);
        }
      }

      // NEW
      // Compute discrete element candidate
      const mesh::PointSetTag ref_topo_id_R = facet_R.pt_set_id();
      const ElemShape shape_R               = ref_topo_id_R.elem_shape();
      const Uint order_R                    = ref_topo_id_R.poly_order();

      mesh::sf::SFTag const facet_sf_type_R = sf_generator(shape_R, order_R);
      const mesh::PointSetTag quad_tag_R    = quad_generator(shape_R, order_R);
      const mesh::PointSetTagExt ext_quad_tag_R(quad_tag_R, P0, local_sub_tag_R.adapt_op_id(),
                                                local_sub_tag_R.local_pos_in_parent());

      const mesh::DiscreteElemKey candidate_R(key_R, facet_sf_type_R, ext_quad_tag_R);
      if (!discrete_elem_found(m_discrete_elems, candidate_R))
      {
        m_discrete_elems.push_back(candidate_R);
      }

    } // If size of this block is 2

  } // Loop over the mesh skeleton
}

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // namespace pdekit

#endif
