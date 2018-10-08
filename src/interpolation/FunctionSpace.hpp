#ifndef PDEKIT_Interpolation_Function_Space_hpp
#define PDEKIT_Interpolation_Function_Space_hpp

#include <type_traits>

#include "interpolation/FunctionSpaceCells.hpp"
#include "interpolation/FunctionSpaceTraces.hpp"

namespace pdekit
{

namespace interpolation
{

template <typename MeshConfig, Uint DIM = MeshConfig::TDIM>
class FunctionSpace
    : public std::conditional<DIM == MeshConfig::TDIM, FunctionSpaceCells<MeshConfig>,
                              FunctionSpaceTraces<MeshConfig, DIM>>::type
{
  private:
  using super_t = typename std::conditional<DIM == MeshConfig::TDIM, FunctionSpaceCells<MeshConfig>,
                                            FunctionSpaceTraces<MeshConfig, DIM>>::type;

  public:
  /// TYPEDEFS:
  using ptr       = std::shared_ptr<FunctionSpace<MeshConfig, DIM>>;
  using const_ptr = std::shared_ptr<FunctionSpace<MeshConfig, DIM> const>;

  using ref_elem_map        = typename super_t::ref_elem_map;
  using discrete_elem_store = typename super_t::discrete_elem_store;
  using sf_generator_fcn    = typename super_t::sf_generator_fcn;
  using quad_generator_fcn  = typename super_t::quad_generator_fcn;

  /// Constructor
  FunctionSpace();

  /// Deleted copy constructor
  FunctionSpace(const FunctionSpace &other) = delete;

  /// Destructor
  ~FunctionSpace();

  /// Deleted assignment operator
  FunctionSpace &operator=(const FunctionSpace &other) = delete;

  /// Set reference values simpulatenously on a pair spaces. The parameter
  /// quad_order_rule selects quadrature based on a pair of element tags
  /// (geo_elem_tag, sol_elem_tag)
  template <typename QRule>
  static void set_reference_fe_values_on_space_pair(
      const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
      const typename result_of::dof_map_t<MeshConfig> &sol_dofs, SFunc sf_type,
      PointSetID quad_type, const QRule &quad_order_rule, FunctionSpace &geo_fe_space,
      FunctionSpace &sol_fe_space,
      std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> &tag_combinations);

  /// Set reference values simpulatenously on a pair spaces. The parameter
  /// quad_order_rule selects quadrature based on a pair of element tags
  /// (geo_elem_tag, sol_elem_tag)

  template <typename QRule>
  static void set_reference_fe_values_on_space_pair(
      mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry_facets,
      const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
      const typename result_of::dof_map_t<MeshConfig> &sol_dofs, SFunc sf_type,
      PointSetID quad_type, const QRule &quad_order_rule, FunctionSpace &geo_fe_space,
      FunctionSpace &sol_fe_space,
      std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> &tag_combinations);

  /// Get all the reference elements in this space
  ref_elem_map const &reference_elements() const;

  /// Get all the reference elements in this space
  discrete_elem_store const &discrete_elements() const;

  /// Return the number of element types in this space
  Uint nb_elem_types() const;

  /// Return a reference element corresponding to certain interpolation type
  /// set
  common::PtrHandle<FEValues> fe(mesh::PointSetTagExt const elem_type);

  /// Return a reference element corresponding to certain interpolation type
  /// set
  const common::PtrHandle<FEValues const> fe(mesh::PointSetTagExt const elem_type) const;

  void print_element_types() const;
  void print_discrete_element_types() const;
};

// ----------------------------------------------------------------------------

/// Free functions
/// Find all pairs of tags on two DOF handlers
template <typename MeshConfig, typename QRule>
void fill_fe_value_pairs(
    const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs,
    const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order,
                                        const PointSetID rt_type)> &sf_generator,
    const std::function<mesh::PointSetTag(const ElemShape shape, const Uint geo_elem_order,
                                          const Uint sol_elem_order)> &quad_generator,
    common::DataMap<mesh::PointSetTagExt, FEValues> &geo_fe_values_map,
    common::DataMap<mesh::PointSetTagExt, FEValues> &sol_fe_values_map,
    std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> &tag_combinations);

/// Find all pairs of tags on the BOUNDARY of two DOF handlers
template <typename MeshConfig, typename QRule>
void fill_fe_value_pairs(
    mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry_facets,
    const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs, SFunc sf_type, PointSetID quad_type,
    const QRule &quad_order_rule,
    common::DataMap<mesh::PointSetTagExt, FEValues> &geo_fe_values_map,
    common::DataMap<mesh::PointSetTagExt, FEValues> &sol_fe_values_map,
    std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> &tag_combinations);

/// Find all pairs of tags on the SKELETON of two DOF handlers
template <typename MeshConfig, typename QRule, Uint FacetDIM = MeshConfig::TDIM - 1>
void fill_fe_value_pairs(
    typename mesh::Tria<MeshConfig> const &cell_topology,
    const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs, SFunc sf_type, PointSetID quad_type,
    const QRule &quad_order_rule,
    common::DataMap<mesh::PointSetTagExt, FEValues> &geo_fe_values_map,
    common::DataMap<mesh::PointSetTagExt, FEValues> &sol_fe_values_map,
    std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> &tag_combinations,
    const bool use_filter = false);

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
FunctionSpace<MeshConfig, DIM>::FunctionSpace() : super_t()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
FunctionSpace<MeshConfig, DIM>::~FunctionSpace()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
template <typename QRule>
void FunctionSpace<MeshConfig, DIM>::set_reference_fe_values_on_space_pair(
    const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs, SFunc sf_type, PointSetID quad_type,
    const QRule &quad_order_rule, FunctionSpace &geo_fe_space, FunctionSpace &sol_fe_space,
    std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> &tag_combinations)
{

  auto sf_generator = [sf_type](const ElemShape shape, const Uint order,
                                const PointSetID rt_type) -> mesh::sf::SFTag {
    return mesh::sf::SFTag(shape, sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [quad_type](const ElemShape shape, const Uint geo_elem_order,
                                    const Uint sol_elem_order) -> mesh::PointSetTag {
    return mesh::PointSetTag(shape, std::max(1u, std::max(2 * geo_elem_order, 2 * sol_elem_order)),
                             quad_type);
  };

  fill_fe_value_pairs<MeshConfig>(geo_dofs, sol_dofs, sf_generator, quad_generator,
                                  geo_fe_space.super_t::m_elements,
                                  sol_fe_space.super_t::m_elements, tag_combinations);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
template <typename QRule>
void FunctionSpace<MeshConfig, DIM>::set_reference_fe_values_on_space_pair(
    mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry_facets,
    const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs, SFunc sf_type, PointSetID quad_type,
    const QRule &quad_order_rule, FunctionSpace &geo_fe_space, FunctionSpace &sol_fe_space,
    std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> &tag_combinations)
{
  fill_fe_value_pairs<MeshConfig, QRule>(bdry_facets, geo_dofs, sol_dofs, sf_type, quad_type,
                                         quad_order_rule, geo_fe_space.super_t::m_elements,
                                         sol_fe_space.super_t::m_elements, tag_combinations);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
typename FunctionSpace<MeshConfig, DIM>::ref_elem_map const &FunctionSpace<
    MeshConfig, DIM>::reference_elements() const
{
  return super_t::m_elements;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
typename FunctionSpace<MeshConfig, DIM>::discrete_elem_store const &FunctionSpace<
    MeshConfig, DIM>::discrete_elements() const
{
  return super_t::m_discrete_elems;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
Uint FunctionSpace<MeshConfig, DIM>::nb_elem_types() const
{
  return super_t::m_elements.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
common::PtrHandle<FEValues> FunctionSpace<MeshConfig, DIM>::fe(mesh::PointSetTagExt const elem_type)
{
  return super_t::m_elements.std_region_data(elem_type);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
const common::PtrHandle<FEValues const> FunctionSpace<MeshConfig, DIM>::fe(
    mesh::PointSetTagExt const elem_type) const
{
  return super_t::m_elements.std_region_data(elem_type);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void FunctionSpace<MeshConfig, DIM>::print_element_types() const
{
  for (common::DataMap<mesh::PointSetTagExt, FEValues>::const_iterator it =
           super_t::m_elements.cbegin();
       it != super_t::m_elements.cend(); ++it)
  {
    std::cout << "Std region " << it.key_value().std_region_tag().as_string()
              << ", local id = " << it.key_value().local_id()
              << ", order = " << it.key_value().key_p_order()
              << ", adapt op id = " << it.key_value().cell_transform_id() << std::endl;
    (*it.data_ptr()).print();
    std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void FunctionSpace<MeshConfig, DIM>::print_discrete_element_types() const
{
  std::cout << "Discrete elements:" << std::endl;
  for (const auto &elem : super_t::m_discrete_elems)
  {
    std::cout << "{" << elem << "}" << std::endl;
  }
}

// ----------------------------------------------------------------------------

class FSpaceTupleHash
    : public std::unary_function<
          std::tuple<mesh::PointSetTagExt, mesh::sf::SFTag, mesh::PointSetTagExt>, std::size_t>
{
  public:
  inline std::size_t operator()(
      const std::tuple<mesh::PointSetTagExt, mesh::sf::SFTag, mesh::PointSetTagExt> &key) const
  {
    const mesh::PointSetTagExt key_part0 = std::get<0>(key);
    const mesh::sf::SFTag key_part1      = std::get<1>(key);
    const mesh::PointSetTagExt key_part2 = std::get<2>(key);

    const std::size_t hash_val_0 =
        key_part0.std_region_tag().store_value() ^ key_part0.local_id() ^ key_part0.key_p_order();
    const std::size_t hash_val_1 = key_part1.store_value();
    const std::size_t hash_val_2 =
        key_part2.std_region_tag().store_value() ^ key_part2.local_id() ^ key_part2.key_p_order();

    return (hash_val_0 << 10) + (hash_val_1 << 5) + hash_val_2;
  }
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void fill_fe_value_pairs(
    const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs,
    const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order,
                                        const PointSetID rt_type)> &sf_generator,
    const std::function<mesh::PointSetTag(const ElemShape shape, const Uint geo_elem_order,
                                          const Uint sol_elem_order)> &quad_generator,
    common::DataMap<mesh::PointSetTagExt, FEValues> &geo_fe_values_map,
    common::DataMap<mesh::PointSetTagExt, FEValues> &sol_fe_values_map,
    std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> &tag_combinations)
{
  geo_fe_values_map.clear();
  sol_fe_values_map.clear();
  tag_combinations.resize(0);

  std::set<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> tmp_tag_pairs;

  mesh::StdPointSet quadrature;

  for (Uint ac = 0; ac < geo_dofs.nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity geo_cell = geo_dofs.active_cell(mesh::ActiveIdx(ac));
    const mesh::MeshEntity sol_cell = sol_dofs.active_cell(mesh::ActiveIdx(ac));

    // std::cout << geo_cell << std::endl;
    // std::cout << sol_cell << std::endl << std::endl;

    const mesh::PointSetTag geo_tag = geo_cell.pt_set_id();
    const mesh::PointSetTag sol_tag = sol_cell.pt_set_id();

    const mesh::PointSetTag quad_tag =
        quad_generator(geo_tag.elem_shape(), geo_tag.poly_order(), sol_tag.poly_order());

    const mesh::PointSetTagExt geo_tag_ext(geo_tag, quad_tag.poly_order(),
                                           mesh::CellTransform::NO_TRANS, 0u);
    common::PtrHandle<FEValues> geo_fe_val_ptr = geo_fe_values_map.std_region_data(geo_tag_ext);

    if (geo_fe_val_ptr.is_null())
    {
      geo_fe_val_ptr = geo_fe_values_map.create(geo_tag_ext);
      quadrature.change_type(quad_tag);

      const mesh::sf::SFTag geo_sf_tag =
          sf_generator(geo_tag.elem_shape(), geo_tag.poly_order(), geo_tag.ref_topology());
      FEValues &geo_fe = (*geo_fe_val_ptr);
      geo_fe.configure(geo_tag, geo_sf_tag);
      geo_fe.fill_Vandermonde(quadrature.get().coordinates(), quadrature.get().weights());
    }

    const mesh::PointSetTagExt sol_tag_ext(sol_tag, quad_tag.poly_order(),
                                           mesh::CellTransform::NO_TRANS, 0u);
    common::PtrHandle<FEValues> sol_fe_val_ptr = sol_fe_values_map.std_region_data(sol_tag_ext);

    if (sol_fe_val_ptr.is_null())
    {
      sol_fe_val_ptr = sol_fe_values_map.create(sol_tag_ext);
      quadrature.change_type(quad_tag);

      const mesh::sf::SFTag sol_sf_tag =
          sf_generator(sol_tag.elem_shape(), sol_tag.poly_order(), sol_tag.ref_topology());
      FEValues &sol_fe = (*sol_fe_val_ptr);
      sol_fe.configure(sol_tag, sol_sf_tag);
      sol_fe.fill_Vandermonde(quadrature.get().coordinates(), quadrature.get().weights());
    } // if

    tmp_tag_pairs.insert(std::make_pair(geo_tag_ext, sol_tag_ext));
  }

  for (auto tag_pair : tmp_tag_pairs)
  {
    tag_combinations.push_back(tag_pair);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename QRule>
void fill_fe_value_pairs(
    mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry_facets,
    const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs, SFunc sf_type, PointSetID quad_type,
    const QRule &quad_order_rule,
    common::DataMap<mesh::PointSetTagExt, FEValues> &geo_fe_values_map,
    common::DataMap<mesh::PointSetTagExt, FEValues> &sol_fe_values_map,
    std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> &tag_combinations)
{
  geo_fe_values_map.clear();
  sol_fe_values_map.clear();
  tag_combinations.resize(0);

  std::set<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> tmp_tag_pairs;

  mesh::StdPointSet quadrature;

  for (Uint ac = 0; ac < bdry_facets.nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity geo_cell = bdry_facets.active_cell(geo_dofs, mesh::ActiveIdx(ac));
    const mesh::MeshEntity sol_cell = bdry_facets.active_cell(sol_dofs, mesh::ActiveIdx(ac));

    // std::cout << geo_cell << std::endl;
    // std::cout << sol_cell << std::endl << std::endl;

    const mesh::PointSetTag geo_tag = geo_cell.pt_set_id();
    const mesh::PointSetTag sol_tag = sol_cell.pt_set_id();

    const Uint quad_order = quad_order_rule(geo_tag, sol_tag);

    const mesh::PointSetTagExt geo_tag_ext(geo_tag, quad_order, mesh::CellTransform::NO_TRANS, 0u);

    common::PtrHandle<FEValues> geo_fe_val_ptr = geo_fe_values_map.std_region_data(geo_tag_ext);

    if (geo_fe_val_ptr.is_null())
    {
      geo_fe_val_ptr = geo_fe_values_map.create(geo_tag_ext);
      const mesh::PointSetTag quad_tag(geo_tag.elem_shape(), quad_order, quad_type);
      quadrature.change_type(quad_tag);

      const mesh::sf::SFTag geo_sf_tag(geo_tag.elem_shape(), sf_type, geo_tag.poly_order(),
                                       ModalBasis::Modal);
      FEValues &geo_fe = (*geo_fe_val_ptr);
      geo_fe.configure(geo_tag, geo_sf_tag);
      geo_fe.fill_Vandermonde(quadrature.get().coordinates(), quadrature.get().weights());
    }

    const mesh::PointSetTagExt sol_tag_ext(sol_tag, quad_order, mesh::CellTransform::NO_TRANS, 0u);
    common::PtrHandle<FEValues> sol_fe_val_ptr = sol_fe_values_map.std_region_data(sol_tag_ext);

    if (sol_fe_val_ptr.is_null())
    {
      sol_fe_val_ptr = sol_fe_values_map.create(sol_tag_ext);
      const mesh::PointSetTag quad_tag(sol_tag.elem_shape(), quad_order, quad_type);
      quadrature.change_type(quad_tag);

      const mesh::sf::SFTag sol_sf_tag(sol_tag.elem_shape(), sf_type, sol_tag.poly_order(),
                                       ModalBasis::Modal);
      FEValues &sol_fe = (*sol_fe_val_ptr);
      sol_fe.configure(sol_tag, sol_sf_tag);
      sol_fe.fill_Vandermonde(quadrature.get().coordinates(), quadrature.get().weights());
    } // if

    tmp_tag_pairs.insert(std::make_pair(geo_tag_ext, sol_tag_ext));
  }

  for (auto tag_pair : tmp_tag_pairs)
  {
    tag_combinations.push_back(tag_pair);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename QRule, Uint FacetDIM>
void fill_fe_value_pairs(
    typename mesh::Tria<MeshConfig> const &cell_topology,
    const typename result_of::dof_map_t<MeshConfig> &geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &sol_dofs, SFunc sf_type, PointSetID quad_type,
    const QRule &quad_order_rule,
    common::DataMap<mesh::PointSetTagExt, FEValues> &geo_fe_values_map,
    common::DataMap<mesh::PointSetTagExt, FEValues> &sol_fe_values_map,
    std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> &tag_combinations,
    const bool use_filter)
{
  geo_fe_values_map.clear();
  sol_fe_values_map.clear();
  tag_combinations.resize(0);

  std::set<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> tmp_tag_pairs;

  mesh::StdPointSet quadrature;

  for (Uint f = 0; f < cell_topology.active_skeleton_size(FacetDIM); ++f)
  {
    const mesh::TraceIncidences facet_block =
        cell_topology.active_skeleton_entry(FacetDIM, mesh::ActiveIdx(f));

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

    const Uint local_id_in_entity_L = facet_block.local_id(local_idx_L);
    const Uint local_id_in_entity_R = facet_block.local_id(local_idx_R);

    // Tags for left geometry and solution facet
    mesh::MeshEntity geo_facet_L = geo_dofs.active_cell(active_cell_id_L);
    geo_facet_L.local_transform(FacetDIM, local_id_in_entity_L);

    mesh::MeshEntity sol_facet_L = sol_dofs.active_cell(active_cell_id_L);
    sol_facet_L.local_transform(FacetDIM, local_id_in_entity_L);

    const mesh::PointSetTag geo_facet_tag_L = geo_facet_L.pt_set_id();
    const mesh::PointSetTag sol_facet_tag_L = sol_facet_L.pt_set_id();

    const Uint quad_order_L = quad_order_rule(geo_facet_tag_L, sol_facet_tag_L);

    // Tags for right geometry and solution facet
    mesh::MeshEntity geo_facet_R = geo_dofs.active_cell(active_cell_id_R);
    geo_facet_R.local_transform(FacetDIM, local_id_in_entity_R);

    mesh::MeshEntity sol_facet_R = sol_dofs.active_cell(active_cell_id_R);
    sol_facet_R.local_transform(FacetDIM, local_id_in_entity_R);

    const mesh::PointSetTag geo_facet_tag_R = geo_facet_R.pt_set_id();
    const mesh::PointSetTag sol_facet_tag_R = sol_facet_R.pt_set_id();

    const Uint quad_order_R = quad_order_rule(geo_facet_tag_R, sol_facet_tag_R);

    const Uint quad_order = std::max(quad_order_L, quad_order_R);

    // Check if pairs [geometry entity, solution entity] on left and right
    // sides of the facet are already present
    const mesh::EntityRealignCode local_p_tag_L = facet_block.permutation(local_idx_L).get().code();

    const mesh::EntityRealignCode local_p_tag_R = facet_block.permutation(local_idx_R).get().code();

    const mesh::PointSetTagExt geo_key_L(geo_facet_L.pt_set_id(), quad_order,
                                         local_p_tag_L.adapt_op_id(),
                                         local_p_tag_L.local_pos_in_parent());

    const mesh::PointSetTagExt sol_key_L(sol_facet_L.pt_set_id(), quad_order,
                                         local_p_tag_L.adapt_op_id(),
                                         local_p_tag_L.local_pos_in_parent());

    common::PtrHandle<FEValues> geo_fe_val_ptr_L = geo_fe_values_map.std_region_data(geo_key_L);

    if (geo_fe_val_ptr_L.is_null())
    {
      geo_fe_val_ptr_L = geo_fe_values_map.create(geo_key_L);

      const mesh::PointSetTag quad_tag_L(geo_facet_tag_L.elem_shape(), quad_order, quad_type);

      quadrature.change_type(quad_tag_L);

      const mesh::sf::SFTag geo_sf_tag_L(geo_facet_tag_L.elem_shape(), sf_type,
                                         geo_facet_tag_L.poly_order(), ModalBasis::Modal);
      FEValues &geo_fe_L = (*geo_fe_val_ptr_L);
      geo_fe_L.configure(geo_facet_tag_L, geo_sf_tag_L);

      const mesh::CellTransform adapt_op_id_L = geo_key_L.cell_transform_id();

      if (adapt_op_id_L == mesh::CellTransform::NO_TRANS)
      {
        geo_fe_L.fill_Vandermonde(quadrature.get().coordinates(local_id_in_entity_L),
                                  quadrature.get().weights(local_id_in_entity_L));
      }
      else
      {
        QuadratureAdaptTransformAlgoFactory::instance_type &quad_adapt_algo_fact =
            QuadratureAdaptTransformAlgoFactory::instance();
        const QuadratureAdaptTransformAlgoFactory::instance_type::const_product_base_ptr
            algo_transform = quad_adapt_algo_fact.create(adapt_op_id_L);

        math::DenseDMat<Real> transf_quad_coords;
        math::DenseDVec<Real> transf_quad_weights;

        algo_transform->compute_transformed_coords(
            quadrature.get().coordinates(local_id_in_entity_L), geo_key_L.local_id(),
            transf_quad_coords);
        algo_transform->compute_transformed_weights(quadrature.get().weights(local_id_in_entity_L),
                                                    geo_key_L.local_id(), transf_quad_weights);

        geo_fe_L.fill_Vandermonde(transf_quad_coords, transf_quad_weights, use_filter);
      }
    }

    common::PtrHandle<FEValues> sol_fe_val_ptr_L = sol_fe_values_map.std_region_data(sol_key_L);

    if (sol_fe_val_ptr_L.is_null())
    {
      sol_fe_val_ptr_L = sol_fe_values_map.create(sol_key_L);

      const mesh::PointSetTag quad_tag_L(sol_facet_tag_L.elem_shape(), quad_order, quad_type);

      quadrature.change_type(quad_tag_L);

      const mesh::sf::SFTag sol_sf_tag_L(sol_facet_tag_L.elem_shape(), sf_type,
                                         sol_facet_tag_L.poly_order(), ModalBasis::Modal);
      FEValues &sol_fe_L = (*sol_fe_val_ptr_L);
      sol_fe_L.configure(sol_facet_tag_L, sol_sf_tag_L);

      const mesh::CellTransform adapt_op_id_L = sol_key_L.cell_transform_id();

      if (adapt_op_id_L == mesh::CellTransform::NO_TRANS)
      {
        sol_fe_L.fill_Vandermonde(quadrature.get().coordinates(local_id_in_entity_L),
                                  quadrature.get().weights(local_id_in_entity_L));
      }
      else
      {
        QuadratureAdaptTransformAlgoFactory::instance_type &quad_adapt_algo_fact =
            QuadratureAdaptTransformAlgoFactory::instance();
        const QuadratureAdaptTransformAlgoFactory::instance_type::const_product_base_ptr
            algo_transform = quad_adapt_algo_fact.create(adapt_op_id_L);

        math::DenseDMat<Real> transf_quad_coords;
        math::DenseDVec<Real> transf_quad_weights;

        algo_transform->compute_transformed_coords(
            quadrature.get().coordinates(local_id_in_entity_L), sol_key_L.local_id(),
            transf_quad_coords);
        algo_transform->compute_transformed_weights(quadrature.get().weights(local_id_in_entity_L),
                                                    sol_key_L.local_id(), transf_quad_weights);

        sol_fe_L.fill_Vandermonde(transf_quad_coords, transf_quad_weights, use_filter);
      }
    }

    tmp_tag_pairs.insert(std::make_pair(geo_key_L, sol_key_L));

    if (facet_block.size() == 2)
    {
      const mesh::PointSetTagExt geo_key_R(geo_facet_R.pt_set_id(), quad_order,
                                           local_p_tag_R.adapt_op_id(),
                                           local_p_tag_R.local_pos_in_parent());

      const mesh::PointSetTagExt sol_key_R(sol_facet_R.pt_set_id(), quad_order,
                                           local_p_tag_R.adapt_op_id(),
                                           local_p_tag_R.local_pos_in_parent());

      common::PtrHandle<FEValues> geo_fe_val_ptr_R = geo_fe_values_map.std_region_data(geo_key_R);

      if (geo_fe_val_ptr_R.is_null())
      {
        geo_fe_val_ptr_R = geo_fe_values_map.create(geo_key_R);

        const mesh::PointSetTag quad_tag_R(geo_facet_tag_R.elem_shape(), quad_order, quad_type);

        quadrature.change_type(quad_tag_R);

        const mesh::sf::SFTag geo_sf_tag_R(geo_facet_tag_R.elem_shape(), sf_type,
                                           geo_facet_tag_R.poly_order(), ModalBasis::Modal);
        FEValues &geo_fe_R = (*geo_fe_val_ptr_R);
        geo_fe_R.configure(geo_facet_tag_R, geo_sf_tag_R);

        const mesh::CellTransform adapt_op_id_R = geo_key_R.cell_transform_id();

        if (adapt_op_id_R == mesh::CellTransform::NO_TRANS)
        {
          geo_fe_R.fill_Vandermonde(quadrature.get().coordinates(local_id_in_entity_R),
                                    quadrature.get().weights(local_id_in_entity_R));
        }
        else
        {
          QuadratureAdaptTransformAlgoFactory::instance_type &quad_adapt_algo_fact =
              QuadratureAdaptTransformAlgoFactory::instance();
          const QuadratureAdaptTransformAlgoFactory::instance_type::const_product_base_ptr
              algo_transform = quad_adapt_algo_fact.create(adapt_op_id_R);

          math::DenseDMat<Real> transf_quad_coords;
          math::DenseDVec<Real> transf_quad_weights;

          algo_transform->compute_transformed_coords(
              quadrature.get().coordinates(local_id_in_entity_R), geo_key_R.local_id(),
              transf_quad_coords);
          algo_transform->compute_transformed_weights(
              quadrature.get().weights(local_id_in_entity_R), geo_key_R.local_id(),
              transf_quad_weights);

          geo_fe_R.fill_Vandermonde(transf_quad_coords, transf_quad_weights, use_filter);
        }
      }

      common::PtrHandle<FEValues> sol_fe_val_ptr_R = sol_fe_values_map.std_region_data(sol_key_R);

      if (sol_fe_val_ptr_R.is_null())
      {
        sol_fe_val_ptr_R = sol_fe_values_map.create(sol_key_R);

        const mesh::PointSetTag quad_tag_R(sol_facet_tag_R.elem_shape(), quad_order, quad_type);

        quadrature.change_type(quad_tag_R);

        const mesh::sf::SFTag sol_sf_tag_R(sol_facet_tag_R.elem_shape(), sf_type,
                                           sol_facet_tag_R.poly_order(), ModalBasis::Modal);
        FEValues &sol_fe_R = (*sol_fe_val_ptr_R);
        sol_fe_R.configure(sol_facet_tag_R, sol_sf_tag_R);

        const mesh::CellTransform adapt_op_id_R = sol_key_R.cell_transform_id();

        if (adapt_op_id_R == mesh::CellTransform::NO_TRANS)
        {
          sol_fe_R.fill_Vandermonde(quadrature.get().coordinates(local_id_in_entity_R),
                                    quadrature.get().weights(local_id_in_entity_R));
        }
        else
        {
          QuadratureAdaptTransformAlgoFactory::instance_type &quad_adapt_algo_fact =
              QuadratureAdaptTransformAlgoFactory::instance();
          const QuadratureAdaptTransformAlgoFactory::instance_type::const_product_base_ptr
              algo_transform = quad_adapt_algo_fact.create(adapt_op_id_R);

          math::DenseDMat<Real> transf_quad_coords;
          math::DenseDVec<Real> transf_quad_weights;

          algo_transform->compute_transformed_coords(
              quadrature.get().coordinates(local_id_in_entity_R), sol_key_R.local_id(),
              transf_quad_coords);
          algo_transform->compute_transformed_weights(
              quadrature.get().weights(local_id_in_entity_R), sol_key_R.local_id(),
              transf_quad_weights);

          sol_fe_R.fill_Vandermonde(transf_quad_coords, transf_quad_weights, use_filter);
        }
      }

      tmp_tag_pairs.insert(std::make_pair(geo_key_R, sol_key_R));

    } // If size of this block is 2

  } // Loop over the mesh skeleton

  for (auto tag_pair : tmp_tag_pairs)
  {
    tag_combinations.push_back(tag_pair);
  }
}

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // namespace pdekit

#endif
