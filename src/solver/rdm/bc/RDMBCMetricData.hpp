#ifndef RDM_Boundary_Condition_Metric_Data_hpp
#define RDM_Boundary_Condition_Metric_Data_hpp

#include "interpolation/FluxSpaceMetric.hpp"

namespace pdekit
{

namespace mesh
{
template <typename MeshConfg>
class Tria;
}

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, Uint BcDim>
struct RDMBCMetricData
{
  using f_space = interpolation::FunctionSpace<MeshConfig, BcDim>;

  /// METHODS

  /// Default constructor
  RDMBCMetricData();

  /// Default destructor
  ~RDMBCMetricData();

  /// Set the function space (reference elements)
  void allocate_cache(const Uint geo_cell_order, const f_space &geo_space, const f_space &sol_space,
                      const Uint max_nb_elem);

  /// Fill geometry cache
  template <typename SolBdryDofIterator>
  void fill_geo_cache(const mesh::Tria<MeshConfig> &tria,
                      const std::vector<common::IteratorRange<SolBdryDofIterator>> &sol_dofs_ranges,
                      const f_space &geo_space);

  /// DATA

  /// Buffer for geometrical data
  interpolation::GeometryCache<MeshConfig::GDIM> m_geo_cache;

  /// Buffer for solution data
  interpolation::SolutionCache m_sol_cache;

  /// Geometry metric
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, BcDim> m_geo_metric;

  /// Solution metric
  interpolation::SolutionSpaceMetric<MeshConfig, BcDim> m_sol_metric;

  /// Flux interpolation metric
  interpolation::FluxSpaceMetric<MeshConfig, Physics, BcDim> m_flux_metric;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
RDMBCMetricData<MeshConfig, Physics, BcDim>::RDMBCMetricData()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
RDMBCMetricData<MeshConfig, Physics, BcDim>::~RDMBCMetricData()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void RDMBCMetricData<MeshConfig, Physics, BcDim>::allocate_cache(const Uint geo_cell_order,
                                                                 const f_space &geo_space,
                                                                 const f_space &sol_space,
                                                                 const Uint max_nb_elem)
{
  m_geo_cache.allocate(geo_space.discrete_elements().cbegin(), geo_space.discrete_elements().cend(),
                       max_nb_elem);
  m_geo_metric.allocate_buffer(geo_space.discrete_elements().cbegin(),
                               geo_space.discrete_elements().cend(), max_nb_elem);

  m_sol_cache.allocate(sol_space.reference_elements().cbegin(),
                       sol_space.reference_elements().cend(), max_nb_elem, Physics::NEQ);
  m_sol_metric.allocate_buffer(sol_space.reference_elements().cbegin(),
                               sol_space.reference_elements().cend(), max_nb_elem, Physics::NEQ);

  m_flux_metric.allocate_buffer(SFunc::Lagrange, geo_cell_order,
                                sol_space.reference_elements().cbegin(),
                                sol_space.reference_elements().cend(), max_nb_elem);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
template <typename SolBdryDofIterator>
void RDMBCMetricData<MeshConfig, Physics, BcDim>::fill_geo_cache(
    const mesh::Tria<MeshConfig> &tria,
    const std::vector<common::IteratorRange<SolBdryDofIterator>> &sol_dofs_ranges,
    const f_space &geo_space)
{
  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  const typename f_space::quad_generator_fcn &geo_quad_gen = geo_space.quad_generator();

  SolBdryDofIterator sol_dof_iterator;

  for (auto const &sol_elem_range : sol_dofs_ranges)
  {
    for (sol_dof_iterator = sol_elem_range.begin(); sol_dof_iterator != sol_elem_range.end();
         ++sol_dof_iterator)
    {
      const mesh::CellGeometry<MeshConfig::GDIM> cell_coords = sol_dof_iterator->cell_geometry();
      const mesh::PointSetTag bdry_facet_tag                 = sol_dof_iterator->geo_pt_set_id();
      const ElemShape elem_shape                             = bdry_facet_tag.elem_shape();
      const Uint elem_order                                  = bdry_facet_tag.poly_order();

      const mesh::PointSetTagExt geo_pt_set_ext(bdry_facet_tag, P0, mesh::CellTransform::NO_TRANS,
                                                0u);
      const mesh::sf::SFTag geo_sf_tag    = sf_generator(elem_shape, elem_order);
      const mesh::PointSetTag quad_pt_set = geo_quad_gen(elem_shape, elem_order);
      const mesh::PointSetTagExt quad_pt_set_ext(quad_pt_set, P0, mesh::CellTransform::NO_TRANS,
                                                 0u);

      const mesh::DiscreteElemKey geo_key(geo_pt_set_ext, geo_sf_tag, quad_pt_set_ext);

      m_geo_cache.push_back_to_buffer(cell_coords, geo_key);
    }
  }
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
