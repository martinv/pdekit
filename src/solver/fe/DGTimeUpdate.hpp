#ifndef PDEKIT_Solver_DG_Time_Update_hpp
#define PDEKIT_Solver_DG_Time_Update_hpp

#include <mutex>

#include "interpolation/GeometryMetric.hpp"
#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/point_set/StdPointSet.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

class DGTimeUpdate
{
  public:
  /// Default constructor
  DGTimeUpdate();

  /// Default destructor
  ~DGTimeUpdate();

  /// Prepare internal data to compute the time step
  template <typename MeshConfig>
  void setup(const mesh::Tria<MeshConfig> &mesh_in,
             const typename result_of::dof_map_t<MeshConfig> &sol_dofs);

  /// Get the time step mesh function, constant version
  const interpolation::ScalarMeshFunction<Real> &nodal_dual_volume() const;

  /// Get the time step mesh function, constant version
  const interpolation::ScalarMeshFunction<Real> &nodal_wave_speed() const;

  /// Compute local time step
  void compute_local_time_step(const Real CFL, interpolation::ScalarMeshFunction<Real> &time_step);

  /// Reset the wave speeds to zero
  void reset_wave_speeds();

  /// Accumulate in the wave speed data: add given value 'ws'
  /// to nodal wave speeds in node 'node_idx'
  void accumulate_nodal_wave_speed(const Uint node_idx, const Real ws);

  /// Accumulate wave speed data from a buffer. This method is thread-safe
  void accumulate_wave_speeds(const std::vector<std::tuple<Uint, Real>> &ws_buffer);

  /// Get the wave speed of node 'node_idx'
  Real wave_speed(const Uint node_idx) const;

  private:
  /// Mutex to make some methods thread-safe
  mutable std::mutex m_mutex;

  /// Volume of each cell in the mesh
  interpolation::ScalarMeshFunction<Real> m_dual_nodal_volume;

  /// Sum of wave speeds per node accumulated by all cells
  /// that touch given node
  interpolation::ScalarMeshFunction<Real> m_wave_speeds_per_node;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DGTimeUpdate::setup(const mesh::Tria<MeshConfig> &mesh_in,
                         const typename result_of::dof_map_t<MeshConfig> &sol_dofs)
{
  std::vector<mesh::DiscreteElemKey> element_types;

  typedef typename result_of::dof_map_t<MeshConfig> dof_map_type;

  for (const typename dof_map_type::const_dof_range_typed &dof_group :
       sol_dofs.all_active_dof_groups())
  {
    typename dof_map_type::const_dof_iterator_typed iter = dof_group.begin();

    const mesh::PointSetTag ref_topo_id = iter->mesh_entity().pt_set_id();
    const ElemShape shape               = ref_topo_id.elem_shape();
    const Uint order                    = ref_topo_id.poly_order();

    mesh::sf::SFTag const sf_type_tag =
        mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);

    mesh::StdPointSet quadrature;
    quadrature.change_type(shape, order, PointSetID::Gauss);

    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      const mesh::PointSetTagExt ref_topo_id_ext(ref_topo_id, P0, mesh::CellTransform::NO_TRANS,
                                                 ie);

      const mesh::PointSetTag qd_pts(shape, order, PointSetID::Gauss);
      const mesh::PointSetTagExt qd_pts_ext(qd_pts, P0, mesh::CellTransform::NO_TRANS, ie);
      const mesh::DiscreteElemKey elem_key(ref_topo_id_ext, sf_type_tag, qd_pts_ext);
      mesh::add_unique_discr_elem_key(element_types, elem_key);
    }
  }

  interpolation::GeometryCache<MeshConfig::GDIM> geo_cache;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric;

  geo_cache.allocate(element_types.cbegin(), element_types.cend(), 1u);
  geo_metric.allocate_buffer(element_types.cbegin(), element_types.cend(), 1u);

  m_dual_nodal_volume.resize(sol_dofs.nb_nodes());
  m_dual_nodal_volume.fill(0.0);

  Real cell_volume = 0.0;

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (const typename dof_map_type::const_dof_range_typed &dof_group :
       sol_dofs.all_active_dof_groups())
  {
    typename dof_map_type::const_dof_iterator_typed iter = dof_group.begin();
    for (; iter != dof_group.end(); ++iter)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_view = iter->tcell();
      const mesh::MeshEntity cell                         = iter->mesh_entity();
      const mesh::PointSetTag std_reg_pts                 = cell.pt_set_id();
      const ElemShape shape                               = std_reg_pts.elem_shape();
      const Uint order                                    = std_reg_pts.poly_order();

      mesh::sf::SFTag const sf_type_tag =
          mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);

      geo_cache.flush();
      geo_metric.empty_buffer();

      const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
          tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

      const mesh::PointSetTagExt std_reg_pts_ext(std_reg_pts, P0, mesh::CellTransform::NO_TRANS, 0);
      const mesh::PointSetTag qd_pts(shape, order, PointSetID::Gauss);
      const mesh::PointSetTagExt qd_pts_ext(qd_pts, P0, mesh::CellTransform::NO_TRANS, 0);
      const mesh::DiscreteElemKey elem_key(std_reg_pts_ext, sf_type_tag, qd_pts_ext);

      geo_cache.push_back_to_buffer(cell_coords, elem_key);
      geo_metric.evaluate(geo_cache, interpolation::RebuildMetricIndex{true});

      const typename interpolation::GeometryMetric<
          MeshConfig::GDIM, MeshConfig::TDIM>::cellwise_metric cell_met_values =
          geo_metric.cellwise_values(0);

      const Uint nb_qd_pts                     = cell_met_values.nb_qd_pts();
      const math::DenseDVec<Real> &w           = cell_met_values.pt_weights();
      const math::DenseConstVecView<Real> jdet = cell_met_values.jdet();

      cell_volume = 0.0;

      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        cell_volume += w[q] * jdet[q];
      }

      /*
      const Real cell_volume_part =
          topo_is_continuous ? (cell_volume / cell.nb_vert()) :
      cell_volume;
      */

      const Real cell_volume_part = cell_volume;
      // topo_is_continuous ? (cell_volume / cell.nb_vert()) :
      // cell_volume;

      for (Uint v = 0; v < cell.nb_vert(); ++v)
      {
        m_dual_nodal_volume[cell.vertex(v)] += cell_volume_part;
      }

    } // Loop over all cells of one type

  } // Loop over all cell type groups in dof handler

  /*
  Real tot_volume = 0.0;
  for (Uint i = 0; i < m_dual_nodal_volume.nb_entries(); ++i)
  {
    tot_volume += m_dual_nodal_volume[i];
  }

  std::cout << "Total volume of all cells in mesh = " << tot_volume <<
  std::endl;
  */

  m_wave_speeds_per_node.resize(sol_dofs.nb_nodes());
  m_wave_speeds_per_node.fill(0.0);
}

// ----------------------------------------------------------------------------

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
