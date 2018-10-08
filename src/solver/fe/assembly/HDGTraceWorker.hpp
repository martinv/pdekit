#ifndef PDEKIT_Solver_HDG_Trace_Worker_hpp
#define PDEKIT_Solver_HDG_Trace_Worker_hpp

/*
#include <memory>

#include "mesh/Mesh.hpp"
*/

#include "common/DataMap.hpp"
#include "interpolation/FEValues.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/containers/DofMap.hpp"
#include "mesh/point_set/StdPointSet.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

namespace detail
{

template <typename MeshConfig>
class HDGTraceWorker
{
  public:
  /// TYPEDEFS
  using tria_t    = typename mesh::Tria<MeshConfig>;
  using dof_map_t = typename pdekit::result_of::dof_map_t<MeshConfig>;

  /// Default constructor
  HDGTraceWorker();

  /// Default destructor
  ~HDGTraceWorker();

  /// Configure the solver for looping over cells
  template <typename QuadOrderRule>
  void setup_traces(const tria_t &input_mesh, const mesh::DofMap<MeshConfig> &sol_dofs,
                    const QuadOrderRule &quad_order_rule);

  private:
  enum
  {
    FACET_DIM = MeshConfig::TDIM - 1
  };

  /// TYPES AND TYPEDEFS
  using fe_values_t = interpolation::FEValues;

  /*
  using Vandermonde_mat_array = std::array<math::DenseDMat<Real>,
  MeshConfig::TDIM>;

  using fe_map_iterator = mesh::StdRegionDataMap<fe_values_t>::iterator;
  using sf_deriv_map_iterator = typename
  mesh::StdRegionDataMap<Vandermonde_mat_array>::iterator;
  */

  // mesh::StdRegionDataMap<fe_values_t> m_geo_facet_fe_map;
  common::DataMap<mesh::PointSetTagExt, fe_values_t> m_geo_facet_L_fe_map;
  common::DataMap<mesh::PointSetTagExt, fe_values_t> m_geo_facet_R_fe_map;

  // For each element type present in the mesh, this map holds
  // one 'FEValues', which has the Vandermonde matrix for shape functions
  // and their derivatives in REFERENCE space
  // mesh::StdRegionDataMap<fe_values_t> m_sol_facet_fe_map;
  common::DataMap<mesh::PointSetTagExt, fe_values_t> m_sol_facet_L_fe_map;
  common::DataMap<mesh::PointSetTagExt, fe_values_t> m_sol_facet_R_fe_map;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
HDGTraceWorker<MeshConfig>::HDGTraceWorker()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
HDGTraceWorker<MeshConfig>::~HDGTraceWorker()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename QuadOrderRule>
void HDGTraceWorker<MeshConfig>::setup_traces(const tria_t &input_mesh,
                                              const mesh::DofMap<MeshConfig> &sol_dofs,
                                              const QuadOrderRule &quad_order_rule)
{
  // Clear map for geometry space reference elements
  m_geo_facet_L_fe_map.clear();
  m_geo_facet_R_fe_map.clear();

  // Clear map for solution space reference elements
  m_sol_facet_L_fe_map.clear();
  m_sol_facet_R_fe_map.clear();

  mesh::StdPointSet facet_quad;

  for (Uint f = 0; f < input_mesh.active_skeleton_size(FACET_DIM); ++f)
  {
    const mesh::TraceIncidences facet =
        input_mesh.active_skeleton_entry(FACET_DIM, mesh::ActiveIdx(f));

    const Uint local_idx_L = 0;
    const Uint local_idx_R = (facet.size() == 2) ? 1 : 0;

    // Get the cell and facet on the left- and right-hand side of interface
    // in GEOMETRY space

    ///*** LEFT****
    const mesh::CellTopologyView<MeshConfig> tcell_L =
        input_mesh.cell(mesh::FlatIdx(facet.cell_id(local_idx_L)));
    const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();

    const mesh::PointSetTag geo_cell_L_tag = tcell_L.pt_set_id();

    const mesh::EntityRealignCode pcode_L = facet.permutation(local_idx_L).get().code();
    const mesh::PointSetTagExt geo_facet_L_key(
        tcell_L.pt_set_id(FACET_DIM, facet.local_id(local_idx_L)), P0, pcode_L.adapt_op_id(),
        pcode_L.local_pos_in_parent());

    ///*** RIGHT ****
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        input_mesh.cell(mesh::FlatIdx(facet.cell_id(local_idx_R)));
    const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

    const mesh::PointSetTag geo_cell_R_tag = tcell_R.pt_set_id();

    const mesh::EntityRealignCode pcode_R = facet.permutation(local_idx_R).get().code();
    const mesh::PointSetTagExt geo_facet_R_key(
        tcell_R.pt_set_id(FACET_DIM, facet.local_id(local_idx_R)), P0, pcode_R.adapt_op_id(),
        pcode_R.local_pos_in_parent());

    // Get the cell and facet on the left- and right-hand side of interface
    // in SOLUTION space

    ///*** LEFT****
    const mesh::MeshEntity sol_cell_L      = sol_dofs.active_cell(active_cell_id_L);
    const mesh::PointSetTag sol_cell_L_tag = sol_cell_L.pt_set_id();
    const mesh::MeshEntity sol_facet_L =
        sol_cell_L.sub_entity(FACET_DIM, facet.local_id(local_idx_L));
    const mesh::PointSetTagExt sol_facet_L_key(sol_facet_L.pt_set_id(), P0, pcode_L.adapt_op_id(),
                                               pcode_L.local_pos_in_parent());

    ///*** RIGHT ****
    const mesh::MeshEntity sol_cell_R      = sol_dofs.active_cell(active_cell_id_R);
    const mesh::PointSetTag sol_cell_R_tag = sol_cell_R.pt_set_id();
    const mesh::MeshEntity sol_facet_R =
        sol_cell_R.sub_entity(FACET_DIM, facet.local_id(local_idx_R));
    const mesh::PointSetTagExt sol_facet_R_key(sol_facet_R.pt_set_id(), P0, pcode_R.adapt_op_id(),
                                               pcode_R.local_pos_in_parent());

    // Compute minimum order for quadrature order
    // 2*(sol_order -1) because of the grad*grad operator

    const Uint required_quad_order_L = quad_order_rule(geo_cell_L_tag, sol_cell_L_tag);
    const Uint required_quad_order_R = quad_order_rule(geo_cell_R_tag, sol_cell_R_tag);
    const Uint required_quad_order   = std::max(required_quad_order_L, required_quad_order_R);

    // Create geometry and solution element values for facet
    // Attach the Vandermonde matrix to the pointer geo_facet_L_fe
    // *** LEFT facet *****
    common::PtrHandle<fe_values_t> geo_facet_L_fe =
        m_geo_facet_L_fe_map.std_region_data(geo_facet_L_key);
    if (geo_facet_L_fe.is_null())
    {
      geo_facet_L_fe = m_geo_facet_L_fe_map.create(geo_facet_L_key);
      const mesh::sf::SFTag geo_facet_sf_tag(geo_cell_L_tag.elem_shape(), SFunc::Lagrange,
                                             geo_cell_L_tag.poly_order(), ModalBasis::Modal);

      (*geo_facet_L_fe).configure(geo_cell_L_tag, geo_facet_sf_tag);
      const mesh::PointSetTag quad_tag =
          mesh::PointSetTag(geo_cell_L_tag.elem_shape(), required_quad_order, PointSetID::Gauss);
      facet_quad.change_type(quad_tag);
      (*geo_facet_L_fe)
          .fill_Vandermonde(facet_quad.get().coordinates(), facet_quad.get().weights());
      (*geo_facet_L_fe).print();
    }

    // *** RIGHT facet *****
    common::PtrHandle<fe_values_t> geo_facet_R_fe =
        m_geo_facet_R_fe_map.std_region_data(geo_facet_R_key);
    if (geo_facet_R_fe.is_null())
    {
      geo_facet_R_fe = m_geo_facet_R_fe_map.create(geo_facet_R_key);
      const mesh::sf::SFTag geo_facet_sf_tag(geo_cell_R_tag.elem_shape(), SFunc::Lagrange,
                                             geo_cell_R_tag.poly_order(), ModalBasis::Modal);

      (*geo_facet_R_fe).configure(geo_cell_R_tag, geo_facet_sf_tag);
      const mesh::PointSetTag quad_tag =
          mesh::PointSetTag(geo_cell_R_tag.elem_shape(), required_quad_order, PointSetID::Gauss);
      facet_quad.change_type(quad_tag);
      (*geo_facet_R_fe)
          .fill_Vandermonde(facet_quad.get().coordinates(), facet_quad.get().weights());
      (*geo_facet_R_fe).print();
    }
    // Create geometry and solution element values for facet
    // Attach the Vandermonde matrix to the pointer sol_facet_L_fe
    // *** LEFT facet *****
    common::PtrHandle<fe_values_t> sol_facet_L_fe =
        m_sol_facet_L_fe_map.std_region_data(sol_facet_L_key);
    if (sol_facet_L_fe.is_null())
    {
      sol_facet_L_fe = m_sol_facet_L_fe_map.create(sol_facet_L_key);
      const mesh::sf::SFTag sol_facet_sf_tag(sol_cell_L_tag.elem_shape(), SFunc::Lagrange,
                                             sol_cell_L_tag.poly_order(), ModalBasis::Modal);

      (*sol_facet_L_fe).configure(sol_cell_L_tag, sol_facet_sf_tag);
      const mesh::PointSetTag quad_tag =
          mesh::PointSetTag(sol_cell_L_tag.elem_shape(), required_quad_order, PointSetID::Gauss);
      facet_quad.change_type(quad_tag);
      (*sol_facet_L_fe)
          .fill_Vandermonde(facet_quad.get().coordinates(), facet_quad.get().weights());
      (*sol_facet_L_fe).print();
    }

    // Create geometry and solution element values for facet
    // Attach the Vandermonde matrix to the pointer sol_facet_R_fe
    // *** RIGHT facet *****
    common::PtrHandle<fe_values_t> sol_facet_R_fe =
        m_sol_facet_R_fe_map.std_region_data(sol_facet_R_key);
    if (sol_facet_R_fe.is_null())
    {
      sol_facet_R_fe = m_sol_facet_R_fe_map.create(sol_facet_R_key);
      const mesh::sf::SFTag sol_facet_sf_tag(sol_cell_R_tag.elem_shape(), SFunc::Lagrange,
                                             sol_cell_R_tag.poly_order(), ModalBasis::Modal);

      (*sol_facet_R_fe).configure(sol_cell_R_tag, sol_facet_sf_tag);
      const mesh::PointSetTag quad_tag =
          mesh::PointSetTag(sol_cell_R_tag.elem_shape(), required_quad_order, PointSetID::Gauss);
      facet_quad.change_type(quad_tag);
      (*sol_facet_R_fe)
          .fill_Vandermonde(facet_quad.get().coordinates(), facet_quad.get().weights());
      (*sol_facet_R_fe).print();
    }
  }
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
