#ifndef PDEKIT_Tools_Mesh_Inspector_hpp
#define PDEKIT_Tools_Mesh_Inspector_hpp

#include <iomanip>

#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/Tria.hpp"

namespace pdekit
{

namespace tools
{

template <typename MeshConfig>
class MeshInspector
{
  public:
  // Default constructor
  MeshInspector();

  // Default destructor
  ~MeshInspector();

  // Compute the volume of the domain covered by mesh by using
  // Green theorem and converting a volume integral (here in 2D):
  // 1/2 * int ( \nabla \cdot (x,y) ) into
  // surface integral 1/2 * \oint ( x*n_x + y*n_y )
  // This requires looping over all mesh edges and serves
  // as detector whether mesh edges are stored correctly
  static Real check_mesh_consistency(const mesh::Tria<MeshConfig> &mesh,
                                     const std::string &dof_handler_name, const Uint quad_order);

  // Generate colors for mesh elements. Each element has a unique color
  void color_mesh(const mesh::Tria<MeshConfig> &mesh,
                  interpolation::VectorMeshFunction<Real> &cell_colors) const;

  // Print all active facets in the mesh on the screen
  void print_active_facets(const mesh::Tria<MeshConfig> &mesh) const;

  // Print info about topology cell
  void print_active_topo_cell_info(const mesh::Tria<MeshConfig> &mesh,
                                   const Uint active_cell_idx) const;

  private:
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
MeshInspector<MeshConfig>::MeshInspector()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
MeshInspector<MeshConfig>::~MeshInspector()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Real MeshInspector<MeshConfig>::check_mesh_consistency(const mesh::Tria<MeshConfig> &mesh,
                                                       const std::string &dof_handler_name,
                                                       const Uint quad_order)
{
  const typename result_of::dof_map_t<MeshConfig> &cell_dofs =
      *(mesh.dof_storage(dof_handler_name));

  std::vector<Real> cell_volume(cell_dofs.nb_active_cells());
  cell_volume.assign(cell_dofs.nb_active_cells(), 0.0);

  interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1> skeleton_space;

  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto quad_generator = [quad_order](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, quad_order, PointSetID::Gauss);
  };

  skeleton_space.set_reference_fe_values_on_skeleton(mesh, cell_dofs, sf_generator, quad_generator);
  // skeleton_space.print_element_types();

  mesh::QuadraturePermutation quad_p_L, quad_p_R;

  /*
  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = cell_dofs.active_cell(c);
    std::cout << "cell [" << cell.idx() << "] = " << cell << std::endl;
  }
  */

  const Uint facet_dim = MeshConfig::TDIM - 1;

  interpolation::GeometryCache<MeshConfig::GDIM> geo_cache_L;
  interpolation::GeometryCache<MeshConfig::GDIM> geo_cache_R;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, facet_dim> geo_metric_L;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, facet_dim> geo_metric_R;

  geo_cache_L.allocate(skeleton_space.discrete_elements().cbegin(),
                       skeleton_space.discrete_elements().cend(), 1);
  geo_cache_R.allocate(skeleton_space.discrete_elements().cbegin(),
                       skeleton_space.discrete_elements().cend(), 1);

  geo_metric_L.allocate_buffer(skeleton_space.discrete_elements().cbegin(),
                               skeleton_space.discrete_elements().cend(), 1);
  geo_metric_R.allocate_buffer(skeleton_space.discrete_elements().cbegin(),
                               skeleton_space.discrete_elements().cend(), 1);

  mesh::adapt::LocalInterpolator loc_interpolator;
  bool mesh_passed_check = true;

  for (Uint f = 0; f < mesh.active_skeleton_size(); ++f)
  {
    const mesh::TraceIncidences facet_block =
        mesh.active_skeleton_entry(facet_dim, mesh::ActiveIdx(f));

    const Uint local_idx_L = 0;
    const Uint local_idx_R = (facet_block.size() == 2) ? 1 : 0;

    // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS OF
    // CELLS, NOT THEIR ACTIVE INDICES !!!
    const mesh::CellTopologyView<MeshConfig> tcell_L =
        mesh.cell(mesh::FlatIdx(facet_block.cell_id(local_idx_L)));
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        mesh.cell(mesh::FlatIdx(facet_block.cell_id(local_idx_R)));

    const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
    const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

    mesh::MeshEntity facet_L = cell_dofs.active_cell(active_cell_id_L);
    facet_L.local_transform(facet_dim, facet_block.local_id(local_idx_L));

    mesh::MeshEntity facet_R = cell_dofs.active_cell(active_cell_id_R);
    facet_R.local_transform(facet_dim, facet_block.local_id(local_idx_R));

    const mesh::EntityRealignCode local_sub_tag_L =
        facet_block.permutation(local_idx_L).get().code();

    const ElemShape facet_shape_L = facet_L.pt_set_id().elem_shape();
    const Uint facet_order_L      = facet_L.pt_set_id().poly_order();
    const mesh::PointSetTagExt pt_set_facet_ext_L(facet_L.pt_set_id(), P0,
                                                  local_sub_tag_L.adapt_op_id(),
                                                  local_sub_tag_L.local_pos_in_parent());
    const mesh::sf::SFTag sf_facet_L     = sf_generator(facet_shape_L, facet_order_L);
    const mesh::PointSetTag quad_facet_L = quad_generator(facet_shape_L, facet_order_L);
    const mesh::PointSetTagExt quad_facet_ext_L(quad_facet_L, P0, local_sub_tag_L.adapt_op_id(),
                                                local_sub_tag_L.local_pos_in_parent());

    const mesh::DiscreteElemKey key_L(pt_set_facet_ext_L, sf_facet_L, quad_facet_ext_L);

    const mesh::EntityRealignCode local_sub_tag_R =
        facet_block.permutation(local_idx_R).get().code();

    const ElemShape facet_shape_R = facet_R.pt_set_id().elem_shape();
    const Uint facet_order_R      = facet_R.pt_set_id().poly_order();
    const mesh::PointSetTagExt pt_set_facet_ext_R(facet_R.pt_set_id(), P0,
                                                  local_sub_tag_R.adapt_op_id(),
                                                  local_sub_tag_R.local_pos_in_parent());
    const mesh::sf::SFTag sf_facet_R     = sf_generator(facet_shape_R, facet_order_R);
    const mesh::PointSetTag quad_facet_R = quad_generator(facet_shape_R, facet_order_R);
    const mesh::PointSetTagExt quad_facet_ext_R(quad_facet_R, P0, local_sub_tag_R.adapt_op_id(),
                                                local_sub_tag_R.local_pos_in_parent());

    const mesh::DiscreteElemKey key_R(pt_set_facet_ext_R, sf_facet_R, quad_facet_ext_R);

    geo_cache_L.flush();
    geo_cache_R.flush();

    const mesh::PointSetTag tcell_facet_L =
        tcell_L.pt_set_id(facet_dim, facet_block.local_id(local_idx_L));

    const math::DenseConstMatView<Real> geo_facet_coords_L = loc_interpolator.transfer_coords(
        tcell_facet_L, facet_L.pt_set_id(),
        tcell_L.coordinates(facet_dim, facet_block.local_id(local_idx_L)));

    geo_cache_L.push_back_to_buffer(geo_facet_coords_L, key_L);
    geo_metric_L.empty_buffer();
    geo_metric_L.evaluate(geo_cache_L, interpolation::RebuildMetricIndex{true});

    const mesh::PointSetTag tcell_facet_R =
        tcell_R.pt_set_id(facet_dim, facet_block.local_id(local_idx_R));

    const math::DenseConstMatView<Real> geo_facet_coords_R = loc_interpolator.transfer_coords(
        tcell_facet_R, facet_R.pt_set_id(),
        tcell_R.coordinates(facet_dim, facet_block.local_id(local_idx_R)));

    geo_cache_R.push_back_to_buffer(geo_facet_coords_R, key_R);
    geo_metric_R.empty_buffer();
    geo_metric_R.evaluate(geo_cache_R, interpolation::RebuildMetricIndex{true});

    // std::cout << key_L << std::endl;
    const typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                                 facet_dim>::cellwise_metric facet_met_data_L =
        geo_metric_L.cellwise_values(0);
    const math::DenseConstMatView<Real> qcoords_L = facet_met_data_L.interpolated_coords();
    // std::cout << qcoords_L << std::endl;
    const math::DenseConstVecView<Real> jdet_L = facet_met_data_L.jdet();
    const math::DenseDVec<Real> &wq_L          = facet_met_data_L.pt_weights();

    // std::cout << facet_met_data_L.interpolated_coords() << std::endl;

    // std::cout << key_R << std::endl;
    const typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                                 facet_dim>::cellwise_metric facet_met_data_R =
        geo_metric_R.cellwise_values(0);
    const math::DenseConstMatView<Real> qcoords_R = facet_met_data_R.interpolated_coords();
    // std::cout << qcoords_R << std::endl;
    const math::DenseConstVecView<Real> jdet_R = facet_met_data_R.jdet();
    const math::DenseDVec<Real> &wq_R          = facet_met_data_R.pt_weights();

    if (qcoords_L.rows() != qcoords_R.rows())
    {
      std::cerr << "MeshInspector: quadratures between active cells " << active_cell_id_L.id()
                << " and " << active_cell_id_R.id() << " differ" << std::endl;
      mesh_passed_check = false;
    }

    const mesh::PointSetTag contour_quad_tag_L =
        mesh::PointSetTag(facet_L.pt_set_id().elem_shape(), quad_order, PointSetID::Gauss);

    const mesh::PointSetTag contour_quad_tag_R =
        mesh::PointSetTag(facet_R.pt_set_id().elem_shape(), quad_order, PointSetID::Gauss);

    quad_p_L.change_type(contour_quad_tag_L, facet_block.local_id(local_idx_L),
                         facet_block.permutation(local_idx_L).get().code());

    quad_p_R.change_type(contour_quad_tag_R, facet_block.local_id(local_idx_R),
                         facet_block.permutation(local_idx_R).get().code());

    math::DenseConstMatView<Real> const normals_L = facet_met_data_L.normals();
    math::DenseConstMatView<Real> const normals_R = facet_met_data_R.normals();

    // std::cout << " ================================================ " <<
    // std::endl;

    for (Uint q = 0; q < qcoords_L.rows(); ++q)
    {
      // (Local) index of left quadrature point
      const Uint qV_L = quad_p_L.get().vertex(q);
      // (Local) index of right quadrature point
      const Uint qV_R = quad_p_R.get().vertex(q);

      /*
      std::cout << "q = " << q << std::endl;
      std::cout << "jL = " << jdet_L[qV_L] << std::endl;
      std::cout << "wL = " << wq_L[qV_L] << std::endl;
      std::cout << "nL = " << normals_L.row(qV_L) << std::endl;
      std::cout << "XL = " << qcoords_L.row(qV_L) << std::endl;

      std::cout << " <<<< " << std::endl;

      std::cout << "jR = " << jdet_R[qV_R] << std::endl;
      std::cout << "wR = " << wq_R[qV_R] << std::endl;
      std::cout << "nR = " << normals_R.row(qV_R) << std::endl;
      std::cout << "XR = " << qcoords_R.row(qV_R) << std::endl;
      */

      const Real one_over_dim = 1. / MeshConfig::GDIM;

      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {

        // Accumulate to the cell volume on the left:
        cell_volume[active_cell_id_L.id()] +=
            one_over_dim * jdet_L[qV_L] * wq_L[qV_L] * (normals_L(qV_L, d) * qcoords_L(qV_L, d));

        if (facet_block.size() == 2)
        {
          // Accumulate to the cell volume on the right:
          cell_volume[active_cell_id_R.id()] +=
              one_over_dim * jdet_R[qV_R] * wq_R[qV_R] * (normals_R(qV_R, d) * qcoords_R(qV_R, d));
        }

        if (std::abs(qcoords_L(qV_L, d) - qcoords_R(qV_R, d)) > 1.e-10)
        {
          std::cerr << "MeshInspector: quadrature coords  between "
                       "active cells "
                    << active_cell_id_L.id() << " and " << active_cell_id_R.id() << " differ\n"
                    << "(qV_L, d) != (qV_R, d) with qV_L = " << qV_L << ", qV_R = " << qV_R
                    << ", d = " << d << std::endl;
          mesh_passed_check = false;
        }
      }
    }
  }

  const std::string verdict = mesh_passed_check ? "passed" : "failed";

  std::cout << "MeshInspector::mesh \"" << mesh.name() << "\" with DofMap \"" << dof_handler_name
            << "\" " << verdict << " sanity check." << std::endl;

  Real tot_volume = 0.0;
  for (Uint i = 0; i < cell_volume.size(); ++i)
  {
    tot_volume += cell_volume[i];
  }

  return tot_volume;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshInspector<MeshConfig>::color_mesh(
    const mesh::Tria<MeshConfig> &mesh, interpolation::VectorMeshFunction<Real> &cell_colors) const
{
  const typename result_of::dof_map_t<MeshConfig> &cell_dofs = mesh.topology().dof_storage();

  cell_colors.resize(1, cell_dofs.nb_active_cells());

  for (Uint ac = 0; ac < cell_dofs.nb_active_cells(); ++ac)
  {
    interpolation::VectorMeshFunction<Real>::entry_type cell_value = cell_colors.value(ac);
    cell_value[0]                                                  = ac;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshInspector<MeshConfig>::print_active_facets(const mesh::Tria<MeshConfig> &mesh) const
{
  const typename result_of::dof_map_t<MeshConfig> &cell_dofs = mesh.topology().dof_storage();

  const Uint FACET_DIM = MeshConfig::TDIM - 1;

  Uint min_node_id = std::numeric_limits<Uint>::max();
  Uint max_node_id = 0;

  std::vector<bool> cell_is_adjacent_to_facet(cell_dofs.nb_active_cells(), false);

  for (Uint f = 0; f < mesh.active_skeleton_size(); ++f)
  {
    // std::cout << "f = " << f << std::endl;

    const mesh::TraceIncidences facet_blk = mesh.active_skeleton_entry(MeshConfig::TDIM - 1, f);

    // Get the permutation sign of the first entity in the block (position
    // 0)
    const mesh::EntityRealignCode permutation_code = facet_blk.permutation(0).get().code();

    // If the entity on position 0 in incidence block needs to be flipped,
    // it is the entity on the rhs of this facet
    const Uint idx_in_block_L = (permutation_code.nb_flips() == 0) ? 0 : 1;

    // If the facet block has size 2, we have 2 entities forming the face
    // (i.e K+ and K-). In that case, idx_R is the 'other' index (choosing
    // from 0 and 1) than idx_L, and can be determined as idx_R = (idx_L +
    // 1) % 2 If the facet block has size 1, then idx_R has to be the same
    // as index left: equal to 0
    const Uint idx_in_block_R = (facet_blk.size() == 2) ? (idx_in_block_L + 1) % 2 : 0;

    // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS OF
    // CELLS, NOT THEIR ACTIVE INDICES !!! IT IS FOR THIS REASON THAT WE
    // FIRST NEED TO RETRIEVE TO TOPOLOGICAL CELLS AND GET THEIR ACTIVE
    // INDICES BEFORE GETTING THE CORRECT CELLS FROM THE DOF HANDLER
    const mesh::CellTopologyView<MeshConfig> tcell_L = mesh.cell(facet_blk.cell_id(idx_in_block_L));
    const mesh::CellTopologyView<MeshConfig> tcell_R = mesh.cell(facet_blk.cell_id(idx_in_block_R));

    const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
    const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

    // Get the cell and facet on the left- and right-hand side of interface
    // in GEOMETRY space
    const mesh::MeshEntity dof_cell_L = cell_dofs.active_cell(active_cell_id_L);
    const mesh::MeshEntity dof_cell_R = cell_dofs.active_cell(active_cell_id_R);

    const mesh::MeshEntity facet_L =
        dof_cell_L.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_L));
    const mesh::MeshEntity facet_R =
        dof_cell_R.sub_entity(FACET_DIM, facet_blk.local_id(idx_in_block_R));

    cell_is_adjacent_to_facet[active_cell_id_L.id()] = true;
    cell_is_adjacent_to_facet[active_cell_id_R.id()] = true;

    if (facet_blk.size() == 1)
    {
      std::cout << "Facet [" << f << "]  "
                << "{" << dof_cell_L.idx() << "/" << active_cell_id_L.id() << "} " << std::setw(15)
                << facet_L << std::endl;

      for (Uint v = 0; v < facet_L.nb_vert(); ++v)
      {
        min_node_id = std::min(min_node_id, facet_L.vertex(v));
        max_node_id = std::max(max_node_id, facet_L.vertex(v));
      }
    }
    else
    {
      std::cout << "Facet [" << f << "]  "
                << "{" << dof_cell_L.idx() << "/" << active_cell_id_L.id() << " | "
                << dof_cell_R.idx() << "/" << active_cell_id_R.id() << "} " << std::setw(11)
                << facet_L << " | " << facet_R << std::endl;

      for (Uint v = 0; v < facet_L.nb_vert(); ++v)
      {
        min_node_id = std::min(min_node_id, facet_L.vertex(v));
        max_node_id = std::max(max_node_id, facet_L.vertex(v));
      }

      for (Uint v = 0; v < facet_R.nb_vert(); ++v)
      {
        min_node_id = std::min(min_node_id, facet_R.vertex(v));
        max_node_id = std::max(max_node_id, facet_R.vertex(v));
      }
    }

  } // Loop over all active facets

  std::cout << "The following active cells were not listed as adjacent to any facet" << std::endl;
  for (Uint i = 0; i < cell_is_adjacent_to_facet.size(); ++i)
  {
    if (!cell_is_adjacent_to_facet[i])
    {
      std::cout << i << " ";
    }
  }
  std::cout << std::endl;
  std::cout << "Minimum vertex number found on facets: " << min_node_id << std::endl;
  std::cout << "Maximum vertex number found on facets: " << max_node_id << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshInspector<MeshConfig>::print_active_topo_cell_info(const mesh::Tria<MeshConfig> &mesh,
                                                            const Uint active_cell_idx) const
{
  const mesh::CellTopologyView<MeshConfig> tcell = mesh.active_cell(active_cell_idx);
  std::cout << "Topology cell " << tcell.linear_pos_idx() << "(absolute position)" << std::endl;
  std::cout << "Active cell idx = " << tcell.active_idx() << std::endl;
  const common::ArrayView<const Uint, _1D, Uint> neighbouring_faces = tcell.incident_facets();

  std::cout << "Storage positions of neighbouring facets = {";
  std::cout << neighbouring_faces[0];
  for (Uint n = 1; n < neighbouring_faces.size(); ++n)
  {
    std::cout << "," << neighbouring_faces[n];
  }
  std::cout << "}" << std::endl;

  std::cout << "Neighbouring cells:" << std::endl;
  for (Uint n = 0; n < neighbouring_faces.size(); ++n)
  {
    // const Uint facet_abs_idx = neighbouring_faces[n];
    const mesh::TraceIncidences facet_blk =
        mesh.skeleton_entry(MeshConfig::TDIM - 1, neighbouring_faces[n]);
    std::cout << "   ";
    for (Uint i = 0; i < facet_blk.size(); ++i)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_neighb = mesh.cell(facet_blk.cell_id(i));
      std::cout << "{ Active(" << tcell_neighb.active_idx() << "), Linear("
                << tcell_neighb.linear_pos_idx() << ") } ";
    }
    std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace tools

} // namespace pdekit

#endif
