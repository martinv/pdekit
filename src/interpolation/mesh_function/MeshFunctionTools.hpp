#ifndef PDEKIT_Interpolation_Mesh_Function_Tools_hpp
#define PDEKIT_Interpolation_Mesh_Function_Tools_hpp

#include <limits>
#include <queue>

#include "common/DataMap.hpp"
#include "interpolation/mesh_function/ScalarMeshFunction.hpp"

#include "common/ArrayView.hpp"
#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "mesh/Tria.hpp"
#include "mesh/std_region/PointSetTagExt.hpp"

namespace pdekit
{

namespace interpolation
{

class MeshFunctionTools
{
  public:
  /// Default constructor
  MeshFunctionTools();

  /// Default destructor
  ~MeshFunctionTools();

  /// Distribute the cell-wise input function to vertices so that
  /// a P1 piecewise-continuous interpolation is created
  /// The values in each P1 vertex are taken as maximum of all elements
  /// incident to this vertex
  template <typename MeshConfig>
  static void interp_from_cells_to_nodes_p1_max(
      mesh::Tria<MeshConfig> &input_mesh,
      typename result_of::dof_map_t<MeshConfig> const &dof_handler,
      interpolation::ScalarMeshFunction<Real> const &f_cells,
      interpolation::ScalarMeshFunction<Real> &f_nodes);

  template <typename MeshConfig>
  static void compute_cell_distance(
      const mesh::Tria<MeshConfig> &input_mesh,
      const common::ArrayView<const mesh::ActiveIdx, _1D, Uint> &level_zero_cells,
      ScalarMeshFunction<Uint> &distances);

  private:
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshFunctionTools::interp_from_cells_to_nodes_p1_max(
    mesh::Tria<MeshConfig> &input_mesh,
    typename result_of::dof_map_t<MeshConfig> const &dof_handler,
    interpolation::ScalarMeshFunction<Real> const &f_cells,
    interpolation::ScalarMeshFunction<Real> &f_nodes)
{
  // -------------------------------------------
  // Move the values from cells to cell vertices
  // -------------------------------------------
  common::BlockArray<std::tuple<Uint, Uint>, Uint> incident_p1_nodes;
  common::BlockArray<Uint, Uint> p1_node_to_cell;

  input_mesh.compute_node_to_cell_connectivity(incident_p1_nodes, p1_node_to_cell);

  f_nodes.resize(dof_handler.nb_nodes());
  f_nodes.fill(0.0);

  // STAGE I
  // Define function values each P1 node of the mesh as maximum over
  // all cellwise function values incident to this node
  for (Uint b = 0; b < incident_p1_nodes.nb_blocks(); ++b)
  {
    const common::ArrayView<const std::tuple<Uint, Uint>, _1D, Uint> inc_p1_nodes =
        incident_p1_nodes.const_block(b);

    Real max_f_value_inc_cells = 0.0;

    // Get the maximum function value from all cells incident
    // to given node
    for (Uint i = 0; i < inc_p1_nodes.size(); ++i)
    {
      const Uint active_cell_id = std::get<0>(inc_p1_nodes[i]);
      const mesh::MeshEntity active_sol_cell =
          dof_handler.active_cell(mesh::ActiveIdx(active_cell_id));

      max_f_value_inc_cells = std::max(max_f_value_inc_cells, f_cells[active_cell_id]);
    }

    for (Uint i = 0; i < inc_p1_nodes.size(); ++i)
    {
      const Uint active_cell_id               = std::get<0>(inc_p1_nodes[i]);
      const Uint local_p1_node_in_active_cell = std::get<1>(inc_p1_nodes[i]);

      const mesh::MeshEntity active_sol_cell =
          dof_handler.active_cell(mesh::ActiveIdx(active_cell_id));

      /*
      f_nodes[active_sol_cell.vertex(local_p1_node_in_active_cell)] =
          std::max(f_nodes[active_sol_cell.vertex(local_p1_node_in_active_cell)],
                   f_cells[active_cell_id]);
      */

      f_nodes[active_sol_cell.vertex(local_p1_node_in_active_cell)] = max_f_value_inc_cells;
    }
  }

  // STAGE II
  // Process hanging nodes: in case there's a facet incident to two refined
  // facets, the cellwise function value on the non-refined side is taken into
  // consideration by the P1 vertices of refined elements (i.e. by their P1
  // hanging nodes)

  std::vector<bool> is_local_hanging_p1_node_L;
  std::vector<bool> is_local_hanging_p1_node_R;

  const Uint FACET_DIM = MeshConfig::TDIM - 1;

  mesh::TraceEntityTuple facet_tuple;

  for (Uint f = 0; f < input_mesh.active_skeleton_size(FACET_DIM); ++f)
  {
    const mesh::TraceIncidences skeleton_block =
        input_mesh.active_skeleton_entry(MeshConfig::TDIM - 1, mesh::ActiveIdx(f));
    if (skeleton_block.size() == 2)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_L =
          input_mesh.cell(mesh::FlatIdx(skeleton_block.cell_id(0)));
      const mesh::CellTopologyView<MeshConfig> tcell_R =
          input_mesh.cell(mesh::FlatIdx(skeleton_block.cell_id(1)));

      const bool facet_is_active = (tcell_L.status() == mesh::EntityStatus::Active) &&
                                   (tcell_R.status() == mesh::EntityStatus::Active);

      if (facet_is_active)
      {
        const mesh::EntityDofRealign perm_L = skeleton_block.permutation(0);
        const mesh::EntityDofRealign perm_R = skeleton_block.permutation(1);

        const mesh::EntityRealignCode pcode_L = perm_L.get().code();
        const mesh::EntityRealignCode pcode_R = perm_R.get().code();

        if (((pcode_L.adapt_op_id() == mesh::CellTransform::NO_TRANS) &&
             (pcode_R.adapt_op_id() != mesh::CellTransform::NO_TRANS)) ||
            ((pcode_L.adapt_op_id() != mesh::CellTransform::NO_TRANS) &&
             (pcode_R.adapt_op_id() == mesh::CellTransform::NO_TRANS)))
        {

          const mesh::MeshEntity sol_cell_L =
              dof_handler.active_cell(mesh::ActiveIdx(tcell_L.active_idx()));
          const mesh::MeshEntity sol_facet_L =
              sol_cell_L.sub_entity(FACET_DIM, skeleton_block.local_id(0));

          const mesh::MeshEntity sol_cell_R =
              dof_handler.active_cell(mesh::ActiveIdx(tcell_R.active_idx()));
          const mesh::MeshEntity sol_facet_R =
              sol_cell_R.sub_entity(FACET_DIM, skeleton_block.local_id(1));

          /*
          std::cout << std::endl;
          std::cout <<
          "***************************************************" <<
          std::endl; std::cout << "Solution cell L = " << sol_cell_L
          << std::endl; std::cout << "Solution facet L = " <<
          sol_facet_L << std::endl; std::cout << "Solution cell R = "
          << sol_cell_R << std::endl; std::cout << "Solution facet R =
          " << sol_facet_R << std::endl;
          */

          const std::tuple<Uint, Uint> rel_refinement_levels = mesh::get_relative_refinement_levels(
              tcell_L.refinement_level(), tcell_R.refinement_level());

          facet_tuple.change_type(std::make_tuple(pcode_L, std::get<LEFT>(rel_refinement_levels),
                                                  pcode_R, std::get<RIGHT>(rel_refinement_levels)));

          const std::vector<std::pair<Uint, Uint>> incident_p1_nodes =
              facet_tuple.get().incident_p1_nodes();

          is_local_hanging_p1_node_L.resize(sol_facet_L.nb_vert());
          is_local_hanging_p1_node_L.assign(sol_facet_L.nb_vert(), true);

          is_local_hanging_p1_node_R.resize(sol_facet_R.nb_vert());
          is_local_hanging_p1_node_R.assign(sol_facet_R.nb_vert(), true);

          // std::cout << "Incident P1 nodes: " << std::endl;

          for (auto inc_p1_node_pair : incident_p1_nodes)
          {
            is_local_hanging_p1_node_L[inc_p1_node_pair.first]  = false;
            is_local_hanging_p1_node_R[inc_p1_node_pair.second] = false;
          }

          for (Uint v = 0; v < sol_facet_L.nb_vert(); ++v)
          {
            if (!sol_facet_L.vert_is_p1(v))
            {
              is_local_hanging_p1_node_L[v] = false;
            }
          }

          for (Uint v = 0; v < sol_facet_R.nb_vert(); ++v)
          {
            if (!sol_facet_R.vert_is_p1(v))
            {
              is_local_hanging_p1_node_R[v] = false;
            }
          }

          /*
          std::cout << "Hanging nodes L: " << std::endl;
          for (Uint i = 0; i < is_local_hanging_p1_node_L.size(); ++i)
          {
            if (is_local_hanging_p1_node_L[i])
            {
              std::cout << i << " ";
            }
          }
          std::cout << std::endl;

          std::cout << "Hanging nodes R: " << std::endl;
          for (Uint i = 0; i < is_local_hanging_p1_node_R.size(); ++i)
          {
            if (is_local_hanging_p1_node_R[i])
            {
              std::cout << i << " ";
            }
          }
          std::cout << std::endl << std::endl;
          */

          // const Real max_f_cell_LR =
          // std::max(f_cells[sol_cell_L.idx()],
          // f_cells[sol_cell_R.idx()]);

          // Case where LEFT entity is refined and only a part of the
          // RIGHT entity is incident to its left counterpart
          if ((pcode_L.adapt_op_id() == mesh::CellTransform::NO_TRANS) &&
              (pcode_R.adapt_op_id() != mesh::CellTransform::NO_TRANS))
          {
            Real face_value_R = 0.0;
            Uint nb_p1_vert   = 0;

            for (Uint v = 0; v < sol_facet_R.nb_vert(); ++v)
            {
              if (sol_facet_R.vert_is_p1(v))
              {
                face_value_R += f_nodes[sol_facet_R.vertex(v)];
                nb_p1_vert++;
              }
            }
            face_value_R /= nb_p1_vert;

            for (Uint v = 0; v < is_local_hanging_p1_node_L.size(); ++v)
            {
              if (is_local_hanging_p1_node_L[v])
              {
                f_nodes[sol_facet_L.vertex(v)] =
                    std::max(f_nodes[sol_facet_L.vertex(v)], face_value_R);
              }
            }
          }
          else
          {
            Real face_value_L = 0.0;
            Uint nb_p1_vert   = 0;

            for (Uint v = 0; v < sol_facet_L.nb_vert(); ++v)
            {
              if (sol_facet_L.vert_is_p1(v))
              {
                face_value_L += f_nodes[sol_facet_L.vertex(v)];
                nb_p1_vert++;
              }
            }
            face_value_L /= nb_p1_vert;

            for (Uint v = 0; v < is_local_hanging_p1_node_R.size(); ++v)
            {
              if (is_local_hanging_p1_node_R[v])
              {
                f_nodes[sol_facet_R.vertex(v)] =
                    std::max(f_nodes[sol_facet_R.vertex(v)], face_value_L);
              }
            }
          } // else
        }   // If one incident facet is adapted and the other is not
      }     // If facet is active
    }       // If skeleton block has size == 2
  }

  // STAGE III
  // Go again through all incident P1 vertices and make sure that they have a
  // common (maximum) value of f
  for (Uint b = 0; b < incident_p1_nodes.nb_blocks(); ++b)
  {
    const common::ArrayView<const std::tuple<Uint, Uint>, _1D, Uint> inc_p1_nodes =
        incident_p1_nodes.const_block(b);

    Real max_artif_visc = 0.0;

    // Get the maximum f from all cells incident
    // to given node
    for (Uint i = 0; i < inc_p1_nodes.size(); ++i)
    {
      const Uint active_cell_id               = std::get<0>(inc_p1_nodes[i]);
      const Uint local_p1_node_in_active_cell = std::get<1>(inc_p1_nodes[i]);

      const mesh::MeshEntity active_sol_cell =
          dof_handler.active_cell(mesh::ActiveIdx(active_cell_id));

      max_artif_visc =
          std::max(max_artif_visc, f_nodes[active_sol_cell.vertex(local_p1_node_in_active_cell)]);
    }

    for (Uint i = 0; i < inc_p1_nodes.size(); ++i)
    {
      const Uint active_cell_id               = std::get<0>(inc_p1_nodes[i]);
      const Uint local_p1_node_in_active_cell = std::get<1>(inc_p1_nodes[i]);

      const mesh::MeshEntity active_sol_cell =
          dof_handler.active_cell(mesh::ActiveIdx(active_cell_id));

      f_nodes[active_sol_cell.vertex(local_p1_node_in_active_cell)] =
          std::max(max_artif_visc, f_nodes[active_sol_cell.vertex(local_p1_node_in_active_cell)]);
    }
  }

  // Position 0: values after interpolation (nb. of HO nodes) x DIM
  // Position 1: interpolation matrix (nb. of HO nodes) x (nb. of P1 nodes)
  // Position 2: values before interpolation (nb. of P1 nodes) x DIM
  typedef std::array<math::DenseDMat<Real>, 3> interp_mat_array_type;

  common::DataMap<mesh::PointSetTagExt, interp_mat_array_type> coord_interp_mat_map;
  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < dof_handler.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = dof_handler.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity active_sol_cell      = dof_handler.active_cell(mesh::ActiveIdx(c));
    const mesh::PointSetTag active_sol_cell_tag = active_sol_cell.pt_set_id();
    const mesh::PointSetTagExt active_sol_cell_tag_ext(active_sol_cell_tag, P0,
                                                       mesh::CellTransform::NO_TRANS, 0u);
    common::PtrHandle<interp_mat_array_type> interp_mat_ptr =
        coord_interp_mat_map.std_region_data(active_sol_cell_tag_ext);

    if (interp_mat_ptr.is_null())
    {
      interp_mat_ptr = coord_interp_mat_map.create(active_sol_cell_tag_ext);
      mesh::StdRegion std_reg_ho(active_sol_cell.pt_set_id());

      interp_mat_array_type &interp_mats = (*interp_mat_ptr);
      interp_mats[0].resize(std_reg_ho.get().nb_nodes(), 1);
      interp_mats[1].resize(std_reg_ho.get().nb_nodes(), std_reg_ho.get().nb_p1_nodes());
      interp_mats[2].resize(std_reg_ho.get().nb_p1_nodes(), 1);

      mesh::sf::SFTag interp_sf_tag(active_sol_cell_tag.elem_shape(), SFunc::Lagrange, P1,
                                    ModalBasis::Modal);

      mesh::sf::ShapeFunction interp_sf(std::make_tuple(active_sol_cell_tag, interp_sf_tag));
      interp_sf.get().compute_ref_values(std_reg_ho.get().coordinates(), interp_mats[1]);
      // std::cout << interp_mats[1] << std::endl;
    }

    interp_mat_array_type &interp_mats = (*interp_mat_ptr);

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), active_sol_cell.pt_set_id(), tcell_view.coordinates());

    Uint pos = 0;
    for (Uint vert = 0; vert < active_sol_cell.nb_vert(); ++vert)
    {
      if (active_sol_cell.vert_is_p1(vert))
      {
        const math::DenseConstVecView<Real> node_coord = cell_coords.row_transpose(vert);

        interp_mats[2](pos, 0) = f_nodes[active_sol_cell.vertex(vert)];
        pos++;
      }
    }

    interp_mats[0] = interp_mats[1] * interp_mats[2];

    /*
    std::cout << "******************************************" << std::endl;
    std::cout << "Function values in P1 nodes: " << std::endl;
    std::cout << interp_mats[2] << std::endl;
    std::cout << "Function values in HO element: " << std::endl;
    std::cout << interp_mats[0] << std::endl;
    */

    for (Uint v = 0; v < active_sol_cell.nb_vert(); ++v)
    {
      if (!active_sol_cell.vert_is_p1(v))
      {
        f_nodes[active_sol_cell.vertex(v)] = interp_mats[0](v, 0);
      }
    }
  }
}
// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshFunctionTools::compute_cell_distance(
    const mesh::Tria<MeshConfig> &input_mesh,
    const common::ArrayView<const mesh::ActiveIdx, _1D, Uint> &level_zero_cells,
    ScalarMeshFunction<Uint> &distances)
{
  const Uint max_distance = std::numeric_limits<Uint>::max();

  distances.resize(input_mesh.nb_active_cells());
  distances.fill(max_distance);

  std::queue<Uint> cells_to_process;
  std::vector<bool> visited;
  visited.resize(input_mesh.nb_active_cells());
  visited.assign(input_mesh.nb_active_cells(), false);

  for (Uint c = 0; c < level_zero_cells.size(); ++c)
  {
    const mesh::ActiveIdx idx = level_zero_cells[c];
    distances[idx.id()]       = 0u;
    cells_to_process.push(idx.id());
  }

  while (!cells_to_process.empty())
  {
    const Uint curr_cell_id = cells_to_process.front();
    cells_to_process.pop();

    const mesh::CellTopologyView<MeshConfig> curr_cell =
        input_mesh.active_cell(mesh::ActiveIdx(curr_cell_id));
    if (!visited[curr_cell.active_idx().id()])
    {
      const std::vector<std::tuple<mesh::CellTopologyView<MeshConfig>, Uint, Uint>> neighbours =
          input_mesh.active_neighbours(curr_cell);

      for (auto neighb_item : neighbours)
      {
        const mesh::CellTopologyView<MeshConfig> neighb_cell = std::get<0>(neighb_item);
        if (!visited[neighb_cell.active_idx().id()])
        {
          const Uint new_dist = distances[curr_cell.active_idx().id()] + 1;
          if (new_dist < distances[neighb_cell.active_idx().id()])
          {
            distances[neighb_cell.active_idx().id()] = new_dist;
          }
          cells_to_process.push(neighb_cell.active_idx().id());
        }
      }
    }

    visited[curr_cell.active_idx().id()] = true;
  }
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
