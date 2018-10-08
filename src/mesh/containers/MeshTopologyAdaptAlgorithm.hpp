#ifndef PDEKIT_Mesh_Containers_MeshTopologyAdaptAlgorithm_hpp
#define PDEKIT_Mesh_Containers_MeshTopologyAdaptAlgorithm_hpp

#include "mesh/adaptation/GeometryAdapter.hpp"
#include "mesh/containers/MeshBoundary.hpp"
#include "mesh/containers/TriaCells.hpp"
#include "mesh/containers/TriaFacets.hpp"
#include "mesh/view/CellTopologyView.hpp"

namespace pdekit
{

namespace mesh
{

namespace internal
{

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class MeshTopologyAdaptAlgorithm
{
  public:
  /// @param mesh_cells          ... mesh containers for all cells in the mesh
  /// @param zero_level_skeleton ... all facets among level-zero cells (i.e.
  /// initial
  ///                                mesh before adaptation
  /// @param active_skeleton     ... container for mesh facets that are
  /// incident
  ///                                to active (leaf) cells
  static void rebuild_active_skeleton(const TriaCells<MeshConfig> &mesh_cells,
                                      const TriaFacets<MeshConfig> &zero_level_skeleton,
                                      TriaFacets<MeshConfig> &active_skeleton);

  /// @param cell_ids          ... numbers (ids) of cells to remove
  /// @param mesh_cells        ... cell container which holds cells that
  /// should be removed
  /// @param skeleton_facets   ... holds facets incident to mesh_cells
  /// @param mesh_boundary_set ... container of all mesh boundaries
  static void remove_cells(TriaCells<MeshConfig> &mesh_cells,
                           const TriaFacets<MeshConfig> &zero_level_skeleton,
                           TriaFacets<MeshConfig> &active_skeleton,
                           MeshBoundarySet<MeshConfig> &mesh_boundary_set,
                           std::vector<CellTransform> const &active_cell_ops);

  /// @brief Adapt the mesh allowing for 2:1 adaptation with hanging nodes
  static void adapt_w_hanging_nodes(TriaCells<MeshConfig> &mesh_cells,
                                    TriaFacets<MeshConfig> &skeleton_facets,
                                    MeshBoundarySet<MeshConfig> &mesh_boundary_set,
                                    std::vector<CellTransform> const &active_cell_ops);

  /// @brief Adapt the mesh using red-green refinement. Only suitable for
  /// triangular meshes in 2D.
  static void adapt_red_green(TriaCells<MeshConfig> &mesh_cells,
                              const TriaFacets<MeshConfig> &zero_level_skeleton,
                              TriaFacets<MeshConfig> &active_skeleton,
                              MeshBoundarySet<MeshConfig> &mesh_boundary_set,
                              std::vector<CellTransform> const &active_cell_ops);

  /// @brief Verify that the cell removal operations can be handled by the
  /// adaptation algorithm.
  /// @param mesh_cells        ... cell container which holds cells that
  /// should be removed
  /// @param skeleton_facets   ... holds facets incident to mesh_cells
  /// @param cell_ids          ... numbers (ids) of cells to remove
  static void prepare_for_cell_removal(TriaCells<MeshConfig> const &mesh_cells,
                                       TriaFacets<MeshConfig> const &skeleton_facets,
                                       std::vector<CellTransform> &adapt_op);

  /// @brief read and eventually modify adaptation operations so that the mesh
  /// changes that they
  ///        prescribe can be handled by the adaptation algorithm
  ///        This is for adaptation which allows hanging nodes with 2:1
  ///        refinement
  /// @param mesh cells      ... set of topological cells in the mesh
  /// @param skeleton facets ... all faces in the mesh
  /// @param adapt_op        ... vector of adaptation operations to apply. May
  /// be modified
  ///                            by this function
  static void prepare_cells_for_hanging_node_adapt(TriaCells<MeshConfig> const &mesh_cells,
                                                   TriaFacets<MeshConfig> const &skeleton_facets,
                                                   std::vector<CellTransform> &adapt_op);

  /// @brief read and eventually modify adaptation operations so that the mesh
  /// changes that they
  ///        prescribe can be handled by the adaptation algorithm
  ///        This is for red-green adaptation on triangular mesh in 2D. No
  ///        hanging nodes allowed.
  /// @param mesh cells      ... set of topological cells in the mesh
  /// @param skeleton facets ... all faces in the mesh
  /// @param adapt_op        ... vector of adaptation operations to apply. May
  /// be modified
  ///                            by this function
  /// @param colors          ... colors of cells in the mesh. Colors determine
  /// type
  ///                            of refinement that will be applied.
  static void generate_red_adapt_ops(const TriaCells<MeshConfig> &mesh_cells,
                                     const TriaFacets<MeshConfig> &skeleton_facets,
                                     std::vector<CellTransform> &adapt_op,
                                     std::vector<Uint> &colors);

  /// @brief read and eventually modify adaptation operations so that the mesh
  /// changes that they
  ///        prescribe can be handled by the adaptation algorithm
  ///        This is for red-green adaptation on triangular mesh in 2D. No
  ///        hanging nodes allowed.
  /// @param mesh cells      ... set of topological cells in the mesh
  /// @param skeleton facets ... all faces in the mesh
  /// @param adapt_op        ... vector of adaptation operations to apply. May
  /// be modified
  ///                            by this function
  /// @param colors          ... colors of cells in the mesh. Colors determine
  /// type
  ///                            of refinement that will be applied.
  static void generate_green_adapt_ops(const TriaCells<MeshConfig> &mesh_cells,
                                       const TriaFacets<MeshConfig> &skeleton_facets,
                                       std::vector<CellTransform> &adapt_op);

  /// @brief generate vector of operations such that all leaf cells which were
  /// obtained
  ///        by green refinement will be removed
  /// @param mesh_cells      ... set of topological cells in the mesh
  /// @param skeleton_facets ... all faces in the mesh
  /// @param adapt_op        ... vector of operations for all active cells
  /// that are defined
  ///                            so that all leaf cells due to green
  ///                            refinement would be removed
  static void generate_ops_to_remove_green_leafs(const TriaCells<MeshConfig> &mesh_cells,
                                                 const TriaFacets<MeshConfig> &skeleton_facets,
                                                 std::vector<CellTransform> &adapt_op);

  /// @brief Return a vector of all topological cells that are neighbours of
  /// given cell
  ///        and are active
  /// @param mesh_cells      ... set of all topological cells in the mesh
  /// @param skeleton_facets ... set of all faces in the mesh
  /// @param tcell           ... a cell whose neighbours we're looking for
  /// @return                ....a vector of tuples, where each tuple stores
  ///                            the following information:
  ///                            <neighbour cell,
  ///                             local facet id in neighbour incident to this
  ///                             cell, local facet id in this cell>
  static const std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> active_neighbours(
      const TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &skeleton_facets,
      const CellTopologyView<MeshConfig> tcell);

  /// @brief add all cells touching given cell that are active to a vector.
  ///        Candidate cells that might be added are listed in TraceIncidences
  ///        parameter. This function does not remove any existing (previously
  ///        added) cells from the vector of adjacent cells.
  /// @param mesh_cells      ... set of all topological cells in the mesh
  /// @param skeleton_facets ... set of all faces in the mesh
  /// @param center_cell     ... a cell whose neighbours we're looking for
  /// @param facet_block     ... holds ids of all incidenc cells that might be
  /// eventually added
  /// @param adjacent_cells  ... vector of active cells adjacent to
  /// 'center_cell'
  ///                            The tuple stores the following information:
  ///                            <neighbour cell,
  ///                             local facet id in neighbour incident to this
  ///                             cell, local facet id in this cell>
  static void add_active_adjacent_cells(
      const TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &skeleton_facets,
      const CellTopologyView<MeshConfig> center_cell, const TraceIncidences &facet_block,
      std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> &adjacent_cells);

  private:
  /// TYPEDEFS
  using const_skeleton_iterator = typename internal::TriaFacets<MeshConfig>::const_iterator;

  /// METHODS

  static void check_mesh_after_topo_change(TriaCells<MeshConfig> &mesh_cells,
                                           TriaFacets<MeshConfig> &skeleton_facets);

  /// Preserve only facets with given status
  /// @param mesh_cells      ... container which holds information the status
  /// of all cells
  /// @param skeleton_facets ... container with cell incidences - knows which
  /// cells share faces
  /// @param filter_status   ... reference entity status. Only facets that
  /// refer exclusively
  ///                            to cells with this status will be preserved.
  ///                            All other facets will be removed
  static void filter_facets(const TriaCells<MeshConfig> &mesh_cells,
                            TriaFacets<MeshConfig> &skeleton_facets,
                            const EntityStatus filter_status);

  static void ensure_local_facet_numbering(TriaCells<MeshConfig> &mesh_cells,
                                           const TriaFacets<MeshConfig> &skeleton_facets);

  static Uint get_local_id_of_child_cell(
      const std::vector<CellTopologyView<MeshConfig>> &all_children,
      const CellTopologyView<MeshConfig> &child_cell);

  static bool check_mesh_h_grading(const TriaCells<MeshConfig> &mesh_cells,
                                   const TriaFacets<MeshConfig> &skeleton_facets,
                                   const std::vector<CellTransform> &adapt_op,
                                   const Uint max_level_diff, const bool verbose = false);

  private:
  static void generate_new_facet_offsets_from_buffer(
      const internal::TriaFacets<MeshConfig> &facet_buffer, const Uint nb_old_cells_total,
      const Uint nb_cells_inserted, std::vector<Uint> &new_cell_facet_offsets,
      std::vector<Uint> &new_incident_facet_ids);
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::rebuild_active_skeleton(
    const TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &zero_level_skeleton,
    TriaFacets<MeshConfig> &active_skeleton)
{
  // [ cell id L, local id L, perm. code L, level L,
  //   cell id R, local id R, perm. code R, level R ]
  typedef std::tuple<Uint, Uint, EntityRealignCode, Uint, Uint, Uint, EntityRealignCode, Uint>
      interior_facet_entry_type;

  // [ cell id L, local id L, perm. code L]
  typedef std::tuple<Uint, Uint, EntityRealignCode> boundary_facet_entry_type;

  constexpr Uint CELL_ID_L   = 0;
  constexpr Uint LOCAL_ID_L  = 1;
  constexpr Uint PERM_CODE_L = 2;
  constexpr Uint LEVEL_L     = 3;

  constexpr Uint CELL_ID_R   = 4;
  constexpr Uint LOCAL_ID_R  = 5;
  constexpr Uint PERM_CODE_R = 6;
  constexpr Uint LEVEL_R     = 7;

  TraceEntityTuple local_trace_topology;

  CellTransform topology_input_transform_L;
  CellTransform topology_input_transform_R;
  std::vector<IncidenceEntry> child_incidences_L;
  std::vector<IncidenceEntry> child_incidences_R;

  std::vector<std::vector<IncidenceEntry>> combined_child_incidences;
  std::vector<std::vector<EntityDofRealign>> combined_child_permutations;

  std::queue<interior_facet_entry_type> interior_facets_queue;
  std::queue<boundary_facet_entry_type> boundary_facets_queue;

  std::vector<bool> cell_processed(mesh_cells.nb_all_cells(), false);

  /*
  std::cout << "Mesh cells before rebuilding edges:" << std::endl;
  std::cout << mesh_cells << std::endl;

  std::cout << "Listing cell children before rebuilding edges:" << std::endl;
  for (Uint c = 0; c < mesh_cells.nb_all_cells(); ++c)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.cell(c);
    if (tcell.nb_children() == 0)
    {
      std::cout << "Topology cell " << tcell.linear_pos_idx() << " does not
  have children"
                << std::endl;
    }
    else
    {
      const std::vector<CellTopologyView<MeshConfig>> children = tcell.children();
      std::cout << "Topology cell " << tcell.linear_pos_idx() << " has
  children"; for (Uint i = 0; i < children.size(); ++i)
      {
        const CellTopologyView<MeshConfig> &child = children[i];
        std::cout << " " << child.linear_pos_idx();
      }
      std::cout << std::endl;
    }
  }
  */

  // Clear the existing contents of the active skeleton
  active_skeleton.clear();

  for (Uint f = 0; f < zero_level_skeleton.nb_all_facets(); ++f)
  {
    const TraceIncidences facet_block = zero_level_skeleton.facet_data(FlatIdx(f));

    /*
    std::cout << "\n************************************" << std::endl;
    std::cout << "Facet block:\n" << facet_block << std::endl;
    */

    // --------------------------------------------------------------
    // BOUNDARY FACET
    // --------------------------------------------------------------
    if (facet_block.size() == 1)
    {
      const CellTopologyView<MeshConfig> tcell0_L =
          mesh_cells.cell(FlatIdx(facet_block.cell_id(LEFT)));
      if (tcell0_L.status() == EntityStatus::Active)
      {
        // std::cout << "Active facet on boundary\n" << facet_block <<
        // std::endl;

        /*
        combined_child_incidences.resize(1);
        combined_child_incidences[0].resize(1);
        combined_child_permutations.resize(1);
        combined_child_permutations[0].resize(1);

        combined_child_incidences[0][0] =
        facet_block.incidence_entry(0);
        combined_child_permutations[0][0] = facet_block.permutation(0);

        active_skeleton.add_facet(combined_child_incidences[0],
        combined_child_permutations[0]);
        */
        active_skeleton.add_facet({facet_block.incidence_entry(0)}, {facet_block.permutation(0)});
      }
      else
      {
        // std::cout << "Processing boudary facet\n" << facet_block <<
        // std::endl;

        const boundary_facet_entry_type init_entry =
            std::make_tuple(facet_block.cell_id(LEFT), facet_block.local_id(LEFT),
                            facet_block.permutation(LEFT).get().code());

        boundary_facets_queue.push(init_entry);

        while (!boundary_facets_queue.empty())
        {
          const boundary_facet_entry_type facet_entry = boundary_facets_queue.front();
          boundary_facets_queue.pop();

          const CellTopologyView<MeshConfig> tcell_L =
              mesh_cells.cell(FlatIdx(std::get<CELL_ID_L>(facet_entry)));

          if (tcell_L.nb_children() == 0)
          {
            /*
            std::cout << "Splitting generated boundary facet [" <<
            tcell_L.linear_pos_idx() << ","
                      << std::get<LOCAL_ID_L>(facet_entry) << ","
                      <<
            std::get<PERM_CODE_L>(facet_entry).as_string() << "]"
            << std::endl;
            */

            combined_child_incidences.resize(1);
            combined_child_incidences[0].resize(1);
            combined_child_permutations.resize(1);
            combined_child_permutations[0].resize(1);

            const Uint local_id_tcell_L = std::get<LOCAL_ID_L>(facet_entry);

            combined_child_incidences[0][0] =
                IncidenceEntry(tcell_L.linear_pos_idx().id(), local_id_tcell_L);

            const StdRegion cell_type_L = tcell_L.std_region();
            const PointSetTag facet_tag =
                (*cell_type_L.get().elem_entity(MeshConfig::TDIM - 1, local_id_tcell_L))
                    .pt_set_id();

            // Make sure that the transformation of leaf boundary
            // face is 'DO_NOTHING'
            const EntityRealignCode old_p_code = std::get<PERM_CODE_L>(facet_entry);

            const EntityRealignCode new_p_code(
                old_p_code.elem_shape(), CellTransform::NO_TRANS, old_p_code.local_pos_in_parent(),
                old_p_code.parent_shape(), old_p_code.nb_flips(), old_p_code.nb_rotations());

            /*
            const EntityPermutation p(
                std::make_pair(facet_tag,
            std::get<PERM_CODE_L>(facet_entry)));
            */
            const EntityDofRealign p(std::make_pair(facet_tag, new_p_code));

            combined_child_permutations[0][0] = p;

            active_skeleton.add_facet(combined_child_incidences[0], combined_child_permutations[0]);

            /*
            if (p.get().code().adapt_op_id() !=
            CellTransform::DO_NOTHING)
            {
              std::cout << "Warning: on facet " <<
            combined_child_incidences[0][0].cell_idx
                        << " (" <<
            combined_child_incidences[0][0].local_id <<
            ")" << std::endl;
            }
            */
          }
          else
          {
            const std::vector<CellTopologyView<MeshConfig>> children_L = tcell_L.children();

            const adapt::CellAdaptOp cell_adapt_op_L = tcell_L.cell_adapt_op();
            const TraceIncidences parent_to_child_on_facet_L =
                cell_adapt_op_L.get().parent_child_incidences_on_facet(
                    std::get<LOCAL_ID_L>(facet_entry));

            for (Uint i = 0; i < parent_to_child_on_facet_L.size(); ++i)
            {
              const boundary_facet_entry_type new_bfacet = std::make_tuple(
                  children_L[parent_to_child_on_facet_L.cell_id(i)].linear_pos_idx().id(),
                  parent_to_child_on_facet_L.local_id(i),
                  parent_to_child_on_facet_L.permutation(i).get().code());

              boundary_facets_queue.push(new_bfacet);
            }
          }

        } // while boundary facet queue is not empty
      }
    }
    // --------------------------------------------------------------
    // INTERIOR FACET
    // --------------------------------------------------------------
    else if (facet_block.size() == 2)
    {
      const CellTopologyView<MeshConfig> tcell0_L =
          mesh_cells.cell(FlatIdx(facet_block.cell_id(LEFT)));
      const CellTopologyView<MeshConfig> tcell0_R =
          mesh_cells.cell(FlatIdx(facet_block.cell_id(RIGHT)));

      if ((tcell0_L.status() == EntityStatus::Active) &&
          (tcell0_R.status() == EntityStatus::Active))
      {
        // std::cout << "Left and right neighbours are active" <<
        // std::endl;

        combined_child_incidences.resize(1);
        combined_child_incidences[0].resize(facet_block.size());
        combined_child_permutations.resize(1);
        combined_child_permutations[0].resize(facet_block.size());

        for (Uint i = 0; i < facet_block.size(); ++i)
        {
          combined_child_incidences[0][i]   = facet_block.incidence_entry(i);
          combined_child_permutations[0][i] = facet_block.permutation(i);
        }

        active_skeleton.add_facet(combined_child_incidences[0], combined_child_permutations[0]);
      }
      else
      {
        /*
        std::cout << "  One or both neighbours are not active" <<
        std::endl; std::cout << "  Status L = " << tcell0_L.status() <<
        std::endl; std::cout << "  Status R = " << tcell0_R.status() <<
        std::endl;
        */

        const interior_facet_entry_type init_entry = std::make_tuple(
            facet_block.cell_id(LEFT), facet_block.local_id(LEFT),
            facet_block.permutation(LEFT).get().code(), 0, facet_block.cell_id(RIGHT),
            facet_block.local_id(RIGHT), facet_block.permutation(RIGHT).get().code(), 0);

        if (!interior_facets_queue.empty())
        {
          std::cerr << "MeshTopologyAdaptAlgorithm::rebuild_active_"
                       "skeleton:: "
                       "error, the \n"
                    << "facet recursion queue is not empty as it should be" << std::endl;
        }

        interior_facets_queue.push(init_entry);

        while (!interior_facets_queue.empty())
        {
          const interior_facet_entry_type facet_entry = interior_facets_queue.front();
          interior_facets_queue.pop();

          /*
          std::cout << "Processing entry {" <<
          std::get<CELL_ID_L>(facet_entry)
          << ","
                    << std::get<LOCAL_ID_L>(facet_entry) << ","
                    << std::get<PERM_CODE_L>(facet_entry).as_string()
          << ","
                    << std::get<LEVEL_L>(facet_entry) << "}\n {"
                    << std::get<CELL_ID_R>(facet_entry) << "," <<
          std::get<LOCAL_ID_R>(facet_entry)
                    << "," <<
          std::get<PERM_CODE_R>(facet_entry).as_string() <<
          ","
                    << std::get<LEVEL_R>(facet_entry) << "}" <<
          std::endl;
          */

          local_trace_topology.change_type(
              std::get<PERM_CODE_L>(facet_entry), std::get<LEVEL_L>(facet_entry),
              std::get<PERM_CODE_R>(facet_entry), std::get<LEVEL_R>(facet_entry));

          const CellTopologyView<MeshConfig> tcell_L =
              mesh_cells.cell(FlatIdx(std::get<CELL_ID_L>(facet_entry)));
          const CellTopologyView<MeshConfig> tcell_R =
              mesh_cells.cell(FlatIdx(std::get<CELL_ID_R>(facet_entry)));

          const adapt::CellAdaptOp cell_adapt_op_L = tcell_L.cell_adapt_op();
          const adapt::CellAdaptOp cell_adapt_op_R = tcell_R.cell_adapt_op();

          /*
          std::cout << "     Cell adapt op L = "
                    <<
          cell_adapt_op_L.get().cell_adapt_op_tag().adapt_op_id()
          << std::endl; std::cout << "     Cell adapt op R = "
                    <<
          cell_adapt_op_R.get().cell_adapt_op_tag().adapt_op_id()
          << std::endl;
          */

          // Push all interior facets of left cell in the queue
          if (!cell_processed[tcell_L.linear_pos_idx().id()])
          {
            if (tcell_L.nb_children() > 0)
            {
              const std::vector<CellTopologyView<MeshConfig>> children_L = tcell_L.children();

              for (Uint f_L = 0; f_L < cell_adapt_op_L.get().nb_internal_child_facets(); ++f_L)
              {
                const TraceIncidences facet = cell_adapt_op_L.get().internal_child_facet(f_L);

                const FlatIdx child_id_L = children_L[facet.cell_id(LEFT)].linear_pos_idx();
                const FlatIdx child_id_R = children_L[facet.cell_id(RIGHT)].linear_pos_idx();

                const interior_facet_entry_type init_entry = std::make_tuple(
                    child_id_L.id(), facet.local_id(LEFT), facet.permutation(LEFT).get().code(), 0,
                    child_id_R.id(), facet.local_id(RIGHT), facet.permutation(RIGHT).get().code(),
                    0);

                interior_facets_queue.push(init_entry);
              }
            }
            cell_processed[tcell_L.linear_pos_idx().id()] = true;
          }

          // Push all interior facets of right cell in the queue
          if (!cell_processed[tcell_R.linear_pos_idx().id()])
          {
            if (tcell_R.nb_children() > 0)
            {
              const std::vector<CellTopologyView<MeshConfig>> children_R = tcell_R.children();

              for (Uint f_R = 0; f_R < cell_adapt_op_R.get().nb_internal_child_facets(); ++f_R)
              {
                const TraceIncidences facet = cell_adapt_op_R.get().internal_child_facet(f_R);

                const FlatIdx child_id_L = children_R[facet.cell_id(LEFT)].linear_pos_idx();
                const FlatIdx child_id_R = children_R[facet.cell_id(RIGHT)].linear_pos_idx();

                const interior_facet_entry_type init_entry = std::make_tuple(
                    child_id_L.id(), facet.local_id(LEFT), facet.permutation(LEFT).get().code(), 0,
                    child_id_R.id(), facet.local_id(RIGHT), facet.permutation(RIGHT).get().code(),
                    0);

                interior_facets_queue.push(init_entry);
              }
            }
            cell_processed[tcell_R.linear_pos_idx().id()] = true;
          }

          const CellTransform cell_transform_L =
              cell_adapt_op_L.get().cell_adapt_op_tag().adapt_op_id();
          const CellTransform cell_transform_R =
              cell_adapt_op_R.get().cell_adapt_op_tag().adapt_op_id();

          // --------------------------------------------------------
          // Generate local topology data on the left side of facet
          // --------------------------------------------------------
          const Uint local_id_in_parent_L = std::get<LOCAL_ID_L>(facet_entry);
          child_incidences_L.resize(0);
          if (cell_transform_L == CellTransform::NO_TRANS)
          {
            topology_input_transform_L = CellTransform::NO_TRANS;
            child_incidences_L.push_back(
                IncidenceEntry(std::get<CELL_ID_L>(facet_entry), local_id_in_parent_L));

            /*
            std::cout << "\n     Left neighbour " <<
            tcell_L.linear_pos_idx()
                      << " does not have children" << std::endl;
            */
          }
          else
          {
            topology_input_transform_L =
                cell_adapt_op_L.get().parent_facet_adapt_op_id(local_id_in_parent_L);

            const std::vector<CellTopologyView<MeshConfig>> children_L = tcell_L.children();
            const TraceIncidences parent_to_child_on_facet_L =
                cell_adapt_op_L.get().parent_child_incidences_on_facet(local_id_in_parent_L);

            child_incidences_L.resize(0);
            for (Uint i = 0; i < parent_to_child_on_facet_L.size(); ++i)
            {
              child_incidences_L.push_back(IncidenceEntry(
                  children_L[parent_to_child_on_facet_L.cell_id(i)].linear_pos_idx().id(),
                  parent_to_child_on_facet_L.local_id(i)));
            }
            /*
            std::cout << "\n     Left neighbour " <<
            tcell_L.linear_pos_idx() << " has children "; for (Uint
            i = 0; i < child_incidences_L.size();
            ++i)
            {
              std::cout << " [" << child_incidences_L[i].cell_idx <<
            ","
                        << child_incidences_L[i].local_id << "]";
            }
            std::cout << std::endl;
            */
          }

          // --------------------------------------------------------
          // Generate local topology data on the right side of facet
          // --------------------------------------------------------
          const Uint local_id_in_parent_R = std::get<LOCAL_ID_R>(facet_entry);
          child_incidences_R.resize(0);
          if (cell_transform_R == CellTransform::NO_TRANS)
          {
            topology_input_transform_R = CellTransform::NO_TRANS;
            child_incidences_R.push_back(
                IncidenceEntry(std::get<CELL_ID_R>(facet_entry), local_id_in_parent_R));

            /*
            std::cout << "\n     Right neighbour " <<
            tcell_R.linear_pos_idx()
                      << " does not have children" << std::endl;
            */
          }
          else
          {
            topology_input_transform_R =
                cell_adapt_op_R.get().parent_facet_adapt_op_id(local_id_in_parent_R);

            const std::vector<CellTopologyView<MeshConfig>> children_R = tcell_R.children();
            const TraceIncidences parent_to_child_on_facet_R =
                cell_adapt_op_R.get().parent_child_incidences_on_facet(local_id_in_parent_R);

            for (Uint i = 0; i < parent_to_child_on_facet_R.size(); ++i)
            {
              child_incidences_R.push_back(IncidenceEntry(
                  children_R[parent_to_child_on_facet_R.cell_id(i)].linear_pos_idx().id(),
                  parent_to_child_on_facet_R.local_id(i)));
            }
            /*
            std::cout << "\n     Right neighbour " <<
            tcell_R.linear_pos_idx()
            << " has children "; for (Uint i = 0; i <
            child_incidences_R.size();
            ++i)
            {
              std::cout << " [" << child_incidences_R[i].cell_idx <<
            ","
                        << child_incidences_R[i].local_id << "]";
            }
            std::cout << std::endl;
            */
          }

          // --------------------------------------------------------
          // Combine left and right data to generate new facets
          // --------------------------------------------------------

          /*
          std::cout << "\n     Local trace topology input:" <<
          std::endl; std::cout << "       [" <<
          child_incidences_L[0].cell_idx << ","
                    << child_incidences_L[0].local_id << "] [" <<
          child_incidences_R[0].cell_idx
                    << "," << child_incidences_R[0].local_id << "]" <<
          std::endl;
          */

          local_trace_topology.get().combine_incidences_on_adjacent_facets(
              topology_input_transform_L, topology_input_transform_R, child_incidences_L,
              child_incidences_R, combined_child_incidences, combined_child_permutations);

          for (Uint i = 0; i < combined_child_incidences.size(); ++i)
          {
            const std::vector<IncidenceEntry> &generated_incidence_pair =
                combined_child_incidences[i];
            const std::vector<EntityDofRealign> &generated_perm_pair =
                combined_child_permutations[i];

            /*
            std::cout << "\n     i = " << i << std::endl;
            std::cout << "     Newly generated:" << std::endl;
            std::cout << "     ";
            for (const auto entry : generated_incidence_pair)
            {
              std::cout << "[" << entry.cell_idx << "," <<
            entry.local_id << "]
            ";
            }
            std::cout << std::endl;
            std::cout << "     ";
            for (const auto entry : generated_perm_pair)
            {
              std::cout << entry.get().code().as_string() << " ";
            }
            std::cout << std::endl;
            */

            // Here we assume that generated_incidence_pair has size
            // 2 ...
            const CellTopologyView<MeshConfig> child_L =
                mesh_cells.cell(FlatIdx(generated_incidence_pair[LEFT].cell_idx));
            const CellTopologyView<MeshConfig> child_R =
                mesh_cells.cell(FlatIdx(generated_incidence_pair[RIGHT].cell_idx));

            if ((child_L.status() == EntityStatus::Active) &&
                (child_R.status() == EntityStatus::Active))
            {
              /*
              std::cout << "     Splitting generated facet pair ["
              << child_L.linear_pos_idx()
                        << " , " <<
              generated_perm_pair[LEFT].get().code().as_string()
              << " ] ["
                        << child_R.linear_pos_idx() << " , "
                        <<
              generated_perm_pair[RIGHT].get().code().as_string()
              << "]" << std::endl;
              */

              active_skeleton.add_facet(generated_incidence_pair, generated_perm_pair);
            }
            else
            {
              const Uint added_level_L = ((generated_perm_pair[LEFT].get().code().adapt_op_id() ==
                                           CellTransform::NO_TRANS) &&
                                          (child_incidences_L.size() > 1))
                                             ? 1
                                             : 0;
              const Uint added_level_R = ((generated_perm_pair[RIGHT].get().code().adapt_op_id() ==
                                           CellTransform::NO_TRANS) &&
                                          (child_incidences_R.size() > 1))
                                             ? 1
                                             : 0;
              /*
              std::cout << "     Added level L = " <<
              added_level_L << std::endl; std::cout << "     Added
              level R = " << added_level_R
              << std::endl;
              */

              const std::tuple<Uint, Uint> new_levels_LR =
                  get_relative_refinement_levels(std::get<LEVEL_L>(facet_entry) + added_level_L,
                                                 std::get<LEVEL_R>(facet_entry) + added_level_R);

              const interior_facet_entry_type new_facet_entry = std::make_tuple(
                  generated_incidence_pair[LEFT].cell_idx, generated_incidence_pair[LEFT].local_id,
                  generated_perm_pair[LEFT].get().code(), std::get<LEFT>(new_levels_LR),
                  generated_incidence_pair[RIGHT].cell_idx,
                  generated_incidence_pair[RIGHT].local_id, generated_perm_pair[RIGHT].get().code(),
                  std::get<RIGHT>(new_levels_LR));

              interior_facets_queue.push(new_facet_entry);
            }
          }
        } // While interior facet queue is not empty
      }

    } // If facet block size == 2
  }   // Loop over zero-level facets

  std::vector<EntityStatus> cell_status;
  cell_status.resize(mesh_cells.nb_all_cells());

  for (Uint i = 0; i < mesh_cells.nb_all_cells(); ++i)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.cell(FlatIdx(i));
    cell_status[tcell.linear_pos_idx().id()] = tcell.status();
  }

  active_skeleton.update_status(cell_status);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::remove_cells(
    TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &zero_level_skeleton,
    TriaFacets<MeshConfig> &active_skeleton, MeshBoundarySet<MeshConfig> &mesh_boundary_set,
    std::vector<CellTransform> const &active_cell_ops)
{

  for (Uint i = 0; i < active_cell_ops.size(); ++i)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ActiveIdx(i));

    if (active_cell_ops[i] == CellTransform::COARSEN)
    {
      mesh_cells.change_status(tcell.linear_pos_idx(), EntityStatus::PendingRemoval);
    }
  }

  std::vector<FlatIdx> new_cell_id;
  mesh_cells.remove_all_leaf_cells_w_status(EntityStatus::PendingRemoval, new_cell_id);
  rebuild_active_skeleton(mesh_cells, zero_level_skeleton, active_skeleton);

  // Re-generate the information about neighbours
  // Make sure that each cells refers only to active neighbours
  mesh_cells.rebuild_neighbour_information(active_skeleton.cbegin(), active_skeleton.cend());
  ensure_local_facet_numbering(mesh_cells, active_skeleton);
  check_mesh_after_topo_change(mesh_cells, active_skeleton);

  // Update boundary data

  for (auto bdry_facets : mesh_boundary_set.all_domains())
  {
    bdry_facets->update();
  }

  // std::cout << mesh_cells << std::endl;

#if 0
  mesh_cells.remove_all_leaf_cells_w_status(EntityStatus::PendingRemoval, new_cell_id);

  // Update status of facets depending on whether they refer to active or inactive cells
  active_skeleton.update_status(mesh_cells.cell_status());

  // Make sure that the local indices of trace facets of each cell
  // are correctly ordered. Each cell should first refer to its facet
  // which is incident to its local face 0, then 1, ...
  ensure_local_facet_numbering(mesh_cells, active_skeleton);
  check_mesh_after_topo_change(mesh_cells, active_skeleton);
#endif
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::adapt_w_hanging_nodes(
    TriaCells<MeshConfig> &mesh_cells, TriaFacets<MeshConfig> &skeleton_facets,
    MeshBoundarySet<MeshConfig> &mesh_boundary_set,
    std::vector<CellTransform> const &active_cell_ops)
{
  std::vector<std::tuple<Uint, CellTransform>> cell_ops_linear;
  cell_ops_linear.resize(mesh_cells.nb_all_cells());

  for (Uint i = 0; i < cell_ops_linear.size(); ++i)
  {
    cell_ops_linear[i] = std::make_tuple(INVALID_CELL_ID, CellTransform::NO_TRANS);
  }

  for (Uint ac = 0; ac < active_cell_ops.size(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ActiveIdx(ac));
    const Uint linear_idx                    = tcell.linear_pos_idx().id();
    cell_ops_linear[linear_idx]              = std::make_tuple(linear_idx, active_cell_ops[ac]);
  }

  // Check the grading
  const Uint max_level_diff = 1;
  check_mesh_h_grading(mesh_cells, skeleton_facets, active_cell_ops, max_level_diff, false);

  // First store the (absolute) positions of all cells that are active now
  // (before adaptation)
  const Uint nb_old_active_cells = mesh_cells.nb_active_cells();
  const Uint nb_old_cells_total  = mesh_cells.nb_all_cells();

  adapt::CellAdaptOp active_cell_adapt_op;
  std::vector<Uint> new_child_cell_ids;

  adapt::CellAdaptOp cell_adapt_op_L;
  adapt::CellAdaptOp cell_adapt_op_R;

  TraceEntityTuple facet_subdomain_tuple;

  CellTransform facet_adapt_op_id_L;
  CellTransform facet_adapt_op_id_R;

  // Permutation code of the left and right cell - cells
  // adjacent to given facet. These must be the cells
  // on the 'parent' level.
  EntityRealignCode pcode_L;
  EntityRealignCode pcode_R;

  EntityDofRealign bdry_facet_permutation;

  // Vector of cell incidences holding all cells adjacent to one facet
  // from the left/right side
  std::vector<IncidenceEntry> subcell_incidences_L;
  std::vector<IncidenceEntry> subcell_incidences_R;

  // Data to hold indexes of newly created facets after refinement of one
  // or more cells adjacent to given 'parent' facet
  std::vector<std::vector<IncidenceEntry>> subcell_incidences_combined;
  std::vector<std::vector<EntityDofRealign>> subcell_facet_permutations;

  // -----------------------------------------------------------------------
  // Record the adaptation operations - changes in numbering of active cells
  // so that the dof handler can perform its own renumbering
  // -----------------------------------------------------------------------

  // const Uint nb_old_all_cells = mesh_cells.nb_all_cells();
  Uint nb_cells_deactivated = 0;
  Uint nb_cells_inserted    = 0;

  // FIRST LOOP OVER CELLS:
  // - GENERATE CHILDREN OF ALL CELLS THAT ARE TO BE REFINED
  // - DON'T SET THE FACET CONNECTIVITY YET (THE INTERIOR FACETS OBTAINED BY
  // ELEMENT
  //   SPLITTING ARE PLACED IN A BUFFER)

  internal::TriaFacets<MeshConfig> facet_buffer;
  facet_buffer.set_dim(MeshConfig::TDIM - 1);

  adapt::GeometryAdapter<MeshConfig> geo_adapter;

  std::vector<IncidenceEntry> tmp_incidences;
  std::vector<EntityDofRealign> tmp_facet_permutations;

  for (Uint cell_pos_idx = 0; cell_pos_idx < cell_ops_linear.size(); ++cell_pos_idx)
  {
    const CellTransform cell_op = std::get<1>(cell_ops_linear[cell_pos_idx]);

    if ((cell_op != CellTransform::NO_TRANS) && (cell_op != CellTransform::COARSEN))
    {
      // ******************************************************
      // 1) Add children to existing cell that is being refined
      // ******************************************************

      const CellTopologyView<MeshConfig> parent_tcell = mesh_cells.cell(FlatIdx(cell_pos_idx));

      const Uint parent_level          = parent_tcell.refinement_level();
      const StdRegion parent_cell_type = parent_tcell.std_region();

      // Set the type of adaptation operation that will be applied to
      // current
      // ('active') cell
      const adapt::CellAdaptOpTag parent_adapt_tag(parent_cell_type.get().pt_set_id().elem_shape(),
                                                   cell_op);
      active_cell_adapt_op.change_type(parent_adapt_tag);

      // Make sure that the cell being split knows about its new
      // adaptation operation type
      mesh_cells.set_cell_adapt_op(FlatIdx(cell_pos_idx), active_cell_adapt_op);

      // Number of child elements after refinement of 'active' cell
      const Uint nb_new_child_elems = active_cell_adapt_op.get().nb_child_elems();

      nb_cells_inserted += nb_new_child_elems;
      nb_cells_deactivated++;

      // Vector of cell numbers of newly created child cells
      new_child_cell_ids.resize(0);

      const common::ArrayView<const math::DenseDMat<Real>, _1D, Uint> child_coords =
          geo_adapter.compute_child_coords(parent_cell_type.get().pt_set_id(), parent_adapt_tag,
                                           parent_tcell.coordinates());

      for (Uint nc = 0; nc < nb_new_child_elems; ++nc)
      {
        // Cell adapt op knows the type of the element, but not its
        // polynomial order We set the child to have the same polynomial
        // order as the parent
        const PointSetTag child_type_tag(active_cell_adapt_op.get().child_elem_shape(nc),
                                         parent_cell_type.get().pt_set_id().poly_order(),
                                         parent_cell_type.get().pt_set_id().ref_topology());

        const StdRegion child_type(child_type_tag);

        const Uint new_child_cell_idx = mesh_cells.add_cell(
            child_type, parent_level + 1, EntityStatus::Active, cell_pos_idx, child_coords[nc]);

        // mesh_cells.set_parent(new_child_cell_idx,
        // active_cell_pos_idx);

        // Register the new cell ids in the linear cell index and
        // collect them in child_cell_ids
        new_child_cell_ids.push_back(new_child_cell_idx);
      }

      // Register the new children as children of active cell (i.e. the
      // 'parent') This simultaneously marks the parent as going through
      // refinement
      mesh_cells.add_child_indexes(cell_pos_idx, new_child_cell_ids,
                                   EntityStatus::PendingRefinement);

      // ************************************************************
      // 2) The cell splitting created new facets that are all placed
      //    IN THE INTERIOR of the parent cell. Add these
      //    facets to a facet buffer
      // ************************************************************

      const Uint nb_new_internal_facets = active_cell_adapt_op.get().nb_internal_child_facets();

      for (Uint f = 0; f < nb_new_internal_facets; ++f)
      {
        // TODO: check 'recreate_local_incidence_patterns'
        const TraceIncidences internal_facet = active_cell_adapt_op.get().internal_child_facet(f);

        const Uint facet_size = internal_facet.size();
        tmp_incidences.resize(facet_size);
        tmp_facet_permutations.resize(facet_size);

        for (Uint i = 0; i < facet_size; ++i)
        {
          tmp_incidences[i].cell_idx = new_child_cell_ids[internal_facet.cell_id(i)];
          tmp_incidences[i].local_id = internal_facet.local_id(i);
          tmp_facet_permutations[i]  = internal_facet.permutation(i);
        }

        /*
        const TraceIncidences new_incidences =
            facet_buffer.insert_facet_block(internal_facet, new_child_cell_ids);
        */
        const TraceIncidences new_incidences =
            facet_buffer.add_facet(tmp_incidences, tmp_facet_permutations);
      }

    } // if (cell_ops[active_cell] > 1)
  }   // End of first pass through all cells

  // The goal of the block of code below is the following:
  // 1) For left facet, fill a vector of incidence pairs that indicate
  //    the relation of the child cell facets to the left (parent) facet.
  //    By 'relation' we mean its position as sub-entity within the parent and
  //    its permutation (rotation, flip)
  // 2) The same is done for the right facet.
  // 3) Use the two arrays filled in 1) and 2) to generate a set of new
  //    facets obtained by the splitting (i.e. match each facet
  //    in 'left' array of incidence pairs with a corresponding entity
  //    in 'right' array
  //    Then add this to the facet buffer.
  // 4) After the facet buffer is filled, join it with the original cells.

  /// TODO: this should be replaced by a loop over ACTIVE facets, instead of
  /// looping over all facets
  const_skeleton_iterator facets_begin(skeleton_facets, 0);
  const_skeleton_iterator facets_end(skeleton_facets, skeleton_facets.nb_all_facets());
  const_skeleton_iterator facet_it;

  for (facet_it = facets_begin; facet_it != facets_end; ++facet_it)
  {
    const TraceIncidences mesh_facet = facet_it->incidences();

    // If the size of the facet block is equal to one, then this is a
    // boundary facet
    if (mesh_facet.size() == 1)
    {

      // Boundary facet ...
      const IncidenceEntry inc_pair = mesh_facet.incidence_entry(0);

      // If this face is not adjacent to active entity, just skip it and
      // move on to the next face
      const CellTopologyView<MeshConfig> tcell = mesh_cells.cell(FlatIdx(inc_pair.cell_idx));

      const bool cell_op_is_refinement = (tcell.status() == EntityStatus::PendingRefinement);

      // -----------------------------------------------------
      // Get the adaptation operation id of the boundary facet
      // -----------------------------------------------------
      if (cell_op_is_refinement)
      {
        const StdRegion cell_type = tcell.std_region();

        // Set the type of adaptation operation that will be applied to
        // the cell
        cell_adapt_op_L.change_type(cell_type.get().pt_set_id().elem_shape(),
                                    std::get<1>(cell_ops_linear[tcell.linear_pos_idx().id()]));

        facet_adapt_op_id_L = cell_adapt_op_L.get().parent_facet_adapt_op_id(inc_pair.local_id);

        // -------------------------------------------------
        // Get the vector of incidences to the facet.
        // This vector should contain the incidence pair
        // of all child cells touching the facet AFTER
        // the local adaptation operation is applied.
        // -------------------------------------------------

        // List of ALL children of the left cell
        const std::vector<CellTopologyView<MeshConfig>> cell_children =
            mesh_cells.cell_children(FlatIdx(tcell.linear_pos_idx()));

        const TraceIncidences child_facet_incidences =
            cell_adapt_op_L.get().parent_child_incidences_on_facet(inc_pair.local_id);

        subcell_incidences_L.resize(child_facet_incidences.size());
        for (Uint i = 0; i < subcell_incidences_L.size(); ++i)
        {
          subcell_incidences_L[i].cell_idx =
              cell_children[child_facet_incidences.cell_id(i)].linear_pos_idx().id();
          subcell_incidences_L[i].local_id = child_facet_incidences.local_id(i);
        }

        for (Uint new_face = 0; new_face < child_facet_incidences.size(); ++new_face)
        {
          // Update the incidence pair of boundary facet so that its
          // cell index is cell index in 'global' cell numbering and
          // not only index withing parent cell
          IncidenceEntry new_inc_pair = child_facet_incidences.incidence_entry(new_face);
          new_inc_pair.cell_idx       = cell_children[new_inc_pair.cell_idx].linear_pos_idx().id();

          // Generate the subdomain tag of the children: it should
          // have DO_NOTHING as adaptation operation and 0 as local
          // position, because it belongs to a newly generated element
          // which does not have a neighbor across this edge

          const EntityDofRealign tmp_permutation = child_facet_incidences.permutation(new_face);
          const EntityRealignCode tmp_code       = tmp_permutation.get().code();
          pcode_L = EntityRealignCode(tmp_code.elem_shape(), CellTransform::NO_TRANS, 0,
                                      tmp_code.parent_shape(), tmp_code.nb_flips(),
                                      tmp_code.nb_rotations());

          bdry_facet_permutation.change_type(tmp_permutation.get().type_id(), pcode_L);
          const TraceIncidences new_incidences =
              facet_buffer.add_facet({new_inc_pair}, {bdry_facet_permutation});
        }

      } // If the adaptation operation on boundary cell is refinement

      // ----------------------------------
      // If it's not refinement, do nothing
      // ----------------------------------

    } // If this is boundary face ...
    // ... else this is interior face
    else
    {
      const IncidenceEntry inc_pair_L = mesh_facet.incidence_entry(0);
      const IncidenceEntry inc_pair_R = mesh_facet.incidence_entry(1);

      const CellTopologyView<MeshConfig> tcell_L = mesh_cells.cell(FlatIdx(inc_pair_L.cell_idx));
      const CellTopologyView<MeshConfig> tcell_R = mesh_cells.cell(FlatIdx(inc_pair_R.cell_idx));

      const bool cell_op_L_is_refinement = (tcell_L.status() == EntityStatus::PendingRefinement);
      const bool cell_op_R_is_refinement = (tcell_R.status() == EntityStatus::PendingRefinement);

      const bool facet_is_not_deactivated = ((tcell_L.status() != EntityStatus::NotActive) &&
                                             (tcell_R.status() != EntityStatus::NotActive));

      if (facet_is_not_deactivated && (cell_op_L_is_refinement || cell_op_R_is_refinement))
      {
        // -------------------------------------------------
        // Get the adaptation operation id of the left facet
        // -------------------------------------------------
        if (cell_op_L_is_refinement)
        {
          const StdRegion left_cell_type = tcell_L.std_region();

          // Set the type of adaptation operation that will be applied
          // to left cell
          cell_adapt_op_L.change_type(left_cell_type.get().pt_set_id().elem_shape(),
                                      std::get<1>(cell_ops_linear[tcell_L.linear_pos_idx().id()]));

          facet_adapt_op_id_L = cell_adapt_op_L.get().parent_facet_adapt_op_id(inc_pair_L.local_id);

          // -------------------------------------------------
          // Get the vector of incidences to the left facet.
          // This vector should contain the incidence pair of
          // all child cells touching left facet AFTER the
          // local adaptation operation is applied.
          // -------------------------------------------------

          // List of ALL children of the left cell
          const std::vector<CellTopologyView<MeshConfig>> cell_children_L =
              mesh_cells.cell_children(FlatIdx(tcell_L.linear_pos_idx()));

          const TraceIncidences child_facet_incidences_L =
              cell_adapt_op_L.get().parent_child_incidences_on_facet(inc_pair_L.local_id);

          subcell_incidences_L.resize(child_facet_incidences_L.size());
          for (Uint i = 0; i < subcell_incidences_L.size(); ++i)
          {
            subcell_incidences_L[i].cell_idx =
                cell_children_L[child_facet_incidences_L.cell_id(i)].linear_pos_idx().id();
            subcell_incidences_L[i].local_id = child_facet_incidences_L.local_id(i);
          }
        }
        // else if (cell_ops[tcell_L.active_idx()] ==
        // CellAdaptOpID::DO_NOTHING)
        else
        {
          facet_adapt_op_id_L = CellTransform::NO_TRANS;
          subcell_incidences_L.resize(1);
          subcell_incidences_L[0] =
              IncidenceEntry(tcell_L.linear_pos_idx().id(), inc_pair_L.local_id);
        }

        // --------------------------------------------------
        // Get the adaptation operation id of the right facet
        // --------------------------------------------------
        if (cell_op_R_is_refinement)
        {
          const StdRegion right_cell_type = tcell_R.std_region();

          cell_adapt_op_R.change_type(right_cell_type.get().pt_set_id().elem_shape(),
                                      std::get<1>(cell_ops_linear[tcell_R.linear_pos_idx().id()]));
          facet_adapt_op_id_R = cell_adapt_op_R.get().parent_facet_adapt_op_id(inc_pair_R.local_id);

          // -------------------------------------------------
          // Get the vector of incidences to the right facet.
          // This vector should contain the incidence pair of
          // all child cells touching right facet AFTER the
          // local adaptation operation is applied.
          // -------------------------------------------------

          // List of ALL children of the right cell
          const std::vector<CellTopologyView<MeshConfig>> cell_children_R =
              mesh_cells.cell_children(FlatIdx(tcell_R.linear_pos_idx()));

          const TraceIncidences child_facet_incidences_R =
              cell_adapt_op_R.get().parent_child_incidences_on_facet(inc_pair_R.local_id);

          subcell_incidences_R.resize(child_facet_incidences_R.size());
          for (Uint i = 0; i < subcell_incidences_R.size(); ++i)
          {
            subcell_incidences_R[i].cell_idx =
                cell_children_R[child_facet_incidences_R.cell_id(i)].linear_pos_idx().id();
            subcell_incidences_R[i].local_id = child_facet_incidences_R.local_id(i);
          }
        }
        // else if (cell_ops[tcell_R.active_idx()] ==
        // CellAdaptOpID::DO_NOTHING)
        else
        {
          facet_adapt_op_id_R = CellTransform::NO_TRANS;
          subcell_incidences_R.resize(1);
          subcell_incidences_R[0] =
              IncidenceEntry(tcell_R.linear_pos_idx().id(), inc_pair_R.local_id);
        }

        // -------------------------------------------------
        // Get the permutation codes of left and right facet
        // -------------------------------------------------
        pcode_L = mesh_facet.permutation(0).get().code();
        pcode_R = mesh_facet.permutation(1).get().code();

        const std::tuple<Uint, Uint> relative_ref_levels =
            get_relative_refinement_levels(tcell_L.refinement_level(), tcell_R.refinement_level());

        facet_subdomain_tuple.change_type(pcode_L, std::get<LEFT>(relative_ref_levels), pcode_R,
                                          std::get<RIGHT>(relative_ref_levels));
        // std::cout << facet_subdomain_tuple.get() << std::endl;
        facet_subdomain_tuple.get().combine_incidences_on_adjacent_facets(
            facet_adapt_op_id_L, facet_adapt_op_id_R, subcell_incidences_L, subcell_incidences_R,
            subcell_incidences_combined, subcell_facet_permutations);

        for (Uint new_face = 0; new_face < subcell_incidences_combined.size(); ++new_face)
        {
          facet_buffer.add_facet(subcell_incidences_combined[new_face],
                                 subcell_facet_permutations[new_face]);
        }
      }
    } // Case when this was internal facet
  }   // Loop over facets

#if 0
  std::vector<Uint> new_cell_facet_offsets(nb_cells_inserted + 1);
  new_cell_facet_offsets.assign(nb_cells_inserted + 1, 0);

  const internal::TriaFacets<MeshConfig> &facet_buffer_const(facet_buffer);

  for (Uint f = 0; f < facet_buffer_const.nb_all_facets(); ++f)
  {
    const TraceIncidences facet_block = facet_buffer_const.facet_data(FlatIdx(f));
    for (Uint i = 0; i < facet_block.size(); ++i)
    {
      if (facet_block.cell_id(i) >= nb_old_cells_total)
      {
        const Uint shifted_cell_idx = facet_block.cell_id(i) - nb_old_cells_total;
        new_cell_facet_offsets[shifted_cell_idx + 1]++;
      }
    }
  }

  for (Uint i = 0; i < (new_cell_facet_offsets.size() - 1); ++i)
  {
    new_cell_facet_offsets[i + 1] += new_cell_facet_offsets[i];
  }

  // Update information about facets incident to new cells
  const Uint nb_new_facet_incidences = new_cell_facet_offsets.back();

  std::vector<Uint> new_incident_facet_ids(nb_new_facet_incidences, INVALID_FACET_ID);
  for (Uint f = 0; f < facet_buffer_const.nb_all_facets(); ++f)
  {
    const TraceIncidences facet_block = facet_buffer_const.facet_data(FlatIdx(f));
    for (Uint i = 0; i < facet_block.size(); ++i)
    {
      if (facet_block.cell_id(i) >= nb_old_cells_total)
      {
        const Uint shifted_cell_idx = facet_block.cell_id(i) - nb_old_cells_total;
        for (Uint j = new_cell_facet_offsets[shifted_cell_idx];
             j < new_cell_facet_offsets[shifted_cell_idx + 1]; ++j)
        {
          if (new_incident_facet_ids[j] == INVALID_FACET_ID)
          {
            new_incident_facet_ids[j] = f;
            break;
          }
        }
      }
    }
  }
#endif

  std::vector<Uint> new_cell_facet_offsets;
  std::vector<Uint> new_incident_facet_ids;
  generate_new_facet_offsets_from_buffer(facet_buffer, nb_old_cells_total, nb_cells_inserted,
                                         new_cell_facet_offsets, new_incident_facet_ids);

  // Update information about neighbours of cells
  mesh_cells.add_new_cell_neighbours(skeleton_facets.nb_all_facets(), new_incident_facet_ids,
                                     new_cell_facet_offsets);

  // Merge the existing skeleton with buffer
  skeleton_facets.merge(facet_buffer);

  mesh_cells.deactivate_nonleafs_and_update_active_cell_positions();
  const Uint nb_new_active_cells = mesh_cells.nb_active_cells();

  std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
  std::cout << "Adaptation (with hanging nodes) stats: " << std::endl;
  std::cout << "   Number of active cells before adaptation: " << nb_old_active_cells << std::endl;
  std::cout << "   Number of deactivated cells:              " << nb_cells_deactivated << std::endl;
  std::cout << "   Number of newly inserted cells:           " << nb_cells_inserted << std::endl;
  std::cout << "   Number of active cells after adaptation:  " << nb_new_active_cells << std::endl;
  std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

  // Update status of facets depending on whether they refer to active or
  // inactive cells
  skeleton_facets.update_status(mesh_cells.cell_status());

  // Make sure that the local indices of trace facets of each cell
  // are correctly ordered. Each cell should first refer to its facet
  // which is incident to its local face 0, then 1, ...
  ensure_local_facet_numbering(mesh_cells, skeleton_facets);

  // Update boundary data
  for (auto bdry_facets : mesh_boundary_set.all_domains())
  {
    bdry_facets->update();
  }

  // Re-generate the information about neighbours
  // Make sure that each cells refers only to active neighbours
  // mesh_cells.rebuild_neighbour_information(skeleton_facets.cbegin(),
  // skeleton_facets.cend());

  check_mesh_after_topo_change(mesh_cells, skeleton_facets);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::adapt_red_green(
    TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &zero_level_skeleton,
    TriaFacets<MeshConfig> &active_skeleton, MeshBoundarySet<MeshConfig> &mesh_boundary_set,
    std::vector<CellTransform> const &active_cell_ops)
{
  // ------------------------------------------------------
  // Phase 1: cache all adaptation operations by tying them
  //          to cell paths. This is necessary as
  //          connectivity is going to change in phase 2
  // ------------------------------------------------------

  std::vector<std::tuple<CellTransform, CellPath>> cell_transf_cache;
  cell_transf_cache.reserve(active_cell_ops.size());
  cell_transf_cache.resize(0);

  std::vector<bool> cell_processed;
  cell_processed.resize(mesh_cells.nb_all_cells());

  // Leave only cells that are currently active as 'not processed'
  cell_processed.assign(mesh_cells.nb_all_cells(), true);
  for (Uint c = 0; c < mesh_cells.nb_all_cells(); ++c)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.cell(FlatIdx(c));
    if (tcell.status() == EntityStatus::Active)
    {
      cell_processed[tcell.linear_pos_idx().id()] = false;
    }
  }

  CellTopologyView<MeshConfig> root_tcell;
  std::vector<Uint> path_entries;

  for (Uint cell_id = 0; cell_id < mesh_cells.nb_active_cells(); ++cell_id)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ActiveIdx(cell_id));

    if (!cell_processed[tcell.linear_pos_idx().id()])
    {
      if (tcell.refinement_level() == 0)
      {
        const CellPath cp(FlatIdx(tcell.linear_pos_idx()), {});
        cell_transf_cache.push_back(std::make_tuple(active_cell_ops[tcell.active_idx().id()], cp));
        cell_processed[tcell.linear_pos_idx().id()] = true;
      }
      else
      {
        const CellTopologyView<MeshConfig> parent_tcell = tcell.parent();
        const adapt::CellAdaptOp parent_adapt_op        = parent_tcell.cell_adapt_op();

        // If the operation applied to parent was 'green' refinement,
        // then store the information about what operation to apply to
        // parent, because green children will be removed in phase 2 If
        // only one green child is scheduled for refinement, then the
        // parent will be refined anyway
        if (CellTransformTraits::is_aniso_refinement(
                parent_adapt_op.get().cell_adapt_op_tag().adapt_op_id()))
        {
          const std::vector<CellTopologyView<MeshConfig>> child_tcells = parent_tcell.children();

          bool refine_parent = false;
          for (const CellTopologyView<MeshConfig> &child_tcell : child_tcells)
          {
            // Mark all children as processed: they will be removed
            // and we will work with their parent in the next phase
            cell_processed[child_tcell.linear_pos_idx().id()] = true;
            if (active_cell_ops[child_tcell.active_idx().id()] != CellTransform::NO_TRANS)
            {
              refine_parent = true;
            }
          }

          mesh_cells.path(parent_tcell, root_tcell, path_entries);
          const CellPath parent_path(FlatIdx(root_tcell.linear_pos_idx()), path_entries);

          const std::tuple<CellTransform, CellPath> parent_cache_entry =
              refine_parent ? std::make_tuple(CellTransform::UNIFORM_REFINE, parent_path)
                            : std::make_tuple(CellTransform::NO_TRANS, parent_path);

          cell_transf_cache.push_back(parent_cache_entry);

        } // If the parent cell went through aniso refinement
        else
        {
          mesh_cells.path(tcell, root_tcell, path_entries);
          const CellPath cp(FlatIdx(root_tcell.linear_pos_idx()), path_entries);
          cell_transf_cache.push_back(
              std::make_tuple(active_cell_ops[tcell.active_idx().id()], cp));
          cell_processed[tcell.linear_pos_idx().id()] = true;
        }

      } // If this cell is not on refinement level 0
    }   // If this cell was not processed yet
  }

  /*
  std::cout << "Debug a: Contents of cell transform cache (after phase 1):" <<
  std::endl; for (auto cache_entry : cell_transf_cache)
  {
    std::cout << std::get<0>(cache_entry) << " --> " <<
  std::get<1>(cache_entry)
  << std::endl;
  }
  */

  // ------------------------------------------------------
  // Phase 2: find all 'green' cells and remove them
  // ------------------------------------------------------

  std::vector<CellTransform> tmp_cell_transf;
  Uint nb_green_cells = 0;

  tmp_cell_transf.resize(mesh_cells.nb_active_cells());
  tmp_cell_transf.assign(mesh_cells.nb_active_cells(), CellTransform::NO_TRANS);

  for (Uint ac = 0; ac < mesh_cells.nb_active_cells(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ActiveIdx(ac));
    const bool check_cell                    = (tcell.refinement_level() > 0) &&
                            (tmp_cell_transf[tcell.active_idx().id()] == CellTransform::NO_TRANS);

    if (check_cell)
    {
      const CellTopologyView<MeshConfig> parent_tcell = tcell.parent();
      const adapt::CellAdaptOp parent_adapt_op        = parent_tcell.cell_adapt_op();

      if (CellTransformTraits::is_aniso_refinement(
              parent_adapt_op.get().cell_adapt_op_tag().adapt_op_id()))
      {
        const std::vector<CellTopologyView<MeshConfig>> children = parent_tcell.children();
        nb_green_cells += children.size();

        for (Uint c = 0; c < children.size(); ++c)
        {
          const CellTopologyView<MeshConfig> child = children[c];
          tmp_cell_transf[child.active_idx().id()] = CellTransform::COARSEN;
        }
      } // If parent was anisotropically refined
    }
  }

  // If there are green cells in the mesh, remove them
  // std::cout << "REMOVAL OF GREEN CELLS START" << std::endl;
  if (nb_green_cells > 0)
  {
    remove_cells(mesh_cells, zero_level_skeleton, active_skeleton, mesh_boundary_set,
                 tmp_cell_transf);
    // std::cout << mesh_cells << std::endl;
  }
  /*
  else
  {
    std::cout << "Removal of green cells: no cells to remove" << std::endl;
  }
  std::cout << "REMOVAL OF GREEN CELLS END" << std::endl;
  */

  // ------------------------------------------------------
  // Phase 3: restore the adaptation operations cached
  //          during phase one
  // ------------------------------------------------------

  tmp_cell_transf.resize(mesh_cells.nb_active_cells());
  tmp_cell_transf.assign(mesh_cells.nb_active_cells(), CellTransform::NO_TRANS);

  for (const std::tuple<CellTransform, CellPath> &cache_entry : cell_transf_cache)
  {
    const CellTopologyView<MeshConfig> leaf_tcell = mesh_cells.path_leaf(std::get<1>(cache_entry));
    if (leaf_tcell.nb_children() > 0)
    {
      std::cerr << "MeshTopologyAdaptAlgorithm::adapt_red_green: cell with "
                   "linear idx "
                << leaf_tcell.linear_pos_idx() << " should have\n"
                << "been a leaf cell, but is not. Error." << std::endl;
    }
    tmp_cell_transf[leaf_tcell.active_idx().id()] = std::get<0>(cache_entry);
  }

  // ------------------------------------------------------
  // Phase 4: perform 'red' adaptation with hanging nodes
  // ------------------------------------------------------
  // prepare_cells_for_hanging_node_adapt(mesh_cells, skeleton_facets,
  // tmp_cell_transf);

  std::vector<Uint> colors;
  generate_red_adapt_ops(mesh_cells, active_skeleton, tmp_cell_transf, colors);

  /*
  std::cout << "Red adaptation operations after check:" << std::endl;
  for (Uint ac = 0; ac < tmp_cell_transf.size(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ac);
    std::cout << "Linear index " << tcell.linear_pos_idx() << ", op = " <<
  tmp_cell_transf[ac]
              << std::endl;
  }
  */

  adapt_w_hanging_nodes(mesh_cells, active_skeleton, mesh_boundary_set, tmp_cell_transf);

  // ------------------------------------------------------
  // Phase 5: add green cells in order to obtain
  //          conforming mesh
  // ------------------------------------------------------
  generate_green_adapt_ops(mesh_cells, active_skeleton, tmp_cell_transf);

  /*
  std::cout << "Green adaptation operations:" << std::endl;
  for (Uint ac = 0; ac < tmp_cell_transf.size(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ac);
    std::cout << "Linear index " << tcell.linear_pos_idx() << ", op = " <<
  tmp_cell_transf[ac]
              << std::endl;
  }
  */

  adapt_w_hanging_nodes(mesh_cells, active_skeleton, mesh_boundary_set, tmp_cell_transf);
  // rebuild_active_skeleton(mesh_cells, zero_level_skeleton,
  // active_skeleton);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::prepare_for_cell_removal(
    TriaCells<MeshConfig> const &mesh_cells, TriaFacets<MeshConfig> const &skeleton_facets,
    std::vector<CellTransform> &adapt_op)
{
  std::vector<bool> cell_will_be_removed(mesh_cells.nb_all_cells());
  cell_will_be_removed.assign(mesh_cells.nb_all_cells(), false);

  bool inadmissible_op_present = false;
  Uint inadmissible_op_cell_id = 0;

  for (Uint ac = 0; ac < adapt_op.size(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ActiveIdx(ac));
    if (adapt_op[ac] == CellTransform::COARSEN)
    {
      cell_will_be_removed[tcell.linear_pos_idx().id()] = true;
    }

    if ((adapt_op[ac] != CellTransform::NO_TRANS) && (adapt_op[ac] != CellTransform::COARSEN))
    {
      inadmissible_op_present = true;
      inadmissible_op_cell_id = tcell.linear_pos_idx().id();
      break;
    }
  }

  if (inadmissible_op_present)
  {
    const CellTopologyView<MeshConfig> inadmissible_op_cell =
        mesh_cells.cell(FlatIdx(inadmissible_op_cell_id));
    std::cerr << "MeshTopologyAdaptAlgorithm::prepare_for_cell_removal: unable "
                 "to proceed with "
                 "cell removal.\n"
              << "Cell with flat idx = " << inadmissible_op_cell.linear_pos_idx()
              << " and active idx " << inadmissible_op_cell.active_idx()
              << " currently scheduled for inadmissible "
              << "operation " << adapt_op[inadmissible_op_cell.active_idx().id()] << "\n"
              << "I won't remove any cells." << std::endl;
    return;
  }

  // Loop over all cells and make sure that if cell is marked for removal,
  // then its siblings (i.e. other cells that have the same parent) are marked
  // for removal as well

  bool level_zero_cell_will_be_removed = false;
  Uint level_zero_cell_id              = 0;

  bool only_leaf_cells_will_be_removed = true;
  Uint non_leaf_cell_id                = 0;

  for (Uint c = 0; c < mesh_cells.nb_all_cells(); ++c)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.cell(FlatIdx(c));
    if (cell_will_be_removed[tcell.linear_pos_idx().id()])
    {
      if (tcell.refinement_level() == 0)
      {
        level_zero_cell_will_be_removed = true;
        level_zero_cell_id              = tcell.linear_pos_idx().id();
        break;
      }
      else
      {
        const CellTopologyView<MeshConfig> tcell_parent = tcell.parent();

        const std::vector<CellTopologyView<MeshConfig>> tcell_children = tcell_parent.children();
        for (auto child : tcell_children)
        {
          if (!cell_will_be_removed[child.linear_pos_idx().id()])
          {
            cell_will_be_removed[child.linear_pos_idx().id()] = true;
            adapt_op[child.active_idx().id()]                 = CellTransform::COARSEN;
          }
        }
      }

      if (cell_will_be_removed[tcell.linear_pos_idx().id()] && (tcell.nb_children() > 0))
      {
        only_leaf_cells_will_be_removed = false;
        non_leaf_cell_id                = c;
        break;
      }
    }
  }

  if (!only_leaf_cells_will_be_removed)
  {
    const CellTopologyView<MeshConfig> non_leaf_cell = mesh_cells.cell(FlatIdx(non_leaf_cell_id));
    std::cerr << "MeshTopologyAdaptAlgorithm::prepare_for_cell_removal: unable "
                 "to remove\n"
              << "cell with flat idx = " << non_leaf_cell.linear_pos_idx() << " and active idx "
              << non_leaf_cell.active_idx() << " currently scheduled for removal.\n"
              << "This is not a leaf-type cell. I won't remove any cells." << std::endl;
    return;
  }

  if (level_zero_cell_will_be_removed)
  {
    const CellTopologyView<MeshConfig> level_zero_cell =
        mesh_cells.cell(FlatIdx(level_zero_cell_id));
    std::cerr << "MeshTopologyAdaptAlgorithm::prepare_for_cell_removal: unable "
                 "to remove\n"
              << "cell with flat idx = " << level_zero_cell.linear_pos_idx() << " and active idx "
              << level_zero_cell.active_idx() << " currently scheduled for removal.\n"
              << "This is zero-level cell and thus non-removable. I won't "
                 "remove any cells."
              << std::endl;
    return;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::prepare_cells_for_hanging_node_adapt(
    TriaCells<MeshConfig> const &mesh_cells, TriaFacets<MeshConfig> const &skeleton_facets,
    std::vector<CellTransform> &adapt_op)
{
#if 0
              const Uint nb_active_cells = mesh_cells.nb_active_cells();
              std::vector<bool> cell_has_been_visited(nb_active_cells, false);

              // Maximum difference between levels of neighbouring cells
              const Uint max_level_diff = 1;
              std::queue<Uint> active_cells_to_check;

              // First insert all cells adjacent to those already marked
              // for refinement into the queue. These adjacent cells
              // will be checked first
              for (Uint ac = 0; ac < nb_active_cells; ++ac)
              {
              if (adapt_op[ac] != CellTransform::DO_NOTHING)
              {
              active_cells_to_check.push(ac);
              }
              }

              while (!active_cells_to_check.empty())
              {
              // Get the index of first cell in queue
              // and REMOVE IT FROM QUEUE
              const Uint ac = active_cells_to_check.front();
              active_cells_to_check.pop();

              // If the cell did not pass check, loop over
              // its neighbours and see what is the difference in their refinement
              // levels

              const CellTopologyView<MeshConfig> center_cell = mesh_cells.active_cell(ac);
              Uint refinement_level_curr_cell = (adapt_op[ac] != CellTransform::DO_NOTHING)
                                ? center_cell.refinement_level() + 1
                                : center_cell.refinement_level();

              const std::vector<CellTopologyView<MeshConfig>> neighbours = active_neighbours(mesh_cells, skeleton_facets, center_cell);

              // Mark the current center as checked
              cell_has_been_visited[ac] = true;

              // Now verify that this is really true by comparing the refinement
              // level of the center cell with neighbours
              for (auto neighb_cell : neighbours)
              {
              const Uint refinement_level_neighb_cell =
              (adapt_op[neighb_cell.active_idx()] != CellTransform::DO_NOTHING)
              ? neighb_cell.refinement_level() + 1
              : neighb_cell.refinement_level();

              const Uint max_level = std::max(refinement_level_curr_cell, refinement_level_neighb_cell);
              const Uint min_level = std::min(refinement_level_curr_cell, refinement_level_neighb_cell);
              const Uint current_level_diff = max_level - min_level;

              // 1) Case when current cell is 'less refined' than neighbour
              if (refinement_level_curr_cell < refinement_level_neighb_cell)
              {
              if (current_level_diff > max_level_diff)
              {
              // If the refinement level difference is bigger than allowed maximum,
              // mark current cell for refinement
              // This will DECREASE the difference in refinement levels by 1
              if ((adapt_op[ac] == CellTransform::DO_NOTHING) &&
              ((current_level_diff - 1) == max_level_diff))
              {
              adapt_op[ac] = CellTransform::UNIFORM_REFINE;
              refinement_level_curr_cell++;
              }
              else
              {
              std::cerr
              << "MeshAdaptSchedule::define_h_adapt_ops::Error, did not manage to configure\n"
              << "cells for refinement: difference between cell refinement levels is too big"
              << std::endl;
              // return;
              }
              }
              }

              if (!cell_has_been_visited[neighb_cell.active_idx()])
              {
              active_cells_to_check.push(neighb_cell.active_idx());
              }

              } // Loop over neighbours of current cell

              } // While the queue is not empty

              // Set the p-adaptation vector size to 0 - only h adaptation is considered
              m_cell_p_order.resize(0);

              const bool mesh_passed_grading_check = check_mesh_h_grading(mesh, adapt_op, 1);

              if (mesh_passed_grading_check)
              {
              std::cout << "MeshTopologyAdaptAlgorithm: mesh passed grading check" << std::endl;
              }
              else
              {
              std::cout << "MeshTopologyAdaptAlgorithm: mesh did not pass grading check" << std::endl;
              }

#else

  const Uint nb_active_cells             = mesh_cells.nb_active_cells();
  Uint nb_cells_scheduled_for_refinement = 0;

  // Maximum difference between levels of neighbouring cells
  const Uint max_level_diff = 1;

  for (Uint i = 0; i < adapt_op.size(); ++i)
  {
    if (adapt_op[i] != CellTransform::NO_TRANS)
    {
      nb_cells_scheduled_for_refinement++;
    }
  }

  std::vector<Uint> new_cell_active_ids_to_refine;
  bool mesh_passed_grading_check = false;

  while (!mesh_passed_grading_check && (nb_cells_scheduled_for_refinement < nb_active_cells))
  {
    for (Uint ac = 0; ac < nb_active_cells; ++ac)
    {
      const CellTopologyView<MeshConfig> center_cell = mesh_cells.active_cell(ActiveIdx(ac));
      const Uint refinement_level_curr_cell          = (adapt_op[ac] != CellTransform::NO_TRANS)
                                                  ? center_cell.refinement_level() + 1
                                                  : center_cell.refinement_level();

      const std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> neighbours =
          active_neighbours(mesh_cells, skeleton_facets, center_cell);

      // Now verify that this is really true by comparing the refinement
      // level of the center cell with neighbours
      for (auto neighb_item : neighbours)
      {
        const CellTopologyView<MeshConfig> neighb_cell = std::get<0>(neighb_item);
        const Uint refinement_level_neighb_cell =
            (adapt_op[neighb_cell.active_idx().id()] != CellTransform::NO_TRANS)
                ? neighb_cell.refinement_level() + 1
                : neighb_cell.refinement_level();

        const Uint max_level = std::max(refinement_level_curr_cell, refinement_level_neighb_cell);
        const Uint min_level = std::min(refinement_level_curr_cell, refinement_level_neighb_cell);
        const Uint current_level_diff = max_level - min_level;

        if (current_level_diff > max_level_diff)
        {
          // 1) Case when current cell is 'less refined' than
          // neighbour
          if (refinement_level_curr_cell < refinement_level_neighb_cell)
          {
            // If the refinement level difference is bigger than
            // allowed maximum, mark current cell for refinement
            // This will DECREASE the difference in refinement
            // levels by 1
            if (adapt_op[ac] == CellTransform::NO_TRANS)
            {
              adapt_op[ac] = CellTransform::UNIFORM_REFINE;
              new_cell_active_ids_to_refine.push_back(ac);
              nb_cells_scheduled_for_refinement++;
            }
          }
          else
          {
            // If the refinement level difference is bigger than
            // allowed maximum, mark current cell for refinement
            // This will DECREASE the difference in refinement
            // levels by 1
            if (adapt_op[neighb_cell.active_idx().id()] == CellTransform::NO_TRANS)
            {
              adapt_op[neighb_cell.active_idx().id()] = CellTransform::UNIFORM_REFINE;
              new_cell_active_ids_to_refine.push_back(neighb_cell.active_idx().id());
              nb_cells_scheduled_for_refinement++;
            }
          }
        }
      } // Loop over neighbours of active cell
    }   // Loop over active cells

    for (Uint j = 0; j < new_cell_active_ids_to_refine.size(); ++j)
    {
      adapt_op[new_cell_active_ids_to_refine[j]] = CellTransform::UNIFORM_REFINE;
    }

    new_cell_active_ids_to_refine.resize(0);
    mesh_passed_grading_check =
        check_mesh_h_grading(mesh_cells, skeleton_facets, adapt_op, max_level_diff);

  } // while
  mesh_passed_grading_check =
      check_mesh_h_grading(mesh_cells, skeleton_facets, adapt_op, max_level_diff);

#endif

  if (mesh_passed_grading_check)
  {
    std::cout << "MeshTopologyAdaptAlgorithm: mesh passed grading check" << std::endl;
  }
  else
  {
    std::cout << "MeshTopologyAdaptAlgorithm: mesh did not pass grading check" << std::endl;
  }

  /*
  if (!mesh_passed_grading_check)
  {
  std::cout << "CHECKME" << std::endl;
  check_mesh_h_grading(mesh_cells, skeleton_facets, adapt_op, max_level_diff,
  true);
  }
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::generate_red_adapt_ops(
    const TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &skeleton_facets,
    std::vector<CellTransform> &adapt_op, std::vector<Uint> &colors)
{
  // ------------------------------------------------------
  // PHASE 1: check that only admissible adaptation
  // operations are required on  mesh cells
  // ------------------------------------------------------
  const Uint nb_active_cells = mesh_cells.nb_active_cells();

  // First check that the cells are scheduled for two types
  // of transformation: either uniform refinement, or no action
  bool cells_have_admissible_cell_transform = true;
  bool green_cells_present                  = false;

  for (Uint i = 0; i < adapt_op.size(); ++i)
  {
    if ((adapt_op[i] != CellTransform::NO_TRANS) && (adapt_op[i] != CellTransform::UNIFORM_REFINE))
    {
      cells_have_admissible_cell_transform = false;
      break;
    }

    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ActiveIdx(i));
    if (tcell.refinement_level() > 0)
    {
      const bool curr_cell_was_obtained_by_aniso_refine = CellTransformTraits::is_aniso_refinement(
          tcell.parent().cell_adapt_op().get().cell_adapt_op_tag().adapt_op_id());

      if (curr_cell_was_obtained_by_aniso_refine)
      {
        green_cells_present = true;
        break;
      }
    }
  }

  if (!cells_have_admissible_cell_transform)
  {
    std::cerr << "MeshTopologyAdaptAlgorithm::prepare_cells_for_red_green_"
                 "adapt: found cells\n"
              << "scheduled for inadmissible transformation (other than "
                 "uniform refinement\n"
              << "or do-nothing operation\n"
              << "Will not prepare cells for red-green adaptation. Quitting." << std::endl;
    return;
  }

  if (green_cells_present)
  {
    std::cerr << "MeshTopologyAdaptAlgorithm:: "
                 "prepare_cells_for_red_green_adapt:\n"
              << "green (anisotropically split) cells in the mesh. Please remove "
                 "green\n"
              << "cells before performing red-green refinement. Quitting." << std::endl;
    return;
  }

  // ------------------------------------------------------
  // PHASE 2: make sure that the mesh can be adapted with
  // hanging nodes, i.e. the 2:1 refinement rule
  // is respected
  // ------------------------------------------------------
  /*
  std::cout << "In generate_red_adapt_ops, before hanging node check is
  performed:" << std::endl; for (Uint ac = 0; ac < adapt_op.size(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ac);
    std::cout << "  Linear [" << tcell.linear_pos_idx() << "] " <<
  adapt_op[ac]
  << std::endl;
  }
  */

  prepare_cells_for_hanging_node_adapt(mesh_cells, skeleton_facets, adapt_op);

  /*
  std::cout << "In generate_red_adapt_ops, after hanging node check is
  performed:" << std::endl; for (Uint ac = 0; ac < adapt_op.size(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ac);
    std::cout << "  Linear [" << tcell.linear_pos_idx() << "] " <<
  adapt_op[ac]
  << std::endl;
  }
  */

  // ------------------------------------------------------
  // PHASE 3: mark all cells whose 2 or more edges
  // will be bisected due to neighbour refinement
  // Such cells must be also isotropically refined
  // ------------------------------------------------------

  colors.resize(nb_active_cells);

  std::vector<bool> cell_is_in_queue(nb_active_cells, false);

  // A queue of cells that still need to be checked
  std::queue<Uint> active_cells_to_check;

  // First insert all cells adjacent to those already marked
  // for refinement into the queue. These adjacent cells
  // will be checked first
  // std::cout << "Initial queue (generate_red_adapt_ops) status -
  // linear/active ids " << std::endl;
  for (Uint ac = 0; ac < nb_active_cells; ++ac)
  {
    if (adapt_op[ac] != CellTransform::NO_TRANS)
    {
      const CellTopologyView<MeshConfig> center_cell = mesh_cells.active_cell(ActiveIdx(ac));
      const std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> neighbours =
          active_neighbours(mesh_cells, skeleton_facets, center_cell);

      for (Uint n = 0; n < neighbours.size(); ++n)
      {
        const CellTopologyView<MeshConfig> neighb = std::get<0>(neighbours[n]);
        const Uint neighb_active_id               = neighb.active_idx().id();
        if ((adapt_op[neighb_active_id] == CellTransform::NO_TRANS) &&
            !cell_is_in_queue[neighb_active_id])
        {
          active_cells_to_check.push(neighb_active_id);
          cell_is_in_queue[neighb_active_id] = true;
          // std::cout << neighb.linear_pos_idx() << "/" <<
          // neighb_active_id << " ";
        }
      }
    }
  }
  // std::cout << std::endl;

  // The purpose of this loop is to parse all cells that are marked
  // for uniform refinement and for each member tcell_N of their immediate
  // neighbour set verify that tcell_N is either 1) also scheduled for uniform
  // refinement or 2) has only one edge that will be bisected due to
  // refinement of exactly
  //    one neighbour of tcell_N and mesh conformity can therefore be
  //    preserved by green-refining tcell_N (at a later stage)

  while (!active_cells_to_check.empty())
  {
    // Get the index of first cell in queue
    // and REMOVE IT FROM QUEUE
    const Uint ac = active_cells_to_check.front();
    active_cells_to_check.pop();

    cell_is_in_queue[ac] = false;

    // If the cell did not pass check, loop over
    // its neighbours and see what is the difference in their refinement
    // levels

    const CellTopologyView<MeshConfig> center_cell = mesh_cells.active_cell(ActiveIdx(ac));

    // If the current cell is supposed to be refined, check if its
    // neighbours need to be refined too

    const std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> neighbours =
        active_neighbours(mesh_cells, skeleton_facets, center_cell);

    // Now verify that this is really true by comparing the refinement
    // level of the center cell with neighbours

    Uint nb_same_level_refined_neighbours   = 0;
    Uint nb_deeper_level_refined_neighbours = 0;

    for (auto neighb_item : neighbours)
    {
      const CellTopologyView<MeshConfig> neighb_cell = std::get<0>(neighb_item);

      const CellTransform neighb_cell_transform = adapt_op[neighb_cell.active_idx().id()];

      CellTransform neighb_cell_parent_transform = CellTransform::NO_TRANS;

      if (neighb_cell.refinement_level() > 0)
      {
        const CellTopologyView<MeshConfig> neighb_cell_parent = neighb_cell.parent();
        neighb_cell_parent_transform =
            neighb_cell_parent.cell_adapt_op().get().cell_adapt_op_tag().adapt_op_id();
      }

      // If the neighbour is going to be refined, it might trigger the
      // refinement of the center cell as well

      /*
      if ((neighb_cell_transform == CellTransform::UNIFORM_REFINE) &&
          (center_cell.refinement_level() <
      neighb_cell.refinement_level()))
      */

      // if (neighb_cell_transform == CellTransform::UNIFORM_REFINE)

      if ((neighb_cell_transform == CellTransform::UNIFORM_REFINE) &&
          (neighb_cell.refinement_level() == center_cell.refinement_level()))
      {
        nb_same_level_refined_neighbours++;
      }
      else if ((neighb_cell_parent_transform == CellTransform::UNIFORM_REFINE) &&
               (neighb_cell.refinement_level() > center_cell.refinement_level()))
      {
        nb_same_level_refined_neighbours++;
      }
      else if ((neighb_cell_transform == CellTransform::UNIFORM_REFINE) &&
               (neighb_cell.refinement_level() > center_cell.refinement_level()))
      {
        nb_deeper_level_refined_neighbours++;
      }
    }

    /*
    std::cout << ">>> Cell with linear id " << center_cell.linear_pos_idx()
    << std::endl; std::cout << "    number of same-level neighbours
    scheduled for refinement: "
              << nb_same_level_refined_neighbours << std::endl;
    std::cout << "    number of deeper-level neighbours scheduled for
    refinement: "
              << nb_deeper_level_refined_neighbours << std::endl;
    */

    if ((nb_same_level_refined_neighbours > 1) || (nb_deeper_level_refined_neighbours > 0))
    {
      adapt_op[ac] = CellTransform::UNIFORM_REFINE;

      for (auto neighb_item : neighbours)
      {
        const CellTopologyView<MeshConfig> neighb_cell = std::get<0>(neighb_item);
        const CellTransform neighb_cell_transform      = adapt_op[neighb_cell.active_idx().id()];

        if ((neighb_cell_transform == CellTransform::NO_TRANS) &&
            !cell_is_in_queue[neighb_cell.active_idx().id()])
        {
          const Uint neighb_active_idx = neighb_cell.active_idx().id();
          active_cells_to_check.push(neighb_active_idx);
          cell_is_in_queue[neighb_active_idx] = true;
        }
      }
    } // If number of refined neighbours > 1

  } // While the queue is not empty

  /*
  std::cout << "In generate_red_adapt_ops, terminating with the following
  setup:" << std::endl; for (Uint ac = 0; ac < adapt_op.size(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ac);
    std::cout << "  Linear [" << tcell.linear_pos_idx() << "] " <<
  adapt_op[ac]
  << std::endl;
  }
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::generate_green_adapt_ops(
    const TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &skeleton_facets,
    std::vector<CellTransform> &adapt_op)
{
  // ------------------------------------------------------
  // This method supposes that we are left with cells which
  // have either no edge scheduled for refinement or
  // exactly one edge must be bisected because of neighbour
  // refinement.
  // To preserve conformity of the mesh, these cells
  // will be anisotropically split (refined as 'green')
  // ------------------------------------------------------

  const Uint nb_active_cells = mesh_cells.nb_active_cells();
  adapt_op.resize(nb_active_cells);
  adapt_op.assign(nb_active_cells, CellTransform::NO_TRANS);

  for (Uint ac = 0; ac < nb_active_cells; ++ac)
  {
    const CellTopologyView<MeshConfig> center_cell = mesh_cells.active_cell(ActiveIdx(ac));
    const std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> neighbours =
        active_neighbours(mesh_cells, skeleton_facets, center_cell);

    for (Uint n = 0; n < neighbours.size(); ++n)
    {
      const CellTopologyView<MeshConfig> neighb = std::get<0>(neighbours[n]);

      if (neighb.refinement_level() > 0)
      {
        const CellTopologyView<MeshConfig> neighb_parent         = neighb.parent();
        const std::tuple<Uint, Uint> rel_levels_center_nb_parent = get_relative_refinement_levels(
            center_cell.refinement_level(), neighb_parent.refinement_level());

        const bool neighb_parent_uniformly_refined =
            neighb_parent.cell_adapt_op().get().cell_adapt_op_tag().adapt_op_id() ==
            CellTransform::UNIFORM_REFINE;

        if ((std::get<0>(rel_levels_center_nb_parent) == 0) &&
            (std::get<1>(rel_levels_center_nb_parent) == 0) && neighb_parent_uniformly_refined)
        {
          const Uint local_aniso_face_id = std::get<2>(neighbours[n]);

          adapt_op[center_cell.active_idx().id()] =
              CellTransformTraits::aniso_refinement_ortho_face(local_aniso_face_id);
          break;
        }
      }
    }

  } // Loop over active cells
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::generate_ops_to_remove_green_leafs(
    const TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &skeleton_facets,
    std::vector<CellTransform> &adapt_op)
{
  adapt_op.resize(mesh_cells.nb_active_cells());
  adapt_op.assign(mesh_cells.nb_active_cells(), CellTransform::NO_TRANS);

  for (Uint ac = 0; ac < mesh_cells.nb_active_cells(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell = mesh_cells.active_cell(ac);
    const bool check_cell = (tcell.refinement_level() > 0) && (tcell.nb_children() == 0) &&
                            (adapt_op[tcell.active_idx()] == CellTransform::NO_TRANS);

    if (check_cell)
    {
      const CellTopologyView<MeshConfig> parent_tcell = tcell.parent();
      const adapt::CellAdaptOp parent_adapt_op        = parent_tcell.cell_adapt_op();

      if (CellTransformTraits::is_aniso_refinement(
              parent_adapt_op.get().cell_adapt_op_tag().adapt_op_id()))
      {
        const std::vector<CellTopologyView<MeshConfig>> children = parent_tcell.children();

        for (Uint c = 0; c < children.size(); ++c)
        {
          const CellTopologyView<MeshConfig> child = children[c];
          adapt_op[child.active_idx()]             = CellTransform::COARSEN;
        }
      } // If parent was anisotropically refined
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> MeshTopologyAdaptAlgorithm<
    MeshConfig>::active_neighbours(const TriaCells<MeshConfig> &mesh_cells,
                                   const TriaFacets<MeshConfig> &skeleton_facets,
                                   const CellTopologyView<MeshConfig> tcell)
{
  const CellTopologyView<MeshConfig> center_cell =
      mesh_cells.active_cell(ActiveIdx(tcell.active_idx()));
  const common::ArrayView<const Uint, _1D, Uint> adjacent_facet_ids = center_cell.incident_facets();

  std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> neighbours;

  for (Uint i = 0; i < adjacent_facet_ids.size(); ++i)
  {
    const TraceIncidences adj_facet = skeleton_facets.facet_data(FlatIdx(adjacent_facet_ids[i]));
    add_active_adjacent_cells(mesh_cells, skeleton_facets, center_cell, adj_facet, neighbours);
  }

  return neighbours;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::add_active_adjacent_cells(
    const TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &skeleton_facets,
    const CellTopologyView<MeshConfig> center_cell, const TraceIncidences &facet_block,
    std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> &adjacent_cells)
{
  bool block_is_adjacent_to_center_cell = false;
  for (Uint i = 0; i < facet_block.size(); ++i)
  {
    if (facet_block.cell_id(i) == center_cell.linear_pos_idx().id())
    {
      block_is_adjacent_to_center_cell = true;
    }
  }

  if (!block_is_adjacent_to_center_cell)
  {
    return;
  }

  // Determine what is the local id of the center cell - this is the facet (or
  // edge) incident to other cells listed in the block
  SUint local_incidence_idx = 0;
  for (Uint i = 0; i < facet_block.size(); ++i)
  {
    if (facet_block.cell_id(i) == center_cell.linear_pos_idx().id())
    {
      local_incidence_idx = facet_block.local_id(i);
    }
  }

  for (Uint i = 0; i < facet_block.size(); ++i)
  {
    if (facet_block.cell_id(i) != center_cell.linear_pos_idx().id())
    {
      const CellTopologyView<MeshConfig> neighbour =
          mesh_cells.cell(FlatIdx(facet_block.cell_id(i)));
      if (neighbour.status() == EntityStatus::Active)
      {
        adjacent_cells.push_back(
            std::make_tuple(neighbour, facet_block.local_id(i), local_incidence_idx));
        return;
      }
      else
      {
        const std::vector<CellTopologyView<MeshConfig>> neighbour_children = neighbour.children();

        for (Uint c = 0; c < neighbour_children.size(); ++c)
        {
          const common::ArrayView<const Uint, _1D, Uint> child_facets =
              neighbour_children[c].incident_facets();

          for (Uint f = 0; f < child_facets.size(); ++f)
          {
            const TraceIncidences facet = skeleton_facets.facet_data(FlatIdx(child_facets[f]));
            if (facet.has_linear_cell_idx(center_cell.linear_pos_idx().id()))
            {
              add_active_adjacent_cells(mesh_cells, skeleton_facets, center_cell, facet,
                                        adjacent_cells);
            }
          }

        } // Loop over neighbour children

      } // else

    } // If this cell is different from center cell

  } // Loop over facet block entries
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::check_mesh_after_topo_change(
    TriaCells<MeshConfig> &mesh_cells, TriaFacets<MeshConfig> &skeleton_facets)
{
  filter_facets(mesh_cells, skeleton_facets, EntityStatus::Active);
  // skeleton_facets.update_status(mesh_cells.cell_status());
  mesh_cells.rebuild_neighbour_information(skeleton_facets.begin(), skeleton_facets.end());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::filter_facets(const TriaCells<MeshConfig> &mesh_cells,
                                                           TriaFacets<MeshConfig> &skeleton_facets,
                                                           const EntityStatus filter_status)
{
  std::vector<FlatIdx> facets_to_remove;
  for (Uint f = 0; f < skeleton_facets.nb_all_facets(); ++f)
  {
    const TraceIncidences facet_block = skeleton_facets.facet_data(FlatIdx(f));

    bool remove_block = false;

    for (Uint i = 0; i < facet_block.size(); ++i)
    {
      const CellTopologyView<MeshConfig> tcell = mesh_cells.cell(FlatIdx(facet_block.cell_id(i)));
      if (tcell.status() != filter_status)
      {
        remove_block = true;
        break;
      }
    }
    if (remove_block)
    {
      facets_to_remove.push_back(FlatIdx(f));
    }
  }

  skeleton_facets.remove_facets(facets_to_remove);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::ensure_local_facet_numbering(
    TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &skeleton_facets)
{
  std::vector<Uint> local_facet_reordering;
  for (Uint c = 0; c < mesh_cells.nb_all_cells(); ++c)
  {
    const CellTopologyView<MeshConfig> cell = mesh_cells.cell(FlatIdx(c));
    const Uint nb_topo_facets = cell.std_region().get().nb_entities(MeshConfig::FACET_DIM);

    // The length of cell_facet_ids can be different from nb_topo_facets.
    // When a non-refined element E1 is incident to a refined element E2,
    // the topological face of E1 incident to E2 consists of several facet
    // blocks, each block representing only a portion of facet seen from E1
    const common::ArrayView<const Uint, _1D, Uint> cell_facet_ids = cell.incident_facets();
    local_facet_reordering.resize(cell_facet_ids.size());

    Uint new_facet_id = 0;
    // Loop over topological faces 0,1,2 ... of each element
    for (Uint topo_f = 0; topo_f < nb_topo_facets; ++topo_f)
    {
      // Go through all facet blocks and pick facets which form the
      // topological face topo_f
      for (Uint f = 0; f < cell_facet_ids.size(); ++f)
      {
        const TraceIncidences facet = skeleton_facets.facet_data(FlatIdx(cell_facet_ids[f]));
        for (Uint i = 0; i < facet.size(); ++i)
        {
          if ((facet.cell_id(i) == cell.linear_pos_idx().id()) && (facet.local_id(i) == topo_f))
          {
            local_facet_reordering[f] = new_facet_id;
            new_facet_id++;
          }
        } // Loop over entries of one facet block
      }   // Loop over facet blocks
    }     // Loop over topological faces

    mesh_cells.reorder_cell_facets(cell, local_facet_reordering);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint MeshTopologyAdaptAlgorithm<MeshConfig>::get_local_id_of_child_cell(
    const std::vector<CellTopologyView<MeshConfig>> &all_children,
    const CellTopologyView<MeshConfig> &child_cell)
{
  for (Uint i = 0; i < all_children.size(); ++i)
  {
    if (all_children[i].linear_pos_idx() == child_cell.linear_pos_idx())
    {
      return i;
    }
  }
  return INVALID_CELL_ID;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
bool MeshTopologyAdaptAlgorithm<MeshConfig>::check_mesh_h_grading(
    const TriaCells<MeshConfig> &mesh_cells, const TriaFacets<MeshConfig> &skeleton_facets,
    const std::vector<CellTransform> &adapt_op, const Uint max_level_diff, const bool verbose)

{
  bool mesh_passed_grading_check = true;

  for (Uint ac = 0; ac < mesh_cells.nb_active_cells(); ++ac)
  {
    const CellTopologyView<MeshConfig> center_cell = mesh_cells.active_cell(ActiveIdx(ac));
    const Uint refinement_level_curr_cell          = (adapt_op[ac] != CellTransform::NO_TRANS)
                                                ? center_cell.refinement_level() + 1
                                                : center_cell.refinement_level();

    const std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> neighbours =
        active_neighbours(mesh_cells, skeleton_facets, center_cell);

    // Now verify that this is really true by comparing the refinement
    // level of the center cell with neighbours
    for (auto neighb_item : neighbours)
    {
      const CellTopologyView<MeshConfig> neighb_cell = std::get<0>(neighb_item);
      const Uint refinement_level_neighb_cell =
          (adapt_op[neighb_cell.active_idx().id()] != CellTransform::NO_TRANS)
              ? neighb_cell.refinement_level() + 1
              : neighb_cell.refinement_level();

      const Uint max_level = std::max(refinement_level_curr_cell, refinement_level_neighb_cell);
      const Uint min_level = std::min(refinement_level_curr_cell, refinement_level_neighb_cell);
      const Uint current_level_diff = max_level - min_level;

      if (current_level_diff > max_level_diff)
      {
        mesh_passed_grading_check = false;
        if (verbose)
        {
          std::cout << "Problematic cells when checking grading: " << center_cell.linear_pos_idx()
                    << " - " << neighb_cell.linear_pos_idx() << std::endl;
        }
      }
    }
  }

  return mesh_passed_grading_check;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshTopologyAdaptAlgorithm<MeshConfig>::generate_new_facet_offsets_from_buffer(
    const internal::TriaFacets<MeshConfig> &facet_buffer, const Uint nb_old_cells_total,
    const Uint nb_cells_inserted, std::vector<Uint> &new_cell_facet_offsets,
    std::vector<Uint> &new_incident_facet_ids)
{
  new_cell_facet_offsets.resize(nb_cells_inserted + 1);
  new_cell_facet_offsets.assign(nb_cells_inserted + 1, 0);

  for (Uint f = 0; f < facet_buffer.nb_all_facets(); ++f)
  {
    const TraceIncidences facet_block = facet_buffer.facet_data(FlatIdx(f));
    for (Uint i = 0; i < facet_block.size(); ++i)
    {
      if (facet_block.cell_id(i) >= nb_old_cells_total)
      {
        const Uint shifted_cell_idx = facet_block.cell_id(i) - nb_old_cells_total;
        new_cell_facet_offsets[shifted_cell_idx + 1]++;
      }
    }
  }

  for (Uint i = 0; i < (new_cell_facet_offsets.size() - 1); ++i)
  {
    new_cell_facet_offsets[i + 1] += new_cell_facet_offsets[i];
  }

  // Update information about facets incident to new cells
  const Uint nb_new_facet_incidences = new_cell_facet_offsets.back();

  new_incident_facet_ids.resize(nb_new_facet_incidences);
  new_incident_facet_ids.assign(nb_new_facet_incidences, INVALID_FACET_ID);
  for (Uint f = 0; f < facet_buffer.nb_all_facets(); ++f)
  {
    const TraceIncidences facet_block = facet_buffer.facet_data(FlatIdx(f));
    for (Uint i = 0; i < facet_block.size(); ++i)
    {
      if (facet_block.cell_id(i) >= nb_old_cells_total)
      {
        const Uint shifted_cell_idx = facet_block.cell_id(i) - nb_old_cells_total;
        for (Uint j = new_cell_facet_offsets[shifted_cell_idx];
             j < new_cell_facet_offsets[shifted_cell_idx + 1]; ++j)
        {
          if (new_incident_facet_ids[j] == INVALID_FACET_ID)
          {
            new_incident_facet_ids[j] = f;
            break;
          }
        }
      }
    }
  }
}
// ----------------------------------------------------------------------------

} // namespace internal

} // namespace mesh

} // namespace pdekit

#endif
