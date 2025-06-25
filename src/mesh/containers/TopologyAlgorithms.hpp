#ifndef PDEKIT_Mesh_Topology_Algorithms_hpp
#define PDEKIT_Mesh_Topology_Algorithms_hpp

#include <algorithm>
#include <ctime>
#include <forward_list>
#include <list>
#include <set>
#include <tuple>

#include "common/IteratorRange.hpp"
#include "mesh/CellBuffer.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/TopologyPredicates.hpp"
#include "mesh/containers/TriaCells.hpp"
#include "mesh/containers/TriaFacets.hpp"
#include "mesh/iterators/TraceTopologyIterator.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"

namespace pdekit
{

namespace mesh
{

class TopologyAlgorithms
{

  public:
  /// Default constructor
  TopologyAlgorithms();

  /// Default destructor
  ~TopologyAlgorithms();

  /// Identify subcells in existing topology. Some subcells (faces, or edges)
  /// are passed in the array 'cells_to_identify' This algorithm tries to
  /// determine to which cells do the edges or faces belong by filling the
  /// variable 'correspondence_vector' correspondence_vector is a vector of
  /// pairs [Cell_id, Local_subcell_id]. This means that if the first pair is
  /// [34214,2], for example, the first subcell in 'cells_to_identify' is in
  /// fact the edge number 2 of cell number 34214
  /// @param topology          ... contains the cells in which we are
  /// searching for the subcells
  /// @param dim               ... dimension of subcells. If dim == 1, we are
  /// searching for edges,
  ///                                                     if dim == 2, we are
  ///                                                     searching for faces
  /// @param cells_to_identify  ... connectivity table of subcells which we
  /// are trying
  ///                                  to match against cells in 'topology'
  /// @param correspondence vector ... output vector with identified subcells

  template <typename MeshConfig>
  static void identify_subcells(
      DofMap<MeshConfig> const &cell_dofs, const Uint subcell_dim,
      CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> const &cells_to_identify,
      std::vector<IncidenceEntry> &correspondence_vector);

  template <typename MeshConfig>
  static void identify_subcells2(
      DofMap<MeshConfig> const &cell_dofs, const Uint subcell_dim,
      CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> const &cells_to_identify,
      std::vector<IncidenceEntry> &correspondence_vector);

  template <typename DofIterT1, typename DofIterT2>
  static void identify_subcells3(const common::IteratorRange<DofIterT1> &dof_range_L,
                                 const common::IteratorRange<DofIterT2> &dof_range_R,
                                 const Uint subcell_dim,
                                 std::vector<IncidenceEntry> &correspondence_vector);

  template <typename MeshConfig>
  static void compute_node_to_cell_connectivity(
      const internal::TriaCells<MeshConfig> &mesh_cells,
      const internal::TriaFacets<MeshConfig> &mesh_facets,
      common::BlockArray<std::tuple<Uint, Uint>, Uint> &node_to_cell_connectivity,
      common::BlockArray<Uint, Uint> &cell_to_p1_node);

  private:
  using TP = TopologyPredicates;

  template <typename DofIterT1, typename DofIterT2>
  static void build_subcells_hash_maps(
      const common::IteratorRange<DofIterT1> &dof_range_long,
      const common::IteratorRange<DofIterT2> &dof_range_short, const Uint subcell_dim,
      std::unordered_map<Uint, std::forward_list<std::tuple<DofIterT1, Uint>>> &cells_L_hash_map,
      std::unordered_map<Uint, std::forward_list<std::tuple<DofIterT2, Uint>>> &cells_S_hash_map);
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TopologyAlgorithms::identify_subcells(
    DofMap<MeshConfig> const &cell_dofs, const Uint subcell_dim,
    CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> const &cells_to_identify,
    std::vector<IncidenceEntry> &correspondence_vector)
{
  if (cells_to_identify.nb_active_cells() == 0)
  {
    std::cout << "TopologyAlgorithms::identify_subcells: no cells in list" << std::endl;
    return;
  }

  using dof_storage_type = DofMap<MeshConfig>;

  Uint max_node_nr = 0;
  for (typename dof_storage_type::const_dof_iterator cell_iter = cell_dofs.begin();
       cell_iter != cell_dofs.end(); ++cell_iter)
  {
    const MeshEntity &cell = *cell_iter;
    max_node_nr            = std::max(TP::max_entity_vertex(cell), max_node_nr);
  }

  /// Keep information about nodes that are also vertices of 'cells to
  /// identify' Later on, we will consider candidate cells to match the
  /// 'cells_to_identify' such that these candidate cells have
  /// at least one vertex common with 'cells_to_identify'
  std::vector<bool> is_on_cells_to_identify(max_node_nr + 1, false);

  for (Uint c = 0; c < cells_to_identify.nb_active_cells(); ++c)
  {
    const MeshEntity cell_to_ident = cells_to_identify.active_cell(ActiveIdx(c));
    for (Uint v = 0; v < cell_to_ident.nb_vert(); ++v)
    {
      is_on_cells_to_identify[cell_to_ident.vertex(v)] = true;
    }
  }

  /// For each node in the mesh, build a list of cells which stores the list
  /// of cells that have (contain) this node Note that we could take the
  /// number of rows in the coordinates matrix to obtain the number of nodes.
  /// However, if the coordinates are not filled yet, this would give 0. To
  /// make the algorithm purely topological, we compute the maximum node
  /// number by checking all cells.
  // std::vector<std::list<Uint> > node_to_cell (
  // geometry.geo_coord().nb_rows()
  // );
  std::vector<std::list<Uint>> node_to_cell(max_node_nr + 1);

  for (typename dof_storage_type::const_dof_iterator cell_iter = cell_dofs.begin();
       cell_iter != cell_dofs.end(); ++cell_iter)
  {
    const MeshEntity &cell = *cell_iter;
    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      if (is_on_cells_to_identify[cell.vertex(n)])
      {
        node_to_cell[cell.vertex(n)].push_back(cell.idx());
      }
    }
  }

  // const typename TopologyType::cell_connectivity_type & cells =
  // topology.cells();

  Uint nb_of_identified_bcells = 0;
  bool found_matching_cell     = false;
  correspondence_vector.resize(cells_to_identify.nb_active_cells());

  EntityDofRealign p;

  for (Uint bc_id = 0; bc_id < cells_to_identify.nb_active_cells(); ++bc_id)
  {
    const MeshEntity boundary_cell = cells_to_identify.active_cell(ActiveIdx(bc_id));

    const Uint min_vertex_in_bcell = TP::min_entity_vertex(boundary_cell);
    found_matching_cell            = false;

    for (std::list<Uint>::const_iterator candidate_cell_iter =
             node_to_cell[min_vertex_in_bcell].begin();
         candidate_cell_iter != node_to_cell[min_vertex_in_bcell].end(); ++candidate_cell_iter)
    {
      const MeshEntity candidate_cell = cell_dofs.active_cell(ActiveIdx(*candidate_cell_iter));

      for (SUint s = 0; s < candidate_cell.nb_sub_elements(subcell_dim); ++s)
      {
        const MeshEntity candidate_subcell = candidate_cell.sub_entity(subcell_dim, s);

        if (TP::entities_match(boundary_cell, candidate_subcell, p) ||
            TP::entities_match_reverse(boundary_cell, candidate_subcell,
                                       p)) /// This should not
                                           /// be needed?
        {
          // std::cout << "\tThe boundary cell " << boundary_cell << "
          // corresponds to [" << candidate_cell.idx()
          //          << "," << s << "]" << std::endl;
          correspondence_vector[nb_of_identified_bcells].cell_idx = candidate_cell.idx();
          correspondence_vector[nb_of_identified_bcells].local_id = s;
          nb_of_identified_bcells++;
          found_matching_cell = true;
          break;

        } // if

        if (found_matching_cell)
          break;

      } // loop over subcells of candidate cell

      if (found_matching_cell)
        break;

    } // loop over candidate cells

    if (!found_matching_cell)
    {
      std::cerr << "TopologyAlgorithms::identify_subcells: boundary cell " << boundary_cell
                << " not identified" << std::endl;
    }

  } // loop over boundary cells

  if (nb_of_identified_bcells < cells_to_identify.nb_active_cells())
  {
    std::cerr << "TopologyAlgorithms::identify_subcells: didn't manage to "
                 "identify all subcells. The provided subcell "
              << "list contains " << cells_to_identify.nb_active_cells()
              << " subcells to be identified, but I managed to identify "
              << "only " << nb_of_identified_bcells << " of them" << std::endl;
  }
  else
  {
    std::cout << "TopologyAlgorithms::identify_subcells: successfully identified "
              << nb_of_identified_bcells << " subcells "
              << "out of " << cells_to_identify.nb_active_cells() << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TopologyAlgorithms::identify_subcells2(
    DofMap<MeshConfig> const &cell_dofs, const Uint subcell_dim,
    CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> const &cells_to_identify,
    std::vector<IncidenceEntry> &correspondence_vector)
{
  if (cells_to_identify.nb_active_cells() == 0)
  {
    std::cout << "TopologyAlgorithms::identify_subcells: no cells in list" << std::endl;
    return;
  }

  using dof_storage_type = DofMap<MeshConfig>;

  Uint max_node_nr = 0;
  for (Uint c = 0; c < cells_to_identify.nb_active_cells(); ++c)
  {
    const MeshEntity cell = cells_to_identify.active_cell(ActiveIdx(c));
    max_node_nr           = std::max(TP::max_entity_vertex(cell), max_node_nr);
  }

  /// Keep information about nodes that are also vertices of 'cells to
  /// identify' Later on, we will consider candidate cells to match the
  /// 'cells_to_identify' such that these candidate cells have
  /// at least one vertex common with 'cells_to_identify'
  std::vector<bool> is_on_cells_to_identify(max_node_nr + 1, false);

  for (Uint c = 0; c < cells_to_identify.nb_active_cells(); ++c)
  {
    const MeshEntity cell_to_ident = cells_to_identify.active_cell(ActiveIdx(c));
    for (Uint v = 0; v < cell_to_ident.nb_vert(); ++v)
    {
      is_on_cells_to_identify[cell_to_ident.vertex(v)] = true;
    }
  }

  /// For each node in the CellBuffer, build a list of cells which stores the
  /// list of cells that have (contain) this node Note that we could take the
  /// number of rows in the coordinates matrix to obtain the number of nodes.
  /// However, if the coordinates are not filled yet, this would give 0. To
  /// make the algorithm purely topological, we compute the maximum node
  /// number by checking all cells.

  std::vector<std::list<Uint>> node_to_cell(max_node_nr + 1);

  for (Uint c = 0; c < cells_to_identify.nb_active_cells(); ++c)
  {
    const MeshEntity cell = cells_to_identify.active_cell(ActiveIdx(c));
    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      if (is_on_cells_to_identify[cell.vertex(n)])
      {
        node_to_cell[cell.vertex(n)].push_back(cell.idx());
      }
    }
  }

  // const typename TopologyType::cell_connectivity_type & cells =
  // topology.cells();

  Uint nb_of_identified_bcells = 0;
  correspondence_vector.resize(cells_to_identify.nb_active_cells());
  for (Uint i = 0; i < correspondence_vector.size(); ++i)
  {
    correspondence_vector[i].cell_idx = INVALID_CELL_ID;
    correspondence_vector[i].local_id = 0;
  }

  EntityDofRealign p;

  /*
  std::cout <<
  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
            << std::endl;
  std::cout << "Number of cells in buffer = " << cells_to_identify.nb_cells()
  << std::endl;

  for (Uint i = 0; i < node_to_cell.size(); ++i)
  {
    const std::list<Uint> &node_list = node_to_cell[i];
    if (!node_list.empty())
    {
      std::cout << "{" << i << "}";
      for (std::list<Uint>::const_iterator it = node_list.cbegin(); it !=
  node_list.cend(); ++it)
      {
        std::cout << " " << *it;
      }
      std::cout << std::endl;
    }
  }


  std::cout <<
  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
            << std::endl;
  */
  for (typename dof_storage_type::const_dof_iterator cell_iter = cell_dofs.begin();
       cell_iter != cell_dofs.end(); ++cell_iter)
  {
    const MeshEntity interior_cell = cell_iter->mesh_entity();

    bool interior_cell_has_vert_on_boundary = false;

    // Go through vertices of this cell and check whether it has at least
    // one vertex that coincides with some vertices of the boundary cells
    // (i.e. 'cells to identify')

    for (Uint v = 0; v < interior_cell.nb_vert(); ++v)
    {
      const Uint interior_cell_vert_nr = interior_cell.vertex(v);
      if ((interior_cell_vert_nr <= max_node_nr) &&
          (is_on_cells_to_identify[interior_cell_vert_nr]))
      {
        interior_cell_has_vert_on_boundary = true;
        break;
      }
    }

    // std::cout << "Cell [" << interior_cell.idx() << "] " <<
    // interior_cell;

    if (interior_cell_has_vert_on_boundary)
    {
      // std::cout << " has vertex on boundary" << std::endl;

      for (SUint s = 0; s < interior_cell.nb_sub_elements(subcell_dim); ++s)
      {
        const MeshEntity interior_subcell          = interior_cell.sub_entity(subcell_dim, s);
        const Uint min_vertex_in_candidate_subcell = TP::min_entity_vertex(interior_subcell);

        // Here we check again whether min_vertex_in_candidate_subcell
        // is still on boundary or not
        if ((min_vertex_in_candidate_subcell <= max_node_nr) &&
            (is_on_cells_to_identify[min_vertex_in_candidate_subcell]))
        {

          EntityDofRealign p;

          // std::cout << "Interior subcell = " << interior_subcell <<
          // std::endl; std::cout << "Min vertex in candidate subcell
          // = " << min_vertex_in_candidate_subcell
          //           << std::endl;
          for (std::list<Uint>::const_iterator candidate_cell_iter =
                   node_to_cell[min_vertex_in_candidate_subcell].begin();
               candidate_cell_iter != node_to_cell[min_vertex_in_candidate_subcell].end();
               ++candidate_cell_iter)
          {
            const MeshEntity boundary_cell =
                cells_to_identify.active_cell(ActiveIdx(*candidate_cell_iter));

            if (TP::entities_match(interior_subcell, boundary_cell, p) ||
                TP::entities_match_reverse(interior_subcell, boundary_cell,
                                           p)) /// This should not
                                               /// be needed?
            {
              /*
              std::cout << "\tThe boundary cell " << boundary_cell
              << " corresponds to [" << interior_cell.idx()
                       << "," << s << "] => " << interior_cell <<
              std::endl;
              */
              // correspondence_vector[nb_of_identified_bcells].cell_idx
              // = interior_cell.idx();
              // correspondence_vector[nb_of_identified_bcells].local_id
              // = s;

              correspondence_vector[boundary_cell.idx()].cell_idx = interior_cell.idx();
              correspondence_vector[boundary_cell.idx()].local_id = s;
              nb_of_identified_bcells++;
            } // if
          }   // Loop over boundary faces with incident minimum vertex
        }

      } // Loop over facets of interior cell which touches boundary

    } // If interior cell touches boundary
      //    else
      //    {
      //      std::cout << " does NOT have vertex on boundary" << std::endl;
      //    }

    /*
    if (!found_matching_cell)
    {
      std::cerr << "TopologyAlgorithms::identify_subcells: boundary cell "
    << boundary_cell
                << " not identified" << std::endl;
    }
    */

  } // loop over interior cells

  if (nb_of_identified_bcells < cells_to_identify.nb_active_cells())
  {
    std::cerr << "TopologyAlgorithms::identify_subcells: didn't manage to "
                 "identify all subcells. The provided subcell "
              << "list contains " << cells_to_identify.nb_active_cells()
              << " subcells to be identified, but I managed to identify "
              << "only " << nb_of_identified_bcells << " of them" << std::endl;
  }
  else
  {
    std::cout << "TopologyAlgorithms::identify_subcells: successfully identified "
              << nb_of_identified_bcells << " subcells "
              << "out of " << cells_to_identify.nb_active_cells() << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <typename DofIterT1, typename DofIterT2>
void TopologyAlgorithms::identify_subcells3(const common::IteratorRange<DofIterT1> &dof_range_L,
                                            const common::IteratorRange<DofIterT2> &dof_range_R,
                                            const Uint subcell_dim,
                                            std::vector<IncidenceEntry> &correspondence_vector)
{
  const DofIterT1 cells_L_begin = dof_range_L.begin();
  const DofIterT1 cells_L_end   = dof_range_L.end();

  const DofIterT2 cells_R_begin = dof_range_R.begin();
  const DofIterT2 cells_R_end   = dof_range_R.end();

  DofIterT1 iter_cells_L = cells_L_begin;
  DofIterT2 iter_cells_R = cells_R_begin;

  Uint nb_cells_L     = 0;
  Uint nb_sub_cells_L = 0;

  Uint nb_cells_R     = 0;
  Uint nb_sub_cells_R = 0;

  for (; iter_cells_L != cells_L_end; ++iter_cells_L)
  {
    nb_cells_L++;
    const mesh::MeshEntity cell_L = iter_cells_L->mesh_entity();
    nb_sub_cells_L += cell_L.nb_sub_elements(subcell_dim);
  }

  for (; iter_cells_R != cells_R_end; ++iter_cells_R)
  {
    nb_cells_R++;
    const mesh::MeshEntity cell_R = iter_cells_R->mesh_entity();
    nb_sub_cells_R += cell_R.nb_sub_elements(subcell_dim);
  }

  // Types of facet lists in left and right container
  using facet_list_L_t = std::forward_list<std::tuple<DofIterT1, Uint>>;
  using facet_list_R_t = std::forward_list<std::tuple<DofIterT2, Uint>>;

  using facet_hash_map_L_t = std::unordered_map<Uint, facet_list_L_t>;
  using facet_hash_map_R_t = std::unordered_map<Uint, facet_list_R_t>;

  // Hash map for left cell range
  facet_hash_map_L_t cells_L_hash_map;
  // Hash map for right cell range
  facet_hash_map_R_t cells_R_hash_map;

  if (nb_sub_cells_L <= nb_sub_cells_R)
  {
    build_subcells_hash_maps(dof_range_R, dof_range_L, subcell_dim, cells_R_hash_map,
                             cells_L_hash_map);
    if (cells_R_hash_map.empty())
    {
      std::cout << "TopologyAlgorithms::identify_subcells3: no intersecting cells \n"
                << "in containers - right hash map empty" << std::endl;
      return;
    }
  }
  else
  {
    build_subcells_hash_maps(dof_range_L, dof_range_R, subcell_dim, cells_L_hash_map,
                             cells_R_hash_map);

    if (cells_L_hash_map.empty())
    {
      std::cout << "TopologyAlgorithms::identify_subcells3: no intersecting cells \n"
                << "in containers - left hash map empty" << std::endl;
      return;
    }
  }

  std::vector<IncidenceEntry> buffer_vector;
  std::vector<Uint> final_position_table;

  EntityDofRealign p;

  Uint nb_matching_sub_cells = 0;

  for (auto it_map_L = cells_L_hash_map.cbegin(); it_map_L != cells_L_hash_map.cend(); ++it_map_L)
  {
    const Uint min_dof_L = it_map_L->first;
    const auto it_map_R  = cells_R_hash_map.find(min_dof_L);

    // Search for matching facets if right list is not empty
    if (it_map_R != cells_R_hash_map.cend())
    {
      const auto &list_L = it_map_L->second;
      const auto &list_R = it_map_R->second;

      auto iter_list_R = list_R.cbegin();

      for (auto iter_list_L = list_L.cbegin(); iter_list_L != list_L.cend(); ++iter_list_L)
      {
        const std::tuple<DofIterT1, Uint> &tuple_L = *iter_list_L;
        DofIterT1 cell_it_L                        = std::get<0>(tuple_L);
        const mesh::MeshEntity cell_L              = cell_it_L->mesh_entity();
        const Uint facet_id_L                      = std::get<1>(tuple_L);
        const mesh::MeshEntity facet_L             = cell_L.sub_entity(subcell_dim, facet_id_L);

        iter_list_R = list_R.cbegin();

        for (; iter_list_R != list_R.cend(); ++iter_list_R)
        {
          const std::tuple<DofIterT2, Uint> &tuple_R = *iter_list_R;
          DofIterT2 cell_it_R                        = std::get<0>(tuple_R);
          const mesh::MeshEntity cell_R              = cell_it_R->mesh_entity();
          const Uint facet_id_R                      = std::get<1>(tuple_R);
          const mesh::MeshEntity facet_R             = cell_R.sub_entity(subcell_dim, facet_id_R);

          if (TP::entities_match(facet_L, facet_R, p) ||
              TP::entities_match_reverse(facet_L, facet_R, p)) /// This should not
                                                               /// be needed?
          {
            buffer_vector.push_back(IncidenceEntry(cell_L.idx(), facet_id_L));
            final_position_table.push_back(facet_R.idx());
            nb_matching_sub_cells++;
          }

        } // Iteration over right list

      } // Iteration over left list

    } // If right list is not empty

  } // Iteration over left map

  correspondence_vector.resize(buffer_vector.size());
  for (Uint i = 0; i < buffer_vector.size(); ++i)
  {
    const Uint target_idx             = final_position_table[i];
    correspondence_vector[target_idx] = buffer_vector[i];
  }

  std::cout << "TopologyAlgorithms::identify_subcells3: successfully identified "
            << nb_matching_sub_cells << " subcells of dimension " << subcell_dim << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TopologyAlgorithms::compute_node_to_cell_connectivity(
    const internal::TriaCells<MeshConfig> &mesh_cells,
    const internal::TriaFacets<MeshConfig> &mesh_facets,
    common::BlockArray<std::tuple<Uint, Uint>, Uint> &node_to_cell_connectivity,
    common::BlockArray<Uint, Uint> &cell_to_p1_node)
{
  // ---------------------------------------------------------------------------
  // PART 1: COMPUTE NODE-TO-CELL CONNECTIVITY - represent each P1 node in
  // the mesh as pair (cell index, local node index in cell) and put all pairs
  // of nodes that are physically located at the same position (i.e. they're
  // incident) into the same block. Each block of node_to_cell_connectivity
  // then stores all incident nodes
  // ---------------------------------------------------------------------------

  // using const_skeleton_iterator = typename
  // internal::TriaFacets<MeshConfig>::const_iterator;

  // This type is supposed to hold a list of pairs (active_cell_id, node_id)
  // such that they are incident (have the same location)
  using node_chain_type = std::forward_list<std::tuple<Uint, Uint>>;

  // The i-th block in node_to_chain_id represents all p1 nodes in the
  // i-th cell. The values in the block are chain ids
  common::BlockArray<Uint, Uint> node_to_chain_id;

  Uint tot_nb_p1_nodes = 0;
  for (Uint c = 0; c < mesh_cells.nb_active_cells(); ++c)
  {
    const CellTopologyView<MeshConfig> topo_cell = mesh_cells.active_cell(ActiveIdx(c));
    const StdRegion std_reg                      = topo_cell.pt_set_id();
    tot_nb_p1_nodes += std_reg.get().nb_p1_nodes();
  }

  // Form initial chains: each chain has exactly
  // one tuple (c,i), where c is an index of active cell
  // and i is an index of one p1 vertex cell in that cell
  // Suppose we have a mesh with 2 triangles and 1 quadrilateral.
  // Then the initial chains should look as follows:
  //
  // First triangle
  // Chain 0: {0, 0} // Element 0, vertex 0
  // Chain 1: {0, 1} // Element 0, vertex 1
  // Chain 2: {0, 2} // Element 0, vertex 2
  //
  // Second triangle:
  // Chain 3: {1, 0} // Element 1, vertex 0
  // Chain 4: {1, 1} // Element 1, vertex 1
  // Chain 5: {1, 2} // Element 1, vertex 2
  //
  // First quadrilateral
  // Chain 6: {2, 0} // Element 2, vertex 0
  // Chain 7: {2, 1} // Element 2, vertex 1
  // Chain 8: {2, 2} // Element 2, vertex 2
  // Chain 9: {2, 3} // Element 2, vertex 3
  //
  // Note: The objective of this function is to merge all chains
  // corresponding to vertices with the same physical location together

  std::vector<Uint> chain_ids;
  Uint current_chain_id = 0;

  std::vector<node_chain_type> node_chains;
  node_chains.resize(tot_nb_p1_nodes);

  for (Uint c = 0; c < mesh_cells.nb_active_cells(); ++c)
  {
    const CellTopologyView<MeshConfig> topo_cell = mesh_cells.active_cell(ActiveIdx(c));
    const StdRegion std_reg                      = topo_cell.pt_set_id();

    const Uint nb_p1_vert_in_cell = std_reg.get().nb_p1_nodes();

    chain_ids.resize(nb_p1_vert_in_cell);
    for (Uint i = 0; i < nb_p1_vert_in_cell; ++i)
    {
      chain_ids[i] = current_chain_id;
      node_chains[current_chain_id].push_front(std::make_tuple(c, i));
      current_chain_id++;
    }
    const common::ArrayView<const Uint, _1D, Uint> block_view(chain_ids.data(), chain_ids.size());

    node_to_chain_id.create_back_block(block_view.size());
    node_to_chain_id.fill_last_block(block_view);
  }

  /*
  std::cout << "Chains: " << std::endl;

  for (Uint c = 0; c < node_chains.size(); ++c)
  {
    std::cout << "Chain {" << c << "}:";

    for (auto chain_elem : node_chains[c])
    {
      std::cout << " [" << std::get<0>(chain_elem) << "," <<
  std::get<1>(chain_elem) << "]";
    }
    std::cout << std::endl;
  }

  std::cout << node_to_chain_id << std::endl;
  */

  node_to_chain_id.reserve(tot_nb_p1_nodes, mesh_cells.nb_active_cells());

  // std::cout << "Nb of active cells = " << mesh_cells.nb_active_cells() <<
  // std::endl; std::cout << "Nb of p1 nodes = " << tot_nb_p1_nodes <<
  // std::endl;

  TraceEntityTuple facet_tuple;

  for (Uint f = 0; f < mesh_facets.nb_active_facets(); ++f)
  {
    const TraceIncidences facet_incidences = mesh_facets.active_facet_data(ActiveIdx(f));
    if (facet_incidences.size() == 2)
    {
      // std::cout << "Processing facet data " << facet_incidences <<
      // std::endl;
      const CellTopologyView<MeshConfig> tcell_L =
          mesh_cells.cell(FlatIdx(facet_incidences.cell_id(LEFT)));
      const CellTopologyView<MeshConfig> tcell_R =
          mesh_cells.cell(FlatIdx(facet_incidences.cell_id(RIGHT)));

      if ((tcell_L.status() == EntityStatus::Active) && (tcell_R.status() == EntityStatus::Active))
      {
        const EntityRealignCode pcode_L = facet_incidences.permutation(LEFT).get().code();
        const EntityRealignCode pcode_R = facet_incidences.permutation(RIGHT).get().code();

        const bool faces_are_conformal = ((pcode_L.adapt_op_id() == CellTransform::NO_TRANS) &&
                                          (pcode_R.adapt_op_id() == CellTransform::NO_TRANS));

        // We can't simply get the refinement levels by directly calling
        // get_relative_refinement_levels
        // In red-green adapted mesh, faces are always conformal, i.e.
        // their relative levels are (0,0), but if one of the neighbours
        // is a 'red' cell and the other one is 'green', then the green
        // neighbour has relative refinement level higher by 1 compared
        // to the red neighbour Therefore we first check whether the
        // faces are conformal: if they are, then the relative levels
        // have to be equal to zero. If the faces are not conformal,
        // then we rely on relative refinement levels
        const std::tuple<Uint, Uint> relative_ref_levels =
            faces_are_conformal ? std::make_tuple(0u, 0u)
                                : get_relative_refinement_levels(tcell_L.refinement_level(),
                                                                 tcell_R.refinement_level());

        facet_tuple.change_type(pcode_L, std::get<LEFT>(relative_ref_levels), pcode_R,
                                std::get<RIGHT>(relative_ref_levels));

        const std::vector<std::pair<Uint, Uint>> &incident_p1_nodes =
            facet_tuple.get().incident_p1_nodes();

        /*
        for (auto p : incident_p1_nodes)
        {
          std::cout << "[" << p.first << "," << p.second << "]" <<
        std::endl;
        }
        std::cout << "-----------------------------------" << std::endl;
        */

        const StdRegion std_reg_L = tcell_L.pt_set_id();
        const StdRegion std_reg_R = tcell_R.pt_set_id();

        const SUint local_id_L = facet_incidences.local_id(LEFT);
        const SUint local_id_R = facet_incidences.local_id(RIGHT);

        const std::shared_ptr<StdRegionEntity const> std_reg_entity_L =
            std_reg_L.get().elem_entity(MeshConfig::TDIM - 1, local_id_L);
        const std::shared_ptr<StdRegionEntity const> std_reg_entity_R =
            std_reg_R.get().elem_entity(MeshConfig::TDIM - 1, local_id_R);

        /*
        for (auto p : incident_p1_nodes)
        {
          std::cout << "[" << (*std_reg_entity_L).vertex(p.first) << ","
                    << (*std_reg_entity_R).vertex(p.second) << "]" <<
        std::endl;
        }
        std::cout << "-----------------------------------" << std::endl;
        */

        // Take each pair of incident p1 nodes and chain them
        for (auto p : incident_p1_nodes)
        {
          const Uint p1_node_id_L = (*std_reg_entity_L).vertex(p.first);
          const Uint p1_node_id_R = (*std_reg_entity_R).vertex(p.second);

          common::ArrayView<const Uint, _1D, Uint> chain_block_L =
              node_to_chain_id.const_block(tcell_L.active_idx().id());
          common::ArrayView<const Uint, _1D, Uint> chain_block_R =
              node_to_chain_id.const_block(tcell_R.active_idx().id());

          /*
          std::cout << "Cell(" << tcell_L.active_idx() << "), node("
          << p1_node_id_L
                    << ") *** Cell(" << tcell_R.active_idx() << "),
          node(" << p1_node_id_R << ")"
                    << std::endl;
          */

          // If the left and right vertices belong to different
          // chains, merge the chains. The chain with higher id will
          // be moved to be merged with the chain with the lower id
          // (the preference is towards chains with lower id by
          // convention)
          const Uint chain_id_min =
              std::min(chain_block_L[p1_node_id_L], chain_block_R[p1_node_id_R]);
          const Uint chain_id_max =
              std::max(chain_block_L[p1_node_id_L], chain_block_R[p1_node_id_R]);

          if (chain_id_min != chain_id_max)
          {
            node_chain_type &chain_max = node_chains[chain_id_max];
            for (auto node : chain_max)
            {
              node_to_chain_id.insert_value_in_block(std::get<0>(node), std::get<1>(node),
                                                     chain_id_min);
            }

            node_chain_type &chain_min = node_chains[chain_id_min];
            chain_min.merge(chain_max);
            // std::sort(chain_min.cbegin(), chain_min.cend()); //
            // can not use std::sort with std::forward_list: list
            // does not have random access iterators
            chain_min.sort();
          }
        }
        // std::cout << "-----------------------------------------" <<
        // std::endl;
      }
    } // If this facet has 2 cells
  }

  Uint nb_blocks_out = 0;

  for (Uint i = 0; i < node_chains.size(); ++i)
  {
    node_chain_type const &chain = node_chains[i];
    if (!chain.empty())
    {
      nb_blocks_out++;

      /*
      for (auto node : chain)
      {
        std::cout << "[" << std::get<0>(node) << "," << std::get<1>(node)
      << "]
      ";
      }
      std::cout << std::endl;
      */
    }
  }

  node_to_cell_connectivity.reserve(tot_nb_p1_nodes, nb_blocks_out);
  node_to_cell_connectivity.resize(0, 0);

  std::vector<std::tuple<Uint, Uint>> buffer;

  for (Uint i = 0; i < node_chains.size(); ++i)
  {
    node_chain_type const &chain = node_chains[i];
    if (!chain.empty())
    {
      buffer.resize(0);
      for (auto node : chain)
      {
        buffer.push_back(node);
      }
      common::ArrayView<const std::tuple<Uint, Uint>, _1D, Uint> buffer_view(buffer.data(),
                                                                             buffer.size());

      node_to_cell_connectivity.create_back_block(buffer_view.size());
      node_to_cell_connectivity.fill_last_block(buffer_view);
    }
  }

  // ---------------------------------------------------------------------------
  // PART 2: COMPUTE THE INVERSE MAPPING - another BlockArray where each block
  // represents all P1 vertices of one ACTIVE cell. The values in the block
  // are position ids of blocks in cell_to_p1_node_array_pos containing the P1
  // verts
  // ---------------------------------------------------------------------------

  // This array knows where in node_to_cell_connectivity can each block
  // corresponding to one p1 vertex of cell be found.
  //
  // Example: if the first block corresponds to triangle, then the first 3
  // entries [10 3 17 ...] would mean that:
  // - the nodes incident to local node 0 of this triangle can be found in
  // block 10 of
  //  node_to_cell_connectivity
  // - the nodes incident to local node 1 of the triangle can be found in
  // block 3 of
  //   cell_to_p1_node_array_pos etc.
  std::unique_ptr<std::vector<Uint>> cell_to_p1_node_array_pos(new std::vector<Uint>());

  // Each entry of this array knows how many p1 nodes does given cell have
  // nb_p1_nodes_in_cell[i] says how many p1 nodes the i-th cell has
  std::unique_ptr<std::vector<Uint>> nb_p1_nodes_in_cell(new std::vector<Uint>());

  cell_to_p1_node_array_pos->resize(node_to_cell_connectivity.size());
  nb_p1_nodes_in_cell->resize(mesh_cells.nb_active_cells());

  std::fill(nb_p1_nodes_in_cell->begin(), nb_p1_nodes_in_cell->end(), 0u);

  // Helper array which says where in in cell_to_p1_node_array_pos is the next
  // free position for given cell
  std::vector<Uint> fill_offset(mesh_cells.nb_active_cells());
  std::fill(fill_offset.begin(), fill_offset.end(), 0u);

  for (Uint b = 0; b < node_to_cell_connectivity.nb_blocks(); ++b)
  {
    const common::ArrayView<const std::tuple<Uint, Uint>, _1D, Uint> incident_verts =
        node_to_cell_connectivity.const_block(b);
    for (Uint i = 0; i < incident_verts.size(); ++i)
    {
      (*nb_p1_nodes_in_cell)[std::get<0>(incident_verts[i])]++;
    }
  }

  // Compute offsets for filling of cell_to_p1_node_array_pos
  fill_offset[0] = 0;
  for (Uint c = 1; c < fill_offset.size(); ++c)
  {
    fill_offset[c] = fill_offset[c - 1] + (*nb_p1_nodes_in_cell)[c];
  }

  // Loop over blocks in node_to_cell_connectivity to fill
  // cell_to_p1_node_array_pos
  for (Uint b = 0; b < node_to_cell_connectivity.nb_blocks(); ++b)
  {
    const common::ArrayView<const std::tuple<Uint, Uint>, _1D, Uint> incident_verts =
        node_to_cell_connectivity.const_block(b);
    for (Uint i = 0; i < incident_verts.size(); ++i)
    {
      const Uint cell_id    = std::get<0>(incident_verts[i]);
      const Uint p1_vert_id = std::get<1>(incident_verts[i]);
      (*cell_to_p1_node_array_pos)[fill_offset[cell_id] + p1_vert_id] = b;
    }
  }

  cell_to_p1_node.build(std::move(cell_to_p1_node_array_pos), std::move(nb_p1_nodes_in_cell));
}

// ----------------------------------------------------------------------------

template <typename DofIterT1, typename DofIterT2>
void TopologyAlgorithms::build_subcells_hash_maps(
    const common::IteratorRange<DofIterT1> &dof_range_long,
    const common::IteratorRange<DofIterT2> &dof_range_short, const Uint subcell_dim,
    std::unordered_map<Uint, std::forward_list<std::tuple<DofIterT1, Uint>>> &cells_long_hash_map,
    std::unordered_map<Uint, std::forward_list<std::tuple<DofIterT2, Uint>>> &cells_short_hash_map)
{
  cells_long_hash_map.clear();
  cells_short_hash_map.clear();

  const DofIterT1 cells_L_begin = dof_range_long.begin();
  const DofIterT1 cells_L_end   = dof_range_long.end();

  const DofIterT2 cells_S_begin = dof_range_short.begin();
  const DofIterT2 cells_S_end   = dof_range_short.end();

  DofIterT2 iter_cells_S = cells_S_begin;

  for (; iter_cells_S != cells_S_end; ++iter_cells_S)
  {
    const mesh::MeshEntity cell = iter_cells_S->mesh_entity();
    const Uint nb_sub_cells     = cell.nb_sub_elements(subcell_dim);

    for (Uint isub = 0; isub < nb_sub_cells; ++isub)
    {
      const mesh::MeshEntity sub_cell = cell.sub_entity(subcell_dim, isub);
      const Uint min_vert_id          = TP::TopologyPredicates::min_entity_vertex(sub_cell);
      cells_short_hash_map[min_vert_id].emplace_front(std::make_tuple(iter_cells_S, isub));
    }
  }

  DofIterT1 iter_cells_L = cells_L_begin;

  for (; iter_cells_L != cells_L_end; ++iter_cells_L)
  {
    const mesh::MeshEntity cell = iter_cells_L->mesh_entity();
    const Uint nb_sub_cells     = cell.nb_sub_elements(subcell_dim);

    for (Uint isub = 0; isub < nb_sub_cells; ++isub)
    {
      const mesh::MeshEntity sub_cell = cell.sub_entity(subcell_dim, isub);
      const Uint min_vert_id          = TP::TopologyPredicates::min_entity_vertex(sub_cell);

      if (cells_short_hash_map.count(min_vert_id))
      {
        cells_long_hash_map[min_vert_id].emplace_front(std::make_tuple(iter_cells_L, isub));
      }
    }
  }
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
