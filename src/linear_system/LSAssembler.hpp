#ifndef PDEKIT_Linear_System_Assembler_hpp
#define PDEKIT_Linear_System_Assembler_hpp

#include <chrono>
#include <ctime>
#include <tuple>
#include <vector>

#include "graph/Graph.hpp"
#include "math/BlockMatrixSparsityPattern.hpp"
#include "mesh/EntityDofRealign.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/containers/DofMap.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"

namespace pdekit
{

namespace ls
{

class LSAssembler
{
  public:
  /// Default constructor
  LSAssembler();

  /// Deleted copy constructor
  LSAssembler(const LSAssembler &other) = delete;

  /// Destructor
  ~LSAssembler();

  /// Deleted assignment operator
  LSAssembler &operator=(const LSAssembler &rhs) = delete;

  /// Perform mock assembly to generate sparsity pattern of given system of
  /// equations on the dof handler. This considers only ONE component of the
  /// system of PDEs In addition, the discretization is assumed to be
  /// CONTINUOUS
  template <typename MeshConfig>
  void build_dof_sparsity_pattern_continuous(const typename result_of::dof_map_t<MeshConfig> &dofs,
                                             graph::Graph<Int> &ls_dofs);

  /// Perform mock assembly to generate sparsity pattern of given system of
  /// equations on the dof handler. This considers only ONE component of the
  /// system of PDEs In addition, the discretization is assumed to be
  /// DISCONTINUOUS
  template <typename MeshConfig>
  void build_dof_sparsity_pattern_discontinuous(
      const mesh::Tria<MeshConfig> &cell_topology,
      const typename result_of::dof_map_t<MeshConfig> &dofs, graph::Graph<Int> &ls_dofs);

  /// Perform mock assembly to find out how big is the matrix fill going to be
  /// Continuous case
  template <typename MeshConfig>
  void build_matrix_sparsity_pattern_continuous(
      const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &dofs,
      graph::Graph<Int> &ls_dofs);

  /// Perform mock assembly to find out how big is the matrix fill going to be
  /// Discontinuous case
  template <typename MeshConfig>
  void build_matrix_sparsity_pattern_discontinuous(
      const Uint NEQ, const mesh::Tria<MeshConfig> &cell_topology,
      const typename result_of::dof_map_t<MeshConfig> &dofs, graph::Graph<Int> &ls_dofs);

  /// Build a block sparsity pattern corresponding to a discontinuous method
  void build_matrix_sparsity_pattern_block_discontinuous(
      const common::BlockArray<Uint, Uint> mesh_dual_graph_crs,
      const std::vector<Uint> &block_sizes, const std::vector<Uint> &cell_reordering,
      math::BlockMatrixSparsityPattern<Int> &coarse_matrix_sparsity, const std::string &name);

  /// Set the mapping between local and global
  /// dof numbering
  void set_local_to_global_dof_map(const Uint NEQ, const mesh::MeshEntity &entity);

  /// Set all matrix entries to given value
  void set_local_mat_entries(const Real value = 0.0);

  /// Set all rhs entries to given values
  void set_local_rhs_entries(const Real value = 0.0);

  /// Insert value to element matrix.
  /// This overwrites previous value in the matrix
  /// @param inode ... local node idx
  /// @param ieq   ... index of equation (in case this is a system of PDEs)
  /// @param jnode ... local node idx
  /// @param jeq   ... index of equation (in case this is a system of PDEs)
  /// @param value ... value that should be inserted
  void insert_in_local_mat(const Uint NEQ, const Uint inode, const Uint ieq, const Uint jnode,
                           const Uint jeq, const Real value);

  /// Add value to element matrix
  /// This accumulates to previous value in the matrix
  /// @param inode ... local node idx
  /// @param ieq   ... index of equation (in case this is a system of PDEs)
  /// @param jnode ... local node idx
  /// @param jeq   ... index of equation (in case this is a system of PDEs)
  /// @param value ... value that should be accumulated
  void add_to_local_mat(const Uint NEQ, const Uint inode, const Uint ieq, const Uint jnode,
                        const Uint jeq, const Real value);

  /// Get one value from the matrix in a tuple together with its global row
  /// and column indices
  /// @param inode  ... local node idx
  /// @param ieq    ... index of equation (in case this is a system of PDEs)
  /// @param jnode  ... local node idx
  /// @param jeq    ... index of equation (in case this is a system of PDEs)
  /// @return tuple [global_dof_id_i, global_dof_id_j, value]
  const std::tuple<Uint, Uint, Real> mat_value(const Uint NEQ, const Uint inode, const Uint ieq,
                                               const Uint jnode, const Uint jeq) const;

  /// Insert value to element right-hand side
  /// This overwrites previous value in the rhs
  /// @param inode ... local node idx
  /// @param ieq   ... index of equation (in case this is a system of PDEs)
  /// @param value ... value that should be inserted
  void insert_in_local_rhs(const Uint NEQ, const Uint inode, const Uint ieq, const Real value);

  /// Add value to element right-hand side
  /// This accumulates to previous value in the matrix
  /// @param inode ... local node idx
  /// @param ieq   ... index of equation (in case this is a system of PDEs)
  /// @param value ... value that should be accumulated
  void add_to_local_rhs(const Uint NEQ, const Uint inode, const Uint ieq, const Real value);

  /// Get one value from the rhs in a tuple together with its global row
  /// index
  /// @param inode ... local node idx
  /// @param ieq   ... index of equation (in case this is a system of PDEs)
  /// @return tuple [global_dof_id_i, value]
  const std::tuple<Uint, Real> rhs_value(const Uint NEQ, const Uint inode, const Uint ieq) const;

  private:
  /// Vector which knows what is the global id of each local DOF
  std::vector<Uint> m_local_to_global_map;

  /// Local element matrix stored row-wise
  std::vector<Real> m_elem_matrix;

  /// Local element RHS
  std::vector<Real> m_elem_rhs;

  /// Number of local degrees of freedom
  Uint m_nb_local_dof;

  /// Number of active cells in mesh
  Uint m_nb_active_cells_in_mesh;

  /// Number of active nodes inmesh
  Uint m_nb_active_nodes_in_mesh;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void LSAssembler::build_dof_sparsity_pattern_continuous(
    const typename result_of::dof_map_t<MeshConfig> &dofs, graph::Graph<Int> &mesh_dofs)
{
  typedef typename result_of::dof_map_t<MeshConfig> cell_dofs_type;

  m_nb_active_cells_in_mesh = dofs.nb_active_cells();
  m_nb_active_nodes_in_mesh = dofs.nb_nodes();

  const Uint ls_size = dofs.nb_nodes();

  if (ls_size != mesh_dofs.nb_vertices())
  {
    std::cerr << "LSAssembler::build_dof_sparsity_pattern_continuous:: graph "
                 "has wrong size\n"
              << " and can not represent dof sparsity. Graph has " << mesh_dofs.nb_vertices()
              << " but DofMap has total " << ls_size << " dofs." << std::endl;
    return;
  }

  mesh_dofs.remove_all_edges();

  Uint max_row_len = 0;

  std::cout << "    LSAssembler: building CONTINUOUS dof sparsity pattern ... " << std::endl;

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::milliseconds elapsed;

  start = std::chrono::high_resolution_clock::now();

  // Minimum and maximum nonzero value in each row to speed up
  // inserting new values
  std::vector<std::pair<Int, Int>> row_min_max(ls_size);
  for (Uint r = 0; r < row_min_max.size(); ++r)
  {
    row_min_max[r].first  = ls_size + 1;
    row_min_max[r].second = 0;
  }

  for (const typename cell_dofs_type::const_dof_range_typed &dof_group :
       dofs.all_active_dof_groups())
  {
    const Uint nb_dof_per_elem = dof_group.begin()->mesh_entity().nb_vert();

    for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_elem; ++i_dof_in_elem)
    {
      for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = dof_group.begin();
           cell_iter != dof_group.end(); ++cell_iter)
      {
        const mesh::MeshEntity solution_elem = cell_iter->mesh_entity();

        for (Uint v = 0; v < solution_elem.nb_vert(); ++v)
        {
          const Int global_row_idx = solution_elem.vertex(v);

          const Int col_idx = solution_elem.vertex(i_dof_in_elem);

          if (col_idx < row_min_max[global_row_idx].first)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].first =
                std::min(row_min_max[global_row_idx].first, col_idx);
          }
          else if (col_idx > row_min_max[global_row_idx].second)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].second =
                std::max(row_min_max[global_row_idx].second, col_idx);
          }
          else // Perform search in row
          {
            mesh_dofs.insert_edge_unique(global_row_idx, col_idx);
            max_row_len = std::max(max_row_len, mesh_dofs.number_adj_vertices(global_row_idx));
          }
        }

      } // Loop over all cells of one group

    } // Loop over local nodes of one element

  } // Lop over all cell groups

  end     = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  const Uint tot_nb_dofs = mesh_dofs.nb_vertices();

  std::cout << "  LSAssembler: Counted nonzero entires in linear system: " << tot_nb_dofs
            << std::endl;
  std::cout << "  LSAssembler: Maximum row length: " << max_row_len << std::endl;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(7);
  std::cout << "LSAssembler: CPU time (inspecting matrix structure) = " << elapsed.count() << " ms"
            << std::endl;

  Real used_storage_size = 0.0;

  for (Uint r = 0; r < mesh_dofs.nb_vertices(); ++r)
  {
    used_storage_size += mesh_dofs.number_adj_vertices(r);
  }
  used_storage_size *= sizeof(Real);

  std::cout << "LSAssembler: used storage size = " << used_storage_size / 1.e6 << " Mb"
            << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void LSAssembler::build_dof_sparsity_pattern_discontinuous(
    const mesh::Tria<MeshConfig> &cell_topology,
    const typename result_of::dof_map_t<MeshConfig> &dofs, graph::Graph<Int> &mesh_dofs)
{
  typedef typename result_of::dof_map_t<MeshConfig> cell_dofs_type;

  m_nb_active_cells_in_mesh = dofs.nb_active_cells();
  m_nb_active_nodes_in_mesh = dofs.nb_nodes();

  const Uint ls_size = dofs.nb_nodes();

  if (ls_size != mesh_dofs.nb_vertices())
  {
    std::cerr << "LSAssembler::build_dof_sparsity_pattern_discontinuous:: "
                 "graph has wrong size\n"
              << " and can not represent dof sparsity. Graph has " << mesh_dofs.nb_vertices()
              << " but DofMap has total " << ls_size << " dofs." << std::endl;
    return;
  }

  mesh_dofs.remove_all_edges();

  Uint max_row_len = 0;

  std::cout << "    LSAssembler: building DISCONTINUOUS dof sparsity pattern ... " << std::endl;

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::milliseconds elapsed;

  mesh::EntityDofRealign sol_facet_permutation_L, sol_facet_permutation_R;

  start = std::chrono::high_resolution_clock::now();

  // Minimum and maximum nonzero value in each row to speed up
  // inserting new values
  std::vector<std::pair<Int, Int>> row_min_max(ls_size);
  for (Uint r = 0; r < row_min_max.size(); ++r)
  {
    row_min_max[r].first  = ls_size + 1;
    row_min_max[r].second = 0;
  }

  // --------------------
  // PROCESSING OF FACETS
  // --------------------

  const Uint nb_facets = cell_topology.active_skeleton_size(MeshConfig::TDIM - 1);

  for (Uint f = 0; f < nb_facets; ++f)
  {
    // std::cout << "f = " << f << std::endl;

    const mesh::TraceIncidences facet_blk =
        cell_topology.active_skeleton_entry(MeshConfig::TDIM - 1, mesh::ActiveIdx(f));

    // Get the permutation sign of the first entity in the block (position
    // 0)
    const mesh::EntityRealignCode permutation_code = facet_blk.permutation(0).get().code();

    // If the entity on position 0 in incidence block needs to be flipped,
    // it is the entity on the rhs of this facet
    const Uint idx_L = (permutation_code.nb_flips() == 0) ? 0 : 1;

    // If the facet block has size 2, we have 2 entities forming the face
    // (i.e K+ and K-). In that case, idx_R is the 'other' index (choosing
    // from 0 and 1) than idx_L, and can be determined as idx_R = (idx_L +
    // 1) % 2 If the facet block has size 1, then idx_R has to be the same
    // as index left: equal to 0
    const Uint idx_R = (facet_blk.size() == 2) ? (idx_L + 1) % 2 : 0;

    const mesh::CellTopologyView<MeshConfig> tcell_L =
        cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_L)));
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_R)));

    const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
    const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

    const mesh::MeshEntity sol_cell_L = dofs.active_cell(active_cell_id_L);
    const mesh::MeshEntity sol_facet_L =
        sol_cell_L.sub_entity(MeshConfig::TDIM - 1, facet_blk.local_id(idx_L));

    const mesh::MeshEntity sol_cell_R = dofs.active_cell(active_cell_id_R);
    const mesh::MeshEntity sol_facet_R =
        sol_cell_R.sub_entity(MeshConfig::TDIM - 1, facet_blk.local_id(idx_R));

    sol_facet_permutation_L.change_type(sol_facet_L.pt_set_id(),
                                        facet_blk.permutation(idx_L).get().code());
    sol_facet_permutation_R.change_type(sol_facet_R.pt_set_id(),
                                        facet_blk.permutation(idx_R).get().code());

    const Uint nb_dof_per_L_facet = sol_facet_L.nb_vert();
    const Uint nb_dof_per_R_facet = sol_facet_R.nb_vert();

    if (idx_L != idx_R)
    {
      // ----------------------------------------------------------------------
      // Entries corresponding to the perturbation of residuals on left
      // facet
      // ----------------------------------------------------------------------
      for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_L_facet; ++i_dof_in_elem)
      {
        const Int col_idx = sol_facet_L.vertex(i_dof_in_elem);

        // Matrix entires corresponding to: dR_L/du_L, where R_L are
        // residual values on Left facet and du_L are solution
        // perturbations on Left facet
        for (Uint v = 0; v < nb_dof_per_L_facet; ++v)
        {
          const Int global_row_idx = sol_facet_L.vertex(v);

          if (col_idx < row_min_max[global_row_idx].first)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].first =
                std::min(row_min_max[global_row_idx].first, col_idx);
          }
          else if (col_idx > row_min_max[global_row_idx].second)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].second =
                std::max(row_min_max[global_row_idx].second, col_idx);
          }
          else // Perform search in row
          {
            mesh_dofs.insert_edge_unique(global_row_idx, col_idx);
            max_row_len = std::max(max_row_len, mesh_dofs.number_adj_vertices(global_row_idx));
          }
        }

        // Matrix entires corresponding to: dR_R/du_L, where R_R are
        // residual values on Right facet and du_L are solution
        // perturbations on Left facet
        for (Uint v = 0; v < nb_dof_per_R_facet; ++v)
        {
          const Int global_row_idx = sol_facet_R.vertex(sol_facet_permutation_R.get().vertex(v));

          if (col_idx < row_min_max[global_row_idx].first)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].first =
                std::min(row_min_max[global_row_idx].first, col_idx);
          }
          else if (col_idx > row_min_max[global_row_idx].second)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].second =
                std::max(row_min_max[global_row_idx].second, col_idx);
          }
          else // Perform search in row
          {
            mesh_dofs.insert_edge_unique(global_row_idx, col_idx);
            max_row_len = std::max(max_row_len, mesh_dofs.number_adj_vertices(global_row_idx));
          }
        }
      } // Loop over DOFs of left facet

      // ----------------------------------------------------------------------
      // Entries corresponding to the perturbation of residuals on right
      // facet
      // ----------------------------------------------------------------------
      for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_R_facet; ++i_dof_in_elem)
      {
        const Int col_idx = sol_facet_R.vertex(sol_facet_permutation_R.get().vertex(i_dof_in_elem));

        // Matrix entires corresponding to: dR_L/du_R, where R_L are
        // residual values on Left facet and du_R are solution
        // perturbations on Right facet
        for (Uint v = 0; v < nb_dof_per_L_facet; ++v)
        {
          const Int global_row_idx = sol_facet_L.vertex(v);

          if (col_idx < row_min_max[global_row_idx].first)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].first =
                std::min(row_min_max[global_row_idx].first, col_idx);
          }
          else if (col_idx > row_min_max[global_row_idx].second)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].second =
                std::max(row_min_max[global_row_idx].second, col_idx);
          }
          else // Perform search in row
          {
            mesh_dofs.insert_edge_unique(global_row_idx, col_idx);
            max_row_len = std::max(max_row_len, mesh_dofs.number_adj_vertices(global_row_idx));
          }
        } // Loop over DOFs of left facet

        // Matrix entires corresponding to: dR_R/du_R, where R_R are
        // residual values on Right facet and du_R are solution
        // perturbations on Right facet
        for (Uint v = 0; v < nb_dof_per_R_facet; ++v)
        {
          const Int global_row_idx = sol_facet_R.vertex(sol_facet_permutation_R.get().vertex(v));

          if (col_idx < row_min_max[global_row_idx].first)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].first =
                std::min(row_min_max[global_row_idx].first, col_idx);
          }
          else if (col_idx > row_min_max[global_row_idx].second)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].second =
                std::max(row_min_max[global_row_idx].second, col_idx);
          }
          else // Perform search in row
          {
            mesh_dofs.insert_edge_unique(global_row_idx, col_idx);
            max_row_len = std::max(max_row_len, mesh_dofs.number_adj_vertices(global_row_idx));
          }
        }
      } // Loop over DOFs of right facet

    } // if ( idx_L != idx_R ) - we are on internal face

  } // Loop over all facets

  // -------------------
  // PROCESSING OF CELLS
  // -------------------

  for (const typename cell_dofs_type::const_dof_range_typed &dof_group :
       dofs.all_active_dof_groups())
  {
    const Uint nb_dof_per_elem = dof_group.begin()->mesh_entity().nb_vert();

    for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_elem; ++i_dof_in_elem)
    {
      for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = dof_group.begin();
           cell_iter != dof_group.end(); ++cell_iter)
      {
        const mesh::MeshEntity solution_elem = cell_iter->mesh_entity();

        for (Uint v = 0; v < solution_elem.nb_vert(); ++v)
        {
          const Int global_row_idx = solution_elem.vertex(v);

          const Int col_idx = solution_elem.vertex(i_dof_in_elem);

          if (col_idx < row_min_max[global_row_idx].first)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].first =
                std::min(row_min_max[global_row_idx].first, col_idx);
          }
          else if (col_idx > row_min_max[global_row_idx].second)
          {
            mesh_dofs.insert_edge(global_row_idx, col_idx);
            row_min_max[global_row_idx].second =
                std::max(row_min_max[global_row_idx].second, col_idx);
          }
          else // Perform search in row
          {
            mesh_dofs.insert_edge_unique(global_row_idx, col_idx);
            max_row_len = std::max(max_row_len, mesh_dofs.number_adj_vertices(global_row_idx));
          }
        }

      } // Loop over all cells of one group

    } // Loop over local nodes of one element

  } // Lop over all cell groups

  end     = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  const Uint tot_nb_dofs = mesh_dofs.nb_vertices();

  std::cout << "  LSAssembler: Counted nonzero entires in linear system: " << tot_nb_dofs
            << std::endl;
  std::cout << "  LSAssembler: Maximum row length: " << max_row_len << std::endl;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(7);
  std::cout << "LSAssembler: CPU time (inspecting matrix structure) = " << elapsed.count() << " ms"
            << std::endl;

  Real used_storage_size = 0.0;

  for (Uint r = 0; r < mesh_dofs.nb_vertices(); ++r)
  {
    used_storage_size += mesh_dofs.number_adj_vertices(r);
  }
  used_storage_size *= sizeof(Real);

  std::cout << "LSAssembler: used storage size = " << used_storage_size / 1.e6 << " Mb"
            << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void LSAssembler::build_matrix_sparsity_pattern_continuous(
    const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &dofs,
    graph::Graph<Int> &ls_dofs)
{
  using cell_dofs_type      = typename result_of::dof_map_t<MeshConfig>;
  using adj_vert_const_iter = graph::Graph<Int>::adj_vertex_const_iterator;

  m_nb_active_cells_in_mesh = dofs.nb_active_cells();
  m_nb_active_nodes_in_mesh = dofs.nb_nodes();

  const Uint ls_size = dofs.nb_nodes() * NEQ;

  if (ls_size != ls_dofs.nb_vertices())
  {
    std::cerr << "LSAssembler::inspect_global_system_structure_continuous:: "
                 "graph has wrong size\n"
              << " and can not represent system dofs. Graph has " << ls_dofs.nb_vertices()
              << " but lin. system has total " << ls_size << " dofs." << std::endl;
    return;
  }

  ls_dofs.remove_all_edges();

  Uint max_row_len = 0;

  std::cout << "    LSAssembler: Counting dofs in (CONTINUOUS FE) linear "
               "system ... "
            << std::endl;
  std::cout << "    LSAssembler: Full matrix size = " << ls_size << " x " << ls_size << std::endl;

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::milliseconds elapsed;

  start = std::chrono::high_resolution_clock::now();

  // Minimum and maximum nonzero value in each row to speed up
  // inserting new values
  std::vector<std::pair<Int, Int>> row_min_max(ls_size);
  for (Uint r = 0; r < row_min_max.size(); ++r)
  {
    row_min_max[r].first  = ls_size + 1;
    row_min_max[r].second = 0;
  }

  std::vector<Uint> old_num_connected_neighbours(ls_size);
  old_num_connected_neighbours.assign(ls_size, 5);

  for (const typename cell_dofs_type::const_dof_range_typed &dof_group :
       dofs.all_active_dof_groups())
  {
    const Uint nb_dof_per_elem = dof_group.begin()->mesh_entity().nb_vert();

    for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_elem; ++i_dof_in_elem)
    {
      for (Uint comp_u = 0; comp_u < NEQ; ++comp_u)
      {
        for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = dof_group.begin();
             cell_iter != dof_group.end(); ++cell_iter)
        {
          const mesh::MeshEntity solution_elem = cell_iter->mesh_entity();

          for (Uint v = 0; v < solution_elem.nb_vert(); ++v)
          {
            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              // For each component
              // When v == i_dof_in_elem, we are computing
              // derivatives of residuals corresponding to node v
              // with respect to components of u associated to v
              // (diagonal blocks) in the system matrix
              const Int global_row_idx = NEQ * solution_elem.vertex(v) + comp;

              const Int col_idx = NEQ * solution_elem.vertex(i_dof_in_elem) + comp_u;

              if (col_idx < row_min_max[global_row_idx].first)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].first = col_idx;
              }
              else if (col_idx > row_min_max[global_row_idx].second)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].second = col_idx;
              }
              else // Perform search in row
              {
                // Old version
                /*
                ls_dofs.insert_edge_unique(global_row_idx,
                col_idx); max_row_len = std::max(max_row_len,
                ls_dofs.number_adj_vertices(global_row_idx));
                */

                // Insert edge to the end of the graph row
                ls_dofs.insert_edge(global_row_idx, col_idx);

                // If too many insertions have been performed,
                // sort the row and remove duplicate elements

                if (ls_dofs.number_adj_vertices(global_row_idx) >
                    2 * old_num_connected_neighbours[global_row_idx])
                {
                  ls_dofs.remove_duplicate_neighbours(global_row_idx);
                  old_num_connected_neighbours[global_row_idx] =
                      ls_dofs.number_adj_vertices(global_row_idx);

                  std::pair<adj_vert_const_iter, adj_vert_const_iter> adj_verts =
                      ls_dofs.adjacent_vertices(global_row_idx);
                  row_min_max[global_row_idx].first  = (*adj_verts.first);
                  row_min_max[global_row_idx].second = (*(--adj_verts.second));
                }
                else
                {
                  row_min_max[global_row_idx].first =
                      std::min(row_min_max[global_row_idx].first, col_idx);
                  row_min_max[global_row_idx].second =
                      std::max(row_min_max[global_row_idx].second, col_idx);
                }
              }
            } // Loop over components
          }

        } // Loop over all cells of one group

      } // Loop over equation components

    } // Loop over local nodes of one element

  } // Lop over all cell groups

  // Final sweep over the graph to remove any remaining duplicate edges
  max_row_len = 0;

  for (Int row_idx = 0; row_idx < ls_dofs.nb_vertices(); ++row_idx)
  {
    ls_dofs.remove_duplicate_neighbours(row_idx);
    std::pair<adj_vert_const_iter, adj_vert_const_iter> adj_verts =
        ls_dofs.adjacent_vertices(row_idx);
    row_min_max[row_idx].first  = (*adj_verts.first);
    row_min_max[row_idx].second = (*(--adj_verts.second));

    max_row_len = std::max(max_row_len, ls_dofs.number_adj_vertices(row_idx));
  }

  end     = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  const Uint tot_nb_dofs = ls_dofs.nb_vertices();

  std::cout << "  LSAssembler: Counted nonzero entires in CONTINUOUS linear "
               "system: "
            << tot_nb_dofs << std::endl;
  std::cout << "  LSAssembler: Maximum row length: " << max_row_len << std::endl;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(7);
  std::cout << "LSAssembler: CPU time (build_matrix_sparsity_pattern_continuous) = "
            << elapsed.count() << " ms" << std::endl;

  Real used_storage_size = 0.0;

  for (Uint r = 0; r < ls_dofs.nb_vertices(); ++r)
  {
    used_storage_size += ls_dofs.number_adj_vertices(r);
  }
  used_storage_size *= sizeof(Real);

  std::cout << "LSAssembler: used storage size = " << used_storage_size / 1.e6 << " Mb"
            << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void LSAssembler::build_matrix_sparsity_pattern_discontinuous(
    const Uint NEQ, const mesh::Tria<MeshConfig> &cell_topology,
    const typename result_of::dof_map_t<MeshConfig> &dofs, graph::Graph<Int> &ls_dofs)
{
  using cell_dofs_type      = typename result_of::dof_map_t<MeshConfig>;
  using adj_vert_const_iter = graph::Graph<Int>::adj_vertex_const_iterator;

  m_nb_active_cells_in_mesh = dofs.nb_active_cells();
  m_nb_active_nodes_in_mesh = dofs.nb_nodes();

  const Uint ls_size = dofs.nb_nodes() * NEQ;

  if (ls_size != ls_dofs.nb_vertices())
  {
    std::cerr << "LSAssembler::inspect_global_system_structure_discontinuous:: "
                 "graph has wrong size\n"
              << " and can not represent system dofs. Graph has " << ls_dofs.nb_vertices()
              << " but lin. system has total " << ls_size << " dofs." << std::endl;
    return;
  }

  ls_dofs.remove_all_edges();

  Uint max_row_len = 0;

  std::cout << "    LSAssembler: Counting dofs in (DISCONTINUOUS FE) linear "
               "system ... "
            << std::endl;
  std::cout << "    LSAssembler: Full matrix size = " << ls_size << " x " << ls_size << std::endl;

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::milliseconds elapsed;

  mesh::EntityDofRealign sol_facet_permutation_L, sol_facet_permutation_R;

  start = std::chrono::high_resolution_clock::now();

  // Minimum and maximum nonzero value in each row to speed up
  // inserting new values
  std::vector<std::pair<Int, Int>> row_min_max(ls_size);
  for (Uint r = 0; r < row_min_max.size(); ++r)
  {
    row_min_max[r].first  = ls_size + 1;
    row_min_max[r].second = 0;
  }

  std::vector<Uint> old_num_connected_neighbours(ls_size);
  old_num_connected_neighbours.assign(ls_size, 5);

  // --------------------
  // PROCESSING OF FACETS
  // --------------------

  // Loop over all facets and positions in the system matrix where the facets
  // will contribute to system Jacobian

  const Uint nb_facets = cell_topology.active_skeleton_size(MeshConfig::TDIM - 1);

  for (Uint f = 0; f < nb_facets; ++f)
  {
    // std::cout << "f = " << f << std::endl;

    const mesh::TraceIncidences facet_blk =
        cell_topology.active_skeleton_entry(MeshConfig::TDIM - 1, mesh::ActiveIdx(f));

    // Get the permutation sign of the first entity in the block (position
    // 0)
    const mesh::EntityRealignCode permutation_code = facet_blk.permutation(0).get().code();

    // If the entity on position 0 in incidence block needs to be flipped,
    // it is the entity on the rhs of this facet
    const Uint idx_L = (permutation_code.nb_flips() == 0) ? 0 : 1;

    // If the facet block has size 2, we have 2 entities forming the face
    // (i.e K+ and K-). In that case, idx_R is the 'other' index (choosing
    // from 0 and 1) than idx_L, and can be determined as idx_R = (idx_L +
    // 1) % 2 If the facet block has size 1, then idx_R has to be the same
    // as index left: equal to 0
    const Uint idx_R = (facet_blk.size() == 2) ? (idx_L + 1) % 2 : 0;

    // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS OF
    // CELLS, NOT THEIR ACTIVE INDICES !!!
    const mesh::CellTopologyView<MeshConfig> tcell_L =
        cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_L)));
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        cell_topology.cell(mesh::FlatIdx(facet_blk.cell_id(idx_R)));

    const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
    const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();

    const mesh::MeshEntity sol_cell_L = dofs.active_cell(active_cell_id_L);
    const mesh::MeshEntity sol_facet_L =
        sol_cell_L.sub_entity(MeshConfig::TDIM - 1, facet_blk.local_id(idx_L));

    const mesh::MeshEntity sol_cell_R = dofs.active_cell(active_cell_id_R);
    const mesh::MeshEntity sol_facet_R =
        sol_cell_R.sub_entity(MeshConfig::TDIM - 1, facet_blk.local_id(idx_R));

    /*
    const mesh::MeshEntity sol_cell_L =
    dofs.active_cell(facet_blk.cell_id(idx_L)); const mesh::MeshEntity
    sol_facet_L = sol_cell_L.sub_entity(MeshConfig::TDIM - 1,
    facet_blk.local_id(idx_L));

    const mesh::MeshEntity sol_cell_R =
    dofs.active_cell(facet_blk.cell_id(idx_R)); const mesh::MeshEntity
    sol_facet_R = sol_cell_R.sub_entity(MeshConfig::TDIM - 1,
    facet_blk.local_id(idx_R));
    */

    sol_facet_permutation_L.change_type(sol_facet_L.pt_set_id(),
                                        facet_blk.permutation(idx_L).get().code());
    sol_facet_permutation_R.change_type(sol_facet_R.pt_set_id(),
                                        facet_blk.permutation(idx_R).get().code());

    const Uint nb_dof_per_L_facet = sol_facet_L.nb_vert();
    const Uint nb_dof_per_R_facet = sol_facet_R.nb_vert();

    if (idx_L != idx_R)
    {
      // ----------------------------------------------------------------------
      // Entries corresponding to the perturbation of residuals on left
      // facet
      // ----------------------------------------------------------------------
      for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_L_facet; ++i_dof_in_elem)
      {
        for (Uint comp_u = 0; comp_u < NEQ; ++comp_u)
        {
          const Int col_idx = NEQ * sol_facet_L.vertex(i_dof_in_elem) + comp_u;

          // Matrix entires corresponding to: dR_L/du_L, where R_L are
          // residual values on Left facet and du_L are solution
          // perturbations on Left facet
          for (Uint v = 0; v < sol_facet_L.nb_vert(); ++v)
          {
            const Uint left_vertex_id = sol_facet_L.vertex(v);

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              const Int global_row_idx = NEQ * left_vertex_id + comp;

              if (col_idx < row_min_max[global_row_idx].first)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].first =
                    std::min(row_min_max[global_row_idx].first, col_idx);
              }
              else if (col_idx > row_min_max[global_row_idx].second)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].second =
                    std::max(row_min_max[global_row_idx].second, col_idx);
              }
              else // Perform search in row
              {
// Old version
#if 1
                ls_dofs.insert_edge_unique(global_row_idx, col_idx);
                max_row_len = std::max(max_row_len, ls_dofs.number_adj_vertices(global_row_idx));
#else

                // Insert edge to the end of the graph row
                ls_dofs.insert_edge(global_row_idx, col_idx);

                // If too many insertions have been performed,
                // sort the row and remove duplicate elements

                if (ls_dofs.number_adj_vertices(global_row_idx) >
                    2 * old_num_connected_neighbours[global_row_idx])
                {
                  ls_dofs.remove_duplicate_neighbours(global_row_idx);
                  old_num_connected_neighbours[global_row_idx] =
                      ls_dofs.number_adj_vertices(global_row_idx);

                  std::pair<adj_vert_const_iter, adj_vert_const_iter> adj_verts =
                      ls_dofs.adjacent_vertices(global_row_idx);
                  row_min_max[global_row_idx].first  = (*adj_verts.first);
                  row_min_max[global_row_idx].second = (*(--adj_verts.second));
                }
                else
                {
                  row_min_max[global_row_idx].first =
                      std::min(row_min_max[global_row_idx].first, col_idx);
                  row_min_max[global_row_idx].second =
                      std::max(row_min_max[global_row_idx].second, col_idx);
                }
#endif
              }
            } // Loop over components
          }

          // Matrix entires corresponding to: dR_R/du_L, where R_R are
          // residual values on Right facet and du_L are solution
          // perturbations on Left facet
          for (Uint v = 0; v < sol_facet_R.nb_vert(); ++v)
          {
            const Uint right_vertex_id =
                sol_facet_R.vertex(sol_facet_permutation_R.get().vertex(v));

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              const Int global_row_idx = NEQ * right_vertex_id + comp;

              if (col_idx < row_min_max[global_row_idx].first)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].first =
                    std::min(row_min_max[global_row_idx].first, col_idx);
              }
              else if (col_idx > row_min_max[global_row_idx].second)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].second =
                    std::max(row_min_max[global_row_idx].second, col_idx);
              }
              else // Perform search in row
              {
// Old version
#if 1
                ls_dofs.insert_edge_unique(global_row_idx, col_idx);
                max_row_len = std::max(max_row_len, ls_dofs.number_adj_vertices(global_row_idx));
#else

                // Insert edge to the end of the graph row
                ls_dofs.insert_edge(global_row_idx, col_idx);

                // If too many insertions have been performed,
                // sort the row and remove duplicate elements

                if (ls_dofs.number_adj_vertices(global_row_idx) >
                    2 * old_num_connected_neighbours[global_row_idx])
                {
                  ls_dofs.remove_duplicate_neighbours(global_row_idx);
                  old_num_connected_neighbours[global_row_idx] =
                      ls_dofs.number_adj_vertices(global_row_idx);

                  std::pair<adj_vert_const_iter, adj_vert_const_iter> adj_verts =
                      ls_dofs.adjacent_vertices(global_row_idx);
                  row_min_max[global_row_idx].first  = (*adj_verts.first);
                  row_min_max[global_row_idx].second = (*(--adj_verts.second));
                }
                else
                {
                  row_min_max[global_row_idx].first =
                      std::min(row_min_max[global_row_idx].first, col_idx);
                  row_min_max[global_row_idx].second =
                      std::max(row_min_max[global_row_idx].second, col_idx);
                }
#endif
              } // Perform search in row
            }
          }
        } // Loop over comp_u
      }   // Loop over DOFs of left facet

      // ----------------------------------------------------------------------
      // Entries corresponding to the perturbation of residuals on right
      // facet
      // ----------------------------------------------------------------------
      for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_R_facet; ++i_dof_in_elem)
      {
        for (Uint comp_u = 0; comp_u < NEQ; ++comp_u)
        {
          const Int col_idx =
              NEQ * sol_facet_R.vertex(sol_facet_permutation_R.get().vertex(i_dof_in_elem)) +
              comp_u;

          // Matrix entires corresponding to: dR_L/du_R, where R_L are
          // residual values on Left facet and du_R are solution
          // perturbations on Right facet
          for (Uint v = 0; v < sol_facet_L.nb_vert(); ++v)
          {
            const Uint left_vertex_id = sol_facet_L.vertex(v);

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              const Int global_row_idx = NEQ * left_vertex_id + comp;

              if (col_idx < row_min_max[global_row_idx].first)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].first =
                    std::min(row_min_max[global_row_idx].first, col_idx);
              }
              else if (col_idx > row_min_max[global_row_idx].second)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].second =
                    std::max(row_min_max[global_row_idx].second, col_idx);
              }
              else // Perform search in row
              {
// Old version
#if 1
                ls_dofs.insert_edge_unique(global_row_idx, col_idx);
                max_row_len = std::max(max_row_len, ls_dofs.number_adj_vertices(global_row_idx));
#else

                // Insert edge to the end of the graph row
                ls_dofs.insert_edge(global_row_idx, col_idx);

                // If too many insertions have been performed,
                // sort the row and remove duplicate elements

                if (ls_dofs.number_adj_vertices(global_row_idx) >
                    2 * old_num_connected_neighbours[global_row_idx])
                {
                  ls_dofs.remove_duplicate_neighbours(global_row_idx);
                  old_num_connected_neighbours[global_row_idx] =
                      ls_dofs.number_adj_vertices(global_row_idx);

                  std::pair<adj_vert_const_iter, adj_vert_const_iter> adj_verts =
                      ls_dofs.adjacent_vertices(global_row_idx);
                  row_min_max[global_row_idx].first  = (*adj_verts.first);
                  row_min_max[global_row_idx].second = (*(--adj_verts.second));
                }
                else
                {
                  row_min_max[global_row_idx].first =
                      std::min(row_min_max[global_row_idx].first, col_idx);
                  row_min_max[global_row_idx].second =
                      std::max(row_min_max[global_row_idx].second, col_idx);
                }
#endif
              }
            }
          } // Loop over DOFs of left facet

          // Matrix entires corresponding to: dR_R/du_R, where R_R are
          // residual values on Right facet and du_R are solution
          // perturbations on Right facet
          for (Uint v = 0; v < sol_facet_R.nb_vert(); ++v)
          {
            const Uint right_vertex_id =
                sol_facet_R.vertex(sol_facet_permutation_R.get().vertex(v));

            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              const Int global_row_idx = NEQ * right_vertex_id + comp;

              if (col_idx < row_min_max[global_row_idx].first)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].first =
                    std::min(row_min_max[global_row_idx].first, col_idx);
              }
              else if (col_idx > row_min_max[global_row_idx].second)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].second =
                    std::max(row_min_max[global_row_idx].second, col_idx);
              }
              else // Perform search in row
              {
// Old version
#if 1
                ls_dofs.insert_edge_unique(global_row_idx, col_idx);
                max_row_len = std::max(max_row_len, ls_dofs.number_adj_vertices(global_row_idx));
#else

                // Insert edge to the end of the graph row
                ls_dofs.insert_edge(global_row_idx, col_idx);

                // If too many insertions have been performed,
                // sort the row and remove duplicate elements

                if (ls_dofs.number_adj_vertices(global_row_idx) >
                    2 * old_num_connected_neighbours[global_row_idx])
                {
                  ls_dofs.remove_duplicate_neighbours(global_row_idx);
                  old_num_connected_neighbours[global_row_idx] =
                      ls_dofs.number_adj_vertices(global_row_idx);

                  std::pair<adj_vert_const_iter, adj_vert_const_iter> adj_verts =
                      ls_dofs.adjacent_vertices(global_row_idx);
                  row_min_max[global_row_idx].first  = (*adj_verts.first);
                  row_min_max[global_row_idx].second = (*(--adj_verts.second));
                }
                else
                {
                  row_min_max[global_row_idx].first =
                      std::min(row_min_max[global_row_idx].first, col_idx);
                  row_min_max[global_row_idx].second =
                      std::max(row_min_max[global_row_idx].second, col_idx);
                }
#endif
              }
            }
          }
        } // Loop over comp_u
      }   // Loop over DOFs of right facet

    } // if ( idx_L != idx_R ) - we are on internal face

  } // Loop over all facets

#if 0
    // Final sweep over the graph to remove any remaining duplicate edges
    max_row_len = 0;

    for (Int row_idx = 0; row_idx < ls_dofs.nb_vertices(); ++row_idx)
    {
      ls_dofs.remove_duplicate_neighbours(row_idx);
      std::pair<adj_vert_const_iter, adj_vert_const_iter> adj_verts =
          ls_dofs.adjacent_vertices(row_idx);
      row_min_max[row_idx].first = (*adj_verts.first);
      row_min_max[row_idx].second = (*(--adj_verts.second));

      max_row_len = std::max(max_row_len, ls_dofs.number_adj_vertices(row_idx));
    }
#endif

  // -------------------
  // PROCESSING OF CELLS
  // -------------------

  for (const typename cell_dofs_type::const_dof_range_typed &dof_group :
       dofs.all_active_dof_groups())
  {
    const Uint nb_dof_per_elem = dof_group.begin()->mesh_entity().nb_vert();

    for (Uint i_dof_in_elem = 0; i_dof_in_elem < nb_dof_per_elem; ++i_dof_in_elem)
    {
      for (Uint comp_u = 0; comp_u < NEQ; ++comp_u)
      {
        for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = dof_group.begin();
             cell_iter != dof_group.end(); ++cell_iter)
        {
          const mesh::MeshEntity solution_elem = cell_iter->mesh_entity();

          for (Uint v = 0; v < solution_elem.nb_vert(); ++v)
          {
            for (Uint comp = 0; comp < NEQ; ++comp)
            {
              // For each component
              // When v == i_dof_in_elem, we are computing
              // derivatives of residuals corresponding to node v
              // with respect to components of u associated to v
              // (diagonal blocks) in the system matrix
              const Int global_row_idx = NEQ * solution_elem.vertex(v) + comp;

              const Int col_idx = NEQ * solution_elem.vertex(i_dof_in_elem) + comp_u;

              if (col_idx < row_min_max[global_row_idx].first)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].first =
                    std::min(row_min_max[global_row_idx].first, col_idx);
              }
              else if (col_idx > row_min_max[global_row_idx].second)
              {
                ls_dofs.insert_edge(global_row_idx, col_idx);
                row_min_max[global_row_idx].second =
                    std::max(row_min_max[global_row_idx].second, col_idx);
              }
              else // Perform search in row
              {
// Old version
#if 0
                ls_dofs.insert_edge_unique(global_row_idx, col_idx);
                max_row_len = std::max(max_row_len, ls_dofs.number_adj_vertices(global_row_idx));
#else

                // Insert edge to the end of the graph row
                ls_dofs.insert_edge(global_row_idx, col_idx);

                // If too many insertions have been performed,
                // sort the row and remove duplicate elements

                if (ls_dofs.number_adj_vertices(global_row_idx) >
                    2 * old_num_connected_neighbours[global_row_idx])
                {
                  ls_dofs.remove_duplicate_neighbours(global_row_idx);
                  old_num_connected_neighbours[global_row_idx] =
                      ls_dofs.number_adj_vertices(global_row_idx);

                  std::pair<adj_vert_const_iter, adj_vert_const_iter> adj_verts =
                      ls_dofs.adjacent_vertices(global_row_idx);
                  row_min_max[global_row_idx].first  = (*adj_verts.first);
                  row_min_max[global_row_idx].second = (*(--adj_verts.second));
                }
                else
                {
                  row_min_max[global_row_idx].first =
                      std::min(row_min_max[global_row_idx].first, col_idx);
                  row_min_max[global_row_idx].second =
                      std::max(row_min_max[global_row_idx].second, col_idx);
                }
#endif
              }
            } // Loop over components
          }

        } // Loop over all cells of one group

      } // Loop over equation components

    } // Loop over local nodes of one element

  } // Lop over all cell groups

  // Final sweep over the graph to remove any remaining duplicate edges
  max_row_len = 0;

  for (Int row_idx = 0; row_idx < ls_dofs.nb_vertices(); ++row_idx)
  {
    ls_dofs.remove_duplicate_neighbours(row_idx);
    std::pair<adj_vert_const_iter, adj_vert_const_iter> adj_verts =
        ls_dofs.adjacent_vertices(row_idx);
    row_min_max[row_idx].first  = (*adj_verts.first);
    row_min_max[row_idx].second = (*(--adj_verts.second));

    max_row_len = std::max(max_row_len, ls_dofs.number_adj_vertices(row_idx));
  }

  end     = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  const Uint tot_nb_dofs = ls_dofs.nb_vertices();

  std::cout << "  LSAssembler: Counted nonzero entires in DISCONTINUOUS linear "
               "system: "
            << tot_nb_dofs << std::endl;
  std::cout << "  LSAssembler: Maximum row length: " << max_row_len << std::endl;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(7);
  std::cout << "LSAssembler: CPU time "
               "(build_matrix_sparsity_pattern_discontinuous) = "
            << elapsed.count() << " ms" << std::endl;

  Real used_storage_size = 0.0;

  for (Uint r = 0; r < ls_dofs.nb_vertices(); ++r)
  {
    used_storage_size += ls_dofs.number_adj_vertices(r);
  }
  used_storage_size *= sizeof(Real);

  std::cout << "LSAssembler: used storage size = " << used_storage_size / 1.e6 << " Mb"
            << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif
