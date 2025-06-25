#include "linear_system/LSAssembler.hpp"

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------

LSAssembler::LSAssembler()
{
  m_local_to_global_map.resize(0);
  m_elem_matrix.resize(0);
  m_elem_rhs.resize(0.0);
  m_nb_local_dof            = 0u;
  m_nb_active_cells_in_mesh = 0u;
  m_nb_active_nodes_in_mesh = 0u;
}

// ----------------------------------------------------------------------------

LSAssembler::~LSAssembler()
{
}

// ----------------------------------------------------------------------------

void LSAssembler::build_matrix_sparsity_pattern_block_discontinuous(
    const common::BlockArray<Uint, Uint> mesh_dual_graph_crs, const std::vector<Uint> &block_sizes,
    const std::vector<Uint> &cell_reordering,
    math::BlockMatrixSparsityPattern<Int> &coarse_matrix_sparsity, const std::string &name)
{
  // Generate the inverse ordering map: for each new block id,
  // the inverse map stores on position [i] the index of the cell before
  // reordering
  std::vector<Uint> inverse_cell_reordering;
  inverse_cell_reordering.resize(cell_reordering.size());

  for (Uint i = 0; i < cell_reordering.size(); ++i)
  {
    inverse_cell_reordering[cell_reordering[i]] = i;
  }

  std::unique_ptr<std::vector<common::Range1D<Int>>> block_row_dof_spans(
      new std::vector<common::Range1D<Int>>());
  block_row_dof_spans->resize(block_sizes.size());

  // Store the width of each block row
  std::unique_ptr<std::vector<Uint>> block_row_span_width(new std::vector<Uint>());
  block_row_span_width->resize(block_sizes.size());

  // The following two loops fill the vector 'block_row_positions' with
  // ranges spanning the DOF ids corresponding to each block row

  // At the same time, the loop counts the number of all blocks in all block
  // rows

  // Pass 1
  std::unique_ptr<std::vector<Uint>> block_row_offsets(new std::vector<Uint>());

  // Assign the block sizes
  block_row_offsets->resize(block_sizes.size() + 1);
  block_row_offsets->assign(block_sizes.size() + 1, 0u);

  Uint tot_num_blocks    = 0;
  Int num_dofs_processed = 0;

  for (Uint block_idx = 0; block_idx < inverse_cell_reordering.size(); ++block_idx)
  {
    const Uint reordered_block_id = inverse_cell_reordering[block_idx];

    const Uint block_size = block_sizes[reordered_block_id];

    (*block_row_dof_spans)[block_idx] =
        common::Range1D<Int>(num_dofs_processed, num_dofs_processed + block_size - 1);
    (*block_row_span_width)[block_idx] = block_size;
    num_dofs_processed += block_size;

    const common::ArrayView<const Uint, _1D, Uint> neighbours =
        mesh_dual_graph_crs.const_block(reordered_block_id);
    (*block_row_offsets)[block_idx + 1] = neighbours.size() + 1;

    tot_num_blocks += neighbours.size() + 1;
  }

  // Convert sizes to row offsets
  for (Uint i = 1; i < block_row_offsets->size(); ++i)
  {
    (*block_row_offsets)[i] += (*block_row_offsets)[i - 1];
  }

  // Pass 2 - fill information about dof ranges and neighbours for each cell

  std::unique_ptr<std::vector<common::Range1D<Int>>> block_rows(
      new std::vector<common::Range1D<Int>>());
  block_rows->resize(tot_num_blocks);

  std::unique_ptr<std::vector<Int>> block_row_ids(new std::vector<Int>());
  block_row_ids->resize(tot_num_blocks);

  for (Uint ac = 0; ac < block_sizes.size(); ++ac)
  {
    const Uint block_idx = cell_reordering[ac];

    const common::ArrayView<const Uint, _1D, Uint> neighbours = mesh_dual_graph_crs.const_block(ac);
    Uint insert_pos                                           = (*block_row_offsets)[block_idx];

    // One entry of the block row must be the cell itself
    (*block_rows)[insert_pos]    = (*block_row_dof_spans)[block_idx];
    (*block_row_ids)[insert_pos] = block_idx;
    insert_pos++;

    // The remaining entries in the block row are entries corresponding to
    // neighbours
    for (Uint n = 0; n < neighbours.size(); ++n)
    {
      const Uint block_neighbour_id = cell_reordering[neighbours[n]];
      (*block_rows)[insert_pos]     = (*block_row_dof_spans)[block_neighbour_id];
      (*block_row_ids)[insert_pos]  = block_neighbour_id;
      insert_pos++;
    }
  }

  // -------------------------

  /*
  std::cout << "Row ranges (unordered):" << std::endl;
  for (Uint ac = 0; ac < cell_source_dofs.nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity active_cell = cell_source_dofs.active_cell(ac);
    std::cout << (*block_row_dof_spans)[active_cell.idx()] << " ";
  }
  std::cout << std::endl;


  std::cout << "Row ranges (ordered):" << std::endl;
  for (Uint ac = 0; ac < cell_source_dofs.nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity active_cell = cell_source_dofs.active_cell(ac);
    const Uint block_idx = cell_reordering[active_cell.idx()];
    std::cout << (*block_row_dof_spans)[block_idx] << " ";
  }
  std::cout << std::endl;
  */

  // -------------------------

  // Store linearly all block-rows (contatenated row-wise)
  std::unique_ptr<common::BlockMultiArray<common::Range1D<Int>, Int>> nonzero_blocks(
      new common::BlockMultiArray<common::Range1D<Int>, Int>());

  std::tuple<std::unique_ptr<std::vector<common::Range1D<Int>>>, std::unique_ptr<std::vector<Int>>>
      tmp;
  std::get<0>(tmp) = std::move(block_rows);
  std::get<1>(tmp) = std::move(block_row_ids);

  nonzero_blocks->build_from_offsets(std::move(tmp), std::move(block_row_offsets));
  coarse_matrix_sparsity.build_sparsity(std::move(nonzero_blocks), std::move(block_row_dof_spans));

  // coarse_matrix_sparsity.print_vtu("block_sparsity.vtu");
  coarse_matrix_sparsity.print_vtu(name);
}

// ----------------------------------------------------------------------------

void LSAssembler::set_local_to_global_dof_map(const Uint NEQ, const mesh::MeshEntity &entity)
{
  m_local_to_global_map.resize(entity.nb_vert() * NEQ);
  for (Uint v = 0; v < entity.nb_vert(); ++v)
  {
    for (Uint eq = 0; eq < NEQ; ++eq)
    {
      m_local_to_global_map[v * NEQ + eq] = entity.vertex(v) * NEQ + eq;
    }
  }

  m_nb_local_dof = entity.nb_vert() * NEQ;
  m_elem_matrix.resize(m_nb_local_dof * m_nb_local_dof);
  m_elem_rhs.resize(m_nb_local_dof);
}

// ----------------------------------------------------------------------------

void LSAssembler::set_local_mat_entries(const Real value)
{
  m_elem_matrix.assign(m_elem_matrix.size(), value);
}

// ----------------------------------------------------------------------------

void LSAssembler::set_local_rhs_entries(const Real value)
{
  m_elem_rhs.assign(m_elem_rhs.size(), value);
}

// ----------------------------------------------------------------------------

void LSAssembler::insert_in_local_mat(const Uint NEQ, const Uint inode, const Uint ieq,
                                      const Uint jnode, const Uint jeq, const Real value)
{
  const Uint row_idx                                = inode * NEQ + ieq;
  const Uint col_idx                                = jnode * NEQ + jeq;
  m_elem_matrix[row_idx * m_nb_local_dof + col_idx] = value;
}

// ----------------------------------------------------------------------------

void LSAssembler::add_to_local_mat(const Uint NEQ, const Uint inode, const Uint ieq,
                                   const Uint jnode, const Uint jeq, const Real value)
{
  const Uint row_idx = inode * NEQ + ieq;
  const Uint col_idx = jnode * NEQ + jeq;
  m_elem_matrix[row_idx * m_nb_local_dof + col_idx] += value;
}

// ----------------------------------------------------------------------------

const std::tuple<Uint, Uint, Real> LSAssembler::mat_value(const Uint NEQ, const Uint inode,
                                                          const Uint ieq, const Uint jnode,
                                                          const Uint jeq) const
{
  const Uint row_idx = inode * NEQ + ieq;
  const Uint col_idx = jnode * NEQ + jeq;
  return std::tuple<Uint, Uint, Real>(row_idx, col_idx,
                                      m_elem_matrix[row_idx * m_nb_local_dof + col_idx]);
}

// ----------------------------------------------------------------------------

void LSAssembler::insert_in_local_rhs(const Uint NEQ, const Uint inode, const Uint ieq,
                                      const Real value)
{
  const Uint row_idx  = inode * NEQ + ieq;
  m_elem_rhs[row_idx] = value;
}

// ----------------------------------------------------------------------------

void LSAssembler::add_to_local_rhs(const Uint NEQ, const Uint inode, const Uint ieq,
                                   const Real value)
{
  const Uint row_idx = inode * NEQ + ieq;
  m_elem_rhs[row_idx] += value;
}

// ----------------------------------------------------------------------------

const std::tuple<Uint, Real> LSAssembler::rhs_value(const Uint NEQ, const Uint inode,
                                                    const Uint ieq) const
{
  const Uint row_idx = inode * NEQ + ieq;
  return std::tuple<Uint, Real>(m_local_to_global_map[row_idx], m_elem_rhs[row_idx]);
}

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit
