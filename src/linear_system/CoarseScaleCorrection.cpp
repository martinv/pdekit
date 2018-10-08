#include "linear_system/CoarseScaleCorrection.hpp"

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------

CoarseScaleCorrection::CoarseScaleCorrection()
    : m_local_restriction_ops(new ls::LocalTransferOps<Real>())
{
}

// ----------------------------------------------------------------------------

CoarseScaleCorrection::~CoarseScaleCorrection()
{
}

// ----------------------------------------------------------------------------

void CoarseScaleCorrection::build_block_matrix_sparsity_pattern(
    const common::BlockArray<Uint, Uint> &mesh_dual_graph_crs,
    const std::vector<Uint> &coarse_block_sizes, const std::vector<Uint> &fine_block_sizes,
    const std::vector<Uint> &cell_reordering,
    std::unique_ptr<ls::LocalTransferOps<Real>> &&restriction_ops,
    std::unique_ptr<ls::LocalTransferOps<Real>> &&prolongation_ops)
{
  m_inverse_cell_reordering.resize(cell_reordering.size());

  for (Uint i = 0; i < cell_reordering.size(); ++i)
  {
    m_inverse_cell_reordering[cell_reordering[i]] = i;
  }

  LSAssembler assembler;

  assembler.build_matrix_sparsity_pattern_block_discontinuous(
      mesh_dual_graph_crs, coarse_block_sizes, cell_reordering, m_coarse_matrix_sparsity,
      "block_sparsity_coarse.vtu");

  assembler.build_matrix_sparsity_pattern_block_discontinuous(
      mesh_dual_graph_crs, fine_block_sizes, cell_reordering, m_fine_matrix_sparsity,
      "block_sparsity_fine.vtu");

  m_local_restriction_ops  = std::move(restriction_ops);
  m_local_prolongation_ops = std::move(prolongation_ops);
}

// ----------------------------------------------------------------------------

void CoarseScaleCorrection::restrict_vec(const TpetraMultiVector<Real> &rhs_fine,
                                         TpetraMultiVector<Real> &rhs_coarse) const
{
  std::vector<Real> workspace_fine, workspace_coarse;

  for (Uint br = 0; br < m_coarse_matrix_sparsity.nb_lines(); ++br)
  {
    const common::Range1D<Int> coarse_rows = m_coarse_matrix_sparsity.sparse_line_limits(br);
    const common::Range1D<Int> fine_rows   = m_fine_matrix_sparsity.sparse_line_limits(br);

    workspace_coarse.resize(coarse_rows.size());
    workspace_fine.resize(fine_rows.size());

    for (Uint r = fine_rows.lbound(), i = 0; r <= fine_rows.ubound(); ++r, ++i)
    {
      workspace_fine[i] = rhs_fine.value(r, 0);
    }

    const math::DenseConstVecView<Real> loc_rhs_fine(workspace_fine.data(), workspace_fine.size());
    math::DenseVecView<Real> loc_rhs_coarse(workspace_coarse.data(), workspace_coarse.size());

    const math::DenseConstMatView<Real> restrict_op = m_local_restriction_ops->const_op(br);

    loc_rhs_coarse = restrict_op * loc_rhs_fine;

    for (Uint r = coarse_rows.lbound(), i = 0; r <= coarse_rows.ubound(); ++r, ++i)
    {
      const Real val = loc_rhs_coarse[i];
      rhs_coarse.insert_value(r, val, 0);
    }
  } // Loop over cells
}

// ----------------------------------------------------------------------------

void CoarseScaleCorrection::prolongate_vec(const TpetraMultiVector<Real> &u_coarse,
                                           TpetraMultiVector<Real> &u_fine) const
{
  std::vector<Real> workspace_fine;
  std::vector<Real> workspace_coarse;
  // std::vector<Real> workspace_op;

  for (Uint br = 0; br < m_coarse_matrix_sparsity.nb_lines(); ++br)
  {
    const common::Range1D<Int> coarse_rows = m_coarse_matrix_sparsity.sparse_line_limits(br);
    const common::Range1D<Int> fine_rows   = m_fine_matrix_sparsity.sparse_line_limits(br);

    workspace_coarse.resize(coarse_rows.size());
    workspace_fine.resize(fine_rows.size());

    for (Uint r = coarse_rows.lbound(), i = 0; r <= coarse_rows.ubound(); ++r, ++i)
    {
      workspace_coarse[i] = u_coarse.value(r, 0);
    }

    const math::DenseConstVecView<Real> loc_u_coarse(workspace_coarse.data(),
                                                     workspace_coarse.size());
    math::DenseVecView<Real> loc_u_fine(workspace_fine.data(), workspace_fine.size());

    /*
    const math::DenseConstMatView<Real> restrict_op =
    m_local_restriction_ops->const_op(br);

    workspace_op.resize(restrict_op.rows() * restrict_op.cols());
    math::DenseMatView<Real> prolongation_op(workspace_op.data(),
    restrict_op.rows(), restrict_op.cols(), restrict_op.rows());

    for (Uint r = 0; r < restrict_op.rows(); ++r)
    {
      for (Uint c = 0; c < restrict_op.cols(); ++c)
      {
        prolongation_op(c, r) = restrict_op(r, c);
      }
    }
    */

    const math::DenseConstMatView<Real> prolongation_op = m_local_prolongation_ops->const_op(br);

    loc_u_fine = prolongation_op * loc_u_coarse;

    for (Uint r = fine_rows.lbound(), i = 0; r <= fine_rows.ubound(); ++r, ++i)
    {
      const Real val = loc_u_fine[i];
      u_fine.insert_value(r, val, 0);
    }
  } // Loop over cells
}

// ----------------------------------------------------------------------------

const math::BlockMatrixSparsityPattern<Int> &CoarseScaleCorrection::coarse_block_sparsity() const
{
  return m_coarse_matrix_sparsity;
}

// ----------------------------------------------------------------------------

const math::BlockMatrixSparsityPattern<Int> &CoarseScaleCorrection::fine_block_sparsity() const
{
  return m_fine_matrix_sparsity;
}

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit
