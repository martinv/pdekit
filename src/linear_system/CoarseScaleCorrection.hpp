#ifndef PDEKIT_Interpolation_Coarse_Scale_Correction_hpp
#define PDEKIT_Interpolation_Coarse_Scale_Correction_hpp

#include "PDEKit_Config.hpp"

#include <ctime>
#include <memory>
#include <unordered_map>

#include "linear_system/LSAssembler.hpp"
#include "linear_system/LocalTransferOps.hpp"
#include "linear_system/TpetraCrsMatrix.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "Tpetra_MultiVector.hpp"
#endif

namespace pdekit
{

namespace ls
{

class CoarseScaleCorrection
{

  public:
  /// Constructor
  CoarseScaleCorrection();

  /// Destructor
  ~CoarseScaleCorrection();

  /// Compute block sparsity pattern of the coarse matrix
  /// @param mesh_dual_graph_crs ... a CRS-style information about neighbours
  /// of each active
  ///                                cell in the mesh
  /// @param cell_source_dofs ... dof handler of the fine-level discretization
  /// @param cell_reordering  ... cell-wise reordering map
  ///                             The degrees of freedom corresponding to cell
  ///                             with active index 'i' are stored in the
  ///                             block row on position cell_reordering[i] in
  ///                             the system matrix
  /// @param nb_eq            ... number of equations (multiplicity of each
  /// DOF)
  /// @return coarse_matrix_sparsity ... sparsity pattern of the coarse-scale
  ///                                    block-sparse matrix

  void build_block_matrix_sparsity_pattern(
      const common::BlockArray<Uint, Uint> &mesh_dual_graph_crs,
      const std::vector<Uint> &coarse_block_sizes, const std::vector<Uint> &fine_block_sizes,
      const std::vector<Uint> &cell_reordering,
      std::unique_ptr<ls::LocalTransferOps<Real>> &&restriction_ops,
      std::unique_ptr<ls::LocalTransferOps<Real>> &&prolongation_ops);

  /// Build an array of matrices for uniform restriction to
  /// lower-order polynomial space
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  void fill_coarse_level_matrix(
      const TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &fine_level_mat,
      TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &coarse_level_mat);

  /// Build the RHS of the restricted system
  void restrict_vec(const TpetraMultiVector<Real> &rhs_fine,
                    TpetraMultiVector<Real> &rhs_coarse) const;

  /// Build the RHS of the original (fine-level) system
  void prolongate_vec(const TpetraMultiVector<Real> &u_coarse,
                      TpetraMultiVector<Real> &u_fine) const;

#if PDEKIT_HAVE_TRILINOS
  /// Build the RHS of the restricted system
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  void restrict_vec(
      const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &u_fine,
      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &u_coarse) const;
#endif

#if PDEKIT_HAVE_TRILINOS
  /// Build the RHS of the original (fine-level) system
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  void prolongate_vec(
      const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &u_coarse,
      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &u_fine) const;
#endif

  /// Return the block sparsity pattern for the coarse-level matrix
  const math::BlockMatrixSparsityPattern<Int> &coarse_block_sparsity() const;

  /// Return the block sparsity pattern for the coarse-level matrix
  const math::BlockMatrixSparsityPattern<Int> &fine_block_sparsity() const;

  private:
  /// TYPES
  struct UnsignedIntPairHasher
  {
    inline std::size_t operator()(const std::pair<Uint, Uint> &key) const
    {
      return key.first ^ key.second;
    }
  };

  /// DATA
  /// Cell reordering - inverse map
  /// Suppose that cell j becomes cell i after reordering
  /// Then m_cell_reordering[i] == j, i.e. m_cell_reordering[i] stores
  /// the original cell index (before reordering)
  std::vector<Uint> m_inverse_cell_reordering;

  /// Sparsity pattern for coarse scale system matrix
  math::BlockMatrixSparsityPattern<Int> m_coarse_matrix_sparsity;

  /// Sparsity pattern for coarse scale system matrix
  math::BlockMatrixSparsityPattern<Int> m_fine_matrix_sparsity;

  /// An array of block transfer operators
  std::unique_ptr<ls::LocalTransferOps<Real>> m_local_restriction_ops;
  std::unique_ptr<ls::LocalTransferOps<Real>> m_local_prolongation_ops;
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void CoarseScaleCorrection::fill_coarse_level_matrix(
    const TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &fine_level_mat,
    TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &coarse_level_mat)
{
  using block_row_type     = math::BlockMatrixSparsityPattern<Int>::sparse_block_line_type;
  using block_row_ids_type = math::BlockMatrixSparsityPattern<Int>::line_block_ids_type;

  std::vector<Real> coarse_block_data;
  std::vector<Real> fine_block_data;
  // std::vector<Real> prolongation_data_R;

  // Loop over block rows, and for each block in row, compute P^T * A * P,
  // where P   ... prolongation operator
  //       A   ... fine level matrix
  //       P^T ... restriction operator
  for (Uint br = 0; br < m_coarse_matrix_sparsity.nb_lines(); ++br)
  {
    const common::Range1D<Int> coarse_rows = m_coarse_matrix_sparsity.sparse_line_limits(br);
    const common::Range1D<Int> fine_rows   = m_fine_matrix_sparsity.sparse_line_limits(br);

    const block_row_type coarse_block_row = m_coarse_matrix_sparsity.sparse_line(br);
    const block_row_type fine_block_row   = m_fine_matrix_sparsity.sparse_line(br);

    const block_row_ids_type coarse_block_ids = m_coarse_matrix_sparsity.sparse_line_block_ids(br);

    // Loop over block columns of given row
    for (Uint bc = 0; bc < coarse_block_row.size(); ++bc)
    {
      const common::Range1D<Int> coarse_cols = coarse_block_row[bc];
      const common::Range1D<Int> fine_cols   = fine_block_row[bc];

      coarse_block_data.resize(coarse_rows.size() * coarse_cols.size());
      math::DenseMatView<Scalar> coarse_block_view(coarse_block_data.data(), coarse_cols.size(),
                                                   coarse_rows.size(), coarse_cols.size());

      fine_block_data.resize(fine_rows.size() * fine_cols.size());
      math::DenseMatView<Scalar> fine_block_view(fine_block_data.data(), fine_cols.size(),
                                                 fine_rows.size(), fine_cols.size());
      fine_level_mat.get_block(fine_rows, fine_cols, fine_block_view);

      // This is P^T in product P^T * A_fine * P
      // The transpose of prolongation operator is restriction operator
      const math::DenseConstMatView<Scalar> restriction_op_L =
          m_local_restriction_ops->const_op(br);

      /*
      // Restriction operator (i.e. the transpose of the last term in P^T
      * A_fine * P) const Int block_col_id = coarse_block_ids[bc]; const
      math::DenseConstMatView<Scalar> restriction_op_R =
          m_local_restriction_ops->const_op(block_col_id);

      // We need to generate P - the transpose of restriction_op_R

      prolongation_data_R.resize(restriction_op_R.rows() *
      restriction_op_R.cols()); math::DenseMatView<Scalar>
      prolongation_op_R(prolongation_data_R.data(),
      restriction_op_R.rows(), restriction_op_R.cols(),
      restriction_op_R.rows());

      for (Uint r = 0; r < restriction_op_R.rows(); ++r)
      {
        for (Uint c = 0; c < restriction_op_R.cols(); ++c)
        {
          prolongation_op_R(c, r) = restriction_op_R(r, c);
        }
      }
      */

      const Int block_col_id = coarse_block_ids[bc];
      const math::DenseConstMatView<Scalar> prolongation_op_R =
          m_local_prolongation_ops->const_op(block_col_id);

      /*
      std::cout << "Multiply: [" << restriction_op_L.rows() << "x" <<
      restriction_op_L.cols()
                << "] [" << fine_block_view.rows() << "x" <<
      fine_block_view.cols() << "] ["
                << prolongation_op_R.rows() << "x" <<
      prolongation_op_R.cols()
      << "]" << std::endl; std::cout << "Store to: [" <<
      coarse_block_view.rows() << "x" << coarse_block_view.cols()
                << "]\n"
                << std::endl;
      */

      coarse_block_view = restriction_op_L * fine_block_view * prolongation_op_R;
      coarse_level_mat.insert_block(coarse_rows, coarse_cols, coarse_block_view);

    } // Loop over blocks of one block-row
  }   // Loop over block rows
}

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void CoarseScaleCorrection::restrict_vec(
    const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &u_fine,
    Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &u_coarse) const
{
  std::vector<Real> workspace_fine, workspace_coarse;
  // Get column 0 of multivector 'u_fine'
  Teuchos::ArrayRCP<const Scalar> u_fine_view = u_fine.getData(0);
  // Get column 0 of multivector 'u_coarse'
  Teuchos::ArrayRCP<Scalar> u_coarse_view = u_coarse.getDataNonConst(0);

  for (Uint br = 0; br < m_coarse_matrix_sparsity.nb_lines(); ++br)
  {
    const common::Range1D<Int> coarse_rows = m_coarse_matrix_sparsity.sparse_line_limits(br);
    const common::Range1D<Int> fine_rows   = m_fine_matrix_sparsity.sparse_line_limits(br);

    workspace_coarse.resize(coarse_rows.size());
    workspace_fine.resize(fine_rows.size());

    for (Uint r = fine_rows.lbound(), i = 0; r <= fine_rows.ubound(); ++r, ++i)
    {
      workspace_fine[i] = u_fine_view[r];
    }

    const math::DenseConstVecView<Real> loc_rhs_fine(workspace_fine.data(), workspace_fine.size());
    math::DenseVecView<Real> loc_u_coarse(workspace_coarse.data(), workspace_coarse.size());

    const math::DenseConstMatView<Real> restrict_op = m_local_restriction_ops->const_op(br);

    loc_u_coarse = restrict_op * loc_rhs_fine;

    for (Uint r = coarse_rows.lbound(), i = 0; r <= coarse_rows.ubound(); ++r, ++i)
    {
      const Real val   = loc_u_coarse[i];
      u_coarse_view[r] = val;
    }
  } // Loop over cells
}
#endif

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void CoarseScaleCorrection::prolongate_vec(
    const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &u_coarse,
    Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &u_fine) const
{
  std::vector<Real> workspace_fine;
  std::vector<Real> workspace_coarse;
  // std::vector<Real> workspace_op;

  // Get column 0 of multivector 'u_coarse'
  Teuchos::ArrayRCP<const Scalar> u_coarse_view = u_coarse.getData(0);
  // Get column 0 of multivector 'u_fine'
  Teuchos::ArrayRCP<Scalar> u_fine_view = u_fine.getDataNonConst(0);

  for (Uint br = 0; br < m_coarse_matrix_sparsity.nb_lines(); ++br)
  {
    const common::Range1D<Int> coarse_rows = m_coarse_matrix_sparsity.sparse_line_limits(br);
    const common::Range1D<Int> fine_rows   = m_fine_matrix_sparsity.sparse_line_limits(br);

    workspace_coarse.resize(coarse_rows.size());
    workspace_fine.resize(fine_rows.size());

    for (Uint r = coarse_rows.lbound(), i = 0; r <= coarse_rows.ubound(); ++r, ++i)
    {
      workspace_coarse[i] = u_coarse_view[r];
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
      u_fine_view[r] = val;
    }
  } // Loop over cells
}
#endif

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif
