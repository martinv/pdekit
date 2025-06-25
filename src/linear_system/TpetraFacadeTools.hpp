#ifndef PDEKIT_Trilinos_Tpetra_Facade_Tools_hpp
#define PDEKIT_Trilinos_Tpetra_Facade_Tools_hpp

#include "PDEKit_Config.hpp"
#include "common/PDEKit.hpp"
#include "graph/Graph.hpp"
#include "math/BlockMatrixSparsityPattern.hpp"
#include "math/DenseMatView.hpp"
#include "math/MatrixSparsityPattern.hpp"

#if PDEKIT_HAVE_TRILINOS

#include "Teuchos_RCP.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"

#endif

// This is an interface class operating on __raw__ Tpetra objects,
// not on PDEKit wrappers.
// Its interface should be exploited in wrappers of Tpetra objects.

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------

class TpetraFacadeTools
{
  public:
#if PDEKIT_HAVE_TRILINOS

  // ----------------------------------------------------------------------------

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  static void create_matrix_sparsity(
      const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
      const Uint nb_entries_per_row,
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &matrix);

  // ----------------------------------------------------------------------------

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  static void create_matrix_sparsity(
      const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
      const graph::Graph<GlobalOrdinal> &matrix_nonzeros,
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &matrix);

  // ----------------------------------------------------------------------------

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  static void create_matrix_sparsity(
      const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
      const math::MatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros,
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &matrix);

  // ----------------------------------------------------------------------------

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  static void create_matrix_sparsity(
      const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
      const math::BlockMatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros,
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &matrix);

  // ----------------------------------------------------------------------------

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  static void get_block(
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &matrix,
      const common::Range1D<GlobalOrdinal> row_range,
      const common::Range1D<GlobalOrdinal> col_range, math::DenseMatView<Scalar> &block_data);

  // ----------------------------------------------------------------------------

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  static void insert_block(
      const common::Range1D<GlobalOrdinal> row_range,
      const common::Range1D<GlobalOrdinal> col_range, const math::DenseMatView<Scalar> &block_data,
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &matrix);

  // ----------------------------------------------------------------------------

#endif
};

// ----------------------------------------------------------------------------
#if PDEKIT_HAVE_TRILINOS
// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraFacadeTools::create_matrix_sparsity(
    const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
    const Uint nb_entries_per_row,
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &matrix)
{
  using trilinos_matrix_type =
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;

  matrix.reset(new trilinos_matrix_type(map, nb_entries_per_row));
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraFacadeTools::create_matrix_sparsity(
    const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
    const graph::Graph<GlobalOrdinal> &matrix_nonzeros,
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &matrix)
{
  using trilinos_matrix_type =
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;

  Teuchos::ArrayRCP<size_t> nb_cols_per_row(matrix_nonzeros.nb_vertices());

  for (size_t i = 0; i < matrix_nonzeros.nb_vertices(); ++i)
  {
    nb_cols_per_row[i] = matrix_nonzeros.number_adj_vertices(i);
  }

  matrix.reset(new trilinos_matrix_type(map, nb_cols_per_row, Tpetra::StaticProfile));

  std::vector<Scalar> zero_values;
  std::vector<GlobalOrdinal> adj_vertices;

  if (!matrix->isFillActive())
  {
    matrix->resumeFill();
  }

  for (size_t r = 0; r < matrix_nonzeros.nb_vertices(); ++r)
  {
    zero_values.resize(matrix_nonzeros.number_adj_vertices(r));
    adj_vertices.resize(matrix_nonzeros.number_adj_vertices(r));

    const std::pair<typename graph::Graph<GlobalOrdinal>::adj_vertex_const_iterator,
                    typename graph::Graph<GlobalOrdinal>::adj_vertex_const_iterator>
        adj_vert_iterators = matrix_nonzeros.adjacent_vertices(r);

    Uint i = 0;
    for (typename graph::Graph<GlobalOrdinal>::adj_vertex_const_iterator it =
             adj_vert_iterators.first;
         it != adj_vert_iterators.second; ++it)
    {
      adj_vertices[i++] = *it;
    }

    const GlobalOrdinal row_idx = static_cast<GlobalOrdinal>(r);

    Teuchos::ArrayView<const GlobalOrdinal> indices_view(adj_vertices.data(), adj_vertices.size());
    Teuchos::ArrayView<const Scalar> values_view(zero_values.data(), zero_values.size());
    matrix->insertGlobalValues(row_idx, indices_view, values_view);
  }

  matrix->fillComplete();
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraFacadeTools::create_matrix_sparsity(
    const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
    const math::MatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros,
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &matrix)
{
  using trilinos_matrix_type =
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;

  Teuchos::ArrayRCP<size_t> nb_cols_per_row(matrix_nonzeros.nb_lines());

  for (size_t i = 0; i < matrix_nonzeros.nb_lines(); ++i)
  {
    const typename math::MatrixSparsityPattern<GlobalOrdinal>::sparse_line_type sparse_line =
        matrix_nonzeros.sparse_line(i);
    nb_cols_per_row[i] = sparse_line.size();
  }

  matrix.reset(new trilinos_matrix_type(map, nb_cols_per_row, Tpetra::StaticProfile));

  std::vector<Scalar> zero_values;
  std::vector<GlobalOrdinal> adj_vertices;

  if (!matrix->isFillActive())
  {
    matrix->resumeFill();
  }

  for (size_t r = 0; r < matrix_nonzeros.nb_lines(); ++r)
  {
    const typename math::MatrixSparsityPattern<GlobalOrdinal>::sparse_line_type sparse_line =
        matrix_nonzeros.sparse_line(r);

    zero_values.resize(sparse_line.size());
    adj_vertices.resize(sparse_line.size());

    for (Uint i = 0; i < sparse_line.size(); ++i)
    {
      adj_vertices[i] = sparse_line[i];
    }

    const GlobalOrdinal row_idx = static_cast<GlobalOrdinal>(r);

    Teuchos::ArrayView<const GlobalOrdinal> indices_view(adj_vertices.data(), adj_vertices.size());
    Teuchos::ArrayView<const Scalar> values_view(zero_values.data(), zero_values.size());
    matrix->insertGlobalValues(row_idx, indices_view, values_view);
  }

  matrix->fillComplete();
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraFacadeTools::create_matrix_sparsity(
    const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map,
    const math::BlockMatrixSparsityPattern<GlobalOrdinal> &matrix_nonzeros,
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &matrix)
{
  using trilinos_matrix_type =
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;

  const common::ArrayShape<_2D, Uint> mat_shape = matrix_nonzeros.size();

  Teuchos::ArrayRCP<size_t> nb_cols_per_row(mat_shape.size(0));

  for (Uint br = 0; br < matrix_nonzeros.nb_lines(); ++br)
  {
    const common::Range1D<GlobalOrdinal> row_limits = matrix_nonzeros.sparse_line_limits(br);
    const common::ArrayView<const common::Range1D<GlobalOrdinal>, _1D, Uint> row_blocks =
        matrix_nonzeros.sparse_line(br);
    Uint line_len = 0;

    // Loop over block columns
    for (Uint bc = 0; bc < row_blocks.size(); ++bc)
    {
      line_len += row_blocks[bc].size();
    }
    for (GlobalOrdinal r = row_limits.lbound(); r <= row_limits.ubound(); ++r)
    {
      nb_cols_per_row[r] = line_len;
    }
  }

  matrix.reset(new trilinos_matrix_type(map, nb_cols_per_row, Tpetra::StaticProfile));

  std::vector<Scalar> zero_values;
  std::vector<GlobalOrdinal> adj_vertices;

  if (!matrix->isFillActive())
  {
    matrix->resumeFill();
  }

  for (Uint br = 0; br < matrix_nonzeros.nb_lines(); ++br)
  {
    const common::Range1D<GlobalOrdinal> row_limits = matrix_nonzeros.sparse_line_limits(br);
    const common::ArrayView<const common::Range1D<GlobalOrdinal>, _1D, Uint> row_blocks =
        matrix_nonzeros.sparse_line(br);

    adj_vertices.resize(0);
    zero_values.resize(0);

    // Loop over block columns
    for (Uint bc = 0; bc < row_blocks.size(); ++bc)
    {
      for (GlobalOrdinal c = row_blocks[bc].lbound(); c <= row_blocks[bc].ubound(); ++c)
      {
        adj_vertices.push_back(c);
        zero_values.push_back(Scalar());
      }
    }

    Teuchos::ArrayView<const GlobalOrdinal> indices_view(adj_vertices.data(), adj_vertices.size());
    Teuchos::ArrayView<const Scalar> values_view(zero_values.data(), zero_values.size());

    for (GlobalOrdinal r = row_limits.lbound(); r <= row_limits.ubound(); ++r)
    {
      const GlobalOrdinal row_idx = static_cast<GlobalOrdinal>(r);
      matrix->insertGlobalValues(row_idx, indices_view, values_view);
    }
  }

  matrix->fillComplete();
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraFacadeTools::get_block(
    const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &matrix,
    const common::Range1D<GlobalOrdinal> row_range, const common::Range1D<GlobalOrdinal> col_range,
    math::DenseMatView<Scalar> &block_data)
{
  // std::cout << "Block size = [" << nb_block_rows << " x " << nb_block_cols
  // <<
  // "]" << std::endl;

  Teuchos::ArrayView<const GlobalOrdinal> col_indices_view;
  Teuchos::ArrayView<const Scalar> values_view;

  const GlobalOrdinal row_offset = row_range.lbound();
  const GlobalOrdinal col_offset = col_range.lbound();

  for (GlobalOrdinal row = 0; row < row_range.size(); ++row)
  {
    // std::cout << "Row = " << row + row_offset << std::endl;

    matrix.getLocalRowView(row + row_offset, col_indices_view, values_view);

    // std::cout << "Column indices = " << col_indices_view << std::endl;
    // std::cout << "Extracted values = " << values_view << std::endl;

    const LocalOrdinal *row_begin      = col_indices_view.getRawPtr();
    const LocalOrdinal *row_end        = row_begin + col_indices_view.size();
    const LocalOrdinal *search_idx_pos = std::lower_bound(row_begin, row_end, col_offset);

    const Uint local_left_offset_in_row = std::distance(row_begin, search_idx_pos);

    // std::cout << "Position of left corner in row = " <<
    // local_left_offset_in_row << std::endl;

    for (Uint col = 0; col < col_range.size(); ++col)
    {
      block_data(row, col) = values_view[local_left_offset_in_row + col];
    }
  }
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraFacadeTools::insert_block(
    const common::Range1D<GlobalOrdinal> row_range, const common::Range1D<GlobalOrdinal> col_range,
    const math::DenseMatView<Scalar> &block_data,
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &matrix)
{
  const bool mat_storage_is_dynamic = (matrix.getProfileType() == Tpetra::DynamicProfile);

  std::vector<GlobalOrdinal> indices(col_range.size());
  std::vector<Scalar> values(col_range.size());

  for (Uint c = 0; c < col_range.size(); ++c)
  {
    indices[c] = col_range.lbound() + c;
  }

  Teuchos::ArrayView<GlobalOrdinal> indices_view(indices.data(), indices.size());
  Teuchos::ArrayView<Scalar> values_view(values.data(), values.size());

  for (Uint r = 0; r < row_range.size(); ++r)
  {
    const GlobalOrdinal row_idx = static_cast<GlobalOrdinal>(row_range.lbound() + r);
    for (Uint c = 0; c < col_range.size(); ++c)
    {
      values_view[c] = block_data(r, c);
    }
    if (mat_storage_is_dynamic)
    {
      matrix.insertGlobalValues(row_idx, indices_view, values_view);
    }
    else // Should be static
    {
      matrix.replaceGlobalValues(row_idx, indices_view, values_view);
    }
  }
}

// ----------------------------------------------------------------------------
#endif
// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif
