#ifndef PDEKIT_Math_Block_Matrix_Sparsity_Pattern_hpp
#define PDEKIT_Math_Block_Matrix_Sparsity_Pattern_hpp

#include <cmath>
#include <fstream>
#include <iostream>

#include "common/ArrayShape.hpp"
#include "common/BlockArray.hpp"
#include "common/BlockMultiArray.hpp"
#include "common/Meta.hpp"
#include "common/Range1D.hpp"
#include "math/Matrix.hpp"
#include "math/MatrixStorageOrder.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename IdxType, bool SO = DefaultMatrixStorageOrder>
class BlockMatrixSparsityPattern
{
  public:
  using sparse_block_line_type = common::ArrayView<const common::Range1D<IdxType>, _1D, Uint>;
  using line_block_ids_type    = common::ArrayView<const IdxType, _1D, Uint>;

  private:
  using block_sparsity_storage_type = common::BlockMultiArray<common::Range1D<IdxType>, IdxType>;
  using block_line_thickness_store_type = std::vector<common::Range1D<IdxType>>;

  public:
  /// Constructor, takes the number of rows and columns
  explicit BlockMatrixSparsityPattern(Uint m = 0, Uint n = 0);

  /// Copy constructor
  BlockMatrixSparsityPattern(const BlockMatrixSparsityPattern &rhs);

  /// Copy constructor from sparsity pattern with another index type
  template <typename OtherIdxType>
  explicit BlockMatrixSparsityPattern(const BlockMatrixSparsityPattern<OtherIdxType, SO> &init);

  /// Destructor
  ~BlockMatrixSparsityPattern();

  /// Assignment operator
  BlockMatrixSparsityPattern &operator=(const BlockMatrixSparsityPattern &rhs);

  /// Resize the matrix
  void resize(const Uint m, const Uint n);

  /// Build the sparsity pattern
  void build_sparsity(
      const std::vector<std::vector<common::Range1D<IdxType>>> &nonzero_block_limits,
      const std::vector<std::vector<IdxType>> &block_ids,
      const std::vector<common::Range1D<IdxType>> &block_line_positions);

  /// Build the sparsity pattern
  void build_sparsity(
      std::unique_ptr<common::BlockMultiArray<common::Range1D<IdxType>, IdxType>> &&nonzero_blocks,
      std::unique_ptr<std::vector<common::Range1D<IdxType>>> &&block_line_positions);

  /// Return the sparsity pattern dimensions
  const common::ArrayShape<_2D, Uint> size() const;

  /// Return number of nonzero values
  Uint nb_nz() const;

  /// Number of lines in sparsity pattern
  Uint nb_lines() const;

  /// Return index position
  /// All nonzero indices are stored as compressed row storage or compressed
  /// column storage, i.e. nonzero entries are concatenated in one large array
  /// going row by row (or column by column)
  /// This function returns the position of pair (row,col) in that storage
  /// array
  /// @param row ... row index
  /// @param col ... col index
  /// @return the position of (row, col) among all indices stored
  std::tuple<Uint, Uint> idx_position(const IdxType row, const IdxType col) const;

  /// Get one sparse line (which can be row or column, depending on the
  /// storage order
  const sparse_block_line_type sparse_line(const IdxType line_id) const;

  /// Get the block ids in one line
  const line_block_ids_type sparse_line_block_ids(const IdxType line_id) const;

  /// Get one sparse line (which can be row or column, depending on the
  /// storage order
  /// @return a tuple containing the first and last index (included) where
  /// line data is stored
  const common::Range1D<IdxType> sparse_line_limits(const IdxType line_id) const;

  /// Return the memory used
  Real mem_used_mb() const;

  /// Print the sparsity pattern
  void print_sparsity() const;

  /// Print the sparsity pattern to svg file
  void print_svg(const std::string &filename) const;

  /// Print the sparsity to a vtu file
  void print_vtu(const std::string &filename) const;

  private:
  /// TYPES

  enum StoredFields
  {
    BLOCK_LIMITS = 0,
    BLOCK_ID     = 1
  };

  struct CompareBlocks
  {

    bool operator()(const common::Range1D<IdxType> &a, const common::Range1D<IdxType> &b) const;
  };

  // -----------------------------------------------------------------------------
  // A helper structure which treats MatrixSparsityPattern as if it
  // represented compressed ROW storage: 'line' is a matrix row and 'index in
  // line' represents a column index in that row
  // -----------------------------------------------------------------------------

  struct CompressedRowStorageIndexer
  {
    static const std::tuple<Uint, Uint> idx_position(
        const std::vector<common::Range1D<IdxType>> &line_positions,
        const common::BlockMultiArray<common::Range1D<IdxType>, IdxType> &line_blocks,
        const IdxType row, const IdxType col, const common::ArrayShape<_2D, Uint> &mat_shape);
  };

  // -----------------------------------------------------------------------------
  // A helper structure which treats MatrixSparsityPattern as if it
  // represented compressed COLUMN storage: 'line' is a matrix column and
  // 'index in line' represents a row index in that column
  // -----------------------------------------------------------------------------

  struct CompressedColStorageIndexer
  {
    static const std::tuple<Uint, Uint> idx_position(
        const std::vector<common::Range1D<IdxType>> &line_positions,
        const common::BlockMultiArray<common::Range1D<IdxType>, IdxType> &line_blocks,
        const IdxType row, const IdxType col, const common::ArrayShape<_2D, Uint> &mat_shape);
  };

  // Set the type indexing either to row-major or column-major
  using indexer_type = typename common::SelectType<SO == RowMajor, CompressedRowStorageIndexer,
                                                   CompressedColStorageIndexer>::type;

  /// METHODS

  /// Compute the number of nonzeros in one block
  inline static IdxType block_volume(const common::Range1D<IdxType> &col_span,
                                     const common::Range1D<IdxType> &row_span)
  {
    return col_span.size() * row_span.size();
  }

  /// DATA
  /// Number of matrix rows and columns
  common::ArrayShape<_2D, Uint> m_shape;

  /// Number of nonzeros
  Uint m_nb_nz;

  /// Sparsity pattern
  /// Each tuple<IdxType, IdxType> in block_sparsity_storage_type
  /// denotes the beginning and end of one block in one line
  /// In other words, each block in block_sparsity_storage_type is a sequence
  /// of pairs [C0_s, C0_e], [C1_s, C1_e], [C2_s, C2_e], [C3_s, C3_e] ...
  /// where each pair [C*_s, C*_e] denotes the beginning and end (inclusive)
  /// of one block
  std::unique_ptr<block_sparsity_storage_type> m_sparse_store;

  /// Array that stores the pair <start index, end index> for each line
  /// In block-wise CRS format, these would be the row index ranges
  std::unique_ptr<block_line_thickness_store_type> m_block_line_positions;
};

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
BlockMatrixSparsityPattern<IdxType, SO>::BlockMatrixSparsityPattern(Uint m, Uint n)
    : m_shape(m, n), m_nb_nz(0u), m_sparse_store(new block_sparsity_storage_type),
      m_block_line_positions(new block_line_thickness_store_type)
{
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
BlockMatrixSparsityPattern<IdxType, SO>::BlockMatrixSparsityPattern(
    const BlockMatrixSparsityPattern &rhs)
    : m_shape(rhs.m_shape), m_nb_nz(rhs.m_nb_nz)
{
  m_sparse_store->resize(rhs->m_sparse_store.size(), rhs->m_sparse_store.nb_blocks());
  (*m_sparse_store) = (*rhs.m_sparse_store);

  m_block_line_positions->resize(rhs->m_line_widths.size());
  (*m_block_line_positions) = (*rhs.m_block_line_positions);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
template <typename OtherIdxType>
BlockMatrixSparsityPattern<IdxType, SO>::BlockMatrixSparsityPattern(
    const BlockMatrixSparsityPattern<OtherIdxType, SO> &init)
    : m_shape(init.m_shape), m_nb_nz(init.m_nb_nz)
{
  m_sparse_store->resize(init->m_sparse_store.size(), init->m_sparse_store.nb_blocks());
  (*m_sparse_store) = (*init.m_sparse_store);

  m_block_line_positions->resize(init->m_line_widths.size());
  (*m_block_line_positions) = (*init.m_line_widths);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
BlockMatrixSparsityPattern<IdxType, SO>::~BlockMatrixSparsityPattern()
{
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
BlockMatrixSparsityPattern<IdxType, SO> &BlockMatrixSparsityPattern<IdxType, SO>::operator=(
    const BlockMatrixSparsityPattern &rhs)
{
  m_shape = rhs.m_shape;
  m_nb_nz = rhs.m_nb_nz;

  m_sparse_store->resize(rhs->m_sparse_store.size(), rhs->m_sparse_store.nb_blocks());
  (*m_sparse_store) = (*rhs.m_sparse_store);

  m_block_line_positions->resize(rhs->m_line_widths.size());
  (*m_block_line_positions) = (*rhs.m_block_line_positions);

  return *this;
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void BlockMatrixSparsityPattern<IdxType, SO>::resize(const Uint m, const Uint n)
{
  m_shape = common::ArrayShape<_2D, Uint>(m, n);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void BlockMatrixSparsityPattern<IdxType, SO>::build_sparsity(
    const std::vector<std::vector<common::Range1D<IdxType>>> &nonzero_block_limits,
    const std::vector<std::vector<IdxType>> &block_ids,
    const std::vector<common::Range1D<IdxType>> &block_line_positions)

{
  Uint nb_block_entries = 0;
  for (std::vector<common::Range1D<IdxType>> const &line : nonzero_block_limits)
  {
    nb_block_entries += line.size();
  }

  m_sparse_store->reserve(nb_block_entries, nonzero_block_limits.size());
  // m_sparsity.resize(0, 0);

  m_nb_nz           = 0u;
  IdxType max_entry = IdxType();

  for (Uint b = 0; b < nonzero_block_limits.size(); ++b)
  {
    const common::Range1D<IdxType> block_line_width = block_line_positions[b];

    const std::vector<common::Range1D<IdxType>> &line = nonzero_block_limits[b];
    const common::ArrayView<const common::Range1D<IdxType>, _1D, Uint> block_limits_view(
        line.data(), line.size());

    const std::vector<IdxType> &line_block_ids = block_ids[b];
    const common::ArrayView<const IdxType, _1D, Uint> block_id_view(line_block_ids.data(),
                                                                    line_block_ids.size());
    m_sparse_store->create_back_block(line.size());
    m_sparse_store->template fill_last_block<BLOCK_LIMITS>(block_limits_view);
    m_sparse_store->template fill_last_block<BLOCK_ID>(block_id_view);

    for (Uint i = 0; i < block_limits_view.size(); ++i)
    {
      max_entry = std::max(max_entry, block_limits_view[i].ubound());
      m_nb_nz += block_volume(block_line_width, block_limits_view[i]);
    }
  }

  m_block_line_positions->resize(block_line_positions.size());
  (*m_block_line_positions) = block_line_positions;

  IdxType max_line_pos = IdxType();
  for (Uint l = 0; l < block_line_positions.size(); ++l)
  {
    max_line_pos = std::max(max_line_pos, block_line_positions[l].ubound());
  }

  m_shape = (SO == RowMajor) ? common::ArrayShape<_2D, Uint>(max_line_pos + 1, max_entry + 1)
                             : common::ArrayShape<_2D, Uint>(max_entry + 1, max_line_pos + 1);

  CompareBlocks comp;
  m_sparse_store->template sort_blocks<BLOCK_LIMITS>(comp);
  // m_sparse_store->template sort_blocks<BLOCK_ID>(std::less<IdxType>());
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void BlockMatrixSparsityPattern<IdxType, SO>::build_sparsity(
    std::unique_ptr<common::BlockMultiArray<common::Range1D<IdxType>, IdxType>> &&nonzero_blocks,
    std::unique_ptr<std::vector<common::Range1D<IdxType>>> &&block_line_positions)
{
  m_sparse_store         = std::move(nonzero_blocks);
  m_block_line_positions = std::move(block_line_positions);

  m_nb_nz           = 0u;
  IdxType max_entry = IdxType();

  for (Uint b = 0; b < m_sparse_store->nb_blocks(); ++b)
  {
    const common::Range1D<IdxType> block_line_width = (*m_block_line_positions)[b];

    const common::ArrayView<const common::Range1D<IdxType>, _1D, Uint> block_view =
        m_sparse_store->template const_block<BLOCK_LIMITS>(b);

    for (Uint i = 0; i < block_view.size(); ++i)
    {
      max_entry = std::max(max_entry, block_view[i].ubound());
      m_nb_nz += block_volume(block_line_width, block_view[i]);
    }
  }

  IdxType max_line_pos = IdxType();
  for (Uint l = 0; l < m_block_line_positions->size(); ++l)
  {
    max_line_pos = std::max(max_line_pos, (*m_block_line_positions)[l].ubound());
  }

  m_shape = (SO == RowMajor) ? common::ArrayShape<_2D, Uint>(max_line_pos + 1, max_entry + 1)
                             : common::ArrayShape<_2D, Uint>(max_entry + 1, max_line_pos + 1);

  CompareBlocks comp;
  m_sparse_store->template sort_blocks<BLOCK_LIMITS>(comp);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
const common::ArrayShape<_2D, Uint> BlockMatrixSparsityPattern<IdxType, SO>::size() const
{
  return m_shape;
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
Uint BlockMatrixSparsityPattern<IdxType, SO>::nb_nz() const
{
  return m_nb_nz;
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
Uint BlockMatrixSparsityPattern<IdxType, SO>::nb_lines() const
{
  return m_sparse_store->nb_blocks();
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
std::tuple<Uint, Uint> BlockMatrixSparsityPattern<IdxType, SO>::idx_position(
    const IdxType row, const IdxType col) const
{
  return indexer_type::idx_position(*m_block_line_positions, *m_sparse_store, row, col, m_shape);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
const typename BlockMatrixSparsityPattern<IdxType, SO>::sparse_block_line_type
BlockMatrixSparsityPattern<IdxType, SO>::sparse_line(const IdxType line_id) const
{
  return m_sparse_store->template const_block<BLOCK_LIMITS>(line_id);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
const typename BlockMatrixSparsityPattern<IdxType, SO>::line_block_ids_type BlockMatrixSparsityPattern<
    IdxType, SO>::sparse_line_block_ids(const IdxType line_id) const
{
  return m_sparse_store->template const_block<BLOCK_ID>(line_id);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
const common::Range1D<IdxType> BlockMatrixSparsityPattern<IdxType, SO>::sparse_line_limits(
    const IdxType line_id) const
{
  // return m_sparse_store->block_limits(line_id);
  return (*m_block_line_positions)[line_id];
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
Real BlockMatrixSparsityPattern<IdxType, SO>::mem_used_mb() const
{
  Real mem_used = 0.0;

  if (m_sparse_store)
  {
    mem_used += m_sparse_store->mem_used_mb();
  }

  if (m_block_line_positions)
  {
    mem_used += 1.e-6 * m_block_line_positions->size() * sizeof(common::Range1D<IdxType>);
  }
  return mem_used;
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void BlockMatrixSparsityPattern<IdxType, SO>::print_sparsity() const
{
  std::cout << "Block matrix sparsity pattern (" << m_shape.size(0) << "," << m_shape.size(1) << ")"
            << std::endl;
  for (Uint b = 0; b < m_sparse_store->nb_blocks(); ++b)
  {
    std::cout << "Sparse line [" << (*m_block_line_positions)[b].lbound() << ","
              << (*m_block_line_positions)[b].ubound() << "] {";

    const common::ArrayView<const common::Range1D<IdxType>, _1D, Uint> block_line =
        m_sparse_store->template const_block<BLOCK_LIMITS>(b);
    for (Uint i = 0; i < block_line.size(); ++i)
    {
      std::cout << " [" << block_line[i].lbound() << "," << block_line[i].ubound() << "]";
    }
    std::cout << " }" << std::endl;
  }
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void BlockMatrixSparsityPattern<IdxType, SO>::print_svg(const std::string &filename) const
{
  std::ofstream out;
  out.open(filename.c_str());

  const Uint m = m_shape.size(0);
  const Uint n = m_shape.size(1);

  const Real border = 0.5;

  out << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" "
         "viewBox=\"0 0 "
      << n + border << " " << m + border
      << " \">\n"
         "<style type=\"text/css\" >\n"
         "     <![CDATA[\n"
         "      rect.pixel {\n"
         "          fill:   #000000;\n"
         "      }\n"
         "    ]]>\n"
         "  </style>\n\n"
         "   <rect width=\""
      << n + border << "\" height=\"" << m + border
      << "\" fill=\"rgb(128, 128, 128)\"/>\n"
         "   <rect x=\""
      << 0.5 * border << "\" y=\"" << 0.5 * border << "\" width=\"" << n << "\" height=\"" << m
      << "\" fill=\"rgb(255, 255, 255)\"/>\n\n";

  for (Uint b = 0; b < m_sparse_store->nb_blocks(); ++b)
  {
    const common::Range1D<IdxType> row_span = (*m_block_line_positions)[b];

    common::ArrayView<const common::Range1D<IdxType>, _1D, Uint> line =
        m_sparse_store->template const_block<BLOCK_LIMITS>(b);
    for (Uint i = 0; i < line.size(); ++i)
    {
      const common::Range1D<IdxType> col_span = line[i];

      for (IdxType r = row_span.lbound(); r <= row_span.ubound(); ++r)
      {
        for (IdxType c = col_span.lbound(); c <= col_span.ubound(); ++c)
        {
          // The 0.05 offset is here to center the actual nonzero
          // entries are squares of size [0.9 x 0.9] . This means that
          // the right and bottom edge of the plot have a gap of
          // width 1.0 - 0.9 = 0.1. Therefore we divide this gap by 2
          // and offset all squares representing nonzero entries by
          // 0.1/2 = 0.05
          out << R"(  <rect class="pixel" x=")" << c + 0.5 * border + 0.05 << "\" y=\""
              << r + 0.5 * border + 0.05 << "\" width=\".9\" height=\".9\"/>\n";
        }
      }
    }
  }
  out << "</svg>" << std::endl;
  out.close();
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void BlockMatrixSparsityPattern<IdxType, SO>::print_vtu(const std::string &filename) const
{
  std::ofstream out;
  out.open(filename.c_str(), std::ios::binary);
  size_t appended_data_offset = 0;

  const Uint nb_blocks           = m_sparse_store->size();
  const uint32_t nb_block_coords = 4 * nb_blocks;

  // out << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
         "byte_order=\"LittleEndian\" "
         "header_type=\"UInt64\">"
      << std::endl;
  out << "  <UnstructuredGrid>" << std::endl;
  // Number of points and cells is equal to number of nonzeros in sparsity
  // pattern
  out << "    <Piece NumberOfPoints=\"" << nb_block_coords << "\" NumberOfCells=\"" << nb_blocks
      << "\">" << std::endl;

  out << "      <PointData>" << std::endl;
  out << "      </PointData>" << std::endl;

  // ---

  out << "      <CellData Scalars=\"scalars\">" << std::endl;
  out << "        <DataArray type=\"UInt32\" Name=\"block_row\" "
         "format=\"appended\" offset=\""
      << appended_data_offset << "\"/>" << std::endl;

  appended_data_offset += nb_blocks * sizeof(uint32_t) + sizeof(uint64_t);

  out << "        <DataArray type=\"UInt32\" Name=\"block_column\" "
         "format=\"appended\" offset=\""
      << appended_data_offset << "\"/>" << std::endl;

  appended_data_offset += nb_blocks * sizeof(uint32_t) + sizeof(uint64_t);

  out << "        <DataArray type=\"UInt32\" Name=\"row\" "
         "format=\"appended\" "
         "offset=\""
      << appended_data_offset << "\"/>" << std::endl;

  appended_data_offset += nb_blocks * sizeof(uint32_t) + sizeof(uint64_t);

  out << "        <DataArray type=\"UInt32\" Name=\"column\" "
         "format=\"appended\" offset=\""
      << appended_data_offset << "\"/>" << std::endl;

  appended_data_offset += nb_blocks * sizeof(uint32_t) + sizeof(uint64_t);

  out << "      </CellData>" << std::endl;

  // ---

  out << "      <Points>" << std::endl;
  out << "        <DataArray type=\"Float64\" Name=\"Nodes\" "
         "NumberOfComponents=\"3\" "
         "format=\"appended\" "
         "offset=\""
      << appended_data_offset << "\"/>" << std::endl;
  out << "      </Points>" << std::endl;

  appended_data_offset += 3 * sizeof(Real) * nb_block_coords + sizeof(uint64_t);

  out << "      <Cells>" << std::endl;

  // Cell connectivity
  out << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
         "format=\"appended\" offset=\""
      << appended_data_offset << "\"/>" << std::endl;
  appended_data_offset += 4 * sizeof(int64_t) * nb_blocks + sizeof(uint64_t);

  // Cell offset data
  out << "        <DataArray type=\"Int64\" Name=\"offsets\" "
         "format=\"appended\" offset=\""
      << appended_data_offset << "\"/>" << std::endl;
  appended_data_offset += sizeof(int64_t) * nb_blocks + sizeof(uint64_t);

  // Cell type data
  out << "        <DataArray type=\"UInt8\" Name=\"types\" "
         "format=\"appended\" "
         "offset=\""
      << appended_data_offset << "\"/>" << std::endl;
  appended_data_offset += sizeof(uint8_t) * nb_blocks + sizeof(uint64_t);

  out << "      </Cells>" << std::endl;
  out << "    </Piece>" << std::endl;
  out << "  </UnstructuredGrid>" << std::endl;
  out << "  <AppendedData encoding=\"raw\">" << std::endl;

  const std::string begin_marker = "                _";
  out << begin_marker;

  // ----------------------------------------
  // Write cell data - block row index
  // ----------------------------------------

  const uint64_t num_block_row_idx_bytes = sizeof(uint32_t) * nb_blocks;
  out.write(reinterpret_cast<const char *>(&num_block_row_idx_bytes),
            sizeof(num_block_row_idx_bytes));

  for (Uint l = 0; l < m_block_line_positions->size(); ++l)
  {
    common::ArrayView<const IdxType, _1D, Uint> line =
        m_sparse_store->template const_block<BLOCK_ID>(l);
    const uint32_t val = static_cast<uint32_t>(l);

    for (Uint i = 0; i < line.size(); ++i)
    {
      out.write(reinterpret_cast<const char *>(&val), sizeof(uint32_t));
    }
  }

  // ----------------------------------------
  // Write cell data - block column index
  // ----------------------------------------

  const uint64_t num_block_col_idx_bytes = sizeof(uint32_t) * nb_blocks;
  out.write(reinterpret_cast<const char *>(&num_block_col_idx_bytes),
            sizeof(num_block_col_idx_bytes));

  for (Uint l = 0; l < m_block_line_positions->size(); ++l)
  {
    common::ArrayView<const IdxType, _1D, Uint> line =
        m_sparse_store->template const_block<BLOCK_ID>(l);

    for (Uint i = 0; i < line.size(); ++i)
    {
      const uint32_t val = static_cast<uint32_t>(line[i]);
      out.write(reinterpret_cast<const char *>(&val), sizeof(uint32_t));
    }
  }

  // ----------------------------------------
  // Write cell data - row index
  // ----------------------------------------

  const uint64_t num_row_cell_data_bytes = sizeof(uint32_t) * nb_blocks;
  out.write(reinterpret_cast<const char *>(&num_row_cell_data_bytes),
            sizeof(num_row_cell_data_bytes));

  for (Uint l = 0; l < m_block_line_positions->size(); ++l)
  {
    const common::Range1D<IdxType> row_span = (*m_block_line_positions)[l];

    common::ArrayView<const common::Range1D<IdxType>, _1D, Uint> line =
        m_sparse_store->template const_block<BLOCK_LIMITS>(l);
    const uint32_t val = static_cast<uint32_t>(row_span.lbound());

    for (Uint i = 0; i < line.size(); ++i)
    {
      out.write(reinterpret_cast<const char *>(&val), sizeof(uint32_t));
    }
  }

  // ----------------------------------------
  // Write cell data - column index
  // ----------------------------------------

  const uint64_t num_col_cell_data_bytes = sizeof(uint32_t) * nb_blocks;
  out.write(reinterpret_cast<const char *>(&num_col_cell_data_bytes),
            sizeof(num_col_cell_data_bytes));

  for (Uint l = 0; l < m_block_line_positions->size(); ++l)
  {
    common::ArrayView<const common::Range1D<IdxType>, _1D, Uint> line =
        m_sparse_store->template const_block<BLOCK_LIMITS>(l);

    for (Uint i = 0; i < line.size(); ++i)
    {
      const common::Range1D<IdxType> col_span = line[i];
      const uint32_t val                      = static_cast<uint32_t>(col_span.lbound());
      out.write(reinterpret_cast<const char *>(&val), sizeof(uint32_t));
    }
  }

  // ----------------------------------------
  // Write nodal coordinates
  // ----------------------------------------
  const uint64_t num_bytes_coord = 3 * sizeof(double) * nb_block_coords;
  out.write(reinterpret_cast<const char *>(&num_bytes_coord), sizeof(num_bytes_coord));
  double coord_buffer[3];

  coord_buffer[X2] = 0.0;

  for (Uint l = 0; l < m_block_line_positions->size(); ++l)
  {
    const common::Range1D<IdxType> row_span = (*m_block_line_positions)[l];

    common::ArrayView<const common::Range1D<IdxType>, _1D, Uint> line =
        m_sparse_store->template const_block<BLOCK_LIMITS>(l);
    for (Uint i = 0; i < line.size(); ++i)
    {
      const common::Range1D<IdxType> col_span = line[i];

      // SW corner
      coord_buffer[X0] = col_span.lbound();
      coord_buffer[X1] = m_shape.size(0) - row_span.ubound() - 1;
      out.write(reinterpret_cast<const char *>(&coord_buffer[0]), 3 * sizeof(double));

      // SE corner
      coord_buffer[X0] = col_span.ubound();
      coord_buffer[X1] = m_shape.size(0) - row_span.ubound() - 1;
      out.write(reinterpret_cast<const char *>(&coord_buffer[0]), 3 * sizeof(double));

      // NE corner
      coord_buffer[X0] = col_span.ubound();
      coord_buffer[X1] = m_shape.size(0) - row_span.lbound() - 1;
      out.write(reinterpret_cast<const char *>(&coord_buffer[0]), 3 * sizeof(double));

      // NW corner
      coord_buffer[X0] = col_span.lbound();
      coord_buffer[X1] = m_shape.size(0) - row_span.lbound() - 1;
      out.write(reinterpret_cast<const char *>(&coord_buffer[0]), 3 * sizeof(double));
    }
  }

  // ----------------------------------------
  // Write connectivity
  // ----------------------------------------

  int64_t cell_id_buffer[4];

  const int64_t num_bytes_connectivity = 4 * sizeof(int64_t) * nb_blocks;
  out.write(reinterpret_cast<const char *>(&num_bytes_connectivity),
            sizeof(num_bytes_connectivity));
  for (Uint b = 0; b < nb_blocks; ++b)
  {
    cell_id_buffer[0] = 4 * b;
    cell_id_buffer[1] = 4 * b + 1;
    cell_id_buffer[2] = 4 * b + 2;
    cell_id_buffer[3] = 4 * b + 3;
    out.write(reinterpret_cast<const char *>(&cell_id_buffer[0]), 4 * sizeof(int64_t));
  }

  // ----------------------------------------
  // Write cell offsets
  // ----------------------------------------
  const int64_t num_bytes_cell_offsets = sizeof(int64_t) * nb_blocks;
  out.write(reinterpret_cast<const char *>(&num_bytes_cell_offsets),
            sizeof(num_bytes_cell_offsets));
  for (size_t b = 0; b < nb_blocks; ++b)
  {
    const int64_t cell_offset = 4 * (b + 1);
    out.write(reinterpret_cast<const char *>(&cell_offset), sizeof(cell_offset));
  }

  // ----------------------------------------
  // Write cell types
  // ----------------------------------------
  const int64_t num_bytes_cell_types = sizeof(uint8_t) * nb_blocks;
  out.write(reinterpret_cast<const char *>(&num_bytes_cell_types), sizeof(num_bytes_cell_types));
  for (size_t b = 0; b < nb_blocks; ++b)
  {
    const uint8_t cell_type = 9;
    out.write(reinterpret_cast<const char *>(&cell_type), sizeof(cell_type));
  }

  out << std::endl;

  out << "  </AppendedData>" << std::endl;
  out << "</VTKFile>" << std::endl;
  out.close();

  out.close();
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
bool BlockMatrixSparsityPattern<IdxType, SO>::CompareBlocks::operator()(
    const common::Range1D<IdxType> &a, const common::Range1D<IdxType> &b) const
{
  // This assumes that (a.lbound() < a.ubound()) && (b.lbound() < b.ubound())
  return a.ubound() < b.lbound();
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
inline const std::tuple<Uint, Uint> BlockMatrixSparsityPattern<IdxType, SO>::
    CompressedRowStorageIndexer::idx_position(
        const std::vector<common::Range1D<IdxType>> &line_positions,
        const common::BlockMultiArray<common::Range1D<IdxType>, IdxType> &line_blocks,
        const IdxType row, const IdxType col, const common::ArrayShape<_2D, Uint> &mat_shape)
{
  /// @TODO: use binary search instead of linear
  /// Note however, that the number of blocks in each block-sparse line is
  /// likely to be small, because in a sparse matrix representing a FE system
  /// the number of blocks on line 'i' is equal to the number of immediate
  /// neighbours of element 'i' + 1 (where +1 is block representing the i-th
  /// element itself)

  std::tuple<Uint, Uint> result(line_positions.size() + 1, line_positions.size() + 1);

  for (Uint idx0 = 0; idx0 < line_positions.size(); ++idx0)
  {
    const common::Range1D<IdxType> &line_pos = line_positions[idx0];
    if (line_pos.in_range(row))
    {
      std::get<0>(result) = idx0;
      break;
    }
  }

  if (std::get<0>(result) < line_positions.size())
  {
    const common::ArrayView<const common::Range1D<IdxType>, _1D, Uint> &line =
        line_blocks.template const_block<BLOCK_LIMITS>(std::get<0>(result));

    for (Uint idx1 = 0; idx1 < line.size(); ++idx1)
    {
      if (line[idx1].in_range(col))
      {
        const common::ArrayView<const IdxType, _1D, Uint> &block_ids =
            line_blocks.template const_block<BLOCK_ID>(std::get<0>(result));

        std::get<1>(result) = block_ids[idx1];
        return result;
      }
    }
  }

  return result;
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
inline const std::tuple<Uint, Uint> BlockMatrixSparsityPattern<IdxType, SO>::
    CompressedColStorageIndexer::idx_position(
        const std::vector<common::Range1D<IdxType>> &line_positions,
        const common::BlockMultiArray<common::Range1D<IdxType>, IdxType> &line_blocks,
        const IdxType row, const IdxType col, const common::ArrayShape<_2D, Uint> &mat_shape)
{
  /// @TODO: use binary search instead of linear
  /// Note however, that the number of blocks in each block-sparse line is
  /// likely to be small, because in a sparse matrix representing a FE system
  /// the number of blocks on line 'i' is equal to the number of immediate
  /// neighbours of element 'i' + 1 (where +1 is block representing the i-th
  /// element itself)

  std::tuple<Uint, Uint> result(line_positions.size() + 1, line_positions.size() + 1);

  for (Uint idx1 = 0; idx1 < line_positions.size(); ++idx1)
  {
    const common::Range1D<IdxType> &line_pos = line_positions[idx1];
    if (line_pos.in_range(col))
    {
      std::get<1>(result) = idx1;
      break;
    }
  }

  if (std::get<1>(result) < line_positions.size())
  {
    const common::ArrayView<const common::Range1D<IdxType>, _1D, Uint> &line =
        line_blocks.const_block<BLOCK_LIMITS>(std::get<1>(result));

    for (Uint idx0 = 0; idx0 < line.size(); ++idx0)
    {
      if (line[idx0].in_range(row))
      {
        const common::ArrayView<const IdxType, _1D, Uint> &block_ids =
            line_blocks.template const_block<BLOCK_ID>(std::get<1>(result));

        std::get<0>(result) = block_ids[idx0];
        return result;
      }
    }
  }

  return result;
}

// -----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif
