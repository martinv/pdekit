#ifndef PDEKIT_Math_Matrix_Sparsity_Pattern_hpp
#define PDEKIT_Math_Matrix_Sparsity_Pattern_hpp

#include <cmath>
#include <iostream>

#include "common/ArrayShape.hpp"
#include "common/BlockArray.hpp"
#include "common/Meta.hpp"
#include "math/Matrix.hpp"
#include "math/MatrixStorageOrder.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

// A structure holding POD and pointers
// to information stored in sparsity pattern
// Needed by external libraries

template <typename IdxType>
struct MatrixSparsityPatternRawData
{
  // Pointer to block offsets
  const IdxType *lines;

  // Pointer to nonzero values
  const IdxType *nonzeros;

  // Number of nonzeros
  Uint nnz;

  // Number of lines (i.e rows or columns)
  Uint nb_lines;
};

// ----------------------------------------------------------------------------

template <typename IdxType, bool SO = DefaultMatrixStorageOrder>
class MatrixSparsityPattern
{
  /// TYPEDEFS

  public:
  using idx_type         = Uint;
  using sparse_line_type = common::ArrayView<const IdxType, _1D, IdxType>;

  private:
  using sparsity_storage_type = common::BlockArray<IdxType, IdxType>;

  public:
  /// Constructor, takes the number of rows and columns
  explicit MatrixSparsityPattern(Uint m = 0, Uint n = 0);

  /// The copy constructor is explicitly defined due to the required
  /// dynamic memory management and in order to enable/facilitate NRV
  /// optimization.
  MatrixSparsityPattern(const MatrixSparsityPattern &rhs);

  template <typename OtherIdxType>
  explicit MatrixSparsityPattern(const MatrixSparsityPattern<OtherIdxType, SO> &init);

  /// Destructor
  ~MatrixSparsityPattern();

  /// Assignment operator
  MatrixSparsityPattern &operator=(const MatrixSparsityPattern &rhs);

  /// Resize the matrix
  void resize(const Uint m, const Uint n);

  /// Build the sparsity pattern
  void build_sparsity(const std::vector<std::vector<IdxType>> &nonzero_coords);

  /// Build the sparsity pattern
  void build_sparsity(std::unique_ptr<common::BlockArray<IdxType, IdxType>> &&nonzero_coords);

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
  Uint idx_position(const IdxType row, const IdxType col) const;

  /// Get one sparse line (which can be row or column, depending on the
  /// storage order
  const sparse_line_type sparse_line(const IdxType line_id) const;

  /// Get one sparse line (which can be row or column, depending on the
  /// storage order
  /// @return a tuple containing the first and last index (included) where
  /// line data is stored
  const common::Range1D<IdxType> sparse_line_limits(const IdxType line_id) const;

  /// Get raw pointers to all data
  /// Needed by some external libraries
  const MatrixSparsityPatternRawData<IdxType> raw_data() const;

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

  // -----------------------------------------------------------------------------
  // A helper structure which treats MatrixSparsityPattern as if it
  // represented compressed ROW storage: 'line' is a matrix row and 'index in
  // line' represents a column index in that row
  // -----------------------------------------------------------------------------

  struct CompressedRowStorageIndexer
  {
    static Uint idx_position(const common::BlockArray<IdxType, IdxType> &spattern,
                             const IdxType row, const IdxType col);
  };

  // -----------------------------------------------------------------------------
  // A helper structure which treats MatrixSparsityPattern as if it
  // represented compressed COLUMN storage: 'line' is a matrix column and
  // 'index in line' represents a row index in that column
  // -----------------------------------------------------------------------------

  struct CompressedColStorageIndexer
  {
    static Uint idx_position(const common::BlockArray<IdxType, IdxType> &spattern,
                             const IdxType row, const IdxType col);
  };

  // Set the type indexing either to row-major or column-major
  using indexer_type = typename common::SelectType<SO == RowMajor, CompressedRowStorageIndexer,
                                                   CompressedColStorageIndexer>::type;

  /// DATA
  /// Number of matrix rows and columns
  common::ArrayShape<_2D, Uint> m_shape;

  /// Sparsity pattern
  std::unique_ptr<sparsity_storage_type> m_sparse_store;
};

// ----------------------------------------------------------------------------

template <typename IdxType, bool SO>
MatrixSparsityPattern<IdxType, SO>::MatrixSparsityPattern(Uint m, Uint n)
    : m_shape(m, n), m_sparse_store(new sparsity_storage_type)
{
}

// ----------------------------------------------------------------------------

template <typename IdxType, bool SO>
MatrixSparsityPattern<IdxType, SO>::MatrixSparsityPattern(const MatrixSparsityPattern &rhs)
    : m_shape(rhs.m_shape)
{
  m_sparse_store->resize(rhs->m_sparse_store.size(), rhs->m_sparse_store.nb_blocks());
  (*m_sparse_store) = (*rhs.m_sparse_store);
}

// ----------------------------------------------------------------------------

template <typename IdxType, bool SO>
template <typename OtherIdxType>
MatrixSparsityPattern<IdxType, SO>::MatrixSparsityPattern(
    const MatrixSparsityPattern<OtherIdxType, SO> &init)
    : m_shape(init.m_shape)
{
  m_sparse_store->resize(init->m_sparse_store.size(), init->m_sparse_store.nb_blocks());
  (*m_sparse_store) = (*init.m_sparse_store);
}

// ----------------------------------------------------------------------------

template <typename IdxType, bool SO>
MatrixSparsityPattern<IdxType, SO>::~MatrixSparsityPattern()
{
}

// ----------------------------------------------------------------------------

template <typename IdxType, bool SO>
MatrixSparsityPattern<IdxType, SO> &MatrixSparsityPattern<IdxType, SO>::operator=(
    const MatrixSparsityPattern<IdxType, SO> &rhs)
{
  m_shape = rhs.m_shape;
  m_sparse_store->resize(rhs->m_sparse_store.size(), rhs->m_sparse_store.nb_blocks());
  (*m_sparse_store) = (*rhs.m_sparse_store);

  return *this;
}

// ----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void MatrixSparsityPattern<IdxType, SO>::resize(const Uint m, const Uint n)
{
  m_shape = common::ArrayShape<_2D, Uint>(m, n);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void MatrixSparsityPattern<IdxType, SO>::build_sparsity(
    const std::vector<std::vector<IdxType>> &nonzero_coords)
{
  Uint nb_entries = 0;
  for (std::vector<IdxType> const &line : nonzero_coords)
  {
    nb_entries += line.size();
  }

  m_sparse_store->reserve(nb_entries, nonzero_coords.size());
  // m_sparsity.resize(0, 0);

  IdxType max_entry = IdxType();
  for (Uint b = 0; b < nonzero_coords.size(); ++b)
  {
    const std::vector<IdxType> &line = nonzero_coords[b];
    const common::ArrayView<const IdxType, _1D, IdxType> block_view(line.data(), line.size());
    m_sparse_store->create_back_block(block_view.size());
    m_sparse_store->fill_last_block(block_view);

    for (Uint i = 0; i < block_view.size(); ++i)
    {
      max_entry = std::max(max_entry, block_view[i]);
    }
  }

  m_shape = (SO == RowMajor) ? common::ArrayShape<_2D, Uint>(nonzero_coords.size(), max_entry + 1)
                             : common::ArrayShape<_2D, Uint>(max_entry + 1, nonzero_coords.size());
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void MatrixSparsityPattern<IdxType, SO>::build_sparsity(
    std::unique_ptr<common::BlockArray<IdxType, IdxType>> &&nonzero_coords)
{
  m_sparse_store = std::move(nonzero_coords);

  IdxType max_entry = IdxType();
  for (Uint b = 0; b < m_sparse_store->nb_blocks(); ++b)
  {
    const common::ArrayView<const IdxType, _1D, IdxType> block_view =
        m_sparse_store->const_block(b);

    for (Uint i = 0; i < block_view.size(); ++i)
    {
      max_entry = std::max(max_entry, block_view[i]);
    }
  }

  m_shape = (SO == RowMajor)
                ? common::ArrayShape<_2D, Uint>(m_sparse_store->nb_blocks(), max_entry + 1)
                : common::ArrayShape<_2D, Uint>(max_entry + 1, m_sparse_store->nb_blocks());
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
const common::ArrayShape<_2D, Uint> MatrixSparsityPattern<IdxType, SO>::size() const
{
  return m_shape;
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
Uint MatrixSparsityPattern<IdxType, SO>::nb_nz() const
{
  return m_sparse_store->size();
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
Uint MatrixSparsityPattern<IdxType, SO>::nb_lines() const
{
  return m_sparse_store->nb_blocks();
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
inline Uint MatrixSparsityPattern<IdxType, SO>::idx_position(const IdxType row,
                                                             const IdxType col) const
{
  return indexer_type::idx_position(*m_sparse_store, row, col);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
const typename MatrixSparsityPattern<IdxType, SO>::sparse_line_type MatrixSparsityPattern<
    IdxType, SO>::sparse_line(const IdxType line_id) const
{
  return m_sparse_store->const_block(line_id);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
inline const common::Range1D<IdxType> MatrixSparsityPattern<IdxType, SO>::sparse_line_limits(
    const IdxType line_id) const
{
  return m_sparse_store->block_limits(line_id);
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void MatrixSparsityPattern<IdxType, SO>::print_sparsity() const
{
  std::cout << "Matrix sparsity pattern (" << m_shape.size(0) << "," << m_shape.size(1) << ")"
            << std::endl;
  for (Uint b = 0; b < m_sparse_store->nb_blocks(); ++b)
  {
    std::cout << m_sparse_store->const_block(b) << std::endl;
  }
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void MatrixSparsityPattern<IdxType, SO>::print_svg(const std::string &filename) const
{
  const Uint m = m_shape.size(0);
  const Uint n = m_shape.size(1);

  const Real border = 0.5;

  std::fstream out;
  out.open(filename.c_str());

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
    common::ArrayView<const IdxType, _1D, IdxType> line = m_sparse_store->const_block(b);
    for (Uint i = 0; i < line.size(); ++i)
    {
      // The 0.05 offset is here to center the actual nonzero entries are
      // squares of size [0.9 x 0.9] . This means that the right and
      // bottom edge of the plot have a gap of width 1.0 - 0.9 = 0.1.
      // Therefore we divide this gap by 2 and offset all squares
      // representing nonzero entries by 0.1/2 = 0.05
      out << R"(  <rect class="pixel" x=")" << line[i] + 0.5 * border + 0.05 << "\" y=\""
          << b + 0.5 * border + 0.05 << "\" width=\".9\" height=\".9\"/>\n";
    }
  }
  out << "</svg>" << std::endl;
  out.close();
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
void MatrixSparsityPattern<IdxType, SO>::print_vtu(const std::string &filename) const
{
  std::ofstream out;
  out.open(filename.c_str(), std::ios::binary);
  size_t appended_data_offset = 0;

  const Uint nb_nonzeros = m_sparse_store->size();

  // out << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
         "byte_order=\"LittleEndian\" "
         "header_type=\"UInt64\">"
      << std::endl;
  out << "  <UnstructuredGrid>" << std::endl;
  // Number of points and cells is equal to number of nonzeros in sparsity
  // pattern
  out << "    <Piece NumberOfPoints=\"" << nb_nonzeros << "\" NumberOfCells=\"" << nb_nonzeros
      << "\">" << std::endl;

  // ---

  out << "      <PointData Scalars=\"nonzeros\">" << std::endl;
  out << "        <DataArray type=\"UInt32\" Name=\"nonzeros\" "
         "format=\"appended\" offset=\""
      << appended_data_offset << "\"/>" << std::endl;
  out << "      </PointData>" << std::endl;

  appended_data_offset += nb_nonzeros * sizeof(uint32_t) + sizeof(uint64_t);

  // ---

  out << "      <CellData>" << std::endl;
  out << "      </CellData>" << std::endl;

  out << "      <Points>" << std::endl;
  out << "        <DataArray type=\"Float64\" Name=\"Nodes\" "
         "NumberOfComponents=\"3\" "
         "format=\"appended\" "
         "offset=\""
      << appended_data_offset << "\"/>" << std::endl;
  out << "      </Points>" << std::endl;

  appended_data_offset += 3 * sizeof(Real) * nb_nonzeros + sizeof(uint64_t);
  std::cout << "appended data offset = " << appended_data_offset << std::endl;

  out << "      <Cells>" << std::endl;

  // Cell connectivity
  out << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
         "format=\"appended\" offset=\""
      << appended_data_offset << "\"/>" << std::endl;
  appended_data_offset += sizeof(int64_t) * nb_nonzeros + sizeof(uint64_t);

  // Cell offset data
  out << "        <DataArray type=\"Int64\" Name=\"offsets\" "
         "format=\"appended\" offset=\""
      << appended_data_offset << "\"/>" << std::endl;
  appended_data_offset += sizeof(int64_t) * nb_nonzeros + sizeof(uint64_t);

  // Cell type data
  out << "        <DataArray type=\"UInt8\" Name=\"types\" "
         "format=\"appended\" "
         "offset=\""
      << appended_data_offset << "\"/>" << std::endl;
  appended_data_offset += sizeof(uint8_t) * nb_nonzeros + sizeof(uint64_t);

  out << "      </Cells>" << std::endl;
  out << "    </Piece>" << std::endl;
  out << "  </UnstructuredGrid>" << std::endl;
  out << "  <AppendedData encoding=\"raw\">" << std::endl;

  const std::string begin_marker = "                _";
  out << begin_marker;

  // ----------------------------------------
  // Write nodal data - nonzero values
  // ----------------------------------------
  const uint64_t num_bytes = sizeof(uint32_t) * nb_nonzeros;
  out.write(reinterpret_cast<const char *>(&num_bytes), sizeof(num_bytes));

  for (Uint b = 0; b < m_sparse_store->nb_blocks(); ++b)
  {
    common::ArrayView<const IdxType, _1D, IdxType> line = m_sparse_store->const_block(b);
    for (Uint i = 0; i < line.size(); ++i)
    {
      const uint32_t val = static_cast<uint32_t>(b + line[i]);
      out.write(reinterpret_cast<const char *>(&val), sizeof(uint32_t));
    }
  }

  // ----------------------------------------
  // Write nodal coordinates
  // ----------------------------------------
  const uint64_t num_bytes_coord = 3 * sizeof(double) * nb_nonzeros;
  out.write(reinterpret_cast<const char *>(&num_bytes_coord), sizeof(num_bytes_coord));
  double coord_buffer[3];

  coord_buffer[X2] = 0.0;
  for (Uint b = 0; b < m_sparse_store->nb_blocks(); ++b)
  {
    coord_buffer[X1]                                    = m_shape.size(0) - b - 1.0;
    common::ArrayView<const IdxType, _1D, IdxType> line = m_sparse_store->const_block(b);
    for (Uint i = 0; i < line.size(); ++i)
    {
      coord_buffer[X0] = line[i];
      out.write(reinterpret_cast<const char *>(&coord_buffer[0]), 3 * sizeof(double));
    }
  }

  // ----------------------------------------
  // Write connectivity
  // ----------------------------------------
  const int64_t num_bytes_connectivity = sizeof(int64_t) * nb_nonzeros;
  out.write(reinterpret_cast<const char *>(&num_bytes_connectivity),
            sizeof(num_bytes_connectivity));
  for (Uint i = 0; i < nb_nonzeros; ++i)
  {
    const int64_t cell_id = i + 1;
    out.write(reinterpret_cast<const char *>(&cell_id), sizeof(cell_id));
  }

  // ----------------------------------------
  // Write cell offsets
  // ----------------------------------------
  const int64_t num_bytes_cell_offsets = sizeof(int64_t) * nb_nonzeros;
  out.write(reinterpret_cast<const char *>(&num_bytes_cell_offsets),
            sizeof(num_bytes_cell_offsets));
  for (size_t i = 0; i < nb_nonzeros; ++i)
  {
    const int64_t cell_offset = i + 1;
    out.write(reinterpret_cast<const char *>(&cell_offset), sizeof(cell_offset));
  }

  // ----------------------------------------
  // Write cell types
  // ----------------------------------------
  const int64_t num_bytes_cell_types = sizeof(uint8_t) * nb_nonzeros;
  out.write(reinterpret_cast<const char *>(&num_bytes_cell_types), sizeof(num_bytes_cell_types));
  for (size_t i = 0; i < nb_nonzeros; ++i)
  {
    const uint8_t cell_type = 1;
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
const MatrixSparsityPatternRawData<IdxType> MatrixSparsityPattern<IdxType, SO>::raw_data() const
{
  std::tuple<const idx_type *, const idx_type *> sparse_store_raw_data = m_sparse_store->raw_data();

  MatrixSparsityPatternRawData<IdxType> raw_data;
  raw_data.lines    = std::get<0>(sparse_store_raw_data);
  raw_data.nonzeros = std::get<1>(sparse_store_raw_data);
  raw_data.nnz      = m_sparse_store->size();
  raw_data.nb_lines = m_sparse_store->nb_blocks();
  return raw_data;
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
Real MatrixSparsityPattern<IdxType, SO>::mem_used_mb() const
{
  Real mem_used = 0.0;

  if (m_sparse_store)
  {
    mem_used += m_sparse_store->mem_used_mb();
  }

  return mem_used;
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
inline Uint MatrixSparsityPattern<IdxType, SO>::CompressedRowStorageIndexer::idx_position(
    const common::BlockArray<IdxType, IdxType> &spattern, const IdxType row, const IdxType col)
{
  common::ArrayView<const IdxType, _1D, IdxType> row_view = spattern.const_block(row);

  const IdxType *row_begin      = &row_view[0];
  const IdxType *row_end        = row_begin + row_view.size();
  const IdxType *search_idx_pos = std::lower_bound(row_begin, row_end, col);

  if ((search_idx_pos != row_end) && !(col < *search_idx_pos))
  {
    return spattern.distance_from_origin(search_idx_pos);
  }
  return spattern.size();
}

// -----------------------------------------------------------------------------

template <typename IdxType, bool SO>
inline Uint MatrixSparsityPattern<IdxType, SO>::CompressedColStorageIndexer::idx_position(
    const common::BlockArray<IdxType, IdxType> &spattern, const IdxType row, const IdxType col)
{
  common::ArrayView<const IdxType, _1D, IdxType> col_view = spattern.const_block(col);

  const IdxType *col_begin      = &col_view[0];
  const IdxType *col_end        = col_begin + col_view.size();
  const IdxType *search_idx_pos = std::lower_bound(col_begin, col_end, row);

  if ((search_idx_pos != col_end) && !(row < *search_idx_pos))
  {
    return spattern.distance_from_origin(search_idx_pos);
  }
  return spattern.size();
}

// -----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif
