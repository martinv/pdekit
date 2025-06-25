#ifndef PDEKIT_Math_Sparse_Dynamic_Matrix_hpp
#define PDEKIT_Math_Sparse_Dynamic_Matrix_hpp

#include <cmath>
#include <iostream>

#include "common/BlockArray.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseVecView.hpp"
#include "math/MathForward.hpp"
#include "math/Matrix.hpp"
#include "math/MatrixSparsityPattern.hpp"
#include "math/MatrixStorageOrder.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

/// A structure holding POD and pointers
/// to information stored in sparse matrix
/// Needed by external libraries

template <typename T, typename IdxType>
struct SparseDMatRawData
{
  /// Default constructor
  SparseDMatRawData() : lines(nullptr), nonzeros(nullptr), nnz(0), nb_lines(0), values(nullptr)
  {
  }

  /// Construct from raw data of sparsity pattern
  /// Pointer to matrix values won't be correctly set at this point
  SparseDMatRawData(const MatrixSparsityPatternRawData<IdxType> &sparsity_raw_data)
      : lines(sparsity_raw_data.lines), nonzeros(sparsity_raw_data.nonzeros),
        nnz(sparsity_raw_data.nnz), nb_lines(sparsity_raw_data.nb_lines), values(nullptr)
  {
  }

  /// Pointer to block offsets
  const IdxType *lines;

  /// Pointer to nonzero values
  const IdxType *nonzeros;

  /// Number of nonzeros
  Uint nnz;

  /// Number of lines (i.e rows or columns)
  Uint nb_lines;

  /// Pointer to the beginning of daata
  const T *values;
};

// ----------------------------------------------------------------------------

template <typename T, bool SO = DefaultMatrixStorageOrder>
class SparseDMat : public SparseMatrix<SparseDMat<T, SO>, SO>
{
  /// TYPEDEFS

  public:
  using idx_type       = Uint;
  using tensor_type    = SparseDMat<T, SO>;
  using composite_type = SparseDMat<T, SO> const &;
  using value_type     = T;

  private:
  using sparsity_pattern_type = math::MatrixSparsityPattern<idx_type, SO>;
  using sparse_line_type      = typename sparsity_pattern_type::sparse_line_type;

  public:
  enum
  {
    is_expression = 0
  };
  enum
  {
    evaluates_fast = 1
  };
  enum
  {
    owns_data = 1
  };

  /// Default constructor
  explicit SparseDMat();

  /// Constructor, takes the number of rows and columns
  explicit SparseDMat(Uint m, Uint n);

  /// The copy constructor is explicitly defined due to the required
  /// dynamic memory management and in order to enable/facilitate NRV
  /// optimization.
  SparseDMat(const SparseDMat &rhs);

  template <typename MT>
  explicit SparseDMat(const SparseMatrix<MT, SO> &init);

  /// Intentionally disabled - SparseDMat cannot be constructed
  /// from initializer list, because the shape of the matrix would be
  /// undefined!
  /// explicit SparseDMat( const TensorInitializer<T>& init );

  /// Destructor
  ~SparseDMat();

  /// Resize the matrix
  void resize(const Uint m, const Uint n);

  /// Build the sparsity pattern
  void build_sparsity(const std::vector<std::vector<idx_type>> &nonzero_coords);

  /// Fill data
  void fill_data(const math::DenseConstVecView<T> &data);

  /// Indexing operator
  /// @param the row and column index of the element to be returned
  T &operator()(const idx_type i, idx_type j);

  /// Indexing operator, const version
  const T &operator()(const idx_type i, const idx_type j) const;

  /// Return line data, const version
  const std::tuple<const sparse_line_type, const math::DenseConstVecView<T>> const_line_data(
      const idx_type line_id) const;

  /// Return line data
  const std::tuple<const sparse_line_type, math::DenseVecView<T>> line_data(const idx_type line_id);

  /// Return the sparsity pattern dimensions
  const common::ArrayShape<_2D, Uint> size() const;

  /// Return number of nonzero values
  const Uint nb_nz() const;

  /// Get raw pointers to all data
  /// Needed by some external libraries
  const SparseDMatRawData<T, idx_type> raw_data() const;

  /// Print the sparsity pattern
  void print() const;

  private:
  /// DATA
  /// Sparsity pattern of the matrix - storage of nonzero entries
  math::MatrixSparsityPattern<idx_type, SO> m_sparsity;

  // Entries in the matrix
  T *m_data;
};

// ----------------------------------------------------------------------------

template <typename T, bool SO>
SparseDMat<T, SO>::SparseDMat()
    : SparseMatrix<SparseDMat<T, SO>, SO>(*this), m_sparsity(0, 0), m_data(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
SparseDMat<T, SO>::SparseDMat(Uint m, Uint n)
    : SparseMatrix<SparseDMat<T, SO>, SO>(*this), m_sparsity(m, n), m_data(nullptr)
{
  // resize(m,n);
  /*
  if (m * n != 0u)
  {
    m_data = new T[m * n];
    m_rows = m;
    m_cols = n;
  }
  else
  {
    m_data = nullptr;
    m_rows = 0u;
    m_cols = 0u;
  }
  */
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
SparseDMat<T, SO>::SparseDMat(const SparseDMat &rhs)
{
  m_sparsity = rhs.m_sparsity;

  const common::ArrayShape<_2D, Uint> shape = rhs.size();

  const Uint nb_nnz = m_sparsity.nb_nz();
  m_data            = new T[nb_nnz];

  for (Uint i = 0; i < nb_nnz; ++i)
  {
    m_data[i] = T();
  }
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
template <typename MT>
SparseDMat<T, SO>::SparseDMat(const SparseMatrix<MT, SO> &init)
{
  /*
  const MT &rhs_expr = init.wrapped_type();

  // resize(rhs_expr.rows(),rhs_expr.cols());

  m_data = new T[rhs_expr.rows() * rhs_expr.cols()];
  m_rows = rhs_expr.rows();
  m_cols = rhs_expr.cols();

  free_assign(*this, init.wrapped_type());
  */
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
SparseDMat<T, SO>::~SparseDMat()
{
  delete[] m_data;
  m_data = nullptr;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void SparseDMat<T, SO>::resize(const Uint m, const Uint n)
{
  m_sparsity.resize(m, n);
  /*
  if (m * n == 0)
  {
    if (m_data)
    {
      delete[] m_data;
      m_data = nullptr;
    }
    m_rows = 0;
    m_cols = 0;
    return;
  }

  if (m_rows * m_cols != m * n)
  {
    if (m_data != nullptr)
    {
      delete[] m_data;
    }
    m_data = new T[m * n];
  }
  m_rows = m;
  m_cols = n;

  fill(T());
  */
}

// -----------------------------------------------------------------------------

template <typename T, bool SO>
void SparseDMat<T, SO>::build_sparsity(const std::vector<std::vector<idx_type>> &nonzero_coords)
{
  const Uint old_nnz = m_sparsity.nb_nz();

  m_sparsity.build_sparsity(nonzero_coords);

  const Uint new_nnz = m_sparsity.nb_nz();

  if (old_nnz != new_nnz)
  {
    delete[] m_data;
    m_data = new T[new_nnz];
  }

  for (Uint i = 0; i < new_nnz; ++i)
  {
    m_data[i] = T();
  }
}

// -----------------------------------------------------------------------------

template <typename T, bool SO>
void SparseDMat<T, SO>::fill_data(const math::DenseConstVecView<T> &data)
{
  if (data.size() != m_sparsity.nb_nz())
  {
    std::cerr << "SparseDMat::fill_data: the input data has wrong size " << data.size()
              << ". I expected " << m_sparsity.nb_nz() << " instead. Not filling data."
              << std::endl;
    return;
  }

  if (m_sparsity.nb_nz() == 0)
  {
    std::cerr << "SparseDMat::fill_data: I don't have any nonzero entries. "
                 "Nothing to fill."
              << std::endl;
    return;
  }

  if (!m_data)
  {
    m_data = new T[m_sparsity.nb_nz()];
  }

  for (Uint i = 0; i < data.size(); ++i)
  {
    m_data[i] = data[i];
  }
}

// -----------------------------------------------------------------------------

template <typename T, bool SO>
inline T &SparseDMat<T, SO>::operator()(const idx_type i, const idx_type j)
{
  return T();
}

// -----------------------------------------------------------------------------

template <typename T, bool SO>
inline const T &SparseDMat<T, SO>::operator()(const Uint i, const Uint j) const
{
  return T();
}

// -----------------------------------------------------------------------------

template <typename T, bool SO>
const std::tuple<const typename SparseDMat<T, SO>::sparse_line_type,
                 const math::DenseConstVecView<T>>
SparseDMat<T, SO>::const_line_data(const idx_type line_id) const
{
  const common::Range1D<idx_type> line_limits = m_sparsity.sparse_line_limits(line_id);
  const math::DenseConstVecView<T> values(m_data + line_limits.lbound(),
                                          line_limits.ubound() - line_limits.lbound() + 1);

  return std::make_tuple(m_sparsity.sparse_line(line_id), values);
}

// -----------------------------------------------------------------------------

template <typename T, bool SO>
const std::tuple<const typename SparseDMat<T, SO>::sparse_line_type, math::DenseVecView<T>> SparseDMat<
    T, SO>::line_data(const idx_type line_id)
{
  const common::Range1D<idx_type> line_limits = m_sparsity.sparse_line_limits(line_id);
  const math::DenseVecView<T> values(m_data + line_limits.lbound(),
                                     line_limits.ubound() - line_limits.lbound() + 1);

  return std::make_tuple(m_sparsity.sparse_line(line_id), values);
}

// -----------------------------------------------------------------------------

template <typename T, bool SO>
const common::ArrayShape<_2D, Uint> SparseDMat<T, SO>::size() const
{
  return m_sparsity.size();
}

// -----------------------------------------------------------------------------

template <typename T, bool SO>
const Uint SparseDMat<T, SO>::nb_nz() const
{
  return m_sparsity.nb_nz();
}

// -----------------------------------------------------------------------------

template <typename T, bool SO>
const SparseDMatRawData<T, typename SparseDMat<T, SO>::idx_type> SparseDMat<T, SO>::raw_data() const
{
  // const std::tuple<idx_type*, idx_type*> sparsity_raw_data =
  // m_sparsity.raw_data();

  SparseDMatRawData<T, idx_type> raw_data(m_sparsity.raw_data());
  raw_data.values = m_data;
  return raw_data;
}

// -----------------------------------------------------------------------------

template <typename T, bool SO>
void SparseDMat<T, SO>::print() const
{
  const common::ArrayShape<_2D, Uint> mat_shape = m_sparsity.size();
  std::cout << "Sparse matrix (" << mat_shape.size(0) << "," << mat_shape.size(0) << ")"
            << std::endl;

  m_sparsity.print_sparsity();

  if (m_data && (m_sparsity.nb_nz() > 0))
  {
    std::cout << "Sparse matrix values:" << std::endl;
    for (Uint i = 0; i < m_sparsity.nb_nz(); ++i)
    {
      std::cout << m_data[i] << " ";
    }
    std::cout << std::endl;
  }
  else
  {
    std::cout << "Sparse matrix does not have any values stored." << std::endl;
  }
}

// -----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif
