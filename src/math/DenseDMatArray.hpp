#ifndef PDEKIT_Math_Dense_Dynamic_Matrix_Array_hpp
#define PDEKIT_Math_Dense_Dynamic_Matrix_Array_hpp

#include "common/ArrayShape.hpp"
#include "common/BlockArray.hpp"
#include "math/DenseConstMatView.hpp"
#include "math/DenseMatView.hpp"

namespace pdekit
{

namespace math
{

template <typename T, bool SO = DefaultMatrixStorageOrder>
class DenseDMatArray
{
  public:
  /// Default constructor
  DenseDMatArray();

  /// Copy constructor
  DenseDMatArray(const DenseDMatArray &other);

  /// Destructor
  ~DenseDMatArray();

  /// Assignment operator
  DenseDMatArray &operator=(const DenseDMatArray &rhs);

  /// Allocate matrix shapes, don't fill matrices
  void allocate(std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> &&mat_shapes);

  /// Allocate matrix shapes and fill matrices
  void allocate(std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> &&mat_shapes,
                std::unique_ptr<common::BlockArray<T, Uint>> &&mat_values);

  /// Reshape the array so that it has the same shape as some other array
  void reshape(const DenseDMatArray &other);

  /// Get one matrix view, non-const version
  DenseMatView<T, SO> mat_view(const Uint mat_id);

  /// Get one matrix view, const version
  DenseConstMatView<T, SO> const_mat_view(const Uint mat_id) const;

  /// Return number of matrices in array
  Uint size() const;

  /// Print information for debugging
  void debug_print() const;

  private:
  /// Local (element-wise) array of local interpolation matrices
  common::BlockArray<T, Uint> m_values;

  /// The shape of each matrix
  std::vector<common::ArrayShape<_2D, SUint>> m_mat_shapes;
};

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMatArray<T, SO>::DenseDMatArray()
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMatArray<T, SO>::DenseDMatArray(const DenseDMatArray &other)
    : m_values(other.m_values), m_mat_shapes(other.m_mat_shapes)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMatArray<T, SO>::~DenseDMatArray()
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMatArray<T, SO> &DenseDMatArray<T, SO>::operator=(const DenseDMatArray &rhs)
{
  m_values = rhs.m_values;

  m_mat_shapes.resize(rhs.m_mat_shapes.size());
  m_mat_shapes = rhs.m_mat_shapes;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void DenseDMatArray<T, SO>::allocate(
    std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> &&mat_shapes)
{
  if (!mat_shapes)
  {
    return;
  }

  std::unique_ptr<std::vector<T>> all_mat_values(new std::vector<T>());
  std::unique_ptr<std::vector<Uint>> all_mat_sizes(new std::vector<Uint>());

  const Uint num_blocks = mat_shapes->size();

  all_mat_sizes->resize(num_blocks);

  Uint num_values = 0;

  for (Uint i = 0; i < num_blocks; ++i)
  {
    const common::ArrayShape<_2D, SUint> one_mat_dims = (*mat_shapes)[i];
    const Uint mat_storage                            = one_mat_dims.size(0) * one_mat_dims.size(1);
    (*all_mat_sizes)[i]                               = mat_storage;
    num_values += mat_storage;
  }

  m_values.resize(num_values, num_blocks);

  all_mat_values->resize(num_values);
  all_mat_values->assign(num_values, T());

  m_values.build(std::move(all_mat_values), std::move(all_mat_sizes));

  std::swap(m_mat_shapes, *mat_shapes);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void DenseDMatArray<T, SO>::allocate(
    std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> &&mat_shapes,
    std::unique_ptr<common::BlockArray<T, Uint>> &&mat_values)
{
  m_values.swap(*mat_values);
  m_mat_shapes.swap(*mat_shapes);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void DenseDMatArray<T, SO>::reshape(const DenseDMatArray &other)
{
  bool shapes_are_same = (m_mat_shapes.size() == other.m_mat_shapes.size());

  if (shapes_are_same)
  {
    for (Uint i = 0; i < m_mat_shapes.size(); ++i)
    {
      if (m_mat_shapes[i] != other.m_mat_shapes[i])
      {
        shapes_are_same = false;
      }
    }
  }
  if (shapes_are_same)
    return;

  m_mat_shapes.resize(other.m_mat_shapes.size());
  std::vector<Uint> new_block_sizes(m_mat_shapes.size());

  for (Uint i = 0; i < m_mat_shapes.size(); ++i)
  {
    m_mat_shapes[i]    = other.m_mat_shapes[i];
    new_block_sizes[i] = m_mat_shapes[i].size(0) * m_mat_shapes[i].size(1);
  }

  m_values.resize_blocks(new_block_sizes);
  m_values.fill(T{});
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseMatView<T, SO> DenseDMatArray<T, SO>::mat_view(const Uint mat_id)
{
  common::ArrayView<T, _1D, Uint> lin_block = m_values.block(mat_id);
  DenseMatView<T, SO> view(&lin_block[0], m_mat_shapes[mat_id].size(1),
                           m_mat_shapes[mat_id].size(0), m_mat_shapes[mat_id].size(1));

  return view;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstMatView<T, SO> DenseDMatArray<T, SO>::const_mat_view(const Uint mat_id) const
{
  const common::ArrayView<const T, _1D, Uint> const_lin_block = m_values.const_block(mat_id);
  DenseConstMatView<T, SO> view(&const_lin_block[0], m_mat_shapes[mat_id].size(1),
                                m_mat_shapes[mat_id].size(0), m_mat_shapes[mat_id].size(1));
  return view;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
Uint DenseDMatArray<T, SO>::size() const
{
  return m_mat_shapes.size();
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void DenseDMatArray<T, SO>::debug_print() const
{
  std::cout << "DenseDMatArray::shapes:" << std::endl;
  for (auto shape : m_mat_shapes)
  {
    std::cout << " [" << shape.size(0) << "," << shape.size(1) << "]";
  }
  std::cout << std::endl;

  std::cout << "DenseDMatArray::values (" << m_values.size() << " values in "
            << m_values.nb_blocks() << " blocks)" << std::endl;
  std::cout << m_values << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif
