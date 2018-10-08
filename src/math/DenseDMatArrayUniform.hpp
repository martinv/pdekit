#ifndef PDEKIT_Math_Dense_Dynamic_Matrix_Array_Uniform_hpp
#define PDEKIT_Math_Dense_Dynamic_Matrix_Array_Uniform_hpp

#include "common/ArrayShape.hpp"
#include "common/BlockArray.hpp"
#include "math/DenseConstMatView.hpp"
#include "math/DenseMatView.hpp"

namespace pdekit
{

namespace math
{

template <typename T, bool SO = DefaultMatrixStorageOrder>
class DenseDMatArrayUniform
{
  public:
  /// Default constructor
  DenseDMatArrayUniform();

  /// Copy constructor
  DenseDMatArrayUniform(const DenseDMatArrayUniform &other);

  /// Destructor
  ~DenseDMatArrayUniform();

  /// Assignment operator
  DenseDMatArrayUniform &operator=(const DenseDMatArrayUniform &rhs);

  /// Allocate data
  void allocate(common::ArrayShape<_2D, SUint> mat_shape, const Uint nb_matrices);

  /// Get one matrix view, non-const version
  DenseMatView<T, SO> mat_view(const Uint mat_id);

  /// Get one matrix view, non-const version
  DenseConstMatView<T, SO> const_mat_view(const Uint mat_id) const;

  /// Print information for debugging
  void debug_print() const;

  private:
  /// Local (element-wise) array of local interpolation matrices
  std::vector<T> m_values;

  /// Number of matrices stored
  Uint m_num_mats;

  /// The shape of each matrix
  common::ArrayShape<_2D, SUint> m_mat_shape;
};

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMatArrayUniform<T, SO>::DenseDMatArrayUniform()
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMatArrayUniform<T, SO>::DenseDMatArrayUniform(const DenseDMatArrayUniform &other)
    : m_values(other.m_values), m_num_mats(0), m_mat_shape(other.m_mat_shape)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMatArrayUniform<T, SO>::~DenseDMatArrayUniform()
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseDMatArrayUniform<T, SO> &DenseDMatArrayUniform<T, SO>::operator=(
    const DenseDMatArrayUniform &rhs)
{
  m_values.resize(rhs.m_values.size());
  m_values = rhs.m_values;

  m_num_mats = rhs.m_num_mats;

  m_mat_shape = rhs.m_mat_shape;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void DenseDMatArrayUniform<T, SO>::allocate(common::ArrayShape<_2D, SUint> mat_shape,
                                            const Uint nb_matrices)
{
  m_mat_shape = mat_shape;
  m_num_mats  = nb_matrices;

  const Uint num_values = nb_matrices * mat_shape.size(0) * mat_shape.size(1);

  m_values.resize(num_values);
  m_values.assign(num_values, T());
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseMatView<T, SO> DenseDMatArrayUniform<T, SO>::mat_view(const Uint mat_id)
{
  const Uint one_mat_size = m_mat_shape.size(0) * m_mat_shape.size(1);
  T *data_ptr             = m_values.data() + one_mat_size * mat_id;
  const DenseMatView<T, SO> view(data_ptr, m_mat_shape.size(1), m_mat_shape.size(0),
                                 m_mat_shape.size(1));
  return view;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
DenseConstMatView<T, SO> DenseDMatArrayUniform<T, SO>::const_mat_view(const Uint mat_id) const
{
  const Uint one_mat_size = m_mat_shape.size(0) * m_mat_shape.size(1);
  const T *data_ptr       = m_values.data() + one_mat_size * mat_id;
  const DenseConstMatView<T, SO> view(data_ptr, m_mat_shape.size(1), m_mat_shape.size(0),
                                      m_mat_shape.size(1));
  return view;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void DenseDMatArrayUniform<T, SO>::debug_print() const
{
  std::cout << "DenseDMatArrayUniform::shapes:" << std::endl;
  std::cout << " [" << m_mat_shape.size(0) << "," << m_mat_shape.size(1) << "]";

  std::cout << std::endl;

  std::cout << "DenseDMatArrayUniform::values (" << m_values.size() << " values in " << m_num_mats
            << " blocks)" << std::endl;
  for (auto val : m_values)
  {
    std::cout << val << " ";
  }
  std::cout << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif
