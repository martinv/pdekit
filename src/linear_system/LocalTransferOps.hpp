#ifndef PDEKIT_Linear_System_Local_Transfer_Ops_hpp
#define PDEKIT_Linear_System_Local_Transfer_Ops_hpp

#include "math/DenseDMatArray.hpp"

namespace pdekit
{

namespace ls
{

template <typename T, bool SO = math::DefaultMatrixStorageOrder>
class LocalTransferOps
{
  public:
  /// Default constructor
  LocalTransferOps();

  /// Copy constructor
  LocalTransferOps(const LocalTransferOps &other);

  /// Destructor
  ~LocalTransferOps();

  /// Assignment operator
  LocalTransferOps &operator=(const LocalTransferOps &rhs);

  /// Allocate member data
  void allocate(std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> &&mat_shapes,
                std::unique_ptr<common::BlockArray<T, Uint>> &&mat_values,
                std::unique_ptr<std::vector<Uint>> &&idx_to_op_map);

  /// Get one matrix view, non-const version
  math::DenseMatView<T, SO> op(const Uint idx);

  /// Get one matrix view, const version
  math::DenseConstMatView<T, SO> const_op(const Uint idx) const;

  /// Return number of matrices in array
  Uint size() const;

  /// Print the data
  void print() const;

  private:
  /// Storage for local operator matrices
  math::DenseDMatArray<T, SO> m_ops;

  /// Mapping index -> local matrix operator
  std::vector<Uint> m_idx_to_op_map;
};

// ----------------------------------------------------------------------------

template <typename T, bool SO>
LocalTransferOps<T, SO>::LocalTransferOps()
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
LocalTransferOps<T, SO>::LocalTransferOps(const LocalTransferOps &other)
    : m_ops(other.m_ops), m_idx_to_op_map(other.m_idx_to_op_map)
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
LocalTransferOps<T, SO>::~LocalTransferOps()
{
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
LocalTransferOps<T, SO> &LocalTransferOps<T, SO>::operator=(const LocalTransferOps &rhs)
{
  m_ops           = rhs.m_ops;
  m_idx_to_op_map = rhs.m_idx_to_op_map;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void LocalTransferOps<T, SO>::allocate(
    std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> &&mat_shapes,
    std::unique_ptr<common::BlockArray<T, Uint>> &&mat_values,
    std::unique_ptr<std::vector<Uint>> &&idx_to_op_map)
{
  m_ops.allocate(std::move(mat_shapes), std::move(mat_values));
  m_idx_to_op_map.swap(*idx_to_op_map);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
math::DenseMatView<T, SO> LocalTransferOps<T, SO>::op(const Uint idx)
{
  return m_ops.mat_view(m_idx_to_op_map[idx]);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
math::DenseConstMatView<T, SO> LocalTransferOps<T, SO>::const_op(const Uint idx) const
{
  return m_ops.const_mat_view(m_idx_to_op_map[idx]);
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
Uint LocalTransferOps<T, SO>::size() const
{
  return m_idx_to_op_map.size();
}

// ----------------------------------------------------------------------------

template <typename T, bool SO>
void LocalTransferOps<T, SO>::print() const
{
  std::cout << "Local transfer ops: " << m_ops.size()
            << " operators; map size = " << m_idx_to_op_map.size() << std::endl;
  for (Uint i = 0; i < this->size(); ++i)
  {
    const math::DenseConstMatView<T, SO> op = this->const_op(i);
    std::cout << "{" << i << "} [" << op.rows() << " x " << op.cols() << "]" << std::endl;
    std::cout << op << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif
