#ifndef PDEKIT_Interpolation_Vector_Mesh_Function_Sub_Expr_hpp
#define PDEKIT_Interpolation_Vector_Mesh_Function_Sub_Expr_hpp

#include "common/Meta.hpp"
#include "interpolation/mesh_function/VectorMeshFunctionBase.hpp"
#include "math/DenseDVec.hpp"
#include "math/binary_ops/proxy_ops/DVecDVecSubExprProxy.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename VectorMeshFunc1, typename VectorMeshFunc2>
class VectorMeshFunctionSubExpr
    : public VectorMeshFunctionBase<VectorMeshFunctionSubExpr<VectorMeshFunc1, VectorMeshFunc2>>
{
  public:
  /// TYPEDEFS
  typedef typename math::DVecDVecSubExprProxy<typename VectorMeshFunc1::const_entry_type,
                                              typename VectorMeshFunc2::const_entry_type,
                                              math::ColumnVector>
      const_entry_type;

  typedef const VectorMeshFunc1 &Lhs;
  typedef const VectorMeshFunc2 &Rhs;

  /// Constructor
  VectorMeshFunctionSubExpr(const VectorMeshFunc1 &lhs, const VectorMeshFunc2 &rhs)
      : m_lhs(lhs), m_rhs(rhs)
  {
  }

  /// Return one column of the data array
  /// @param index of the row
  /// @return Slice vector representing one row
  const const_entry_type const_value(const Uint i) const
  {
    return const_entry_type(m_lhs.const_value(i), m_rhs.const_value(i));
  }

  /// Get number of rows
  Uint nb_fields() const
  {
    return m_lhs.nb_fields();
  }

  /// Get number of columns
  Uint nb_entries() const
  {
    return m_lhs.nb_entries();
  }

  private:
  Lhs m_lhs;

  Rhs m_rhs;
};

// ----------------------------------------------------------------------------

/// Operator for Subition of two vector mesh functions
template <typename VectorMeshFunc1, typename VectorMeshFunc2>
inline const VectorMeshFunctionSubExpr<VectorMeshFunc1, VectorMeshFunc2> operator-(
    const VectorMeshFunctionBase<VectorMeshFunc1> &lhs,
    const VectorMeshFunctionBase<VectorMeshFunc2> &rhs)
{
  return VectorMeshFunctionSubExpr<VectorMeshFunc1, VectorMeshFunc2>(lhs.wrapped_type(),
                                                                     rhs.wrapped_type());
}

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // Namespace pdekit

#endif
