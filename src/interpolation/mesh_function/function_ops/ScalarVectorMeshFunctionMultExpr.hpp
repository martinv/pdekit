#ifndef PDEKIT_Interpolation_Scalar_Vector_Mesh_Function_Mult_Expr_hpp
#define PDEKIT_Interpolation_Scalar_Vector_Mesh_Function_Mult_Expr_hpp

#include "common/Meta.hpp"
#include "interpolation/mesh_function/VectorMeshFunctionBase.hpp"
#include "math/binary_ops/proxy_ops/ScalarDVecMultExprProxy.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename ST, typename VectorMeshFunc>
class ScalarVectorMeshFunctionMultExpr
    : public VectorMeshFunctionBase<ScalarVectorMeshFunctionMultExpr<ST, VectorMeshFunc>>
{
  public:
  /// TYPEDEFS
  typedef typename math::ScalarDVecMultExprProxy<ST, typename VectorMeshFunc::const_entry_type,
                                                 math::ColumnVector>
      const_entry_type;

  typedef const ST Lhs;
  typedef const VectorMeshFunc &Rhs;

  /// Constructor
  ScalarVectorMeshFunctionMultExpr(const ST lhs, const VectorMeshFunc &rhs) : m_lhs(lhs), m_rhs(rhs)
  {
  }

  /// Return one column of the data array
  /// @param index of the row
  /// @return Slice vector representing one row
  const const_entry_type const_value(const Uint i) const
  {
    return const_entry_type(m_lhs, m_rhs.const_value(i));
  }

  /// Get number of rows
  Uint nb_fields() const
  {
    return 1u;
  }

  /// Get number of columns
  Uint nb_entries() const
  {
    return m_rhs.nb_entries();
  }

  private:
  Lhs m_lhs;

  Rhs m_rhs;
};

// ----------------------------------------------------------------------------

/// Operator for multiplication of a scalar mesh function by a scalar of type
/// 'Real'

template <typename VectorMeshFunc>
inline const ScalarVectorMeshFunctionMultExpr<Real, VectorMeshFunc> operator*(
    const Real factor, const VectorMeshFunctionBase<VectorMeshFunc> &rhs)
{
  return ScalarVectorMeshFunctionMultExpr<Real, VectorMeshFunc>(factor, rhs.wrapped_type());
}

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // Namespace pdekit

#endif
