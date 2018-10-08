#ifndef PDEKIT_Interpolation_Scalar_Mesh_Function_Vector_Mesh_Function_Mult_Expr_hpp
#define PDEKIT_Interpolation_Scalar_Mesh_Function_Vector_Mesh_Function_Mult_Expr_hpp

#include "common/Meta.hpp"
#include "interpolation/mesh_function/ScalarMeshFunctionBase.hpp"
#include "interpolation/mesh_function/VectorMeshFunctionBase.hpp"
#include "math/VectorTransposeFlag.hpp"
#include "math/binary_ops/proxy_ops/ScalarDVecMultExprProxy.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename ScalMeshFunc, typename VectMeshFunc>
class ScalarMeshFVectorMeshFMultExpr
    : public VectorMeshFunctionBase<ScalarMeshFVectorMeshFMultExpr<ScalMeshFunc, VectMeshFunc>>
{
  public:
  /// TYPEDEFS
  typedef typename math::ScalarDVecMultExprProxy<typename ScalMeshFunc::value_type,
                                                 typename VectMeshFunc::const_entry_type,
                                                 math::ColumnVector>
      const_entry_type;

  typedef const ScalMeshFunc &Lhs;
  typedef const VectMeshFunc &Rhs;

  /// Constructor
  ScalarMeshFVectorMeshFMultExpr(const ScalMeshFunc &lhs, const VectMeshFunc &rhs)
      : m_lhs(lhs), m_rhs(rhs)
  {
  }

  /// Return one column of the data array
  /// @param index of the row
  /// @return Slice vector representing one row
  const const_entry_type const_value(const Uint i) const
  {
    return const_entry_type(m_lhs[i], m_rhs.const_value(i));
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

template <typename ScalMeshFunc, typename VectMeshFunc>
inline const ScalarMeshFVectorMeshFMultExpr<ScalMeshFunc, VectMeshFunc> operator*(
    const ScalarMeshFunctionBase<ScalMeshFunc> &lhs,
    const VectorMeshFunctionBase<VectMeshFunc> &rhs)
{
  return ScalarMeshFVectorMeshFMultExpr<ScalMeshFunc, VectMeshFunc>(lhs.wrapped_type(),
                                                                    rhs.wrapped_type());
}

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // Namespace pdekit

#endif
