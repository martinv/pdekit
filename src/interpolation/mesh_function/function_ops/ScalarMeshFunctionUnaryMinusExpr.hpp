#ifndef PDEKIT_Interpolation_Scalar_Mesh_Function_Unary_Minus_Expr_hpp
#define PDEKIT_Interpolation_Scalar_Mesh_Function_Unary_Minus_Expr_hpp

#include "interpolation/mesh_function/ScalarMeshFunctionBase.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename ScalMeshFunc>
class ScalarMeshFunctionUnaryMinusExpr
    : public ScalarMeshFunctionBase<ScalarMeshFunctionUnaryMinusExpr<ScalMeshFunc>>
{
  public:
  /// TYPEDEFS
  typedef typename ScalMeshFunc::const_entry_type const_entry_type;
  typedef typename ScalMeshFunc::entry_type entry_type;

  typedef const ScalMeshFunc &Rhs;

  /// Constructor
  ScalarMeshFunctionUnaryMinusExpr(const ScalMeshFunc &rhs) : m_rhs(rhs)
  {
  }

  /// Return one column of the data array
  /// @param index of the row
  /// @return Slice vector representing one row
  const const_entry_type operator[](const Uint i) const
  {
    return -m_rhs[i];
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
  Rhs m_rhs;
};

// ----------------------------------------------------------------------------

/// Operator for multiplication of a scalar mesh function by a scalar of type
/// 'Real'

template <typename ScalMeshFunc>
inline const ScalarMeshFunctionUnaryMinusExpr<ScalMeshFunc> operator-(
    const ScalarMeshFunctionBase<ScalMeshFunc> &scalar_mesh_function)
{
  return ScalarMeshFunctionUnaryMinusExpr<ScalMeshFunc>(scalar_mesh_function.wrapped_type());
}

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // Namespace pdekit

#endif
