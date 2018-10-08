#ifndef PDEKIT_Interpolation_Scalar_Scalar_Mesh_Function_Mult_Expr_hpp
#define PDEKIT_Interpolation_Scalar_Scalar_Mesh_Function_Mult_Expr_hpp

#include "common/Meta.hpp"
#include "interpolation/mesh_function/ScalarMeshFunctionBase.hpp"
#include "math/traits/MultTrait.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename ST, typename ScalMeshFunc>
class ScalarScalarMeshFunctionMultExpr
    : public ScalarMeshFunctionBase<ScalarScalarMeshFunctionMultExpr<ST, ScalMeshFunc>>
{
  public:
  /// TYPEDEFS
  typedef typename math::MultTrait<ST, typename ScalMeshFunc::entry_type>::type entry_type;

  typedef entry_type const const_entry_type;

  typedef const ST Lhs;
  typedef const ScalMeshFunc &Rhs;

  /// Constructor
  ScalarScalarMeshFunctionMultExpr(const ST lhs, const ScalMeshFunc &rhs) : m_lhs(lhs), m_rhs(rhs)
  {
  }

  /// Return one column of the data array
  /// @param index of the row
  /// @return Slice vector representing one row
  const const_entry_type operator[](const Uint i) const
  {
    return m_lhs * m_rhs[i];
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

template <typename ScalMeshFunc>
inline const ScalarScalarMeshFunctionMultExpr<Real, ScalMeshFunc> operator*(
    const Real factor, const ScalarMeshFunctionBase<ScalMeshFunc> &rhs)
{
  return ScalarScalarMeshFunctionMultExpr<Real, ScalMeshFunc>(factor, rhs.wrapped_type());
}

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // Namespace pdekit

#endif
