#ifndef PDEKIT_Interpolation_Scalar_Mesh_Function_Sub_Expr_hpp
#define PDEKIT_Interpolation_Scalar_Mesh_Function_Sub_Expr_hpp

#include "common/Meta.hpp"
#include "interpolation/mesh_function/ScalarMeshFunctionBase.hpp"
#include "math/traits/SubTrait.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename ScalMeshFunc1, typename ScalMeshFunc2>
class ScalarMeshFunctionSubExpr
    : public ScalarMeshFunctionBase<ScalarMeshFunctionSubExpr<ScalMeshFunc1, ScalMeshFunc2>>
{
  public:
  /// TYPEDEFS
  typedef typename math::SubTrait<typename ScalMeshFunc1::entry_type,
                                  typename ScalMeshFunc2::entry_type>::type entry_type;

  typedef entry_type const const_entry_type;

  typedef const ScalMeshFunc1 &Lhs;
  typedef const ScalMeshFunc2 &Rhs;

  /// Constructor
  ScalarMeshFunctionSubExpr(const ScalMeshFunc1 &lhs, const ScalMeshFunc2 &rhs)
      : m_lhs(lhs), m_rhs(rhs)
  {
  }

  /// Return one column of the data array
  /// @param index of the row
  /// @return Slice vector representing one row
  const const_entry_type operator[](const Uint i) const
  {
    return m_lhs[i] - m_rhs[i];
  }

  /// Get number of rows
  Uint nb_fields() const
  {
    return 1u;
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

/// Operator for subtraction of two scalar mesh functions
template <typename ScalMeshFunc1, typename ScalMeshFunc2>
inline const ScalarMeshFunctionSubExpr<ScalMeshFunc1, ScalMeshFunc2> operator-(
    const ScalarMeshFunctionBase<ScalMeshFunc1> &lhs,
    const ScalarMeshFunctionBase<ScalMeshFunc2> &rhs)
{
  return ScalarMeshFunctionSubExpr<ScalMeshFunc1, ScalMeshFunc2>(lhs.wrapped_type(),
                                                                 rhs.wrapped_type());
}

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // Namespace pdekit

#endif
