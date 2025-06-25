#ifndef PDEKIT_Interpolation_Vector_Mesh_Function_Assign_hpp
#define PDEKIT_Interpolation_Vector_Mesh_Function_Assign_hpp

#include "interpolation/mesh_function/VectorMeshFunctionBase.hpp"

namespace pdekit
{

namespace interpolation
{

namespace detail
{

struct VectorMeshFunctionAssign
{
  // --------------------------------------------------------------------------

  template <typename VF1, typename VF2>
  inline static void assign(VectorMeshFunctionBase<VF1> &expr_lhs,
                            const VectorMeshFunctionBase<VF2> &expr_rhs)
  {
    VF1 &lhs       = expr_lhs.wrapped_type();
    const VF2 &rhs = expr_rhs.wrapped_type();

    for (Uint n = 0; n < lhs.nb_entries(); ++n)
    {
      lhs.value(n) = rhs.const_value(n);
    }
  }

  // --------------------------------------------------------------------------
};

} // namespace detail

} // namespace interpolation

} // namespace pdekit

#endif
