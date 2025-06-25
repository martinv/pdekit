#ifndef PDEKIT_Interpolation_Mesh_Function_Norm_hpp
#define PDEKIT_Interpolation_Mesh_Function_Norm_hpp

#include "interpolation/mesh_function/ScalarMeshFunctionBase.hpp"
#include "interpolation/mesh_function/VectorMeshFunctionBase.hpp"
#include "math/DenseDVec.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename ScalMeshFunc, typename T, bool TF>
void norm_L1(const ScalarMeshFunctionBase<ScalMeshFunc> &mesh_func_expr,
             math::DenseDVec<T, TF> &norm)
{
  const ScalMeshFunc &mesh_func = mesh_func_expr.wrapped_type();

  if (norm.size() != 1)
  {
    norm.resize(1);
  }
  norm[0] = T();

  for (Uint i = 0; i < mesh_func.nb_entries(); ++i)
  {
    norm[0] += std::abs(mesh_func[i]);
  }

  norm[0] /= mesh_func.nb_entries();
}

// ----------------------------------------------------------------------------

template <typename ScalMeshFunc, typename T, bool TF>
void norm_L2(const ScalarMeshFunctionBase<ScalMeshFunc> &mesh_func_expr,
             math::DenseDVec<T, TF> &norm)
{
  const ScalMeshFunc &mesh_func = mesh_func_expr.wrapped_type();

  if (norm.size() != 1)
  {
    norm.resize(1);
  }
  norm[0] = T();

  for (Uint i = 0; i < mesh_func.nb_entries(); ++i)
  {
    norm[0] += (mesh_func[i] * mesh_func[i]);
  }
  norm[0] = std::sqrt(norm[0] / mesh_func.nb_entries());
}

// ----------------------------------------------------------------------------

template <typename VectMeshFunc, typename T, bool TF>
void norm_L1(const VectorMeshFunctionBase<VectMeshFunc> &mesh_func_expr,
             math::DenseDVec<T, TF> &norm)
{
  const VectMeshFunc &mesh_func = mesh_func_expr.wrapped_type();
  typedef typename VectMeshFunc::const_entry_type const_entry_type;

  const Uint nb_fields = mesh_func.nb_fields();

  if (norm.size() != nb_fields)
  {
    norm.resize(nb_fields);
  }
  norm.fill(T());

  for (Uint i = 0; i < mesh_func.nb_entries(); ++i)
  {
    const const_entry_type entry = mesh_func.const_value(i);

    for (Uint f = 0; f < nb_fields; ++f)
    {
      norm[f] += std::abs(entry[f]);
    }
  }

  for (Uint f = 0; f < nb_fields; ++f)
  {
    norm[f] /= mesh_func.nb_entries();
  }
}

// ----------------------------------------------------------------------------

template <typename VectMeshFunc, typename T, bool TF>
void norm_L2(const VectorMeshFunctionBase<VectMeshFunc> &mesh_func_expr,
             math::DenseDVec<T, TF> &norm)
{
  const VectMeshFunc &mesh_func = mesh_func_expr.wrapped_type();
  typedef typename VectMeshFunc::const_entry_type const_entry_type;

  const Uint nb_fields = mesh_func.nb_fields();

  if (norm.size() != nb_fields)
  {
    norm.resize(nb_fields);
  }
  norm.fill(T());

  for (Uint i = 0; i < mesh_func.nb_entries(); ++i)
  {
    const const_entry_type entry = mesh_func.const_value(i);

    for (Uint f = 0; f < nb_fields; ++f)
    {
      norm[f] += (entry[f] * entry[f]);
    }
  }

  for (Uint f = 0; f < nb_fields; ++f)
  {
    norm[f] = std::sqrt(norm[f] / mesh_func.nb_entries());
  }
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
