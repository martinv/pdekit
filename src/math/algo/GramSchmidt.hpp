#ifndef PDEKIT_Math_Algo_Gram_Schmidt_hpp
#define PDEKIT_Math_Algo_Gram_Schmidt_hpp

#include "common/PDEKit.hpp"
#include <iostream>
#include <vector>

namespace pdekit
{

namespace math
{

namespace algo
{

// ----------------------------------------------------------------------------

template <typename VectorType, typename ScalarProduct>
void GramSchmidt(const ScalarProduct &scal_prod, const std::vector<VectorType> &vec_in,
                 std::vector<VectorType> &vec_out)
{
  vec_out.resize(vec_in.size());
  for (Uint i = 0; i < vec_in.size(); ++i)
  {
    vec_out[i].resize(vec_in[i].size());
  }

  vec_out[0] = vec_in[0];

  // Storage for the norm squares of every output vector
  std::vector<Real> norm2_out(vec_in.size());

  VectorType tmp_vec;
  // ScalarProduct scal_prod;

  for (Uint i = 1; i < vec_out.size(); ++i)
  {
    norm2_out[i - 1] = scal_prod(vec_out[i - 1], vec_out[i - 1]);

    if (norm2_out[i - 1] < 1.e-13)
    {
      std::cerr << "math::algo::GramSchmidt:: norm of output vector nr. " << i - 1 << " is zero"
                << std::endl;
      std::cerr << "                          will not orthogonalize vectors" << std::endl;
      return;
    }

    vec_out[i] = vec_in[i];
    for (Uint j = 0; j < i; ++j)
    {
      const Real prod = scal_prod(vec_in[i], vec_out[j]);

      tmp_vec = prod / norm2_out[j] * vec_out[j];
      vec_out[i] -= tmp_vec;
    }
  }
}

// ----------------------------------------------------------------------------

} // namespace algo

} // namespace math

} // namespace pdekit

#endif
