/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gram_schmidt_test
#include <boost/test/unit_test.hpp>

/// PDEKIT headers
#include "math/DenseDVec.hpp"
#include "math/algo/GramSchmidt.hpp"

using namespace pdekit;

class ScalarProduct
{
  public:
  Real operator()(const math::DenseDVec<Real> &a, const math::DenseDVec<Real> &b) const
  {
    if (a.size() != b.size())
    {
      return 0.0;
    }

    if (a.size() == 0)
    {
      return 0.0;
    }

    Real result = a[0] * b[0];
    for (Uint i = 1; i < a.size(); ++i)
    {
      result += a[i] * b[i];
    }
    return result;
  }
};

BOOST_AUTO_TEST_CASE(gram_schmidt_utest)
{
  std::vector<math::DenseDVec<Real>> vec_in;
  std::vector<math::DenseDVec<Real>> vec_out;

  vec_in.resize(2);

  math::DenseDVec<Real> &vec0 = vec_in[0];
  vec0.resize(2);
  vec0[0] = 1.0;
  vec0[1] = 2.0;

  math::DenseDVec<Real> &vec1 = vec_in[1];
  vec1.resize(2);
  vec1[0] = 3.0;
  vec1[1] = 4.0;

  // math::algo::GramSchmidt<math::DenseDVec<Real>, ScalarProduct>(vec_in,
  // vec_out);

  ScalarProduct scal_prod;
  math::algo::GramSchmidt(scal_prod, vec_in, vec_out);

  for (const auto &v : vec_out)
  {
    std::cout << v << std::endl;
  }
}
