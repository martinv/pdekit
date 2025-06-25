#ifndef PDEKIT_Math_Tensor_Rank_hpp
#define PDEKIT_Math_Tensor_Rank_hpp

namespace pdekit
{

namespace math
{

enum TensorRankValueList
{
  tensor_rank_0 = 0,
  tensor_rank_1 = 1,
  tensor_rank_2 = 2
};

struct RankValueScalar
{
  enum
  {
    value = tensor_rank_0
  };
  enum
  {
    nb_dims = 0u
  };
};

struct RankValueVector
{
  enum
  {
    value = tensor_rank_1
  };
  enum
  {
    nb_dims = 1u
  };
};

struct RankValueTensor2
{
  enum
  {
    value = tensor_rank_2
  };
  enum
  {
    nb_dims = 2u
  };
};

// This class should be specialized for particular types
// like DynamicVector  (rank = 1), StaticMatrix (rank = 2),
// MatrixView (rank = 2) etc.

template <typename Tensor>
struct TensorRank
{
};

} // namespace math

} // namespace pdekit

#endif
