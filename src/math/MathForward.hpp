#ifndef PDEKIT_Math_Math_Types_hpp
#define PDEKIT_Math_Math_Types_hpp

#include "common/PDEKit.hpp"
#include "math/MathMeta.hpp"
#include "math/VectorTransposeFlag.hpp"

namespace pdekit
{

namespace math
{

/// Forward declarations

template <typename T>
class ScalarConstant;

template <typename T, bool TF>
class DenseDVec;

template <typename T, Uint N, bool TF>
class DenseSVec;

template <typename T, bool TF>
class DenseConstVecView;

template <typename T, bool TF>
class DenseVecView;

template <typename T, bool SO>
class DenseDMat;

template <typename T, Uint M, Uint N, bool SO>
class DenseSMat;

template <typename T, bool SO>
class DenseConstMatView;

template <typename T, bool SO>
class DenseMatView;

template <typename T, bool SO>
class SparseDMat;

} // namespace math

} // namespace pdekit

#endif
