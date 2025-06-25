#ifndef PDEKIT_Math_Matrix_Determinant_hpp
#define PDEKIT_Math_Matrix_Determinant_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace math
{

/// ============================================================================
/// MATRIX DETERMINANT COMPUTATION
/// ============================================================================

template <typename T, Uint M>
class MatrixDeterminant
{
  public:
  /// TYPEDEFS:

  // Compute the determinant
  // This method is intentionally not implemented - one of the specializations
  // of MatrixDeterminant<T,SIZE> which follow below should be always called
  // inline static T compute(const T* const data, const Uint m);

  private:
};

/// ===============================
/// Specialization for case SIZE=1:
template <typename T>
class MatrixDeterminant<T, 1>
{
  public:
  /// TYPEDEFS:

  // Compute the determinant
  inline static T compute(const T *const data, const Uint m = 1)
  {
    return data[0];
  }

  private:
};

/// ===============================
/// Specialization for case SIZE=2:
template <typename T>
class MatrixDeterminant<T, 2>
{
  public:
  /// TYPEDEFS:

  // Compute the determinant
  inline static T compute(const T *const data, const Uint m = 2)
  {
    return data[0] * data[3] - data[1] * data[2];
  }

  private:
};

/// ===============================
/// Specialization for case SIZE=3:
template <typename T>
class MatrixDeterminant<T, 3>
{
  public:
  /// TYPEDEFS:

  // Compute the determinant
  inline static T compute(const T *const data, const Uint m = 3)
  {
    return data[0] * (data[4] * data[8] - data[5] * data[7]) -
           data[1] * (data[3] * data[8] - data[5] * data[6]) +
           data[2] * (data[3] * data[7] - data[4] * data[6]);
  }

  private:
};

/// ===============================
/// Specialization for case SIZE=4:
template <typename T>
class MatrixDeterminant<T, 4>
{
  public:
  /// TYPEDEFS:

  // Compute the determinant
  inline static T compute(const T *const data, const Uint m = 4)
  {
    T d00 = data[0];
    T d01 = data[1];
    T d02 = data[2];
    T d03 = data[3];
    T d04 = data[4];
    T d05 = data[5];
    T d06 = data[6];
    T d07 = data[7];

    T d10d15_d14d11 = data[10] * data[15] - data[14] * data[11];
    T d09d15_d13d11 = data[9] * data[15] - data[13] * data[11];
    T d09d14_d13d10 = data[9] * data[14] - data[13] * data[10];
    T d08d15_d12d11 = data[8] * data[15] - data[12] * data[11];
    T d08d13_d12d09 = data[8] * data[13] - data[12] * data[9];
    T d08d14_d12d10 = data[8] * data[14] - data[12] * data[10];

    return d00 * (d05 * (d10d15_d14d11)-d06 * (d09d15_d13d11) + d07 * (d09d14_d13d10)) -
           d01 * (d04 * (d10d15_d14d11)-d06 * (d08d15_d12d11) + d07 * (d08d14_d12d10)) +
           d02 * (d04 * (d09d15_d13d11)-d05 * (d08d15_d12d11) + d07 * (d08d13_d12d09)) -
           d03 * (d04 * (d09d14_d13d10)-d05 * (d08d14_d12d10) + d06 * (d08d13_d12d09));
  }

  private:
};

/// ===============================
/// Specialization for case SIZE=DYNAMIC:
template <typename T>
class DynamicMatrixDeterminant
{
  public:
  /// TYPEDEFS:

  // Compute the determinant
  inline static T compute(const T *const data, const Uint m)
  {

    switch (m)
    {
      case 1:
      {
        return data[0];
        break;
      }
      case 2:
      {
        return MatrixDeterminant<T, 2>::compute(data);
        break;
      }
      case 3:
      {
        return MatrixDeterminant<T, 3>::compute(data);
        break;
      }
      case 4:
      {
        return MatrixDeterminant<T, 4>::compute(data);
        break;
      }
      default:
        return T();
    };
  }

  private:
};

} // Namespace math

} // Namespace pdekit

#endif // Matrix_Determinant_hpp
