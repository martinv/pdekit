#ifndef PDEKIT_Math_Matrix_hpp
#define PDEKIT_Math_Matrix_hpp

namespace pdekit
{

namespace math
{

/// ---------------------------------------------------------------------------
/// Template class to mark a class 'MT' as Matrix
/// MT ... matrix type
/// SO ... storage order
/// ---------------------------------------------------------------------------

template <typename MT, bool SO>
class Matrix
{
  public:
  /// Typedef for the matrix type wrapped by this class
  using MatrixType = MT;

  /// Return this converted to wrapped type
  inline MatrixType &wrapped_type()
  {
    return static_cast<MatrixType &>(*this);
  }

  /// Return this instance converted to wrapped type, const version
  inline const MatrixType &wrapped_type() const
  {
    return static_cast<const MatrixType &>(*this);
  }
};

/// ---------------------------------------------------------------------------
/// Template class to mark a class 'MT' as dense matrix
/// MT ... Matrix type
/// SO ... storage order (row-major or column major)
/// Currently, a dense Matrix can be the class StaticMatrix, DynamicMatrix,
/// MatrixBlock or any other class that represents a mathematical expression
/// whose result is a dense Matrix
/// (for example, a sum DynamicMatrix + DynamicMatrix, or StaticMatrix +
/// MatrixBlock or the product constant * MatrixBlock )
/// ---------------------------------------------------------------------------

template <typename MT, bool SO>
struct DenseMatrix : public Matrix<MT, SO>
{
};

/// ---------------------------------------------------------------------------
/// Template class to mark a class 'MT' as sparse matrix
/// MT ... Matrix type
/// SO ... storage order (row-major or column major)
/// ---------------------------------------------------------------------------

template <typename MT, bool SO>
struct SparseMatrix : public Matrix<MT, SO>
{
};

} // namespace math

} // namespace pdekit

#endif
