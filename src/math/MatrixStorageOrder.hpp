#ifndef PDEKIT_Math_Matrix_Storage_Order_hpp
#define PDEKIT_Math_Matrix_Storage_Order_hpp

namespace pdekit
{

namespace math
{

/// This variable is used a 'storage order (SO)' template parameter for matrices
const bool RowMajor = false;

/// This variable is used a 'storage order (SO)' template parameter for matrices
const bool ColumnMajor = true;

/// Extra variable to denote default storage order for matrices
const bool DefaultMatrixStorageOrder = RowMajor;

} // namespace math

} // namespace pdekit

#endif
