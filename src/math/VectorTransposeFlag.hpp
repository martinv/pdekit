#ifndef PDEKIT_Math_Vector_Transpose_Flag_hpp
#define PDEKIT_Math_Vector_Transpose_Flag_hpp

namespace pdekit
{

namespace math
{

/// This variable is used a 'transpose flag' template parameters for vectors
const bool ColumnVector = false;

/// This variable is used a 'transpose flag' template parameters for vectors
const bool RowVector = true;

/// Extra variable to denote default transpose flag for vectors
const bool DefaultVectorTransposeFlag = ColumnVector;

} // namespace math

} // namespace pdekit

#endif
