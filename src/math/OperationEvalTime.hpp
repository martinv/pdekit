#ifndef PDEKIT_Math_Operation_Eval_Time_hpp
#define PDEKIT_Math_Operation_Eval_Time_hpp

namespace pdekit
{

namespace math
{

enum OperationEvalComplexity
{
  eval_complexity_constant  = 0,
  eval_complexity_linear    = 1,
  eval_complexity_quadratic = 2,
  eval_complexity_cubic     = 3
};

} // namespace math

} // namespace pdekit

#endif
