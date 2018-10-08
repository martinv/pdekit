#ifndef PDEKIT_Tools_Loop_hpp
#define PDEKIT_Tools_Loop_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace tools
{

template <typename FUNC, Uint N>
class Loop0
{
  public:
  static inline void run()
  {
    Loop0<FUNC, N - 1>::run();
    FUNC::template exec<N - 1>();
  }
};

template <typename FUNC>
class Loop0<FUNC, 0>
{

  public:
  static inline void run()
  {
  }
};

//####################################################

template <typename FUNC, Uint N>
class Loop1
{

  public:
  template <typename P1>
  static inline void run(P1 &p1)
  {
    Loop1<FUNC, N - 1>::run(p1);
    FUNC::template exec<N - 1>(p1);
  }
};

// Instiantiate class to stop recursion:
template <typename FUNC>
class Loop1<FUNC, 0>
{

  public:
  template <typename P1>
  static inline void run(P1 &p1)
  {
  }
};

} // Namespace tools

} // Namespace pdekit

#endif
