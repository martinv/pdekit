#ifndef PDEKIT_Common_IOStream_Utils_hpp
#define PDEKIT_Common_IOStream_Utils_hpp

#include <iostream>

namespace pdekit
{

// ----------------------------------------------------------------------------

// Output std::pair to output stream

template <typename T1, typename T2>
std::ostream &operator<<(std::ostream &os, std::pair<T1, T2> const &key)
{
  os << "(" << key.first << "," << key.second << ")";
  return os;
}

// ----------------------------------------------------------------------------

} // Namespace pdekit

#endif
