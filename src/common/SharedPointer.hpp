#ifndef Shared_Pointer_hpp
#define Shared_Pointer_hpp

namespace pdekit
{

namespace common
{

/// Helper structure to avoid the shared pointer from deleting the
/// instance of the singleton

struct NullDeleter
{
  inline void operator()(void const *) const
  {
  }
};

} // namespace common

} // namespace pdekit

#endif
