#ifndef PDEKIT_Common_PtrHandle_hpp
#define PDEKIT_Common_PtrHandle_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

template <typename T>
struct DefaultMemoryManagement
{
  typedef T *raw_ptr_type;

  // Default: do nothing
  static void release_memory(raw_ptr_type &ptr)
  {
    ptr = nullptr;
  }
};

// ----------------------------------------------------------------------------

template <typename T>
struct DeleteMemoryOnDestruction
{
  typedef T *raw_ptr_type;

  // Release the memory on destruction:
  static void release_memory(raw_ptr_type &ptr)
  {
    delete ptr;
    ptr = nullptr;
  }
};

// ----------------------------------------------------------------------------

template <typename T, typename Policy = DefaultMemoryManagement<T>>
class PtrHandle
{
  public:
  /// Empty constructor
  PtrHandle();

  /// Constructor with a parameter - naked pointer
  PtrHandle(T *ptr);

  /// Copy constructor
  PtrHandle(const PtrHandle &other);

  /// Copy constructor - construct from another pointer
  template <typename T2, typename Policy2>
  PtrHandle(const PtrHandle<T2, Policy2> &other);

  /// Destructor
  ~PtrHandle();

  /// Overloaded assignement operator
  const PtrHandle &operator=(const PtrHandle &other);

  /// Overloaded assignement operator
  const PtrHandle &operator=(const T *other);

  /// Reset the wrapped pointer to some other naked pointer
  void reset(T *ptr);

  /// Reset the wrapped pointer to some other wrapped pointer
  void reset(const PtrHandle &other);

  /// Overloaded equality operator
  bool operator==(const PtrHandle &other);

  /// Overloading of "!="
  bool operator!=(const PtrHandle &other);

  /// Overloading of "<"
  bool operator<(const PtrHandle &other);

  /// Overloading of ">"
  bool operator>(const PtrHandle &other);

  /// Check if the wrapped pointer is null
  bool is_null() const;

  /// Check if the wrapped pointer is different from null
  bool is_not_null() const;

  /// Overloading of "*"
  /// @return reference to the object pointed
  T &operator*() const;

  /// Get the naked underlying pointer
  // T* get() const;

  private:
  /// FRIENDS

  template <typename T2, typename Policy2>
  friend class PtrHandle;

  /// The actual raw pointer to data
  T *m_ptr;
};

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
PtrHandle<T, Policy>::PtrHandle() : m_ptr(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
PtrHandle<T, Policy>::PtrHandle(T *ptr) : m_ptr(ptr)
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
PtrHandle<T, Policy>::PtrHandle(const PtrHandle &other) : m_ptr(other.m_ptr)
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
template <typename T2, typename Policy2>
PtrHandle<T, Policy>::PtrHandle(const PtrHandle<T2, Policy2> &other) : m_ptr(other.m_ptr)
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
PtrHandle<T, Policy>::~PtrHandle()
{
  Policy::release_memory(m_ptr);
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
inline const PtrHandle<T, Policy> &PtrHandle<T, Policy>::operator=(const PtrHandle &other)
{
  reset(other);
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
inline const PtrHandle<T, Policy> &PtrHandle<T, Policy>::operator=(const T *other)
{
  reset(other);
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
void PtrHandle<T, Policy>::reset(T *ptr)
{
  m_ptr = ptr;
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
void PtrHandle<T, Policy>::reset(const PtrHandle &other)
{
  m_ptr = other.m_ptr;
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
inline bool PtrHandle<T, Policy>::operator==(const PtrHandle &other)
{
  return (m_ptr == other.m_ptr);
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
inline bool PtrHandle<T, Policy>::operator!=(const PtrHandle &other)
{
  return (m_ptr == other.m_ptr);
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
inline bool PtrHandle<T, Policy>::operator<(const PtrHandle &other)
{
  return (m_ptr < other.m_ptr);
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
inline bool PtrHandle<T, Policy>::operator>(const PtrHandle &other)
{
  return (m_ptr > other.m_ptr);
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
inline bool PtrHandle<T, Policy>::is_null() const
{
  return m_ptr == nullptr;
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
inline bool PtrHandle<T, Policy>::is_not_null() const
{
  return (m_ptr != nullptr);
}

// ----------------------------------------------------------------------------

template <typename T, typename Policy>
inline T &PtrHandle<T, Policy>::operator*() const
{
  return *m_ptr;
}

// ----------------------------------------------------------------------------

/*
template<typename T, typename Policy>
T* PtrHandle<T,Policy>::get() const
{
  return m_ptr;
}
*/

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
