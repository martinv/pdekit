// A singleton defined by using CRTP
// The original code was found here:
// http://enki-tech.blogspot.be/2012/08/c11-generic-singleton.html
// Usage:
//
// class MyClass: public Singleton<MyClass>
// {
//   private:
//      friend class Singleton<MyClass>;
//
//      // Make sure the class is not copyable!
//         MyClass(const MyClass& other);
//         MyClass& operator=(const MyClass& other);
//
//      // You can define custom constructor
//      // MyClass(Real param1, Uint param2);
//
// };
//
//    int main()
//    {
//      Real p1;
//      Uint p2;
//      MyClass& my_singleton = MyClass::instance(p1,p2);
//
//    };

#ifndef PDEKIT_common_CRTP_Singleton_Shell_hh
#define PDEKIT_common_CRTP_Singleton_Shell_hh

#include <cstdlib>
#include <exception>
#include <iostream>

#include "common/SharedPointer.hpp"

namespace pdekit
{

namespace common
{

template <class T>
class CRTP_Singleton
{
  public:
  typedef T instance_type;
  typedef T &reference_type;

  /// Get the instance of the singleton
  template <typename... Args>
  static T &instance(Args... args);

  protected:
  // Plain constructor
  CRTP_Singleton();

  // Private copy constructor:
  CRTP_Singleton(const CRTP_Singleton &source);

  // Private destructor:
  // Client holding pointer to Singleton cannot delete it accidentally
  ~CRTP_Singleton();

  // Private assignement operator:
  CRTP_Singleton &operator=(const CRTP_Singleton &other);

  // Create one instance of a singleton
  template <typename... Args>
  static void create(Args... args);

  static void on_dead_reference();

  /// Data

  static T *m_instance;

  static bool m_destroyed;
};

/// ============================================================================

/// Initialization of static variables:
template <typename T>
T *CRTP_Singleton<T>::m_instance = 0;

template <typename T>
bool CRTP_Singleton<T>::m_destroyed = false;

/// ======================================================================================

/// Member functions:

template <typename T>
CRTP_Singleton<T>::CRTP_Singleton()
{
}

/// ============================================================================

template <typename T>
CRTP_Singleton<T>::CRTP_Singleton(const CRTP_Singleton<T> &source)
{
}

/// ============================================================================

template <typename T>
CRTP_Singleton<T>::~CRTP_Singleton()
{
  m_instance  = 0;
  m_destroyed = true;
}

/// ============================================================================

template <typename T>
template <typename... Args>
void CRTP_Singleton<T>::create(Args... args)
{
  static T the_instance(std::forward<Args>(args)...);
  m_instance = &the_instance;
}

/// ============================================================================

template <typename T>
void CRTP_Singleton<T>::on_dead_reference()
{
  // throw std::runtime_error("Dead reference of SingletonShell detected");
  std::cerr << "Dead reference of SingletonShell detected";
}

/// ============================================================================

template <typename T>
template <typename... Args>
T &CRTP_Singleton<T>::instance(Args... args)
{
  if (!m_instance)
  {
    // Check for dead reference:
    if (m_destroyed)
    {
      on_dead_reference();
    }
    else
    {
      create(std::forward<Args>(args)...);
      // Initializer::setup_instance(*m_instance);
    }
  }
  return *m_instance;
}

} // Namespace common

} // Namespace pdekit

#endif // Singleton_Shell_hh
