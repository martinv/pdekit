#ifndef PDEKIT_common_Singleton_Shell_hh
#define PDEKIT_common_Singleton_Shell_hh

#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory> // For shared_ptr

#include "common/SharedPointer.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

template <typename T>
class DefaultSingletonSetup
{
  public:
  static void setup_instance(T &instance)
  {
  }
};

// ----------------------------------------------------------------------------

template <typename T, typename Initializer = DefaultSingletonSetup<T>>
class Singleton
{

  public:
  /// TYPEDEFS:
  typedef T instance_type;
  typedef T &reference_type;

  /// Define a constructor which can fill the underlying instance of T:
  template <typename... Args>
  static T &instance(Args &&... args);

  /// Return the only instance of this class:
  // static T& instance();

  static std::shared_ptr<T> shared_ptr_to_instance();

  void print();

  private:
  /// Make T a singleton:

  // Private constructor:
  //    SingletonShell() { }

  /// Private copy constructor:
  Singleton(const Singleton &source);

  /// Private destructor:
  /// Client holding pointer to SingletonShell cannot delete it accidentally
  ~Singleton();

  /// Private assignement operator:
  Singleton &operator=(const Singleton &other);

  // static void create();

  /// Create an instance
  /// @param ... ...Args ... parameter pack which should be passed to
  ///                        instance upon creation
  /// @note  Args&& ... is rvalue reference to parameter pack to forward
  ///                   the parameters to the instance
  template <typename... Args>
  static void create(Args &&... args);

  static void on_dead_reference();

  /// Data

  static T *m_instance;

  static bool m_destroyed;
};

// ----------------------------------------------------------------------------==========

// Initialization of static variables:
template <typename T, typename Initializer>
T *Singleton<T, Initializer>::m_instance = 0;

template <typename T, typename Initializer>
bool Singleton<T, Initializer>::m_destroyed = false;

// ----------------------------------------------------------------------------==========

// Member functions:

template <typename T, typename Initializer>
Singleton<T, Initializer>::Singleton(const Singleton<T, Initializer> &source)
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Initializer>
Singleton<T, Initializer>::~Singleton()
{
  m_instance  = 0;
  m_destroyed = true;
}

// ----------------------------------------------------------------------------

// template<typename T, typename Initializer>
// void Singleton<T,Initializer>::create() {
//  static T the_instance;
//  m_instance = & the_instance;
//}

// ----------------------------------------------------------------------------

template <typename T, typename Initializer>
template <typename... Args>
void Singleton<T, Initializer>::create(Args &&... args)
{
  static T the_instance(std::forward<Args>(args)...);
  m_instance = &the_instance;
}

// ----------------------------------------------------------------------------

template <typename T, typename Initializer>
void Singleton<T, Initializer>::on_dead_reference()
{
  // throw std::runtime_error("Dead reference of SingletonShell detected");
  std::cerr << "Dead reference of SingletonShell detected";
}

// ----------------------------------------------------------------------------

// template<typename T,typename Initializer>
// T& Singleton<T,Initializer>::instance()
//{
//  if (!m_instance)
//  {
//    // Check for dead reference:
//    if(m_destroyed) {
//      on_dead_reference();
//    }
//    else
//    {
//      create();
//      Initializer::setup_instance(*m_instance);
//    }
//  }
//  return *m_instance;
//}

// ----------------------------------------------------------------------------

template <typename T, typename Initializer>
template <typename... Args>
T &Singleton<T, Initializer>::instance(Args &&... args)
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
      Initializer::setup_instance(*m_instance);
    }
  }
  return *m_instance;
}

// ----------------------------------------------------------------------------

template <typename T, typename Initializer>
std::shared_ptr<T> Singleton<T, Initializer>::shared_ptr_to_instance()
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
      create();
      Initializer::setup_instance(*m_instance);
    }
  }

  std::shared_ptr<T> ptr_to_instance(m_instance, NullDeleter());
  return ptr_to_instance;
}

// ----------------------------------------------------------------------------

template <typename T, typename Initializer>
void Singleton<T, Initializer>::print()
{
  std::cout << "This is Singleton shell" << std::endl;
}

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit

#endif // Singleton_Shell_hh
