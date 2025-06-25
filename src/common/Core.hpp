#ifndef PDEKIT_Common_Environment_hpp
#define PDEKIT_Common_Environment_hpp

#include "common/PDEKit.hpp"
#include "common/Singleton.hpp"

namespace pdekit
{

namespace common
{

class CoreImpl
{
  public:
  /// Default constructor
  /// @param argc         ... number of arguments
  /// @param char* argv[] ... argument values
  CoreImpl(int argc, char *argv[]);

  /// Default destructor
  ~CoreImpl();

  /// Initialize the environment
  void init(int argc, char *argv[]);

  /// Terminate the environment
  void finalize();

  /// Number of threads available
  Uint nb_threads() const;

  /// Maximum number of threads on this hardware
  Uint hardware_concurrency() const;

  /// Set the number of threads
  void set_nb_threads(const Uint nb_threads);

  private:
  /// Command-line arguments count
  int m_argc;
  /// Command-line arguments values
  char **m_argv;

  /// Number of hardware cores
  Uint m_nb_threads;
};

typedef Singleton<CoreImpl> Core;

} // namespace common

} // namespace pdekit

#endif
