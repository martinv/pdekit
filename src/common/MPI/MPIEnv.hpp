#ifndef PDEKIT_Common_MPI_MPI_Env_hpp
#define PDEKIT_Common_MPI_MPI_Env_hpp

#include "common/MPI/MPITypes.hpp"
#include "common/PDEKit.hpp"
#include "common/Singleton.hpp"

namespace pdekit
{

namespace common
{

namespace mpi
{

class MPIEnvImpl
{
  public:
  /// Default constructor
  MPIEnvImpl();

  /// Constructor which takes arguments to initialize the parallel
  /// environment
  MPIEnvImpl(int argc, char *argv[]);

  /// Constructor, supply additional information about threading level
  MPIEnvImpl(int argc, char *argv[], int required);

  /// Default destructor
  ~MPIEnvImpl();

  /// Initialize the parallel environment after construction
  void init(int argc, char *argv[]);

  /// Initialize the parallel environment after construction
  /// Provide information about required threading level
  void init(int argc, char *argv[], int required);

  /// Finalize the parallel environment
  void finalize();

  /// Return true if the parallel environment is initialized
  bool is_initialized() const;

  /// Return true if the parallel environment is finalized
  bool is_finalized() const;

  /// Return the communicator
  Communicator comm() const;

  /// Get the size of the communicator
  Uint comm_size() const;

  /// Get the rank of this process
  Uint rank() const;

  private:
  /// Communicator
  Communicator m_comm;

  /// Thread level
  int m_thread_level;
};

typedef Singleton<MPIEnvImpl> MPIEnv;

} // namespace mpi

} // namespace common

} // namespace pdekit

#endif
