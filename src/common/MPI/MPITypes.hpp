#ifndef PDEKIT_Common_MPI_MPI_Types_hpp
#define PDEKIT_Common_MPI_MPI_Types_hpp

#include "mpi.h"

namespace pdekit
{

namespace common
{

namespace mpi
{

/// communicator
using Communicator = MPI_Comm;

/// operation (mostly for reduce and all_reduce)
using Operation = MPI_Op;

/// datatype
using Datatype = MPI_Datatype;

/// Wrapper to translate between MPI types and PDEKit types

template <typename T>
class MPITypes
{
  public:
  static Datatype type;
};

} // namespace mpi

} // namespace common

} // namespace pdekit

#endif
