#include "common/MPI/MPITypes.hpp"
#include "common/PDEKit.hpp"

namespace pdekit
{

namespace common
{

namespace mpi
{

// Default MPI type - this should never be used
template <typename T>
Datatype MPITypes<T>::type = MPI_DATATYPE_NULL;

template <>
Datatype MPITypes<Char>::type = MPI_UNSIGNED_CHAR;

template <>
Datatype MPITypes<SUint>::type = MPI_SHORT;

template <>
Datatype MPITypes<Int>::type = MPI_INT;

template <>
Datatype MPITypes<Uint>::type = MPI_UNSIGNED;

template <>
Datatype MPITypes<Lint>::type = MPI_LONG_INT;

template <>
Datatype MPITypes<LUint>::type = MPI_UNSIGNED_LONG;

// #ifdef PDEKIT_HAVE_MPI_LONG_LONG
// template <>
// Datatype MPITypes<LLint>::type = MPI_LONG_LONG_INT;
// template <>
// Datatype MPITypes<LLUint>::type = MPI_UNSIGNED_LONG_LONG;
// #endif

template <>
Datatype MPITypes<Float>::type = MPI_FLOAT;

template <>
Datatype MPITypes<Real>::type = MPI_DOUBLE;

} // namespace mpi

} // namespace common

} // namespace pdekit
