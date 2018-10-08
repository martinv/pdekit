#ifndef PDEKIT_Common_MPI_MPI_RECV_HPP
#define PDEKIT_Common_MPI_MPI_RECV_HPP

#include <vector>

#include "common/MPI/MPITypes.hpp"
#include "common/PDEKit.hpp"

namespace pdekit
{

namespace common
{

namespace mpi
{

namespace internal
{

template <typename T>
class Recv;

template <typename T>
class Recv<std::vector<T>>
{
  public:
  static int execute(std::vector<T> &data, int source, int tag, MPI_Comm communicator)
  {
    int size = 0;
    MPI_Status status;
    MPI_Recv(&size, 1, MPI_INT, source, tag, communicator, &status);
    data.resize(size);
    return MPI_Recv(data.data(), size, MPITypes<T>::type, source, tag, communicator, &status);
  }
};

} // namespace internal

// Free wrapper function around MPI send - calls internal implementation

template <typename T>
int recv(T &data, int source, int tag, MPI_Comm communicator)
{
  return internal::Recv<T>::execute(data, source, tag, communicator);
}

} // namespace mpi

} // namespace common

} // namespace pdekit

#endif
