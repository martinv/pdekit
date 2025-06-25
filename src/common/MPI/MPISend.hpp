#ifndef PDEKIT_Common_MPI_MPI_SEND_HPP
#define PDEKIT_Common_MPI_MPI_SEND_HPP

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
class Send;

template <typename T>
class Send<std::vector<T>>
{
  public:
  static int execute(const std::vector<T> &data, int destination, int tag, MPI_Comm communicator)
  {
    int size = data.size();
    MPI_Send(&size, 1, MPI_INT, destination, tag, communicator);
    return MPI_Send(data.data(), size, MPITypes<T>::type, destination, tag, communicator);
  }
};

} // namespace internal

// Free wrapper function around MPI send - calls internal implementation

template <typename T>
int send(const T &data, int destination, int tag, MPI_Comm communicator)
{
  return internal::Send<T>::execute(data, destination, tag, communicator);
}

} // namespace mpi

} // namespace common

} // namespace pdekit

#endif
