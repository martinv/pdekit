#include "common/MPI/MPIEnv.hpp"

namespace pdekit
{

namespace common
{

namespace mpi
{

// ----------------------------------------------------------------------------

MPIEnvImpl::MPIEnvImpl() : m_comm(nullptr), m_thread_level(MPI_THREAD_SINGLE)
{
}

// ----------------------------------------------------------------------------

MPIEnvImpl::MPIEnvImpl(int argc, char *argv[])
{
  m_comm = nullptr;
  init(argc, argv);
}

// ----------------------------------------------------------------------------

MPIEnvImpl::MPIEnvImpl(int argc, char *argv[], int required)
{
  m_comm = nullptr;
  init(argc, argv, required);
}

// ----------------------------------------------------------------------------

MPIEnvImpl::~MPIEnvImpl()
{
}

// ----------------------------------------------------------------------------

void MPIEnvImpl::init(int argc, char *argv[])
{
  if (is_finalized())
  {
    std::cerr << "Should not initialize MPI environment after call to "
                 "'finalize()'"
              << std::endl;
    return;
  }

  if (!is_initialized()) // then initialize
  {
    MPI_Init(&argc, &argv);
    // MPI_CHECK_RESULT(MPI_Init,(&argc,&argv));
  }

  m_comm = MPI_COMM_WORLD;
}

// ----------------------------------------------------------------------------

void MPIEnvImpl::init(int argc, char *argv[], int required)
{
  if (is_finalized())
  {
    std::cerr << "Should not initialize MPI environment after call to "
                 "'finalize()'"
              << std::endl;
    return;
  }

  if (!is_initialized()) // then initialize
  {
    MPI_Init_thread(&argc, &argv, required, &m_thread_level);
  }

  m_comm = MPI_COMM_WORLD;
}

// ----------------------------------------------------------------------------

void MPIEnvImpl::finalize()
{
  if (is_initialized() && !is_finalized()) // then finalize
  {
    MPI_Finalize();
  }

  m_comm = nullptr;
}

// ----------------------------------------------------------------------------

bool MPIEnvImpl::is_initialized() const
{
  int is_initialized = 0;
  MPI_Initialized(&is_initialized);
  return bool(is_initialized);

  // int is_initialized = 0;
  // MPI_CHECK_RESULT(MPI_Initialized,(&is_initialized));
  // return bool(is_initialized);
}

// ----------------------------------------------------------------------------

bool MPIEnvImpl::is_finalized() const
{
  int is_finalized = 0;
  MPI_Finalized(&is_finalized);
  return bool(is_finalized);
}

// ----------------------------------------------------------------------------

Communicator MPIEnvImpl::comm() const
{
  return m_comm;
}

// ----------------------------------------------------------------------------

Uint MPIEnvImpl::comm_size() const
{
  Int comm_size;
  MPI_Comm_size(m_comm, &comm_size);
  return comm_size;
}

// ----------------------------------------------------------------------------

Uint MPIEnvImpl::rank() const
{
  Int rank;
  MPI_Comm_rank(m_comm, &rank);
  return rank;
}

// ----------------------------------------------------------------------------

} // namespace mpi

} // namespace common

} // namespace pdekit
