#include <string>
#include <thread>

#include "common/Core.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

CoreImpl::CoreImpl(int argc, char *argv[]) : m_nb_threads(1u)
{
  m_argc = argc;
  m_argv = argv;

  std::string arguments;

  for (int i = 1; i < (argc - 1); ++i)
  {
    arguments = arguments + argv[i] + " ";
  }

  if (argc > 1)
  {
    arguments = arguments + argv[argc - 1];
  }

  std::cout << "PDEKIT: Initiating environment with arguments [ " << arguments << " ]" << std::endl;
}

// ----------------------------------------------------------------------------

CoreImpl::~CoreImpl()
{
  finalize();
}

// ----------------------------------------------------------------------------

void CoreImpl::init(int argc, char *argv[])
{
  m_argc = argc;
  m_argv = argv;
}

// ----------------------------------------------------------------------------

void CoreImpl::finalize()
{
  std::cout << "PDEKIT: Terminating environment" << std::endl;
}

// ----------------------------------------------------------------------------

Uint CoreImpl::nb_threads() const
{
  return m_nb_threads;
}

// ----------------------------------------------------------------------------

Uint CoreImpl::hardware_concurrency() const
{
  return std::thread::hardware_concurrency();
}

// ----------------------------------------------------------------------------

void CoreImpl::set_nb_threads(const Uint nb_threads)
{
  m_nb_threads = nb_threads;
}

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit
