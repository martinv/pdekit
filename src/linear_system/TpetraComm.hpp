#ifndef PDEKIT_Trilinos_Tpetra_Comm_hpp
#define PDEKIT_Trilinos_Tpetra_Comm_hpp

#include <iostream>

#include "PDEKit_Config.hpp"
#include "common/MPI/MPITypes.hpp"

#if PDEKIT_HAVE_TRILINOS

#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_RCP.hpp"

#endif

namespace pdekit
{

namespace ls
{

namespace internal
{
class TpetraInternalAccess;
}

template <typename Ordinal>
class TpetraComm
{
#if PDEKIT_HAVE_TRILINOS
  private:
  using trilinos_comm_type = Teuchos::Comm<Ordinal>;
#endif

  public:
  /// Default constructor
  TpetraComm();

/// Construct from pointer to the native (wrapped) communicator
#if PDEKIT_HAVE_TRILINOS
  TpetraComm(const Teuchos::RCP<const Teuchos::Comm<Ordinal>> &comm);
#endif

  /// Build Tpetra communicator from MPI communicator
  TpetraComm(const MPI_Comm communicator);

  /// Copy constructor
  TpetraComm(const TpetraComm &rhs);

  /// Assignment operator
  TpetraComm &operator=(const TpetraComm &rhs);

  /// Destructor
  ~TpetraComm();

  private:
  /// FRIENDS
  friend class internal::TpetraInternalAccess;

#if PDEKIT_HAVE_TRILINOS
  Teuchos::RCP<const Teuchos::Comm<Ordinal>> get_communicator() const;

  /// Trilinos communicator object
  /// Takes MPI communicator in constructor
  Teuchos::RCP<const Teuchos::Comm<Ordinal>> m_comm;
#endif
};

// ----------------------------------------------------------------------------

template <typename Ordinal>
TpetraComm<Ordinal>::TpetraComm()
{
#if PDEKIT_HAVE_TRILINOS
  m_comm = Teuchos::RCP<const Teuchos::Comm<Ordinal>>(nullptr);
#endif
}

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Ordinal>
TpetraComm<Ordinal>::TpetraComm(const Teuchos::RCP<const Teuchos::Comm<Ordinal>> &comm)
{
  m_comm = comm;
}
#endif

// ----------------------------------------------------------------------------

template <typename Ordinal>
TpetraComm<Ordinal>::TpetraComm(const MPI_Comm communicator)
{
#if PDEKIT_HAVE_TRILINOS
  m_comm = Teuchos::RCP<const Teuchos::Comm<Ordinal>>(new Teuchos::MpiComm<Ordinal>(communicator));
#endif
}

// ----------------------------------------------------------------------------

template <typename Ordinal>
TpetraComm<Ordinal>::TpetraComm(const TpetraComm &rhs)
{
#if PDEKIT_HAVE_TRILINOS
  m_comm = rhs.m_comm;
#endif
}

// ----------------------------------------------------------------------------

template <typename Ordinal>
TpetraComm<Ordinal> &TpetraComm<Ordinal>::operator=(const TpetraComm &rhs)
{
#if PDEKIT_HAVE_TRILINOS
  m_comm = rhs.m_comm;
#endif
  return *this;
}

// ----------------------------------------------------------------------------

template <typename Ordinal>
TpetraComm<Ordinal>::~TpetraComm()
{
}

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal>> TpetraComm<Ordinal>::get_communicator() const
{
  return m_comm;
}
#endif

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif
