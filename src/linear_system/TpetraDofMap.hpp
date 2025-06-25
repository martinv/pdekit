#ifndef PDEKIT_Trilinos_Tpetra_Dof_Map_hpp
#define PDEKIT_Trilinos_Tpetra_Dof_Map_hpp

#include <iostream>

#include "PDEKit_Config.hpp"
#include "common/MPI/MPIEnv.hpp"
#include "common/PDEKit.hpp"
#include "graph/Graph.hpp"
#include "linear_system/TpetraComm.hpp"
#include "linear_system/TpetraInternalAccess.hpp"
#include "linear_system/TpetraTraits.hpp"
#include "math/DenseConstVecView.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "Tpetra_Map.hpp"
#endif

namespace pdekit
{

namespace ls
{

namespace internal
{
class TpetraInternalAccess;
}

template <typename LocalOrdinal  = TpetraDefaultTraits::LocalOrdinal,
          typename GlobalOrdinal = TpetraDefaultTraits::GlobalOrdinal,
          typename Node          = TpetraDefaultTraits::Node>
class TpetraDofMap
{
  private:
#if PDEKIT_HAVE_TRILINOS
  using trilinos_map_type = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
#endif

  public:
  /// Default constructor
  TpetraDofMap();

/// Construct from a Teuchos RCP to Tpetra map
#if PDEKIT_HAVE_TRILINOS
  TpetraDofMap(const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map_rcp);
#endif

  /// Construct
  /// @param num_global_elements ... total number of elements in map
  /// @param index_base          ... offset for map entries
  /// @param tpetra_comm         ... communicator
  TpetraDofMap(size_t num_global_elements, GlobalOrdinal index_base,
               const TpetraComm<GlobalOrdinal> &tpetra_comm);

  /// Copy constructor
  TpetraDofMap(const TpetraDofMap &rhs);

  /// Assignment operator
  TpetraDofMap &operator=(const TpetraDofMap &rhs);

  /// Destructor
  ~TpetraDofMap();

  /// Reset the map
  void reset(size_t num_global_elements, GlobalOrdinal index_base,
             const TpetraComm<GlobalOrdinal> &tpetra_comm);

  /// Get the global number of elements
  Uint global_num_elements() const;

  private:
#if PDEKIT_HAVE_TRILINOS
  /// FRIENDS
  friend class internal::TpetraInternalAccess;

  Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> get_map() const;

  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> tpetra_map_type;

  Teuchos::RCP<tpetra_map_type const> m_map;
#endif
};

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::TpetraDofMap() : m_map(nullptr)
{
}
#else
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::TpetraDofMap()
{
}
#endif

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::TpetraDofMap(
    const Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> &map_rcp)
    : m_map(map_rcp)
{
}
#endif

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::TpetraDofMap(
    size_t num_global_elements, GlobalOrdinal index_base,
    const TpetraComm<GlobalOrdinal> &tpetra_comm)
    : m_map(new tpetra_map_type(num_global_elements, index_base,
                                internal::TpetraInternalAccess::get_tpetra_comm(tpetra_comm)))
{
}
#else
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::TpetraDofMap(
    size_t num_global_elements, GlobalOrdinal index_base,
    const TpetraComm<GlobalOrdinal> &tpetra_comm)
{
}
#endif

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::TpetraDofMap(const TpetraDofMap &rhs)
    : m_map(rhs.m_map)
{
}
#else
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::TpetraDofMap(const TpetraDofMap &rhs)
{
}
#endif

// ----------------------------------------------------------------------------

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>
    &TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::operator=(const TpetraDofMap &rhs)
{
#if PDEKIT_HAVE_TRILINOS
  m_map = rhs.m_map;
  return *this;
#endif
}

// ----------------------------------------------------------------------------

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::~TpetraDofMap()
{
}

// ----------------------------------------------------------------------------

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::reset(
    size_t num_global_elements, GlobalOrdinal index_base,
    const TpetraComm<GlobalOrdinal> &tpetra_comm)
{
#if PDEKIT_HAVE_TRILINOS
  m_map.reset(new tpetra_map_type(num_global_elements, index_base,
                                  internal::TpetraInternalAccess::get_tpetra_comm(tpetra_comm)));
#endif
}

// ----------------------------------------------------------------------------

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Uint TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::global_num_elements() const
{
#if PDEKIT_HAVE_TRILINOS
  return static_cast<Uint>(m_map->getGlobalNumElements());
#else
  return 0;
#endif
}

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> const> TpetraDofMap<
    LocalOrdinal, GlobalOrdinal, Node>::get_map() const
{
  return m_map;
}
#endif

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif
