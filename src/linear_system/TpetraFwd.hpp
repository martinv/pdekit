#ifndef PDEKIT_Linear_System_Tpetra_Fwd_hpp
#define PDEKIT_Linear_System_Tpetra_Fwd_hpp

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------
// Forward declarations

template <typename Ordinal>
class TpetraComm;

template <class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraDofMap;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
class TpetraMultiVector;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
class TpetraCrsMatrix;

// ----------------------------------------------------------------------------
} // namespace ls

} // namespace pdekit

#endif
