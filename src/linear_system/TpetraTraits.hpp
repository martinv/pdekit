#ifndef PDEKIT_Trilinos_Tpetra_Traits_hpp
#define PDEKIT_Trilinos_Tpetra_Traits_hpp

#include "PDEKit_Config.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "Tpetra_ConfigDefs.hpp"
#endif

namespace pdekit
{

namespace ls
{

struct VoidTpetraNodeType
{
};

struct TpetraDefaultTraits
{
#if PDEKIT_HAVE_TRILINOS
  using Scalar              = Tpetra::Details::DefaultTypes::scalar_type;
  using LocalOrdinal        = Tpetra::Details::DefaultTypes::local_ordinal_type;
  using GlobalOrdinal       = Tpetra::Details::DefaultTypes::global_ordinal_type;
  using Node                = Tpetra::Details::DefaultTypes::node_type;
  static const bool classic = Node::classic;
#else
  using Scalar              = Real;
  using LocalOrdinal        = Int;
  using GlobalOrdinal       = Int;
  using Node                = VoidTpetraNodeType;
  static const bool classic = false;
#endif
};

} // namespace ls

} // namespace pdekit

#endif
