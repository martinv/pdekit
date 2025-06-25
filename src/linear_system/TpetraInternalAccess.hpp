#ifndef PDEKIT_Linear_System_Tpetra_Internal_Access_hpp
#define PDEKIT_Linear_System_Tpetra_Internal_Access_hpp

#include "PDEKit_Config.hpp"
#include "linear_system/TpetraFwd.hpp"

#include <iostream>

#if PDEKIT_HAVE_TRILINOS
#include "Tpetra_Operator.hpp"
#endif

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
class TrilinosPC;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
class IfpackPC;

// ----------------------------------------------------------------------------

namespace internal
{

class TpetraInternalAccess
{
  public:
#if PDEKIT_HAVE_TRILINOS
  template <typename Ordinal>
  static Teuchos::RCP<typename TpetraComm<Ordinal>::trilinos_comm_type const> get_tpetra_comm(
      const TpetraComm<Ordinal> &tpetra_comm)
  {
    return tpetra_comm.get_communicator();
  }

  template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  static Teuchos::RCP<
      typename TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node>::trilinos_map_type const>
  get_dof_map(const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> &tpetra_map)
  {
    return tpetra_map.get_map();
  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  static Teuchos::RCP<typename TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                                               classic>::trilinos_matrix_type const>
  get_matrix_data(
      const TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &tpetra_crs_matrix)
  {
    return tpetra_crs_matrix.get_matrix();
  }

  // Get internal data representation from Tpetra multivector
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  static Teuchos::RCP<typename TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                                                 classic>::trilinos_multivector_type>
  get_vector_data(
      ls::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &tpetra_multivector)
  {
    return tpetra_multivector.get_vector();
  }

  // Get internal data representation from Tpetra multivector, const version
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  static Teuchos::RCP<typename TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                                                 classic>::trilinos_multivector_type const>
  get_vector_data(ls::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> const
                      &tpetra_multivector)
  {
    return tpetra_multivector.get_vector();
  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            const bool classic>
  static Teuchos::RCP<
      typename TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::tpetra_op_type const>
  get_preconditioner(
      TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> const &trilinos_prec)
  {
    return trilinos_prec.get_preconditioner_operator();
  }

#endif // PDEKIT_HAVE_TRILINOS
};

} // namespace internal

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_Linear_System_Tpetra_Internal_Access_hpp
