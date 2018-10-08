#ifndef PDEKIT_Linear_System_AMG_PC_hpp
#define PDEKIT_Linear_System_AMG_PC_hpp

#include <iostream>

#include "linear_system/TpetraInternalAccess.hpp"
#include "linear_system/preconditioner/TrilinosPC.hpp"

#if PDEKIT_HAVE_TRILINOS

#include "MueLu.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_TrilinosSmoother.hpp"

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------

namespace internal
{
class TpetraInternalAccess;
}

// ----------------------------------------------------------------------------

template <typename Scalar        = Tpetra::Details::DefaultTypes::scalar_type,
          typename LocalOrdinal  = Tpetra::Details::DefaultTypes::local_ordinal_type,
          typename GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
          typename Node          = Tpetra::Details::DefaultTypes::node_type,
          const bool classic     = Node::classic>
class AMGPC : public TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>
{
  public:
  /// TYPEDEFSpublic:
  typedef TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> matrix_type;

  /// Constructor
  AMGPC();

  /// Destructor
  ~AMGPC();

  /// Create the preconditioner
  void create(const std::string &prec_type, const std::shared_ptr<matrix_type> &A);

  /// Create the preconditioner and pass block structure for creation
  void create(const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
              std::unique_ptr<std::vector<common::Range1D<LocalOrdinal>>> &&blocks);

  /// Create the preconditioner and pass extra data for coarse-scale
  /// correction
  void create(const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
              const common::BlockArray<Uint, Uint> &mesh_dual_graph_crs,
              const std::vector<Uint> &coarse_block_sizes,
              const std::vector<Uint> &fine_block_sizes, const std::vector<Uint> &cell_reordering,
              std::unique_ptr<ls::LocalTransferOps<Real>> &&restriction_ops,
              std::unique_ptr<ls::LocalTransferOps<Real>> &&prolongation_ops);

  /// Prepare the internal data of the preconditioner
  void initialize(const std::shared_ptr<matrix_type> &A);

  /// Compute the preconditioner
  void compute(const std::shared_ptr<matrix_type> &A);

  /// Apply the matrix as linear operator to vector
  /// Computes Y := beta*Y + alpha*Op(A)*X with Op(A) = A or Op(A) =
  /// transpose(A)
  void apply(const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &X,
             TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &Y,
             const bool transp = false, Scalar alpha = 1.0, Scalar beta = 0.0) const;

  /// Return the estimate of condition number
  Real cond_estimate(const TpetraCondestType cdest_type) const;

  private:
  /// Typedef 'inherited' from base class
  typedef typename TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::tpetra_op_type
      tpetra_op_type;

  /// FRIENDS
  friend class internal::TpetraInternalAccess;

  /// Get a const RCP pointer to the preconditioner
  /// The pointer is of type 'pointer to const Tpetra::Operator', which
  /// is a base class of MueLu preconditioner
  virtual Teuchos::RCP<tpetra_op_type const> get_preconditioner_operator() const;

  /// Set the default parameters
  void set_default_parameters();

  /// Recompute the preconditioner
  /// The reference param_list is not const, because this method calls a
  /// Trilinos method which only takes non-const reference ...
  void compute_preconditioner(
      const Teuchos::RCP<
          const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &A,
      Teuchos::ParameterList &param_list);

  /// MueLu preconditioner
  Teuchos::RCP<tpetra_op_type> m_muelu_prec;

  /// List of preconditioner parameters
  Teuchos::ParameterList m_param_list;
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::AMGPC()
    : TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>()
{
  set_default_parameters();
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::~AMGPC()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create(
    const std::string &prec_type, const std::shared_ptr<matrix_type> &A)
{
  typedef internal::TpetraInternalAccess iaccess;

  // --------------------------------------------
  // Verbosity configuration
  // --------------------------------------------

  Teuchos::oblackholestream blackHole;

  const bool print_output = false;
  std::ostream &out       = (print_output) ? std::cout : blackHole;

  // ---------------------------------
  // Prepare preconditioner parameters
  // ---------------------------------
  set_default_parameters();

  // Create timers to show how long it takes for AMGPreconditioner to do
  // various operations.
  Teuchos::RCP<Teuchos::Time> init_timer =
      Teuchos::TimeMonitor::getNewCounter("AMGPreconditioner::Preconditioner::initialize");
  Teuchos::RCP<Teuchos::Time> compute_timer =
      Teuchos::TimeMonitor::getNewCounter("AMGPreconditioner::Preconditioner::compute");

  /*
  if (!print_output)
  {
    Teuchos::TimeMonitor::disableTimer("AMGPreconditioner::Preconditioner::initialize");
    Teuchos::TimeMonitor::disableTimer("AMGPreconditioner::Preconditioner::compute");
  }
  */

  out << "Creating MueLu preconditioner" << std::endl << "-- Configuring" << std::endl;

  // --------------------------------------------
  // Create the preconditioner and set parameters
  // --------------------------------------------
  // RCP<Matrix> A = TpetraCrs_To_XpetraMatrix<SC,LO,GO,NO>(inA);

  compute_preconditioner(iaccess::get_matrix_data(*A), m_param_list);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create(
    const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
    std::unique_ptr<std::vector<common::Range1D<LocalOrdinal>>> &&blocks)
{
  this->create(prec_type, A);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create(
    const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
    const common::BlockArray<Uint, Uint> &mesh_dual_graph_crs,
    const std::vector<Uint> &coarse_block_sizes, const std::vector<Uint> &fine_block_sizes,
    const std::vector<Uint> &cell_reordering,
    std::unique_ptr<ls::LocalTransferOps<Real>> &&restriction_ops,
    std::unique_ptr<ls::LocalTransferOps<Real>> &&prolongation_ops)
{
  this->create(prec_type, A);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::initialize(
    const std::shared_ptr<matrix_type> &A)
{
  using iaccess = internal::TpetraInternalAccess;
  compute_preconditioner(iaccess::get_matrix_data(*A), m_param_list);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::compute(
    const std::shared_ptr<matrix_type> &A)
{
  using iaccess = internal::TpetraInternalAccess;
  compute_preconditioner(iaccess::get_matrix_data(*A), m_param_list);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::apply(
    const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &X,
    TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &Y, const bool transp,
    Scalar alpha, Scalar beta) const
{
  const Teuchos::ETransp mode = (transp == true) ? Teuchos::TRANS : Teuchos::NO_TRANS;

  m_muelu_prec->apply(*internal::TpetraInternalAccess::get_vector_data(X),
                      *internal::TpetraInternalAccess::get_vector_data(Y), mode, alpha, beta);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Real AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::cond_estimate(
    const TpetraCondestType cdest_type) const
{
  return 0.0;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Teuchos::RCP<
    typename AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::tpetra_op_type const>
AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::get_preconditioner_operator() const
{
  // std::cout << "GETTING MUELU PRECONDITIONER OPERATOR" << std::endl;
  return m_muelu_prec;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::set_default_parameters()
{
  m_param_list.set("verbosity", "none"); // none, low, medium, high, extreme
  m_param_list.set("max levels", 10);
  m_param_list.set("coarse: max size", 10);
  m_param_list.set("multigrid algorithm", "pg"); // sa, unsmoothed, pg, emin
  m_param_list.set("problem: symmetric",
                   false); // is the underlying problem symmetric?
  m_param_list.set("number of equations",
                   1); // number of equations at each grid node
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void AMGPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::compute_preconditioner(
    const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>>
        &A,
    Teuchos::ParameterList &param_list)
{
  // --------------------------------------------
  // Create the preconditioner and set parameters
  // --------------------------------------------

  m_muelu_prec.reset();
  /*
  // First const-cast the original pointer 'A' to a pointer to non-const
  object, because
  // that's what MueLu::CreateTpetraPreconditioner expects
  const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal,
  Node, classic>> pA_non_const = Teuchos::rcp_const_cast<
          Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node,
  classic>>(A);

  // Second, generate a pointer to a Tpetra::operator (i.e. parent class of
  Tpetra::CrsMatrix)
  // from pA_non_const
  const Teuchos::RCP<tpetra_op_type> p_opA_non_const = pA_non_const;
  m_muelu_prec = MueLu::CreateTpetraPreconditioner(p_opA_non_const,
  param_list);
  */

  // To compile with newer versions of Trilinos:
  const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>>
      pA_non_const = Teuchos::rcp_const_cast<
          Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>>(A);
  /*
  m_muelu_prec = MueLu::CreateTpetraPreconditioner<Scalar, LocalOrdinal,
  GlobalOrdinal, Node>( pA_non_const, param_list);
  */
  // Second, generate a pointer to a Tpetra::operator (i.e. parent class of
  // Tpetra::CrsMatrix) from pA_non_const
  const Teuchos::RCP<tpetra_op_type> p_opA_non_const = pA_non_const;
  m_muelu_prec = MueLu::CreateTpetraPreconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
      p_opA_non_const, param_list);
}

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_HAVE_TRILINOS

#endif
