#ifndef PDEKIT_Linear_System_Tpetra_Block_Diag_PC_hpp
#define PDEKIT_Linear_System_Tpetra_Block_Diag_PC_hpp

#include <iostream>

#include "linear_system/preconditioner/TpetraBlockDiagPCImpl.hpp"

#if PDEKIT_HAVE_TRILINOS

#include "Ifpack2_Factory.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_oblackholestream.hpp"

namespace pdekit
{

namespace ls
{

template <typename Scalar        = Tpetra::Details::DefaultTypes::scalar_type,
          typename LocalOrdinal  = Tpetra::Details::DefaultTypes::local_ordinal_type,
          typename GlobalOrdinal = Tpetra::Details::DefaultTypes::global_ordinal_type,
          typename Node          = Tpetra::Details::DefaultTypes::node_type,
          const bool classic     = Node::classic>
class TpetraBlockDiagPC : public TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>
{
  public:
  /// TYPEDEFSpublic:
  typedef TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> matrix_type;

  /// Constructor
  TpetraBlockDiagPC();

  /// Destructor
  ~TpetraBlockDiagPC() override;

  /// Create the preconditioner
  void create(const std::string &prec_type, const std::shared_ptr<matrix_type> &A) override;

  /// Create the preconditioner and pass block structure for creation
  void create(const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
              std::unique_ptr<std::vector<common::Range1D<LocalOrdinal>>> &&blocks) override;

  /// Create the preconditioner and pass extra data for coarse-scale
  /// correction
  void create(const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
              const common::BlockArray<Uint, Uint> &mesh_dual_graph_crs,
              const std::vector<Uint> &coarse_block_sizes,
              const std::vector<Uint> &fine_block_sizes, const std::vector<Uint> &cell_reordering,
              std::unique_ptr<ls::LocalTransferOps<Real>> &&restriction_ops,
              std::unique_ptr<ls::LocalTransferOps<Real>> &&prolongation_ops) override;

  /// Prepare the internal data of the preconditioner
  void initialize(const std::shared_ptr<matrix_type> &A) override;

  /// Compute the preconditioner
  void compute(const std::shared_ptr<matrix_type> &A) override;

  /// Apply the matrix as linear operator to vector
  /// Computes Y := beta*Y + alpha*Op(A)*X with Op(A) = A or Op(A) =
  /// transpose(A)
  void apply(const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &X,
             TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &Y,
             const bool transp = false, Scalar alpha = 1.0, Scalar beta = 0.0) const override;

  /// Return the estimate of condition number
  Real cond_estimate(const TpetraCondestType cdest_type) const override;

  /// Print the contents of the preconditioner
  void print() const;

  private:
  /// Typedef 'inherited' from base class
  using tpetra_op_type =
      typename TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::tpetra_op_type;

  using block_precond_impl_type =
      typename detail::TpetraBlockDiagPCImpl<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  /// FRIENDS
  friend class internal::TpetraInternalAccess;

  /// Get a const RCP pointer to the preconditioner
  Teuchos::RCP<tpetra_op_type const> get_preconditioner_operator() const override;

  /// Set default parameters for the preconditioner of type ILU
  void set_default_parameters_block_diag(Teuchos::ParameterList &param_list) const;

  /// Create block-diagonal preconditioner preconditioner
  /// @param param-list ... settings for the preconditioner
  /// @param A          ... matrix that should be preconditioned
  void create_block_diag_preconditioner(
      const Teuchos::ParameterList &param_list,
      const Teuchos::RCP<
          const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &A);

  /// The actual preconditioner implementation
  Teuchos::RCP<block_precond_impl_type> m_block_diag_prec;
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraBlockDiagPC()
    : TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::~TpetraBlockDiagPC()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create(
    const std::string &prec_type, const std::shared_ptr<matrix_type> &A)
{

  using iaccess = internal::TpetraInternalAccess;

  Teuchos::ParameterList param_list;

  set_default_parameters_block_diag(param_list);
  create_block_diag_preconditioner(param_list, iaccess::get_matrix_data(*A));
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create(
    const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
    std::unique_ptr<std::vector<common::Range1D<LocalOrdinal>>> &&blocks)
{
  using iaccess = internal::TpetraInternalAccess;

  Teuchos::ParameterList param_list;

  set_default_parameters_block_diag(param_list);
  create_block_diag_preconditioner(param_list, iaccess::get_matrix_data(*A));

  m_block_diag_prec->initialize_structure(iaccess::get_matrix_data(*A), std::move(blocks));
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create(
    const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
    const common::BlockArray<Uint, Uint> &mesh_dual_graph_crs,
    const std::vector<Uint> &coarse_block_sizes, const std::vector<Uint> &fine_block_sizes,
    const std::vector<Uint> &cell_reordering,
    std::unique_ptr<ls::LocalTransferOps<Real>> &&restriction_ops,
    std::unique_ptr<ls::LocalTransferOps<Real>> &&prolongation_ops)
{
  using iaccess = internal::TpetraInternalAccess;

  std::unique_ptr<std::vector<common::Range1D<LocalOrdinal>>> block_ranges(
      new std::vector<common::Range1D<LocalOrdinal>>());

  block_ranges->resize(fine_block_sizes.size());

  LocalOrdinal block_offset = LocalOrdinal();
  for (Uint b = 0; b < fine_block_sizes.size(); ++b)
  {
    (*block_ranges)[b] =
        common::Range1D<LocalOrdinal>(block_offset, block_offset + fine_block_sizes[b] - 1);
    block_offset += fine_block_sizes[b];
  }

  Teuchos::ParameterList param_list;

  set_default_parameters_block_diag(param_list);
  create_block_diag_preconditioner(param_list, iaccess::get_matrix_data(*A));

  m_block_diag_prec->initialize_structure(iaccess::get_matrix_data(*A), std::move(block_ranges));
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::initialize(
    const std::shared_ptr<matrix_type> &A)
{
  m_block_diag_prec->initialize();
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::compute(
    const std::shared_ptr<matrix_type> &A)
{
  m_block_diag_prec->compute_impl(A);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::apply(
    const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &X,
    TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &Y, const bool transp,
    Scalar alpha, Scalar beta) const
{
  const Teuchos::ETransp mode = (transp == true) ? Teuchos::TRANS : Teuchos::NO_TRANS;

  m_block_diag_prec->apply(*internal::TpetraInternalAccess::get_vector_data(X),
                           *internal::TpetraInternalAccess::get_vector_data(Y), mode, alpha, beta);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Real TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::cond_estimate(
    const TpetraCondestType cdest_type) const
{
  return 0.0;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::print() const
{
  m_block_diag_prec->print_impl();
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Teuchos::RCP<typename TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node,
                                        classic>::tpetra_op_type const>
TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::get_preconditioner_operator()
    const
{
  return m_block_diag_prec;
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
    set_default_parameters_block_diag(Teuchos::ParameterList &param_list) const
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraBlockDiagPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
    create_block_diag_preconditioner(
        const Teuchos::ParameterList &param_list,
        const Teuchos::RCP<
            const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &A)
{
  // --------------------------------------------
  // Verbosity configuration
  // --------------------------------------------

  Teuchos::oblackholestream blackHole;

  const bool print_output = false;
  std::ostream &out       = (print_output) ? std::cout : blackHole;

  // -----------------------------
  // Preconditioner initialization
  // -----------------------------

  // Create timers to show how long it takes for Ifpack2 to do various
  // operations.
  Teuchos::RCP<Teuchos::Time> init_timer =
      Teuchos::TimeMonitor::getNewCounter("BlockDiagonal::Preconditioner::initialize");
  Teuchos::RCP<Teuchos::Time> compute_timer =
      Teuchos::TimeMonitor::getNewCounter("BlockDiagonal::Preconditioner::compute");

  /*
  if (!print_output)
  {
    Teuchos::TimeMonitor::disableTimer("Ifpack2::Preconditioner::initialize");
    Teuchos::TimeMonitor::disableTimer("Ifpack2::Preconditioner::compute");
  }
  */

  out << "Creating BlockDiagonal preconditioner" << std::endl << "-- Configuring" << std::endl;

  // --------------------------------------------
  // Create the preconditioner and set parameters
  // --------------------------------------------

  // This doesn't actually _compute_ the preconditioner. It just sets up
  // the specific type of preconditioner and its associated parameters
  // (which depend on the type).

  m_block_diag_prec.reset();
  m_block_diag_prec = Teuchos::RCP<block_precond_impl_type>(
      new block_precond_impl_type(A->getGlobalNumRows(), A->getComm()));
  m_block_diag_prec->setParameters(param_list);

  // Initialize the preconditioner.  If the sparse matrix's _structure_
  // changes, you have to call initialize() and compute() again, in that
  // sequence, before you may use ("apply") the preconditioner.

  out << "-- Initializing" << std::endl;
  {
    Teuchos::TimeMonitor mon(*init_timer);
    m_block_diag_prec->initialize();
  }
}

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_HAVE_TRILINOS

#endif
