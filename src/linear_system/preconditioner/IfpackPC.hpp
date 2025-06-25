#ifndef PDEKIT_Linear_System_Ifpack_PC_hpp
#define PDEKIT_Linear_System_Ifpack_PC_hpp

#include <iostream>

#include "linear_system/TpetraInternalAccess.hpp"
#include "linear_system/preconditioner/TrilinosPC.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "Ifpack2_CondestType.hpp"
#include "Ifpack2_Factory.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_oblackholestream.hpp"
#endif

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

template <typename Scalar        = TpetraDefaultTraits::Scalar,
          typename LocalOrdinal  = TpetraDefaultTraits::LocalOrdinal,
          typename GlobalOrdinal = TpetraDefaultTraits::GlobalOrdinal,
          typename Node          = TpetraDefaultTraits::Node,
          const bool classic     = TpetraDefaultTraits::classic>
class IfpackPC : public TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>
{
  public:
  /// TYPEDEFSpublic:
  typedef TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> matrix_type;

  /// Constructor
  IfpackPC();

  /// Destructor
  ~IfpackPC() override;

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

  private:
#if PDEKIT_HAVE_TRILINOS
  /// Typedef 'inherited' from base class
  using tpetra_op_type =
      typename TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::tpetra_op_type;

  // Ifpack2's generic Preconditioner interface.
  // An Ifpack2::Preconditioner "is a" Tpetra::Operator.
  using ifpack_prec_type =
      typename Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  /// FRIENDS
  friend class internal::TpetraInternalAccess;

  /// Get a const RCP pointer to the preconditioner
  /// The pointer is of type 'pointer to const Tpetra::Operator', which
  /// is a base class of Ifpack preconditioner
  Teuchos::RCP<tpetra_op_type const> get_preconditioner_operator() const override;

  /// Set default parameters for the preconditioner of type ILU
  void set_default_parameters_ilut(Teuchos::ParameterList &param_list) const;

  /// Set default parameters for relaxation preconditioner (Jacobi,
  /// Gauss-Seidel)
  void set_default_parameters_relaxation(Teuchos::ParameterList &param_list) const;

  /// Create ILUT preconditioner
  /// @param param-list ... settings for the preconditioner
  /// @param A          ... matrix that should be preconditioned
  void create_ilut_preconditioner(
      const Teuchos::ParameterList &param_list,
      const Teuchos::RCP<
          const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &A);

  /// Create relaxation (Jacobi, Gauss-Seidel) preconditioner
  /// @param param-list ... settings for the preconditioner
  /// @param A          ... matrix that should be preconditioned
  void create_relaxation_preconditioner(
      const Teuchos::ParameterList &param_list,
      const Teuchos::RCP<
          const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> &A);

  /// Ifpack preconditioner
  Teuchos::RCP<ifpack_prec_type> m_ifpack_prec;
#endif
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::IfpackPC()
    : TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::~IfpackPC()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create(
    const std::string &prec_type, const std::shared_ptr<matrix_type> &A)
{
#if PDEKIT_HAVE_TRILINOS
  // Set up the preconditioner of the given type.
  // The preconditioner type parameter can have the following values:
  //
  //    "DIAGONAL"
  //    "RELAXATION"
  //    "CHEBYSHEV"
  //    "ILUT"
  //    "RILUK"
  //    "RBILUK"

  typedef internal::TpetraInternalAccess iaccess;

  Teuchos::ParameterList ifpack_param_list;

  if (prec_type == "RELAXATION")
  {
    set_default_parameters_relaxation(ifpack_param_list);
    create_relaxation_preconditioner(ifpack_param_list, iaccess::get_matrix_data(*A));
  }
  else if (prec_type == "ILUT")
  {
    set_default_parameters_ilut(ifpack_param_list);
    create_ilut_preconditioner(ifpack_param_list, iaccess::get_matrix_data(*A));
  }
  else
  {
    std::cerr << "IfpackPreconditioner::create::don't know how to create "
              << "preconditioner of type '" << prec_type << "'" << std::endl;
  }
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create(
    const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
    std::unique_ptr<std::vector<common::Range1D<LocalOrdinal>>> &&blocks)
{
  this->create(prec_type, A);
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create(
    const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
    const common::BlockArray<Uint, Uint> &mesh_dual_graph_crs,
    const std::vector<Uint> &coarse_block_sizes, const std::vector<Uint> &fine_block_sizes,
    const std::vector<Uint> &cell_reordering,
    std::unique_ptr<ls::LocalTransferOps<Real>> &&restriction_ops,
    std::unique_ptr<ls::LocalTransferOps<Real>> &&prolongation_ops)
{
#if PDEKIT_HAVE_TRILINOS
  this->create(prec_type, A);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::initialize(
    const std::shared_ptr<matrix_type> &A)
{
#if PDEKIT_HAVE_TRILINOS
  m_ifpack_prec->initialize();
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::compute(
    const std::shared_ptr<matrix_type> &A)
{
#if PDEKIT_HAVE_TRILINOS
  m_ifpack_prec->compute();
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::apply(
    const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &X,
    TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &Y, const bool transp,
    Scalar alpha, Scalar beta) const
{
#if PDEKIT_HAVE_TRILINOS
  const Teuchos::ETransp mode = (transp == true) ? Teuchos::TRANS : Teuchos::NO_TRANS;

  m_ifpack_prec->apply(*internal::TpetraInternalAccess::get_vector_data(X),
                       *internal::TpetraInternalAccess::get_vector_data(Y), mode, alpha, beta);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Real IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::cond_estimate(
    const TpetraCondestType cdest_type) const
{
#if PDEKIT_HAVE_TRILINOS
  /*
  if (m_ifpack_prec.is_null())
  {
    return 0.0;
  }

  if (cdest_type == TpetraCondestType::Cheap)
  {
    return m_ifpack_prec-computeCondEst(Ifpack2::Cheap);
  }
  else if (cdest_type == TpetraCondestType::CG)
  {
    return m_ifpack_prec->computeCondEst(Ifpack2::CG);
  }
  else if (cdest_type == TpetraCondestType::GMRES)
  {
    return m_ifpack_prec->computeCondEst(Ifpack2::GMRES);
  }
  */
  return 0.0;
#else
  return 0.0;
#endif
}

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Teuchos::RCP<
    typename IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::tpetra_op_type const>
IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::get_preconditioner_operator() const
{
  return m_ifpack_prec;
}
#endif

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::set_default_parameters_ilut(
    Teuchos::ParameterList &param_list) const
{
  // Ifpack2 expects arguments of type 'double' here, regardless of
  // the scalar or magnitude types of the entries of the sparse
  // matrix.
  const double fill_level = 1.0;
  const double drop_tol   = 1.e-9;
  // const double absThreshold = 0.1; //  What does this do???

  param_list.set("fact: ilut level-of-fill", fill_level);
  param_list.set("fact: drop tolerance", drop_tol);
}
#endif

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node,
              classic>::set_default_parameters_relaxation(Teuchos::ParameterList &param_list) const
{
  // Relaxation types:
  // Jacobi
  // Gauss-Seidel
  // Symmetric Gauss-Seidel

  const std::string relaxation_type = "Gauss-Seidel";
  const int nb_sweeps               = 10;
  const double damping_factor       = 0.5;
  // Used for Gauss-Seidel:
  const bool backward_mode = false;

  param_list.set("relaxation: type", relaxation_type);
  param_list.set("relaxation: sweeps", nb_sweeps);
  param_list.set("relaxation: damping factor", damping_factor);
  param_list.set("relaxation: backward mode", backward_mode);
}
#endif

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create_ilut_preconditioner(
    const Teuchos::ParameterList &ifpack_param_list,
    const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>>
        &A)
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
      Teuchos::TimeMonitor::getNewCounter("Ifpack2::Preconditioner::initialize");
  Teuchos::RCP<Teuchos::Time> compute_timer =
      Teuchos::TimeMonitor::getNewCounter("Ifpack2::Preconditioner::compute");

  /*
  if (!print_output)
  {
    Teuchos::TimeMonitor::disableTimer("Ifpack2::Preconditioner::initialize");
    Teuchos::TimeMonitor::disableTimer("Ifpack2::Preconditioner::compute");
  }
  */

  out << "Creating ILUT preconditioner" << std::endl << "-- Configuring" << std::endl;

  // --------------------------------------------
  // Create the preconditioner and set parameters
  // --------------------------------------------

  // This doesn't actually _compute_ the preconditioner. It just sets up
  // the specific type of preconditioner and its associated parameters
  // (which depend on the type).

  Ifpack2::Factory factory;

  m_ifpack_prec.reset();
  m_ifpack_prec = factory.create("ILUT", A);

  m_ifpack_prec->setParameters(ifpack_param_list);

  // Initialize the preconditioner.  If the sparse matrix's _structure_
  // changes, you have to call initialize() and compute() again, in that
  // sequence, before you may use ("apply") the preconditioner.

  out << "-- Initializing" << std::endl;
  {
    Teuchos::TimeMonitor mon(*init_timer);
    m_ifpack_prec->initialize();
  }
}
#endif

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void IfpackPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::create_relaxation_preconditioner(
    const Teuchos::ParameterList &ifpack_param_list,
    const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>>
        &A)
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
      Teuchos::TimeMonitor::getNewCounter("Ifpack2::Preconditioner::initialize");
  Teuchos::RCP<Teuchos::Time> compute_timer =
      Teuchos::TimeMonitor::getNewCounter("Ifpack2::Preconditioner::compute");

  /*
  if (!print_output)
  {
    Teuchos::TimeMonitor::disableTimer("Ifpack2::Preconditioner::initialize");
    Teuchos::TimeMonitor::disableTimer("Ifpack2::Preconditioner::compute");
  }
  */

  out << "Creating RELAXATION preconditioner" << std::endl << "-- Configuring" << std::endl;

  // --------------------------------------------
  // Create the preconditioner and set parameters
  // --------------------------------------------

  // This doesn't actually _compute_ the preconditioner. It just sets up
  // the specific type of preconditioner and its associated parameters
  // (which depend on the type).

  Ifpack2::Factory factory;

  m_ifpack_prec.reset();
  m_ifpack_prec = factory.create("RELAXATION", A);

  m_ifpack_prec->setParameters(ifpack_param_list);

  // Initialize the preconditioner.  If the sparse matrix's _structure_
  // changes, you have to call initialize() and compute() again, in that
  // sequence, before you may use ("apply") the preconditioner.

  out << "-- Initializing" << std::endl;
  {
    Teuchos::TimeMonitor mon(*init_timer);
    m_ifpack_prec->initialize();
  }
}
#endif

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif
