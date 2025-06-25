#ifndef PDEKIT_Linear_System_Trilinos_PC_hpp
#define PDEKIT_Linear_System_Trilinos_PC_hpp

#include "PDEKit_Config.hpp"
#include "common/PDEKit.hpp"
#include "linear_system/LocalTransferOps.hpp"
#include "linear_system/TpetraCrsMatrix.hpp"
#include "linear_system/TpetraMultiVector.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "Kokkos_DefaultNode.hpp"
#endif

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------

enum class TpetraCondestType
{
  Cheap = 0,
  CG    = 1,
  GMRES = 2
};

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
class TrilinosPC
{
  public:
  /// TYPEDEFS
  using matrix_type      = TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;
  using multivector_type = TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>;

  /// Constructor
  TrilinosPC();

  /// Destructor
  virtual ~TrilinosPC();

  /// Create the preconditioner
  virtual void create(const std::string &prec_type, const std::shared_ptr<matrix_type> &A) = 0;

  /// Create the preconditioner and pass block structure for creation
  virtual void create(const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
                      std::unique_ptr<std::vector<common::Range1D<LocalOrdinal>>> &&blocks) = 0;

  /// Create the preconditioner and pass extra data for coarse-scale
  /// correction
  virtual void create(const std::string &prec_type, const std::shared_ptr<matrix_type> &A,
                      const common::BlockArray<Uint, Uint> &mesh_dual_graph_crs,
                      const std::vector<Uint> &coarse_block_sizes,
                      const std::vector<Uint> &fine_block_sizes,
                      const std::vector<Uint> &cell_reordering,
                      std::unique_ptr<ls::LocalTransferOps<Real>> &&restriction_ops,
                      std::unique_ptr<ls::LocalTransferOps<Real>> &&prolongation_ops) = 0;

  /// Prepare the internal data of the preconditioner
  virtual void initialize(const std::shared_ptr<matrix_type> &A) = 0;

  /// Compute the preconditioner
  virtual void compute(const std::shared_ptr<matrix_type> &A) = 0;

  /// Apply the preconditioner as linear operator to vector
  /// Computes Y := beta*Y + alpha*Op(A)*X with Op(A) = A or Op(A) =
  /// transpose(A)
  virtual void apply(const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &X,
                     TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &Y,
                     const bool transp = false, Scalar alpha = 1.0, Scalar beta = 0.0) const = 0;

  /// Return the condition number estimate
  virtual Real cond_estimate(const TpetraCondestType cdest_type) const = 0;

  protected:
#if PDEKIT_HAVE_TRILINOS
  using tpetra_op_type = Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  /// FRIENDS
  friend class internal::TpetraInternalAccess;

  virtual Teuchos::RCP<tpetra_op_type const> get_preconditioner_operator() const = 0;
#endif
};

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TrilinosPC()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TrilinosPC<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::~TrilinosPC()
{
}

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif
