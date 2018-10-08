#ifndef PDEKIT_Trilinos_Tpetra_Multi_Vector_hpp
#define PDEKIT_Trilinos_Tpetra_Multi_Vector_hpp

#include <iostream>
#include <mutex>

#include "PDEKit_Config.hpp"
#include "common/PDEKit.hpp"

#if PDEKIT_HAVE_TRILINOS
#include "Tpetra_MultiVector.hpp"
#endif

#include "linear_system/TpetraDofMap.hpp"
#include "linear_system/TpetraTraits.hpp"

namespace pdekit
{

namespace ls
{

namespace internal
{
class TpetraInternalAccess;
}

// ----------------------------------------------------------------------------
//             Forward declarations to enable friend functions
// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
class TpetraMultiVector;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
std::ostream &operator<<(
    std::ostream &os,
    const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &tv);

// ----------------------------------------------------------------------------

template <typename Scalar        = TpetraDefaultTraits::Scalar,
          typename LocalOrdinal  = TpetraDefaultTraits::LocalOrdinal,
          typename GlobalOrdinal = TpetraDefaultTraits::GlobalOrdinal,
          typename Node          = TpetraDefaultTraits::Node,
          const bool classic     = TpetraDefaultTraits::classic>
class TpetraMultiVector
{
  public:
#if PDEKIT_HAVE_TRILINOS
  using trilinos_multivector_type = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using scalar_type               = typename trilinos_multivector_type::scalar_type;
  using local_ordinal_type        = typename trilinos_multivector_type::local_ordinal_type;
  using global_ordinal_type       = typename trilinos_multivector_type::global_ordinal_type;
  using node_type                 = typename trilinos_multivector_type::node_type;
#endif

  /// Default constructor - do nothing
  TpetraMultiVector();

  /// Constructor - a map has to be supplied
  TpetraMultiVector(const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> &trilinos_map,
                    const LocalOrdinal nb_vectors = 1);

  /// Destructor
  ~TpetraMultiVector();

  /// Initialize from trilinos map, set number of vectors
  void init(const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> &trilinos_map,
            const LocalOrdinal nb_vectors = 1);

  /// Element access operator
  /// @param i ... index IN vector
  /// @param j ... index OF vector in case multiple vectors are present
  void insert_value(const Uint i, const Scalar value, const Uint j = 0);

  /// Element access operator, const version
  /// @param i ... index IN vector
  /// @param j ... index OF vector in case multiple vectors are present
  void add_value(const Uint i, const Scalar value, const Uint j = 0);

  /// Add values to the j-th vector in multivector
  /// @param j_vec ... index of the vector in multivector
  /// @param buffer ... vector of values to add
  void add_values(const Uint j_vec, const std::vector<std::tuple<Uint, Real>> &buffer);

  /// Get one value
  const Scalar value(const Uint i, const Uint j = 0) const;

  /// Return the length of one vector
  Uint size() const;

  /// Return the total number of vectors stored in the multi-vector
  Uint num_vectors() const;

  /// Fill the whole vector with given value
  void fill(const Scalar value);

  /// Assign from another vector
  void assign(const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &src);

  /// Print the values of the vector
  friend std::ostream &operator<<<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>(
      std::ostream &os,
      const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &tv);

  private:
#if PDEKIT_HAVE_TRILINOS
  friend class internal::TpetraInternalAccess;

  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> get_vector()
      const;

  /// The actual trilinos vector
  Teuchos::RCP<trilinos_multivector_type> m_vector;

  /// Mutex to make some methods thread-safe
  mutable std::mutex m_mutex;
#endif
};

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraMultiVector()
    : m_vector(nullptr)
{
}
#else
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraMultiVector()
{
}
#endif

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraMultiVector(
    const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> &trilinos_map,
    const LocalOrdinal nb_vectors)
    : m_vector(new Tpetra::MultiVector<Scalar, LocalOrdinal>(
          internal::TpetraInternalAccess::get_dof_map(trilinos_map), nb_vectors))
{
}
#else
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::TpetraMultiVector(
    const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> &trilinos_map,
    const LocalOrdinal nb_vectors)
{
}
#endif

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::~TpetraMultiVector()
{
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::init(
    const TpetraDofMap<LocalOrdinal, GlobalOrdinal, Node> &trilinos_map,
    const LocalOrdinal nb_vectors)
{
#if PDEKIT_HAVE_TRILINOS
  m_vector = Teuchos::RCP<trilinos_multivector_type>(new Tpetra::MultiVector<Scalar, LocalOrdinal>(
      internal::TpetraInternalAccess::get_dof_map(trilinos_map), nb_vectors));
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::insert_value(
    const Uint i, const Scalar value, const Uint j)
{
#if PDEKIT_HAVE_TRILINOS
  /*
  // This is LOCAL access!
  Teuchos::ArrayRCP<Scalar> vec = m_vector->getDataNonConst(j);
  vec[i] = value;
  */
  m_vector->replaceGlobalValue(i, j, value);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::add_value(
    const Uint i, const Scalar value, const Uint j)
{
#if PDEKIT_HAVE_TRILINOS
  /*
  // This is LOCAL access!
  Teuchos::ArrayRCP<Scalar> vec = m_vector->getDataNonConst(j);
  vec[i] += value;
  */
  m_vector->sumIntoGlobalValue(i, j, value);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::add_values(
    const Uint j_vec, const std::vector<std::tuple<Uint, Real>> &buffer)
{
#if PDEKIT_HAVE_TRILINOS
  std::lock_guard<std::mutex> lock(m_mutex);

  // This is LOCAL access!
  Teuchos::ArrayRCP<Scalar> vec = m_vector->getDataNonConst(j_vec);

  for (Uint i = 0; i < buffer.size(); ++i)
  {
    const Uint local_idx = std::get<0>(buffer[i]);
    const Scalar value   = static_cast<Scalar>(std::get<1>(buffer[i]));
    vec[local_idx] += value;
  }
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
const Scalar TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::value(
    const Uint i, const Uint j) const
{
#if PDEKIT_HAVE_TRILINOS
  // This is LOCAL access!
  Teuchos::ArrayRCP<Scalar const> const_vec = m_vector->getData(j);
  return const_vec[i];
#else
  return Scalar();
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Uint TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::size() const
{
#if PDEKIT_HAVE_TRILINOS
  return static_cast<Uint>(m_vector->getGlobalLength());
#else
  return 0;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Uint TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::num_vectors() const
{
#if PDEKIT_HAVE_TRILINOS
  return static_cast<Uint>(m_vector->getNumVectors());
#else
  return 0;
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::fill(const Scalar value)
{
#if PDEKIT_HAVE_TRILINOS
  m_vector->putScalar(value);
#endif
}

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::assign(
    const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &src)
{
#if PDEKIT_HAVE_TRILINOS
  m_vector->assign(*src.m_vector);
#endif
}

// ----------------------------------------------------------------------------

#if PDEKIT_HAVE_TRILINOS
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>> TpetraMultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::get_vector() const
{
  return m_vector;
}
#endif

// ----------------------------------------------------------------------------

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node,
          const bool classic>
std::ostream &operator<<(
    std::ostream &os,
    const TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> &tv)
{
#if PDEKIT_HAVE_TRILINOS
  os << (*tv.m_vector);
#endif
  return os;
}

// ----------------------------------------------------------------------------

// Defer instantiation of Real TpetraMultiVector to the cpp file
extern template class TpetraMultiVector<Real>;

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_Trilinos_Sparse_Matrix_hpp
