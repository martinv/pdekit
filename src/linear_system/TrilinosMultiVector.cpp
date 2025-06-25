#include "linear_system/TrilinosMultiVector.hpp"

#if PDEKIT_HAVE_TRILINOS

namespace pdekit
{

namespace ls
{

// ----------------------------------------------------------------------------

TrilinosMultiVector::TrilinosMultiVector(const Epetra_Map &trilinos_map, const Uint nb_vectors)
    : m_vector(new Epetra_MultiVector(trilinos_map, static_cast<int>(nb_vectors)))
{
}

// ----------------------------------------------------------------------------

TrilinosMultiVector::~TrilinosMultiVector()
{
}

// ----------------------------------------------------------------------------

Uint TrilinosMultiVector::size() const
{
  return static_cast<Uint>(m_vector->GlobalLength());
}

// ----------------------------------------------------------------------------

Uint TrilinosMultiVector::num_vectors() const
{
  return static_cast<Uint>(m_vector->NumVectors());
}

// ----------------------------------------------------------------------------

void TrilinosMultiVector::fill(const Real value)
{
  m_vector->PutScalar(value);
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const TrilinosMultiVector &tv)
{
  os << (*tv.m_vector);
  return os;
}

// ----------------------------------------------------------------------------

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_HAVE_TRILINOS
