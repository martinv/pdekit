#include "linear_system/TrilinosCrsMatrix.hpp"

#if PDEKIT_HAVE_TRILINOS

#include "Ifpack_Utils.h"

#include "common/MPI/MPIEnv.hpp"

namespace pdekit
{

namespace ls
{

// ============================================================================

TrilinosCrsMatrix::TrilinosCrsMatrix(const Uint global_nb_rows)
    : m_comm(common::mpi::MPIEnv::instance().comm())
{
  m_map.reset(new Epetra_Map(static_cast<Int>(global_nb_rows), 0, m_comm));
  // m_epetra_matrix.reset( new Epetra_CrsMatrix(Copy,*m_map,global_nb_rows)
  // );
  m_epetra_matrix.reset(new Epetra_CrsMatrix(Copy, *m_map, 8));
}

// ============================================================================

TrilinosCrsMatrix::~TrilinosCrsMatrix()
{
}

// ============================================================================

Uint TrilinosCrsMatrix::nb_global_elem() const
{
  return static_cast<Uint>(m_map->NumGlobalElements());
}

// ============================================================================

void TrilinosCrsMatrix::lock_structure()
{
  m_epetra_matrix->FillComplete();
}

// ============================================================================

void TrilinosCrsMatrix::insert_values_in_row(const Uint row_idx, const std::vector<Real> &values,
                                             const std::vector<Int> &indices)
{
  m_epetra_matrix->InsertGlobalValues(row_idx, indices.size(), values.data(), indices.data());
}

// ============================================================================

void TrilinosCrsMatrix::add_values_to_row(const Uint row_idx, const std::vector<Real> &values,
                                          const std::vector<Int> &indices)
{
  m_epetra_matrix->SumIntoGlobalValues(row_idx, indices.size(), values.data(), indices.data());
}

// ============================================================================

const std::shared_ptr<Epetra_Map> TrilinosCrsMatrix::map() const
{
  return m_map;
}

// ============================================================================

const Epetra_Map &TrilinosCrsMatrix::domain_map() const
{
  return m_epetra_matrix->OperatorDomainMap();
}

// ============================================================================

void TrilinosCrsMatrix::fill(const Real value)
{
  m_epetra_matrix->PutScalar(value);
}

// ============================================================================

void TrilinosCrsMatrix::print_structure_to_file(const std::string &filename) const
{
  Ifpack_PrintSparsity(*m_epetra_matrix, filename.c_str());
}

// ============================================================================

std::ostream &operator<<(std::ostream &os, const TrilinosCrsMatrix &tril_crs_mat)
{
  os << *tril_crs_mat.m_epetra_matrix;
  return os;
}

// ============================================================================

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_HAVE_TRILINOS
