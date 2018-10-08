#ifndef PDEKIT_Trilinos_Sparse_Matrix_hpp
#define PDEKIT_Trilinos_Sparse_Matrix_hpp

#include "PDEKit_Config.hpp"
#if PDEKIT_HAVE_TRILINOS

#include <iostream>
#include <memory>
#include <vector>

#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_RCP.hpp"

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace ls
{

class TrilinosCrsMatrix
{
  public:
  /// Default constructor
  TrilinosCrsMatrix(const Uint global_nb_rows);

  /// Destructor
  ~TrilinosCrsMatrix();

  /// Get the number of global elements in the matrix
  Uint nb_global_elem() const;

  /// Lock the structure of the matrix
  void lock_structure();

  /// Insert values in one row
  void insert_values_in_row(const Uint row_idx, const std::vector<Real> &values,
                            const std::vector<Int> &indices);

  /// Accumulate values in one row
  void add_values_to_row(const Uint row_idx, const std::vector<Real> &values,
                         const std::vector<Int> &indices);

  /// Get the underlying map
  const std::shared_ptr<Epetra_Map> map() const;

  /// Get the domain map of this matrix
  /// @note This should be used to construct the vectors
  ///       for linear systems involving this matrix!
  const Epetra_Map &domain_map() const;

  /// Assign the same scalar value to all entries in the sparse matrix
  void fill(const Real value);

  /// Print the structure of the matrix to file
  void print_structure_to_file(const std::string &filename) const;

  /// Print info
  void print_info() const
  {
#if PDEKIT_HAVE_TRILINOS
    std::cout << "LSS with trilinos support" << std::endl;
#else
    std::cout << "LSS without trilinos support" << std::endl;
#endif
  }

  /// Print the values of the matrix
  friend std::ostream &operator<<(std::ostream &os, const TrilinosCrsMatrix &tril_crs_mat);

  private:
  friend class LSTrilinos;

  /// Trilinos communicator object
  /// Takes MPI communicator in constructor
  Epetra_MpiComm m_comm;

  /// Map - needed to create the matrix
  std::shared_ptr<Epetra_Map> m_map;

  /// Sparse distributed matrix from Trilinos
  Teuchos::RCP<Epetra_CrsMatrix> m_epetra_matrix;
};

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_HAVE_TRILINOS

#endif // PDEKIT_Trilinos_Sparse_Matrix_hpp
