#ifndef PDEKIT_Trilinos_Multi_Vector_hpp
#define PDEKIT_Trilinos_Multi_Vector_hpp

#include "PDEKit_Config.hpp"
#if PDEKIT_HAVE_TRILINOS

#include <iostream>
#include <memory>
#include <vector>

#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace ls
{

class TrilinosMultiVector
{
  public:
  /// Constructor - a map has to be supplied
  TrilinosMultiVector(const Epetra_Map &trilinos_map, const Uint nb_vectors = 1);

  /// Destructor
  ~TrilinosMultiVector();

  /// Element access operator
  /// @param i ... index IN vector
  /// @param j ... index OF vector in case multiple vectors are present
  inline Real &operator()(const Uint i, const Uint j = 0)
  {
    // Epetra_Vector* & vec = (*m_vector)(0);
    // return (*vec)[static_cast<Int>(i)];
    // return (*m_vector)[static_cast<Int>(i)];

    double *&vec = (*m_vector)[j];
    return vec[static_cast<Int>(i)];
  }

  /// Element access operator, const version
  /// @param i ... index IN vector
  /// @param j ... index OF vector in case multiple vectors are present
  inline Real const &operator()(const Uint i, const Uint j = 0) const
  {
    double *&vec = (*m_vector)[j];
    return vec[static_cast<Int>(i)];
    // return (*m_vector)[static_cast<Int>(i)];
  }

  /// Return the length of one vector
  Uint size() const;

  /// Return the total number of vectors stored in the multi-vector
  Uint num_vectors() const;

  /// Fill the whole vector with given value
  void fill(const Real value);

  /// Print the values of the vector
  friend std::ostream &operator<<(std::ostream &os, const TrilinosMultiVector &tv);

  private:
  friend class LSTrilinos;

  /// The actual trilinos vector
  Teuchos::RCP<Epetra_MultiVector> m_vector;
};

} // namespace ls

} // namespace pdekit

#endif // PDEKIT_HAVE_TRILINOS

#endif // PDEKIT_Trilinos_Sparse_Matrix_hpp
