#ifndef PDEKIT_Solver_Time_BDF1_hpp
#define PDEKIT_Solver_Time_BDF1_hpp

#include "linear_system/LSTpetra.hpp"

namespace pdekit
{

namespace solver
{

namespace time
{

// ----------------------------------------------------------------------------

template <typename SpaceDiscretization>
class BDF1
{
  public:
  /// Default constructor
  BDF1();

  /// Delete copy constructor
  BDF1(const BDF1 &other) = delete;

  /// Default destructor
  ~BDF1();

  /// Delete assignment operator
  BDF1 &operator=(const BDF1 &other) = delete;

  private:
  /// Matrices and vectors for the linear system
  std::shared_ptr<ls::TpetraCrsMatrix<Real>> m_mat;
  std::shared_ptr<ls::TpetraMultiVector<Real>> m_du;
  std::shared_ptr<ls::TpetraMultiVector<Real>> m_rhs;

  /// Tpetra linear system
  ls::LSTpetra<Real> m_lin_system;
};

// ----------------------------------------------------------------------------

template <typename SpaceDiscretization>
BDF1<SpaceDiscretization>::BDF1()
    : m_mat(std::make_shared < ls::TpetraCrsMatrix<Real>(new ls::TpetraCrsMatrix<Real>())),
      m_du(std::make_shared < ls::TpetraMultiVector<Real>(new ls::TpetraMultiVector<Real>())),
      m_rhs(std::make_shared < ls::TpetraMultiVector<Real>(new ls::TpetraMultiVector<Real>()))
{
}

// ----------------------------------------------------------------------------

template <typename SpaceDiscretization>
BDF1<SpaceDiscretization>::BDF1()
{
}

// ----------------------------------------------------------------------------

} // namespace time

} // namespace solver

} // namespace pdekit
