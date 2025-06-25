#ifndef PDEKIT_Physics_PhysModelT_hpp
#define PDEKIT_Physics_PhysModelT_hpp

#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"
#include "physics/PhysicsClass.hpp"

namespace pdekit
{

namespace physics
{

// ----------------------------------------------------------------------------

template <typename PhysProperties>
class PhysModelT
{
  public:
  /// ENUMERATIONS AND TYPEDEFS

  enum
  {
    DIM = PhysProperties::DIM
  };
  enum
  {
    NEQ = PhysProperties::NEQ
  };

  typedef PhysProperties Properties;

  typedef typename PhysProperties::CoordV CoordV;     /// type of geometry
                                                      /// coordinates matrix
  typedef typename PhysProperties::SolM SolM;         /// type of solution variables
                                                      /// matrix
  typedef typename PhysProperties::SolV SolV;         /// type of solution variables
                                                      /// vector
  typedef typename PhysProperties::FluxV FluxV;       /// type of flux variables
                                                      /// vector
  typedef typename PhysProperties::SolGradV SolGradV; /// type of solution
                                                      /// gradient vector
  typedef typename PhysProperties::SolGradM SolGradM; /// type of solution
                                                      /// gradient matrix

  typedef math::DenseSMat<Real, NEQ, DIM> FM; /// Flux matrix - each column
                                              /// is one flux vector

  typedef math::DenseSMat<Real, NEQ, NEQ> JM; /// Flux Jacobian matrix

  // Constructor
  PhysModelT();

  // Destructor
  ~PhysModelT();

  protected:
  private:
};

// ----------------------------------------------------------------------------

template <typename PhysProperties>
PhysModelT<PhysProperties>::PhysModelT()
{
}

// ----------------------------------------------------------------------------

template <typename PhysProperties>
PhysModelT<PhysProperties>::~PhysModelT()
{
}

// ----------------------------------------------------------------------------

} // namespace physics

} // namespace pdekit

#endif
