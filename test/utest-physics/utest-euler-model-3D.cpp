/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE euler_model_3D_utest
#include <boost/test/unit_test.hpp>

#include "common/PDEKit.hpp"
#include "math/unary_ops/MatrixNorm.hpp"
#include "physics/euler/Euler3DCons.hpp"

using namespace pdekit;
using namespace pdekit::physics;
using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE(physics_navier_stokes_2d_utest)
{

  Euler3DCons model;

  math::DenseSMat<Real, 1, _3D> coords;
  Euler3DCons::SolM u;
  Euler3DCons::SolGradM grad;
  Euler3DCons::Properties props;

  std::array<Euler3DCons::JM, Euler3DCons::DIM> jacobians;
  Euler3DCons::FluxV res;

  coords(0, X) = 1.0;
  coords(0, Y) = 2.0;
  coords(0, Z) = 3.0;

  u.resize(1, Euler3DCons::NEQ);

  u(0, 0) = 1.0; // rho
  u(0, 1) = 1.0; // rho * u
  u(0, 2) = 2.0; // rho * v
  u(0, 3) = 1.5; // rho * w
  u(0, 4) = 5.0; // e

  // Gradient of u
  math::DenseSVec<Real, Euler3DCons::DIM> direction;
  direction[X] = 1.;
  direction[Y] = 2.;
  direction[Z] = -3;

  for (Uint eq = 0; eq < Euler3DCons::NEQ; ++eq)
  {
    grad(eq, X) = direction[X];
    grad(eq, Y) = direction[Y];
    grad(eq, Z) = direction[Z];
  }

  std::cout << "=======================================================" << std::endl;

  std::cout << "Computing properties at state:" << std::endl;
  std::cout.setf(std::ios::fixed);
  std::cout.precision(2);
  std::cout << u << std::endl;

  model.compute_properties(coords.const_row_transp(0), u.const_row_transp(0), grad, props);

  model.residual(props, jacobians, res);

  std::cout << "Jacobians:" << std::endl;
  std::cout << "dF1/du = " << std::endl;
  std::cout << jacobians[X] << std::endl;
  std::cout << "dF2/du = " << std::endl;
  std::cout << jacobians[Y] << std::endl;
  std::cout << "dF3/du = " << std::endl;
  std::cout << jacobians[Z] << std::endl;

  std::cout << "Direction [nx,ny,nz] in which the linear combination" << std::endl;
  std::cout << "dF1/du * nx + dF2/du * ny + dF3/du * nz will be computed" << std::endl;

  std::cout << props.grad_vars << std::endl;

  std::cout << "The residual dF/du * nx + dF/du * ny + dF/du * nz = " << std::endl;
  std::cout << res << std::endl;

  std::cout << "=======================================================" << std::endl;

  std::cout << "Trying to decompose the linear combination L:" << std::endl;

  Euler3DCons::JM R, L, D, Product, Reference, Diff;

  Reference =
      direction[X] * jacobians[X] + direction[Y] * jacobians[Y] + direction[Z] * jacobians[Z];

  R.fill(0.0);
  L.fill(0.0);
  D.fill(0.0);

  model.flux_jacobian_eigen_structure(props, direction, R, L, D);

  std::cout << "Matrix of right eigenvectors:" << std::endl;
  std::cout << R << std::endl;
  std::cout << "Matrix of left eigenvectors:" << std::endl;
  std::cout << L << std::endl;
  std::cout << "Matrix of eigenvalues:" << std::endl;
  std::cout << D << std::endl;

  Product = R * D * L;

  std::cout << "The product R * D * L = " << std::endl;
  std::cout << Product << std::endl;

  std::cout << "Reference value:" << std::endl;
  std::cout << Reference << std::endl;

  Diff = Product - Reference;

  std::cout << "The norm of (reference - computed value) of dF1/du*nx + "
               "dF2/du*ny + dF3/du*nz = "
            << math::norm_max(Diff) << std::endl;
  BOOST_CHECK_LT(math::norm_max(Diff), 1.e-12);

  std::cout << "=======================================================" << std::endl;

  //  props.print();

  //  model.compute_properties();
}
