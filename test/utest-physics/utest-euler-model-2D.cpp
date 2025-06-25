/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE euler_model_2D_utest
#include <boost/test/unit_test.hpp>

#include "common/PDEKit.hpp"
#include "math/unary_ops/MatrixNorm.hpp"
#include "physics/euler/Euler2DCons.hpp"

using namespace pdekit;
using namespace pdekit::physics;
using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE(physics_navier_stokes_2d_utest)
{

  Euler2DCons model;

  math::DenseSMat<Real, 1, _2D> coords;
  Euler2DCons::SolM u;
  Euler2DCons::SolGradM grad;
  Euler2DCons::Properties props;

  std::array<Euler2DCons::JM, Euler2DCons::DIM> jacobians;
  Euler2DCons::FluxV res;

  coords(0, X) = 1.0;
  coords(0, Y) = 2.0;

  u.resize(1, Euler2DCons::NEQ);

  u(0, 0) = 1.0; // rho
  u(0, 1) = 1.0; // rho * u
  u(0, 2) = 2.0; // rho * v
  u(0, 3) = 5.0; // e

  // Gradient of u
  math::DenseSVec<Real, Euler2DCons::DIM> direction;
  direction[X] = 7.;
  direction[Y] = 2.;

  for (Uint eq = 0; eq < Euler2DCons::NEQ; ++eq)
  {
    grad(eq, X) = direction[X];
    grad(eq, Y) = direction[Y];
  }

  std::cout << "=======================================================" << std::endl;

  std::cout << "Computing properties at state:" << std::endl;
  std::cout.setf(std::ios::fixed);
  std::cout.precision(2);
  std::cout << u << std::endl;

  model.compute_properties(coords.const_row_transp(0), u.const_row_transp(0), grad, props);

  model.residual(props, jacobians, res);

  std::cout << "Jacobians:" << std::endl;
  std::cout << "dF/du = " << std::endl;
  std::cout << jacobians[X] << std::endl;
  std::cout << "dG/du = " << std::endl;
  std::cout << jacobians[Y] << std::endl;

  std::cout << "Direction [nx,ny] in which the linear combination" << std::endl;
  std::cout << "dF/du * nx + dG/du * ny will be computed" << std::endl;

  std::cout << props.grad_vars << std::endl;

  std::cout << "The residual dF/du * nx + dF/du * ny = " << std::endl;
  std::cout << res << std::endl;

  std::cout << "=======================================================" << std::endl;

  std::cout << "Trying to decompose the linear combination L:" << std::endl;

  Euler2DCons::JM R, L, D, Product, Reference, Diff;

  Reference = direction[X] * jacobians[X] + direction[Y] * jacobians[Y];

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

  std::cout << "The norm of (reference - computed value) of dF/du*nx + dG/du*ny = "
            << math::norm_max(Diff) << std::endl;
  BOOST_CHECK_LT(norm_max(Diff), 1.e-12);

  std::cout << "=======================================================" << std::endl;

  //  props.print();

  //  model.compute_properties();
}
