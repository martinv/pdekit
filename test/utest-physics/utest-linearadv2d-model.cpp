#include "common/PDEKit.hpp"

#include "physics/scalar/RotationAdvection2D.hpp"

using namespace pdekit;
using namespace pdekit::physics;

int main()
{

  RotationAdvection2D<AroundOrigin2D> model;

  math::DenseSMat<Real, 1, _2D> coords;
  RotationAdvection2D<AroundOrigin2D>::SolM u;
  RotationAdvection2D<AroundOrigin2D>::SolGradM grad;
  RotationAdvection2D<AroundOrigin2D>::Properties props;

  coords(0, X) = 1.0;
  coords(0, Y) = 2.0;

  u.resize(1, 1);

  u(0, 0) = 3.0;

  model.compute_properties(coords.const_row_transp(0), u.const_row_transp(0), grad, props);

  props.print();

  //  model.compute_properties();

  return 0;
}
