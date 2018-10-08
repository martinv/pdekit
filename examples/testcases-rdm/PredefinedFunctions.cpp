#include "examples/testcases-rdm/PredefinedFunctions.hpp"

namespace pdekit
{

// ----------------------------------------------------------------------------

Real Zero::value(const math::DenseConstVecView<Real> &point_coord,
                 const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                 const Uint component)
{
  return 0.0;
}

// ----------------------------------------------------------------------------

Real BurgersInlet2D::value(
    const math::DenseConstVecView<Real> &point_coord,
    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution, const Uint component)
{
  return 1.5 - 2.0 * point_coord[X];
}

// ----------------------------------------------------------------------------

Real CosineHat2D::value(const math::DenseConstVecView<Real> &point_coord,
                        const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                        const Uint component)
{
  if ((point_coord[X] >= -1.4) && (point_coord[X] <= -0.6))
  {
    const Real cosvalue = std::cos(2.5 * math::pi * (point_coord[X] + 1));
    return std::pow(0.5 * (cosvalue + 1), 4);
  }
  else
  {
    return 0.0;
  }
}

// ----------------------------------------------------------------------------

Real DiffusionBC2D::value(const math::DenseConstVecView<Real> &point_coord,
                          const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                          const Uint component)
{
  if (component == 0)
  {
    /*
    return std::sin(math::pi * (point_coord[X0] + 0.5)) *
           std::sin(math::pi * (point_coord[X1] + 0.5));
    */
    return point_coord[X0] + point_coord[X1];
  }
  else
  {
    return solution[component];
  }
}

// ----------------------------------------------------------------------------

Real SineWave2D::value(const math::DenseConstVecView<Real> &point_coord,
                       const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                       const Uint component)
{
  // sin^2(5*pi*x)
  const Real sine_value = std::sin(5. * math::pi * point_coord[X]);
  return sine_value * sine_value;
}

// ----------------------------------------------------------------------------

Real InletJump2D::value(const math::DenseConstVecView<Real> &point_coord,
                        const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                        const Uint component)
{
  if ((point_coord[X] >= 0.3) && (point_coord[X] <= 0.8) && (point_coord[Y] < 1.e-14))
  {
    return -5.0;
  }
  else
  {
    return 0.0;
  }
}

// ----------------------------------------------------------------------------

Real CosineHat3D::value(const math::DenseConstVecView<Real> &point_coord,
                        const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                        const Uint component)
{
  const Real dist = (point_coord[X] + 1) * (point_coord[X] + 1) + (point_coord[Z] * point_coord[Z]);

  if (dist < 4. / 25.)
  {
    const Real cosvalue = std::cos(2.5 * math::pi * std::sqrt(dist));
    return std::pow(0.5 * (cosvalue + 1), 4);
  }
  else
  {
    return 0.0;
  }
}

// ----------------------------------------------------------------------------

Real RiemannFansInlet2D::value(
    const math::DenseConstVecView<Real> &point_coord,
    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution, const Uint component)
{
  Real value;

  if (point_coord[X] <= 0.5)
  {
    if (component == 0)
    {
      value = 1.0;
    }
    else if (component == 1)
    {
      value = 0.0;
    }
    else if (component == 2)
    {
      value = 2.83972;
    }
    else
    {
      value = 6.532;
    }
  }
  else
  {
    if (component == 0)
    {
      value = 0.5;
    }
    else if (component == 1)
    {
      value = 0.0;
    }
    else if (component == 2)
    {
      value = 1.67332;
    }
    else
    {
      value = 3.425;
    }
  }

  return value;
}

// ----------------------------------------------------------------------------

Real RiemannFansInlet3D::value(
    const math::DenseConstVecView<Real> &point_coord,
    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution, const Uint component)
{
  Real value;

  if (point_coord[X] <= 0.5)
  {
    if (component == 0)
    {
      value = 1.0;
    }
    else if (component == 1)
    {
      value = 0.0;
    }
    else if (component == 2)
    {
      value = 2.83972;
    }
    else if (component == 3)
    {
      value = 0.0;
    }
    else
    {
      value = 6.532;
    }
  }
  else
  {
    if (component == 0)
    {
      value = 0.5;
    }
    else if (component == 1)
    {
      value = 0.0;
    }
    else if (component == 2)
    {
      value = 1.67332;
    }
    else if (component == 3)
    {
      value = 0.0;
    }
    else
    {
      value = 3.425;
    }
  }

  return value;
}

// ----------------------------------------------------------------------------

Real LinearAdvDiff2DSolution::value(
    const math::DenseConstVecView<Real> &point_coord,
    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution, const Uint component)
{
  const Real ax = 0.0;
  const Real ay = 1.0;
  const Real nu = 0.01;

  const Real eta = ay * point_coord[Y] - ax * point_coord[X];
  const Real xi  = ax * point_coord[X] + ay * point_coord[Y];

  const Real cosvalue = std::cos(2. * math::pi * eta);
  const Real expvalue =
      std::exp(xi * (1. - std::sqrt(1. + 16. * math::pi * math::pi * nu * nu)) / (2. * nu));

  return -cosvalue * expvalue;
}

// ----------------------------------------------------------------------------

Real Ringleb2DSolution::value(
    const math::DenseConstVecView<Real> &point_coord,
    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution, const Uint component)
{
  const Real gamma = 1.4;
  const Real a     = Newton_isotach_solve(point_coord);
  const Real rho   = std::pow(a, 2. / (gamma - 1.0));

  const Real J = 1. / a + 1. / (3. * a * a * a) + 1. / (5. * a * a * a * a * a) -
                 0.5 * std::log((1. + a) / (1. - a));
  const Real V2 = 2. * (1. - a * a) / (gamma - 1);
  const Real V  = std::sqrt(V2);

  // Use expression for x(,V) to solve for k:
  const Real k2 = 2. / (1. / V2 - 2. * rho * point_coord[X0] + J * rho);

  // We assume only positive k
  const Real k = std::sqrt(k2);

  // We know k = 1/psi = V/sin(theta), where theta ... hodograph
  // transformation angle Due to roundoff errors, it can happen that sin_theta
  // is slightly bigger than 1
  const Real sin_theta = std::min(1.0, V / k);
  const Real cos_theta = std::sqrt(1. - sin_theta * sin_theta);

  const Real u_sign = (point_coord[X1] >= 0.0) ? -1.0 : 1.0;

  const Real u = u_sign * V * cos_theta;
  const Real v = -V * sin_theta;

  if (component == 0)
  {
    return rho;
  }
  else if (component == 1)
  {
    return rho * u;
  }
  else if (component == 2)
  {
    return rho * v;
  }
  else
  {
    const Real pressure = 1. / gamma * std::pow(a, 2. * gamma / (gamma - 1));
    // p = (gamma-1)*(e - 0.5 * rho * (u^2+v^2)
    const Real e = pressure / (gamma - 1) + 0.5 * rho * (u * u + v * v);
    return e;
  }
}

Real Ringleb2DSolution::F_eval(const math::DenseConstVecView<Real> &point_coord, const Real a)
{
  const Real gamma = 1.4;
  // Note that the domain of log((1+a)/(1-a)) is the interval -1 < a < 1
  // a is speed of sound -> negative values don't make physical sense -> 0 < a
  // <
  // 1
  const Real J = 1. / a + 1. / (3. * a * a * a) + 1. / (5. * a * a * a * a * a) -
                 0.5 * std::log((1. + a) / (1. - a));
  const Real rho = std::pow(a, (2. / (gamma - 1)));
  const Real V2  = 2. * (1. - a * a) / (gamma - 1);

  const Real result = (point_coord[X0] - 0.5 * J) * (point_coord[X0] - 0.5 * J) +
                      point_coord[X1] * point_coord[X1] - 1. / (4. * rho * rho * V2 * V2);

  return result;
}

Real Ringleb2DSolution::Newton_isotach_solve(const math::DenseConstVecView<Real> &point_coord)
{
  // First find a good initial solve
  Real a_test        = 0.1;
  const Real da_test = 0.01;

  Real err_min = std::abs(F_eval(point_coord, a_test));

  Real a   = a_test;
  Real err = err_min;

  while (a_test < 1.0) // See remark in function F_eval about admissible values of a_test
  {
    const Real err_test = std::abs(F_eval(point_coord, a_test));
    if (err_test < err)
    {
      a   = a_test;
      err = err_test;
    }
    a_test += da_test;
  }

  // std::cout << "Initial a = " << a << std::endl;
  // std::cout << "Initial error = " << err << std::endl;

  const Real tol        = 1.e-14;
  const size_t max_iter = 20;
  err                   = 2. * tol;

  // Real a = 0.5;
  // Real err = 2. * tol;
  size_t iter = 0;

  while ((err > tol) && (iter < max_iter))
  {
    const Real da   = std::max(1.e-10, 1.e-7 * std::abs(a));
    const Real dFda = (F_eval(point_coord, a + da) - F_eval(point_coord, a - da)) / (2. * da);

    // std::cout << "a[" << iter << "] = " << a << ", err = " << err <<
    // std::endl;;

    a   = a - F_eval(point_coord, a) / dFda;
    err = std::abs(F_eval(point_coord, a));
    iter++;
  }
  // std::cout << "a[final] = " << a << ", err = " << err << std::endl;;

  return a;
}

// ----------------------------------------------------------------------------

} // namespace pdekit
