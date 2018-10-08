#ifndef Predefined_Functions_hpp
#define Predefined_Functions_hpp

#include <cmath>
#include <iostream>

#include "common/PDEKit.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/MathConstants.hpp"

namespace pdekit
{

// ----------------------------------------------------------------------------
// Functions that define various initial values / boundary conditions
// for the RDS test cases
// ----------------------------------------------------------------------------

class Zero
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);
};

// ----------------------------------------------------------------------------

class BurgersInlet2D
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);
};

// ----------------------------------------------------------------------------

class CosineHat2D
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);
};

// ----------------------------------------------------------------------------

class DiffusionBC2D
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);
};

// ----------------------------------------------------------------------------

class SineWave2D
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);
};

// ----------------------------------------------------------------------------

class InletJump2D
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);
};

// ----------------------------------------------------------------------------

class CosineHat3D
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);
};

// ----------------------------------------------------------------------------

class RiemannFansInlet2D
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);
};

// ----------------------------------------------------------------------------

class RiemannFansInlet3D
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);
};

// ----------------------------------------------------------------------------

class Rotation2DSolution
{
  public:
  template <typename CoordVector>
  static Real value(const CoordVector &point_coord)
  {
    const Real dist = std::sqrt(point_coord[X] * point_coord[X] + point_coord[Y] * point_coord[Y]);

    if ((dist >= 0.6) && (dist <= 1.4))
    {
      const Real cosvalue = std::cos(math::pi * (dist - 1.0) / 0.4);
      return std::pow(0.5 * (cosvalue + 1), 4);
    }

    return 0.0;
  }
};

// ----------------------------------------------------------------------------

class SineWave2DSolution
{
  public:
  template <typename CoordVector>
  static Real value(const CoordVector &point_coord)
  {
    const Real sine_value = std::sin(5. * math::pi * point_coord[X]);
    return sine_value * sine_value;
  }
};

// ----------------------------------------------------------------------------

class LinearAdvDiff2DSolution
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);
};

// ----------------------------------------------------------------------------

class Ringleb2DSolution
{
  public:
  static Real value(const math::DenseConstVecView<Real> &point_coord,
                    const interpolation::VectorMeshFunction<Real>::const_entry_type &solution,
                    const Uint component);

  private:
  static Real F_eval(const math::DenseConstVecView<Real> &point_coord, const Real a);

  static Real Newton_isotach_solve(const math::DenseConstVecView<Real> &point_coord);
};

// ----------------------------------------------------------------------------

} // namespace pdekit

#endif
