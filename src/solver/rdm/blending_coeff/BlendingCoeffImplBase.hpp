#ifndef PDEKIT_Solver_RDM_Blending_Coeff_Implementation_Base_hpp
#define PDEKIT_Solver_RDM_Blending_Coeff_Implementation_Base_hpp

#include <string>

#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/DenseDVec.hpp"
#include "math/MathConstants.hpp"
#include "mesh/MeshConfig.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

namespace detail
{

template <typename MeshConfig>
class BlendingCoeffImplBase
{
  public:
  /// Default constructor
  BlendingCoeffImplBase();

  /// Default destructor
  virtual ~BlendingCoeffImplBase();

  /// Resize internal data
  virtual void setup(typename result_of::dof_map_t<MeshConfig> const &sol_dofs) = 0;

  /// Evaluate the blending coefficient
  virtual void calculate_impl(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                              const interpolation::VectorMeshFunction<Real> &u,
                              interpolation::ScalarMeshFunction<Real> &function) = 0;

  /// Clear all internal data
  virtual void clear() = 0;

  /// Set a parameter
  virtual void set_param(const std::string &param_name, const Real param_value) = 0;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
BlendingCoeffImplBase<MeshConfig>::BlendingCoeffImplBase()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
BlendingCoeffImplBase<MeshConfig>::~BlendingCoeffImplBase()
{
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
