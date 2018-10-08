#ifndef PDEKIT_Solver_RDM_Blending_Coeff_hpp
#define PDEKIT_Solver_RDM_Blending_Coeff_hpp

#include <string>

#include "interpolation/mesh_function/MeshFunctionTools.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "physics/euler/Euler3DCons.hpp"
#include "solver/rdm/blending_coeff/ArtViscIndicatorCoeffImpl.hpp"
#include "solver/rdm/blending_coeff/BonanniBlendingCoeffImpl.hpp"
#include "solver/rdm/blending_coeff/LFBlendingCoeffImpl.hpp"
#include "solver/rdm/blending_coeff/WongJansenBlendingCoeffImpl.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig>
class RDBlendingCoeff
{
  public:
  /// Default constructor
  RDBlendingCoeff();

  /// Copy constructor is deleted
  RDBlendingCoeff(const RDBlendingCoeff &other_coeff) = delete;

  /// Default destructor
  ~RDBlendingCoeff();

  /// Assignment operator is deleted
  RDBlendingCoeff &operator=(const RDBlendingCoeff &other_coeff) = delete;

  /// Resize internal data
  void setup(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
             const std::string &coeff_type);

  /// Evaluate the blending coefficient
  void calculate(mesh::Tria<MeshConfig> &input_mesh,
                 typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                 const interpolation::VectorMeshFunction<Real> &u,
                 interpolation::ScalarMeshFunction<Real> &function);

  /// Evaluate the blending coefficient
  void calculate_in_cells(mesh::Tria<MeshConfig> &input_mesh,
                          typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                          const interpolation::VectorMeshFunction<Real> &u,
                          interpolation::ScalarMeshFunction<Real> &function);

  /// Clear all internal data
  void clear();

  /// Set a parameter
  void set_param(const std::string &param_name, const Real param_value);

  private:
  /// TYPES
  enum CoeffType : unsigned short
  {
    Undefined = 0,
    LF_Blend  = 1,          // Blending coeff. for Lax-Friedrichs scheme of Abgrall
                            // and coworkers
    ViscosityIndicator = 2, // Indicator based on artificial viscosity paper
                            // by Persson & Peraire
    WongJansen = 3,         // Indicator from the paper of Wong & Jansen: Residual
                            // distribution finite element method for
                            // convection-dominated problems
    Bonanni2D = 4,
    Bonanni3D = 5
  };

  /// DATA

  /// Information about the algorithm that is used to compute
  CoeffType m_coeff_type;

  std::unique_ptr<detail::BlendingCoeffImplBase<MeshConfig>> m_blend_coeff_impl;

  /// Cellwise values of blending coefficient
  interpolation::ScalarMeshFunction<Real> m_cell_blend_values;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
RDBlendingCoeff<MeshConfig>::RDBlendingCoeff()
    : m_coeff_type(CoeffType::Undefined), m_cell_blend_values("", "cell_blend_values")
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
RDBlendingCoeff<MeshConfig>::~RDBlendingCoeff()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void RDBlendingCoeff<MeshConfig>::setup(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                                        const std::string &coeff_type)
{
  if (coeff_type == "LF_Blend")
  {
    detail::LFBlendingCoeffImpl<MeshConfig> *blend_coeff_impl_ptr =
        new detail::LFBlendingCoeffImpl<MeshConfig>();
    m_blend_coeff_impl =
        std::move(std::unique_ptr<detail::LFBlendingCoeffImpl<MeshConfig>>(blend_coeff_impl_ptr));
    m_blend_coeff_impl->setup(sol_dofs);

    m_coeff_type = CoeffType::LF_Blend;
  }
  else if (coeff_type == "ViscosityIndicator")
  {
    detail::ArtViscIndicatorCoeffImpl<MeshConfig> *blend_coeff_impl_ptr =
        new detail::ArtViscIndicatorCoeffImpl<MeshConfig>();
    m_blend_coeff_impl = std::move(
        std::unique_ptr<detail::ArtViscIndicatorCoeffImpl<MeshConfig>>(blend_coeff_impl_ptr));
    m_blend_coeff_impl->setup(sol_dofs);

    m_coeff_type = CoeffType::ViscosityIndicator;
  }
  else if (coeff_type == "WongJansen")
  {
    detail::WongJansenBlendingCoeffImpl<MeshConfig> *blend_coeff_impl_ptr =
        new detail::WongJansenBlendingCoeffImpl<MeshConfig>();
    m_blend_coeff_impl = std::move(
        std::unique_ptr<detail::WongJansenBlendingCoeffImpl<MeshConfig>>(blend_coeff_impl_ptr));
    m_blend_coeff_impl->setup(sol_dofs);

    m_coeff_type = CoeffType::WongJansen;
  }
  else if (coeff_type == "Bonanni2D")
  {
    detail::BonanniBlendingCoeffImpl<MeshConfig, physics::Euler2DCons> *blend_coeff_impl_ptr =
        new detail::BonanniBlendingCoeffImpl<MeshConfig, physics::Euler2DCons>();
    m_blend_coeff_impl = std::move(
        std::unique_ptr<detail::BonanniBlendingCoeffImpl<MeshConfig, physics::Euler2DCons>>(
            blend_coeff_impl_ptr));
    m_blend_coeff_impl->setup(sol_dofs);

    m_coeff_type = CoeffType::Bonanni2D;
  }
  else if (coeff_type == "Bonanni3D")
  {
    detail::BonanniBlendingCoeffImpl<MeshConfig, physics::Euler3DCons> *blend_coeff_impl_ptr =
        new detail::BonanniBlendingCoeffImpl<MeshConfig, physics::Euler3DCons>();
    m_blend_coeff_impl = std::move(
        std::unique_ptr<detail::BonanniBlendingCoeffImpl<MeshConfig, physics::Euler3DCons>>(
            blend_coeff_impl_ptr));
    m_blend_coeff_impl->setup(sol_dofs);

    m_coeff_type = CoeffType::Bonanni3D;
  }

  m_cell_blend_values.resize(sol_dofs.nb_active_cells());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void RDBlendingCoeff<MeshConfig>::calculate(
    mesh::Tria<MeshConfig> &input_mesh, typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &u,
    interpolation::ScalarMeshFunction<Real> &function)
{
  m_cell_blend_values.resize(sol_dofs.nb_active_cells());
  m_blend_coeff_impl->calculate_impl(sol_dofs, u, m_cell_blend_values);

  /*
  for (Uint ac = 0; ac < m_cell_blend_values.nb_entries(); ++ac)
  {
    function[ac] = m_cell_blend_values[ac];
  }
  */

  interpolation::MeshFunctionTools::interp_from_cells_to_nodes_p1_max(
      input_mesh, sol_dofs, m_cell_blend_values, function);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void RDBlendingCoeff<MeshConfig>::calculate_in_cells(
    mesh::Tria<MeshConfig> &input_mesh, typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &u,
    interpolation::ScalarMeshFunction<Real> &function)
{
  m_cell_blend_values.resize(sol_dofs.nb_active_cells());
  m_blend_coeff_impl->calculate_impl(sol_dofs, u, function);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void RDBlendingCoeff<MeshConfig>::clear()
{
  m_coeff_type = CoeffType::Undefined;
  m_blend_coeff_impl->clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void RDBlendingCoeff<MeshConfig>::set_param(const std::string &param_name, const Real param_value)
{
  m_blend_coeff_impl->set_param(param_name, param_value);
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
