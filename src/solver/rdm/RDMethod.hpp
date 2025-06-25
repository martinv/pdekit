#ifndef PDEKIT_PG_RD_Method_hpp
#define PDEKIT_PG_RD_Method_hpp

#include <memory>

#include "interpolation/FunctionSpace.hpp"
#include "solver/SolverBase.hpp"
#include "solver/rdm/RDTimeUpdate.hpp"
#include "solver/rdm/bc/RDMBCManager.hpp"

namespace pdekit
{

namespace mesh
{
template <typename MeshConfig>
class Mesh;
}

namespace solver
{

namespace rdm
{

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
class RDMBCBase;

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
class RDMethod : public SolverBase
{
  public:
  /// TYPEDEFS
  using mesh_type     = typename mesh::Tria<MeshConfig>;
  using mesh_boundary = typename result_of::mesh_boundary_set_t<MeshConfig>;
  using cell_dofs     = typename result_of::dof_map_t<MeshConfig>;

  using f_space_cells  = interpolation::FunctionSpace<MeshConfig>;
  using f_space_facets = interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1>;
  using bc_base_type   = RDMBCBase<MeshConfig, Physics, Physics::DIM - 1>;

  /// Default constructor
  RDMethod();

  /// Destructor
  ~RDMethod() override;

  /// Set the data:
  /// a) Topology:
  ///   - cell topology
  ///   - mesh boundaries
  /// b) Geometrical support:
  ///   - geometry dofs
  /// c) Solution support:
  ///   - solution dofs
  void configure_mesh_data(const std::shared_ptr<const mesh_type> &mesh,
                           const common::PtrHandle<const cell_dofs> geo_dofs,
                           const common::PtrHandle<const cell_dofs> sol_dofs);

  /// Prepare the method - allocate data etc.
  virtual void initialize_work_data(const SFunc sf_type, const PointSetID quad_type,
                                    const Uint quadrature_order, const Uint nb_threads,
                                    const Uint nb_blocks) = 0;

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  virtual void compute_node_reordering(const mesh_type &mesh, const cell_dofs &dofs,
                                       std::vector<Int> &reordering) = 0;

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  virtual void compute_cell_reordering(const mesh_type &mesh, const cell_dofs &dofs,
                                       std::vector<Int> &reordering) = 0;

  /// Set the computation of blending coefficient
  virtual void set_blending_coeff(
      const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff);

  /// Set the smoothness/discontinuity sensor
  virtual void set_artificial_viscosity(
      const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity);

  /// Add a new boundary condition to all levels
  /// @param bc_type_name   ... type of this boundary condition (e.g.
  /// WeakSubInlet)
  /// @param condition_name ... (user defined) name of the boundary condition
  /// @param boundary_name  ... name of the mesh boundary segment where the bc
  /// should be applied
  virtual std::shared_ptr<bc_base_type> add_boundary_condition(const std::string &bc_type_name,
                                                               const std::string &condition_name,
                                                               const std::string &boundary_name);

  /// Return existing boundary condition
  /// @param condition_name ... (user defined) name of the boundary condition
  virtual std::shared_ptr<bc_base_type> boundary_condition(
      const std::string &condition_name) const = 0;

  /// Remove a boundary condition
  virtual void remove_boundary_condition(const std::string &condition_name) = 0;

  /// Remove all boundary conditions
  virtual void remove_all_boundary_conditions() = 0;

  /// Compute L1 norm of the residuals (to monitor convergence)
  virtual void compute_residual_norm(math::DenseDVec<Real> &norm_L2) const = 0;

  protected:
  /// METHODS
  std::shared_ptr<f_space_facets> setup_boundary_facet_sol_space(
      const std::string &boundary_name) const;

  /// DATA
  SFunc m_sf_type;
  PointSetID m_quad_type;
  Uint m_quad_order;

  /// Handles to mesh data: support and solution
  std::shared_ptr<const mesh_type> m_tria;
  common::PtrHandle<const cell_dofs> m_geo_dofs;
  common::PtrHandle<const cell_dofs> m_sol_dofs;

  /// Manager to create boundary conditions
  RDMBCManager<MeshConfig, Physics> m_bc_manager;

  /// Pointers to the functions holding update coefficient,
  /// blending coefficient and artificial viscosity
  RDTimeUpdate m_rd_time_update;
  interpolation::ScalarMeshFunction<Real>::ptr m_blending_coeff;
  interpolation::ScalarMeshFunction<Real>::ptr m_art_viscosity;

  private:
  /// Handles to approximation spaces
  typename f_space_cells::ptr m_geo_cell_space;
  typename f_space_cells::ptr m_sol_cell_space;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
RDMethod<MeshConfig, Physics>::RDMethod()
    : SolverBase(), m_quad_order(0), m_geo_cell_space(nullptr), m_sol_cell_space(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
RDMethod<MeshConfig, Physics>::~RDMethod()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDMethod<MeshConfig, Physics>::configure_mesh_data(
    const std::shared_ptr<const mesh_type> &mesh, const common::PtrHandle<const cell_dofs> geo_dofs,
    const common::PtrHandle<const cell_dofs> sol_dofs)
{
  m_tria     = mesh;
  m_geo_dofs = geo_dofs;
  m_sol_dofs = sol_dofs;

  // Setup the time update object
  m_rd_time_update.setup<MeshConfig>(*mesh, *sol_dofs);

  auto time_step = scal_function(SolverScalFn::time_step);
  time_step->resize((*sol_dofs).nb_nodes());
  time_step->fill(0.0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDMethod<MeshConfig, Physics>::initialize_work_data(const SFunc sf_type,
                                                         const PointSetID quad_type,
                                                         const Uint quadrature_order,
                                                         const Uint nb_threads,
                                                         const Uint nb_blocks)
{
  m_sf_type    = sf_type;
  m_quad_type  = quad_type;
  m_quad_order = quadrature_order;

  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [quadrature_order, quad_type](const ElemShape shape,
                                                      const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order, quad_type);
  };

  /*
  const SFGenerator sf_generator{sf_type};
  const QuadGenerator quad_generator{quad_type, quadrature_order};
  */

  m_geo_cell_space = std::make_shared<f_space_cells>();
  m_geo_cell_space->set_reference_fe_values((*m_geo_dofs).as_range(), sf_generator, quad_generator);

  m_sol_cell_space = std::make_shared<f_space_cells>();
  m_sol_cell_space->set_reference_fe_values((*m_sol_dofs).as_range(), sf_generator, quad_generator);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDMethod<MeshConfig, Physics>::set_blending_coeff(
    const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff)
{
  m_blending_coeff = blending_coeff;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void RDMethod<MeshConfig, Physics>::set_artificial_viscosity(
    const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity)
{
  m_art_viscosity = artificial_viscosity;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
std::shared_ptr<typename RDMethod<MeshConfig, Physics>::bc_base_type> RDMethod<
    MeshConfig, Physics>::add_boundary_condition(const std::string &bc_type_name,
                                                 const std::string &condition_name,
                                                 const std::string &boundary_name)
{
  std::shared_ptr<bc_base_type> new_bc = m_bc_manager.create_bc(bc_type_name, condition_name);

  new_bc->configure_mesh_data(m_tria, m_geo_dofs, m_sol_dofs, boundary_name);
  new_bc->configure_spaces(m_sf_type, m_quad_type, m_quad_order);
  new_bc->set_solution(vec_function(SolverVecFn::solution));
  new_bc->set_residuals(vec_function(SolverVecFn::residuals));

  return new_bc;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
std::shared_ptr<typename RDMethod<MeshConfig, Physics>::f_space_facets> RDMethod<
    MeshConfig, Physics>::setup_boundary_facet_sol_space(const std::string &boundary_name) const
{
  mesh::MeshBoundarySet<MeshConfig> const &mesh_boundaries = m_tria->all_boundaries();
  mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry =
      *(mesh_boundaries.domain(boundary_name));

  typename f_space_facets::ptr bdry_space = std::make_shared<f_space_facets>();

  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto bdry_quad_generator = [this](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, m_quad_order + 1, m_quad_type);
  };

  bdry_space->set_reference_fe_values(bdry, *m_sol_dofs, sf_generator, bdry_quad_generator);

  return bdry_space;
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
