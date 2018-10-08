#ifndef PDEKIT_DG_Solver_hpp
#define PDEKIT_DG_Solver_hpp

#include <memory>

#include "interpolation/FunctionSpace.hpp"
#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "mesh/Tria.hpp"
#include "mesh/containers/DofMap.hpp"
#include "solver/fe/DGTimeUpdate.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
class DGSolver
{
  public:
  /// TYPEDEFS

  /// TYPEDEFS
  using mesh_type     = typename mesh::Tria<MeshConfig>;
  using mesh_boundary = typename result_of::mesh_boundary_set_t<MeshConfig>;
  using cell_dofs     = typename result_of::dof_map_t<MeshConfig>;

  using f_space = interpolation::FunctionSpace<MeshConfig>;
  // using bc_base_type = RDMBCBase<MeshConfig, Physics, Physics::DIM - 1>;

  /// Constructor
  DGSolver();

  /// Copy constructor: deleted
  DGSolver(const DGSolver &rhs) = delete;

  /// Assignement operator: deleted
  DGSolver &operator=(const DGSolver &rhs) = delete;

  /// Destructor
  ~DGSolver();

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

  /// Configure underlying function space for geometry and solution
  virtual void build_approximation_spaces(const SFunc sf_type, const PointSetID quad_type,
                                          const Uint quadrature_order) = 0;

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  virtual void compute_node_reordering(const mesh_type &mesh, const cell_dofs &dofs,
                                       std::vector<Int> &reordering) = 0;

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  virtual void compute_cell_reordering(const mesh_type &mesh, const cell_dofs &dofs,
                                       std::vector<Int> &reordering) = 0;

  /// Prepare the method - allocate data etc.
  virtual void initialize_work_data(const Uint nb_threads, const Uint nb_blocks) = 0;

  /// Return pointer to the underlying geometry space
  typename f_space::ptr geo_space() const;

  /// Return pointer to the underlying solution space
  typename f_space::ptr sol_space() const;

  /// Set the solution function
  void set_solution(const interpolation::VectorMeshFunction<Real>::ptr &solution);

  interpolation::VectorMeshFunction<Real>::ptr solution() const;

  /// Set the (nodal) residuals function
  void set_residuals(const interpolation::VectorMeshFunction<Real>::ptr &nodal_residuals);

  /// Get the (nodal) residual function - element residuals
  interpolation::VectorMeshFunction<Real>::ptr residuals() const;

  /// Get the (nodal) residual function - element residuals
  interpolation::ScalarMeshFunction<Real> const &time_step() const;

  protected:
  /// DATA
  /// Handles to approximation spaces

  typename f_space::ptr m_geo_cell_space;
  typename f_space::ptr m_sol_cell_space;

  SFunc m_sf_type;
  PointSetID m_quad_type;
  Uint m_quad_order;

  /// Handles to mesh data: support and solution
  std::shared_ptr<const mesh_type> m_mesh;
  common::PtrHandle<const cell_dofs> m_geo_dofs;
  common::PtrHandle<const cell_dofs> m_sol_dofs;

  interpolation::VectorMeshFunction<Real>::ptr m_solution;
  interpolation::VectorMeshFunction<Real>::ptr m_residuals;
  interpolation::ScalarMeshFunction<Real> m_timestep;

  /// Pointers to the functions holding update coefficient,
  /// blending coefficient and artificial viscosity
  DGTimeUpdate m_dg_time_update;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
DGSolver<MeshConfig, Physics>::DGSolver() : m_timestep("", "time_step")
{
  m_timestep.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
DGSolver<MeshConfig, Physics>::~DGSolver()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void DGSolver<MeshConfig, Physics>::configure_mesh_data(
    const std::shared_ptr<const mesh_type> &mesh, const common::PtrHandle<const cell_dofs> geo_dofs,
    const common::PtrHandle<const cell_dofs> sol_dofs)
{
  m_mesh     = mesh;
  m_geo_dofs = geo_dofs;
  m_sol_dofs = sol_dofs;

  // Setup the time update object
  m_dg_time_update.setup<MeshConfig>(*mesh, *sol_dofs);

  // SolverBase::m_timestep.resize((*sol_dofs).nb_nodes());
  // SolverBase::m_timestep.fill(0.0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
typename DGSolver<MeshConfig, Physics>::f_space::ptr DGSolver<MeshConfig, Physics>::geo_space()
    const
{
  return m_geo_cell_space;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
typename DGSolver<MeshConfig, Physics>::f_space::ptr DGSolver<MeshConfig, Physics>::sol_space()
    const
{
  return m_sol_cell_space;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void DGSolver<MeshConfig, Physics>::set_solution(
    const interpolation::VectorMeshFunction<Real>::ptr &solution)
{
  m_solution = solution;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
interpolation::VectorMeshFunction<Real>::ptr DGSolver<MeshConfig, Physics>::solution() const
{
  return m_solution;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
void DGSolver<MeshConfig, Physics>::set_residuals(
    const interpolation::VectorMeshFunction<Real>::ptr &nodal_residuals)
{
  m_residuals = nodal_residuals;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
interpolation::VectorMeshFunction<Real>::ptr DGSolver<MeshConfig, Physics>::residuals() const
{
  return m_residuals;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics>
interpolation::ScalarMeshFunction<Real> const &DGSolver<MeshConfig, Physics>::time_step() const
{
  return m_timestep;
}

// ----------------------------------------------------------------------------

} // namespace fe

} // namespace solver

} // namespace pdekit

// ----------------------------------------------------------------------------

#endif
