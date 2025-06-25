#ifndef PDEKIT_PG_RD_Method_Implicit_hpp
#define PDEKIT_PG_RD_Method_Implicit_hpp

#include <map>

#include "linear_system/LSTpetra.hpp"
#include "solver/rdm/ImplicitPGRDFullJacobianImpl.hpp"
#include "solver/rdm/RDMethod.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, typename SchemeTraits>
class PGRDMethodImplicit : public RDMethod<MeshConfig, Physics>
{

  public:
  /// TYPEDEFS

  using f_space_cells = interpolation::FunctionSpace<MeshConfig>;
  using bc_base_type  = RDMBCBase<MeshConfig, Physics, Physics::DIM - 1>;

  /// Constructor
  PGRDMethodImplicit();

  /// Destructor
  ~PGRDMethodImplicit() override;

  /// Flag saying whether solver is continuous or not
  bool is_continuous() const override;

  /// Flag saying whether solver is explicit
  bool is_explicit() const override;

  /// Prepare the method - allocate data etc.
  void initialize_work_data(const SFunc sf_type, const PointSetID quad_type,
                            const Uint quadrature_order, const Uint nb_threads,
                            const Uint nb_blocks) override;

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  void compute_node_reordering(typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
                               typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs,
                               std::vector<Int> &reordering) override;

  /// Compute ordering of dofs in that will lead to better-conditioned linear
  /// system
  void compute_cell_reordering(typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
                               typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs,
                               std::vector<Int> &reordering) override;

  /// Configure a vector function field
  void set_vec_function(const SolverVecFn vec_fn_id,
                        const interpolation::VectorMeshFunction<Real>::ptr &vec_fn) override;

  void set_blending_coeff(
      const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff) override;

  void set_artificial_viscosity(
      const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity) override;

  /// Add a new boundary condition
  /// @param bc_type_name   ... type of this boundary condition (e.g.
  /// WeakSubInlet)
  /// @param condition_name ... (user defined) name of the boundary condition
  /// @param boundary_name  ... name of the mesh boundary segment where the bc
  /// should be applied
  std::shared_ptr<bc_base_type> add_boundary_condition(const std::string &bc_type_name,
                                                       const std::string &condition_name,
                                                       const std::string &boundary_name) override;

  /// Get existing boundary condition
  /// @param condition_name ... (user defined) name of the boundary condition
  std::shared_ptr<bc_base_type> boundary_condition(
      const std::string &condition_name) const override;

  /// Remove a boundary condition with given name
  void remove_boundary_condition(const std::string &condition_name) override;

  /// Remove all boundary conditions
  void remove_all_boundary_conditions() override;

  /// Assemble both lhs and rhs simultaneously
  void assemble_lhs_and_rhs(const Real CFL) override;

  /// Assemble only rhs (assuming that left-hand side) is already
  /// assembled
  void assemble_rhs() override;

  /// Solve the underlying system
  void solve(std::vector<SolverOption> const &solver_options = {}) override;

  /// Compute L1 norm of the residuals (to monitor convergence)
  void compute_residual_norm(math::DenseDVec<Real> &norm_L1) const override;

  protected:
  /// TYPEDEFS
  using rd_method_base = RDMethod<MeshConfig, Physics>;

  using geo_cache_type  = interpolation::GeometryCache<MeshConfig::GDIM>;
  using geo_metric_type = interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>;

  using sol_cache_type  = interpolation::SolutionCache;
  using sol_metric_type = interpolation::SolutionSpaceMetric<MeshConfig, MeshConfig::GDIM>;

  using flux_metric_type = interpolation::FluxSpaceMetric<MeshConfig, Physics, MeshConfig::GDIM>;

  using bc_metric_data_type = RDMBCMetricData<MeshConfig, Physics, MeshConfig::TDIM - 1>;

  /// METHODS
  /// Find which degrees of freedom are on Dirichlet boundary
  void detect_Dirichlet_dofs();

  /// TYPES
  using BCScratchDataType =
      typename detail::ImplicitPGRDFullJacobianImpl<MeshConfig, Physics,
                                                    SchemeTraits>::BCScratchDataType;
  /// DATA

  /// Vector determining which DOF is on Dirichlet boundary
  std::vector<bool> m_is_Dirichlet_node;

  /// Data (geometry and solution cache and metric)
  /// for weak boundary conditions
  std::map<std::string, BCScratchDataType *> m_weak_bc_data;

  /// Matrices and vectors for the linear system
  std::shared_ptr<ls::TpetraCrsMatrix<Real>> m_mat;
  std::shared_ptr<ls::TpetraMultiVector<Real>> m_rhs;
  std::shared_ptr<ls::TpetraMultiVector<Real>> m_du;

  std::shared_ptr<ls::LSTpetra<Real>> m_lin_system;
  std::shared_ptr<ls::TrilinosPC<Real>> m_preconditioner;

  ls::SparseSolverType m_sparse_solver_type;

  typename detail::ImplicitPGRDFullJacobianImpl<MeshConfig, Physics, SchemeTraits> m_impl;
  // typename detail::ImplicitPGRDLUSGSImpl<MeshConfig, Physics, SchemeTraits>
  // m_impl;

  Uint m_iter_counter;
};

// ----------------------------------------------------------------------------
// METHODS OF CONTINUOUS (PG) RD METHOD - IMPLICIT
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::PGRDMethodImplicit()
    : RDMethod<MeshConfig, Physics>(), m_mat(std::make_shared<ls::TpetraCrsMatrix<Real>>()),
      m_rhs(std::make_shared<ls::TpetraMultiVector<Real>>()),
      m_du(std::make_shared<ls::TpetraMultiVector<Real>>()),
      m_sparse_solver_type(ls::SparseSolverType::eIterative), m_iter_counter(0)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::~PGRDMethodImplicit()
{
  for (typename std::map<std::string, BCScratchDataType *>::iterator it = m_weak_bc_data.begin();
       it != m_weak_bc_data.end(); ++it)
  {
    delete it->second;
  }
  m_weak_bc_data.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
bool PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::is_continuous() const
{
  return true;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
bool PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::is_explicit() const
{
  return false;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::initialize_work_data(
    const SFunc sf_type, const PointSetID quad_type, const Uint quadrature_order,
    const Uint nb_threads, const Uint nb_blocks)
{
  rd_method_base::initialize_work_data(sf_type, quad_type, quadrature_order, nb_threads, nb_blocks);

#if 0
  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [quadrature_order, quad_type](const ElemShape shape,
                                                      const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order, quad_type);
  };

  rd_method_base::m_geo_cell_space = std::make_shared<interpolation::FunctionSpace<MeshConfig>>();

  rd_method_base::m_geo_cell_space->set_reference_fe_values(
    (*rd_method_base::m_geo_dofs).as_range(), sf_generator, quad_generator);

  rd_method_base::m_sol_cell_space = std::make_shared<interpolation::FunctionSpace<MeshConfig>>();

  rd_method_base::m_sol_cell_space->set_reference_fe_values(
    (*rd_method_base::m_sol_dofs).as_range(), sf_generator, quad_generator);
#endif

  m_impl.configure_cell_spaces((*rd_method_base::m_tria), (*rd_method_base::m_sol_dofs), *m_mat,
                               *m_rhs, *m_du, nb_threads, nb_blocks, rd_method_base::m_sf_type,
                               rd_method_base::m_quad_type, rd_method_base::m_quad_order);

  m_lin_system = std::make_shared<ls::LSTpetra<Real>>(m_sparse_solver_type);
  m_lin_system->initialize_solver(m_mat, m_rhs, m_du, false);

  if (m_sparse_solver_type == ls::SparseSolverType::eIterative)
  {
    m_preconditioner = std::make_shared<ls::IfpackPC<Real>>();
    // m_preconditioner->create("ILUT",
    // ls::internal::TpetraInternalAccess::get_matrix_data(*m_mat));
    m_preconditioner->create("ILUT", m_mat);

    m_lin_system->connect_preconditioner(m_preconditioner);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::compute_node_reordering(
    typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
    typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs, std::vector<Int> &reordering)
{
  m_impl.compute_node_reordering_impl(topology, dofs, reordering);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::compute_cell_reordering(
    typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
    typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs, std::vector<Int> &reordering)
{
  m_impl.compute_cell_reordering_impl(topology, dofs, reordering);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::set_vec_function(
    const SolverVecFn vec_fn_id, const interpolation::VectorMeshFunction<Real>::ptr &vec_fn)
{
  rd_method_base::set_vec_function(vec_fn_id, vec_fn);

  if (SolverVecFn::solution == vec_fn_id)
  {
    // Make sure that the boundary conditions are updated as well regarding
    // the active solution field
    using bc_iterator = typename RDMBCManager<MeshConfig, Physics>::const_iterator;
    for (bc_iterator bc_it = rd_method_base::m_bc_manager.cbegin();
         bc_it != rd_method_base::m_bc_manager.cend(); ++bc_it)
    {
      bc_it->second->set_solution(this->vec_function(SolverVecFn::solution));
    }
  }
  else if (SolverVecFn::residuals == vec_fn_id)
  {
    // Make sure that the boundary conditions are updated as well regarding
    // the active residual field
    using bc_iterator = typename RDMBCManager<MeshConfig, Physics>::const_iterator;
    for (bc_iterator bc_it = rd_method_base::m_bc_manager.cbegin();
         bc_it != rd_method_base::m_bc_manager.cend(); ++bc_it)
    {
      bc_it->second->set_residuals(this->vec_function(SolverVecFn::residuals));
    }
  }
  else if (SolverVecFn::sources == vec_fn_id)
  {
    m_impl.set_sources(vec_fn);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::set_blending_coeff(
    const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff)
{
  rd_method_base::set_blending_coeff(blending_coeff);
  m_impl.set_blending_coeff(blending_coeff);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::set_artificial_viscosity(
    const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity)
{
  rd_method_base::set_artificial_viscosity(artificial_viscosity);
  m_impl.set_artificial_viscosity(artificial_viscosity);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
std::shared_ptr<typename PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::bc_base_type>
PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::add_boundary_condition(
    const std::string &bc_type_name, const std::string &condition_name,
    const std::string &boundary_name)
{
  auto new_bc = rd_method_base::add_boundary_condition(bc_type_name, condition_name, boundary_name);

  RDMBCType const rdm_bc_type = new_bc->bc_type();

  // If this is a weak boundary condition, we need to allocate
  // geometry/solution cache and metric data
  if (rdm_bc_type == BC_TYPE_WEAK)
  {
    /*
    mesh::MeshBoundarySet<MeshConfig> const &mesh_boundaries =
      rd_method_base::m_tria->all_boundaries();
    mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &sol_bdry =
      *(mesh_boundaries.domain(boundary_name));

    typename interpolation::FunctionSpace<MeshConfig, Physics::DIM - 1>::ptr
    bdry_sol_space =
    std::make_shared<interpolation::FunctionSpace<MeshConfig, Physics::DIM -
    1>>();

    auto sf_generator = [=](const ElemShape shape, const Uint order) {
      return mesh::sf::SFTag(shape, SFunc::Lagrange, order,
    ModalBasis::Modal);
    };

    // 1) Get the facet quadrature generator
    const typename f_space_cells::quad_generator_fcn &sol_cell_quad_gen =
      rd_method_base::m_sol_cell_space->quad_generator();

    // 2) Capture facet quadrature generator and use it to get 'higher'
    quadrature rule auto bdry_sol_quad_generator = [sol_cell_quad_gen](const
    ElemShape shape, const Uint elem_order) { const mesh::PointSetTag tmp =
    sol_cell_quad_gen(shape, elem_order); return
    mesh::PointSetTag(tmp.elem_shape(), tmp.poly_order() + 1,
    tmp.ref_topology());
    };

    bdry_sol_space->set_reference_fe_values(sol_bdry,
    *rd_method_base::m_sol_dofs, sf_generator, bdry_sol_quad_generator);
    */

    mesh::MeshBoundarySet<MeshConfig> const &mesh_boundaries =
        rd_method_base::m_tria->all_boundaries();
    mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &sol_bdry =
        *(mesh_boundaries.domain(boundary_name));

    const std::shared_ptr<typename rd_method_base::f_space_facets> bdry_sol_space =
        rd_method_base::setup_boundary_facet_sol_space(boundary_name);

    BCScratchDataType *new_bc_data = new BCScratchDataType();
    const Uint nb_bdry_cells       = sol_bdry.nb_active_cells();

    new_bc_data->residual_cache.allocate(bdry_sol_space->reference_elements().cbegin(),
                                         bdry_sol_space->reference_elements().cend(), nb_bdry_cells,
                                         Physics::NEQ);
    new_bc_data->residual_cache_perturbed.allocate(bdry_sol_space->reference_elements().cbegin(),
                                                   bdry_sol_space->reference_elements().cend(),
                                                   nb_bdry_cells, Physics::NEQ);

    m_weak_bc_data.insert(std::pair<std::string, BCScratchDataType *>(condition_name, new_bc_data));
  }

  // Update information about which nodes are on Dirichlet boundaries
  detect_Dirichlet_dofs();

  return new_bc;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
std::shared_ptr<typename PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::bc_base_type>
PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::boundary_condition(
    const std::string &condition_name) const
{
  return rd_method_base::m_bc_manager.get_bc(condition_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::remove_boundary_condition(
    const std::string &condition_name)
{
  rd_method_base::m_bc_manager.remove_bc(condition_name);
  m_weak_bc_data.erase(condition_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::remove_all_boundary_conditions()
{
  rd_method_base::m_bc_manager.clear_all_bcs();
  m_weak_bc_data.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::assemble_lhs_and_rhs(const Real CFL)
{
  m_impl.assemble_lhs_and_rhs_impl(*rd_method_base::m_tria, *rd_method_base::m_sol_dofs,
                                   m_is_Dirichlet_node, rd_method_base::m_bc_manager,
                                   m_weak_bc_data,
                                   *rd_method_base::vec_function(SolverVecFn::solution), *m_mat,
                                   *m_rhs, rd_method_base::m_rd_time_update, rd_method_base::m_cfl);
  rd_method_base::m_cfl = CFL;

  // rd_method_base::m_rd_time_update.compute_local_time_step(CFL,
  // rd_method_base::m_timestep);

  m_impl.update_residuals(*m_rhs, *rd_method_base::vec_function(SolverVecFn::residuals));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::assemble_rhs()
{
  m_impl.assemble_rhs_impl(*rd_method_base::m_tria, *rd_method_base::m_sol_dofs,
                           m_is_Dirichlet_node, rd_method_base::m_bc_manager, m_weak_bc_data,
                           *rd_method_base::vec_function(SolverVecFn::solution), *m_rhs,
                           rd_method_base::m_rd_time_update);

  m_impl.update_residuals(*m_rhs, *rd_method_base::vec_function(SolverVecFn::residuals));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::solve(
    std::vector<SolverOption> const &solver_options)
{
  bool recompute_preconditioner = false;
  for (auto opt : solver_options)
  {
    if (opt == SolverOption::RecomputePreconditioner)
    {
      recompute_preconditioner = true;
    }
  }
  m_impl.solve(*rd_method_base::m_sol_dofs, m_mat, m_rhs, m_du, *m_lin_system,
               *rd_method_base::vec_function(SolverVecFn::solution), recompute_preconditioner);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::compute_residual_norm(
    math::DenseDVec<Real> &norm_L1) const
{
  m_impl.compute_residual_norm(*rd_method_base::m_sol_dofs, *m_rhs, norm_L1);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>::detect_Dirichlet_dofs()
{
  using bc_manager_type = RDMBCManager<MeshConfig, Physics>;

  typename result_of::dof_map_t<MeshConfig> const &sol_cell_connectivity =
      *rd_method_base::m_sol_dofs;
  typename result_of::mesh_boundary_set_t<MeshConfig> const &mesh_boundaries =
      rd_method_base::m_tria->all_boundaries();

  m_is_Dirichlet_node.resize(sol_cell_connectivity.nb_nodes() * Physics::NEQ);
  m_is_Dirichlet_node.assign(m_is_Dirichlet_node.size(), false);

  for (typename bc_manager_type::const_iterator it = rd_method_base::m_bc_manager.cbegin();
       it != rd_method_base::m_bc_manager.cend(); ++it)
  {
    if (it->second->bc_type() == BC_TYPE_STRONG)
    {
      // std::cout << "BC " << it->second->name() << " is Dirichlet BC" <<
      // std::endl;

      const std::string domain_name = it->second->domain_name();

      const mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> &bdry =
          *(mesh_boundaries.domain(domain_name));

      for (Uint c = 0; c < bdry.nb_active_cells(); ++c)
      {
        /*
        const mesh::IncidencePair bcell_id = bdry.bdry_cell_id(c);
        mesh::MeshEntity sol_cell =
        sol_cell_connectivity.active_cell(bcell_id.cell_idx);
        sol_cell.local_transform(MeshConfig::TDIM - 1,
        bcell_id.local_id);
        */

        const mesh::MeshEntity sol_cell =
            bdry.active_cell(sol_cell_connectivity, mesh::ActiveIdx(c));

        for (Uint v = 0; v < sol_cell.nb_vert(); ++v)
        {
          for (Uint eq = 0; eq < Physics::NEQ; ++eq)
          {
            m_is_Dirichlet_node[Physics::NEQ * sol_cell.vertex(v) + eq] = true;
          }
        }
      }

      // Apply the boundary condition
      // it->second->apply(time_update);
    } // If this is a strong BC
  }   // Loop over all boundary conditions
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
