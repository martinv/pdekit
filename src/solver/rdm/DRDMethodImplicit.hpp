#ifndef PDEKIT_PG_D_RD_Method_Implicit_hpp
#define PDEKIT_PG_D_RD_Method_Implicit_hpp

#include <map>

#include "linear_system/LSTpetra.hpp"
#include "solver/rdm/ImplicitDRDFullJacobianImpl.hpp"
#include "solver/rdm/RDMethod.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
class DRDMethodImplicit : public RDMethod<MeshConfig, Physics>
{
  private:
  using cell_scheme_type  = typename internal::CellSchemeSelector<Physics, CellScheme>::type;
  using facet_scheme_type = typename internal::FacetSchemeSelector<Physics, FacetScheme>::type;

  public:
  /// TYPEDEFS

  using f_space_cells  = interpolation::FunctionSpace<MeshConfig>;
  using f_space_facets = interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1>;
  using bc_base_type   = RDMBCBase<MeshConfig, Physics, Physics::DIM - 1>;

  enum
  {
    is_var_beta_type_rdm = cell_scheme_type::is_var_beta_type_rdm
  };
  enum
  {
    needs_volume_derivatives_on_trace = facet_scheme_type::needs_volume_derivatives_on_trace
  };

  /// Constructor
  DRDMethodImplicit();

  /// Destructor
  ~DRDMethodImplicit() override;

  /// Flag saying whether solver is continuous or not
  bool is_continuous() const override;

  /// Flag saying whether solver is explicit
  bool is_explicit() const override;

  /// Configure underlying function space for geometry
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

  /// Add a new boundary condition to all levels
  /// @param bc_type_name   ... type of this boundary condition (e.g.
  /// WeakSubInlet)
  /// @param condition_name ... (user defined) name of the boundary condition
  /// @param boundary_name  ... name of the mesh boundary segment where the bc
  /// should be applied
  std::shared_ptr<bc_base_type> add_boundary_condition(const std::string &bc_type_name,
                                                       const std::string &condition_name,
                                                       const std::string &boundary_name) override;

  /// Return existing boundary condition
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
  void compute_residual_norm(math::DenseDVec<Real> &norm_L2) const override;

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
      typename detail::ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme,
                                                   FacetScheme>::BCScratchDataType;

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

  typename detail::ImplicitDRDFullJacobianImpl<MeshConfig, Physics, CellScheme, FacetScheme> m_impl;
};

// ----------------------------------------------------------------------------
// METHODS OF DISCONTINUOUS RD METHOD - IMPLICIT
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::DRDMethodImplicit()
    : RDMethod<MeshConfig, Physics>(), m_mat(std::make_shared<ls::TpetraCrsMatrix<Real>>()),
      m_rhs(std::make_shared<ls::TpetraMultiVector<Real>>()),
      m_du(std::make_shared<ls::TpetraMultiVector<Real>>()),
      m_sparse_solver_type(ls::SparseSolverType::eIterative)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::~DRDMethodImplicit()
{
  for (typename std::map<std::string, BCScratchDataType *>::iterator it = m_weak_bc_data.begin();
       it != m_weak_bc_data.end(); ++it)
  {
    delete it->second;
  }
  m_weak_bc_data.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
bool DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::is_continuous() const
{
  return false;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
bool DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::is_explicit() const
{
  return false;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::initialize_work_data(
    const SFunc sf_type, const PointSetID quad_type, const Uint quadrature_order,
    const Uint nb_threads, const Uint nb_blocks)
{
  rd_method_base::initialize_work_data(sf_type, quad_type, quadrature_order, nb_threads, nb_blocks);

  // Configure the actual discretization
  m_impl.configure_cell_spaces((*rd_method_base::m_tria), (*rd_method_base::m_sol_dofs), *m_mat,
                               *m_rhs, *m_du, nb_threads, nb_blocks, rd_method_base::m_sf_type,
                               rd_method_base::m_quad_type, rd_method_base::m_quad_order);

  m_impl.configure_facet_spaces((*rd_method_base::m_tria), (*rd_method_base::m_sol_dofs),
                                rd_method_base::m_sf_type, rd_method_base::m_quad_type,
                                rd_method_base::m_quad_order);

  m_lin_system = std::make_shared<ls::LSTpetra<Real>>(m_sparse_solver_type);
  m_lin_system->initialize_solver(m_mat, m_rhs, m_du, false);

  if (m_sparse_solver_type == ls::SparseSolverType::eIterative)
  {
#if 1
    m_preconditioner = std::make_shared<ls::IfpackPC<Real>>();
    // m_preconditioner->create("ILUT",
    // ls::internal::TpetraInternalAccess::get_matrix_data(*m_mat));
    m_preconditioner->create("ILUT", m_mat);

    m_lin_system->connect_preconditioner(m_preconditioner);
#else

    common::BlockArray<Uint, Uint> mesh_dual_graph_crs;
    (*rd_method_base::m_tria).build_dual_graph_undirected(mesh_dual_graph_crs);

    const mesh::DofMap<MeshConfig> &sol_dofs = (*rd_method_base::m_sol_dofs);

    std::vector<Uint> block_sizes_coarse(sol_dofs.nb_active_cells());
    std::vector<Uint> block_sizes_fine(sol_dofs.nb_active_cells());

    mesh::StdRegion std_reg;
    for (Uint ac = 0; ac < sol_dofs.nb_active_cells(); ++ac)
    {
      const mesh::MeshEntity cell          = sol_dofs.active_cell(ac);
      const mesh::PointSetTag original_tag = cell.std_region_id();

      const mesh::PointSetTag std_reg_tag(original_tag.elem_shape(), P1,
                                          original_tag.ref_topology());
      std_reg.change_type(std_reg_tag);
      block_sizes_fine[ac]   = Physics::NEQ * cell.nb_vert();
      block_sizes_coarse[ac] = Physics::NEQ * std_reg.get().nb_nodes();
    }

    std::unique_ptr<ls::LocalTransferOps<Real>> restriction_ops(new ls::LocalTransferOps<Real>());
    interpolation::CoarseScaleCorrectionOpBuilder::build_restriction_ops<MeshConfig>(
        Physics::NEQ, sol_dofs, sol_dofs, sol_dofs.cell_reordering(), P1, *restriction_ops);

    std::unique_ptr<ls::LocalTransferOps<Real>> prolongation_ops(new ls::LocalTransferOps<Real>());
    interpolation::CoarseScaleCorrectionOpBuilder::build_prolongation_ops<MeshConfig>(
        Physics::NEQ, sol_dofs, sol_dofs, sol_dofs.cell_reordering(), P1, *prolongation_ops);

    m_preconditioner = std::make_shared<ls::TpetraCoarseScaleCorrPC<Real>>();

    m_preconditioner->create("", m_mat, mesh_dual_graph_crs, block_sizes_coarse, block_sizes_fine,
                             sol_dofs.cell_reordering(), std::move(restriction_ops),
                             std::move(prolongation_ops));

    m_lin_system->connect_preconditioner(m_preconditioner);

/*
std::unique_ptr<std::vector<common::Range1D<Int>>> block_ranges(
    new std::vector<common::Range1D<Int>>());
const mesh::DofMap<MeshConfig> &sol_dofs = (*rd_method_base::m_sol_dofs);

block_ranges->resize(sol_dofs.nb_active_cells());
const std::vector<Uint> &cell_ordering = sol_dofs.cell_reordering();

std::vector<Uint> inverse_cell_ordering(cell_ordering.size());
for (Uint i = 0; i < cell_ordering.size(); ++i)
{
  inverse_cell_ordering[cell_ordering[i]] = i;
}

Uint block_offset = 0;
Uint b = 0;
for (auto idx : inverse_cell_ordering)
{
  const mesh::MeshEntity active_cell = sol_dofs.active_cell(idx);
  (*block_ranges)[b] = common::Range1D<Int>(
      block_offset, block_offset + active_cell.nb_vert() * Physics::NEQ - 1);
  block_offset += active_cell.nb_vert() * Physics::NEQ;
  b++;
}

m_preconditioner = std::make_shared<ls::TpetraBlockDiagPreconditioner<Real>>();
m_preconditioner->create("", m_mat, std::move(block_ranges));

m_lin_system->connect_preconditioner(m_preconditioner);
*/
#endif
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::compute_node_reordering(
    typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
    typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs, std::vector<Int> &reordering)
{
  m_impl.compute_node_reordering_impl(topology, dofs, reordering);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::compute_cell_reordering(
    typename RDMethod<MeshConfig, Physics>::mesh_type const &topology,
    typename RDMethod<MeshConfig, Physics>::cell_dofs const &dofs, std::vector<Int> &reordering)
{
  m_impl.compute_cell_reordering_impl(topology, dofs, reordering);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::set_vec_function(
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

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::set_blending_coeff(
    const interpolation::ScalarMeshFunction<Real>::ptr &blending_coeff)
{
  rd_method_base::set_blending_coeff(blending_coeff);
  m_impl.set_blending_coeff(blending_coeff);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::set_artificial_viscosity(
    const interpolation::ScalarMeshFunction<Real>::ptr &artificial_viscosity)
{
  rd_method_base::set_artificial_viscosity(artificial_viscosity);
  m_impl.set_artificial_viscosity(artificial_viscosity);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
std::shared_ptr<
    typename DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::bc_base_type>
DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::add_boundary_condition(
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

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
std::shared_ptr<
    typename DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::bc_base_type>
DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::boundary_condition(
    const std::string &condition_name) const
{
  return rd_method_base::m_bc_manager.get_bc(condition_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::remove_boundary_condition(
    const std::string &condition_name)
{
  rd_method_base::m_bc_manager.remove_bc(condition_name);
  m_weak_bc_data.erase(condition_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme,
                       FacetScheme>::remove_all_boundary_conditions()
{
  rd_method_base::m_bc_manager.clear_all_bcs();
  m_weak_bc_data.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::assemble_lhs_and_rhs(
    const Real CFL)
{
  m_impl.assemble_lhs_and_rhs_impl(*rd_method_base::m_tria, *rd_method_base::m_sol_dofs,
                                   m_is_Dirichlet_node, rd_method_base::m_bc_manager,
                                   m_weak_bc_data,
                                   *rd_method_base::vec_function(SolverVecFn::solution), *m_mat,
                                   *m_rhs, *rd_method_base::vec_function(SolverVecFn::residuals),
                                   rd_method_base::m_rd_time_update, rd_method_base::m_cfl);
  rd_method_base::m_cfl = CFL;

  // rd_method_base::m_rd_time_update.compute_local_time_step(CFL,
  // rd_method_base::m_timestep);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::assemble_rhs()
{
  m_impl.assemble_rhs_impl(*rd_method_base::m_tria, *rd_method_base::m_sol_dofs,
                           m_is_Dirichlet_node, rd_method_base::m_bc_manager, m_weak_bc_data,
                           *rd_method_base::vec_function(SolverVecFn::solution), *m_rhs,
                           *rd_method_base::vec_function(SolverVecFn::residuals),
                           rd_method_base::m_rd_time_update);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::solve(
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

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::compute_residual_norm(
    math::DenseDVec<Real> &norm_L2) const
{
  m_impl.compute_residual_norm(*rd_method_base::m_sol_dofs, *m_rhs, norm_L2);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename CellScheme, typename FacetScheme>
void DRDMethodImplicit<MeshConfig, Physics, CellScheme, FacetScheme>::detect_Dirichlet_dofs()
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
