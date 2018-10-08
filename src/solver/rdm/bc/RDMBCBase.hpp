#ifndef RDM_Boundary_Condition_Base_hpp
#define RDM_Boundary_Condition_Base_hpp

#include "interpolation/FluxSpaceMetric.hpp"
#include "mesh/Tria.hpp"
#include "solver/rdm/RDTimeUpdate.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

// Forward declarations

template <typename MeshConfig, typename Physics, Uint BcDIM>
struct RDMBCMetricData;

enum RDMBCType
{
  BC_TYPE_STRONG = 0,
  BC_TYPE_WEAK   = 1
};

// RDMBCBase::Base class to implement boundary conditions
// MeshConfig ... specifies the type of mesh (topological and geometrical
// dimension) that this
//                BC can be used on
// Physics    ... physical model of given boundary condition
// BcDIM      ... space dimension of entity to which the boundary condition is
// applied Example: RDMBCBase<Cart3D, Euler3D, _2D> would be a boundary
// condition
//          applied to a surface (2D entity) of a 3D mesh

template <typename MeshConfig, typename Physics, Uint BcDIM = Physics::DIM - 1>
class RDMBCBase
{
  public:
  enum
  {
    DIM = BcDIM
  };

  /// TYPEDEFS
  using mesh_config   = MeshConfig;
  using physics_type  = Physics;
  using mesh_type     = typename mesh::Tria<MeshConfig>;
  using mesh_boundary = typename result_of::mesh_boundary_set_t<MeshConfig>;
  using cell_dofs     = typename result_of::dof_map_t<MeshConfig>;

  using f_space = interpolation::FunctionSpace<MeshConfig, BcDIM>;

  protected:
  /// TYPEDEFS

  using geo_cell_metric_type =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             BcDIM>::cellwise_metric;
  using sol_cell_metric_type =
      typename interpolation::SolutionSpaceMetric<MeshConfig, BcDIM>::cellwise_metric;
  using flux_cell_metric_type =
      typename interpolation::FluxSpaceMetric<MeshConfig, Physics, BcDIM>::cellwise_metric;

  public:
  static std::string type_name()
  {
    return "RDMBCBase";
  }

  /// Default constructor
  /// Needed when we want to create a factory of boundary conditions
  /// Factory creates objects with 0 parameters
  RDMBCBase();

  /// Constructor
  RDMBCBase(const std::string &name);

  /// Destructor
  virtual ~RDMBCBase();

  /// Set the name of this bc
  virtual void set_name(const std::string &bc_name);

  /// Return the name of this bc
  virtual const std::string &name() const;

  /// Return the name of the active domain
  virtual const std::string domain_name() const;

  /// Return the type of the boundary condition
  virtual RDMBCType bc_type() const = 0;

  virtual void configure_mesh_data(const std::shared_ptr<const mesh_type> &mesh,
                                   const common::PtrHandle<const cell_dofs> geo_dofs,
                                   const common::PtrHandle<const cell_dofs> sol_dofs,
                                   const std::string &domain_name);

  /// Configure the solution function
  virtual void configure_spaces(const SFunc sf_type, const PointSetID quad_type,
                                const Uint quadrature_order);

  /// Configure the residual function
  virtual void set_solution(const interpolation::VectorMeshFunction<Real>::ptr solution);

  /// Configure the residual function
  virtual void set_residuals(const interpolation::VectorMeshFunction<Real>::ptr residuals);

  /// Set an analytical expression if needed for the bc
  virtual void set_expression(Real (*expr_ptr)(
      const math::DenseConstVecView<Real> &,
      const interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint));

  /// Set reference values if needed by BC
  virtual void set_reference_state(const math::DenseDVec<Real> &state);

  /// Set a parameter for the boundary conditions
  virtual void set_parameter(const std::string &param_name, const Real value);

  /// Return the dimension of the active domain
  virtual Uint dimension() const;

  /// Apply the boundary condition
  virtual void apply(RDTimeUpdate &time_update) = 0;

  /// Apply the boundary condition given cached data and metric
  /*
  virtual void compute_residuals(
      RDMBCMetricData<MeshConfig, Physics, BcDIM> const &data,
      typename mesh::BoundaryFacets<MeshConfig>::const_dof_iterator_typed
  const &sol_iter_begin, typename
  mesh::BoundaryFacets<MeshConfig>::const_dof_iterator_typed const
  &sol_iter_end, RDTimeUpdate &time_update, interpolation::SolutionCache
  &res_cache);
  */

  virtual void global_lsys_fin_diff_jacobian_entries(
      RDTimeUpdate &time_update, interpolation::SolutionCache &res_cache,
      interpolation::SolutionCache &res_cache_perturbed,
      std::vector<std::tuple<Uint, Uint, Real>> &mat_buffer,
      std::vector<std::tuple<Uint, Real>> &rhs_buffer);

  virtual void global_lsys_rhs_entries(RDTimeUpdate &time_update,
                                       interpolation::SolutionCache &res_cache,
                                       std::vector<std::tuple<Uint, Real>> &rhs_buffer);

  /// Print the value of stored parameters
  virtual void print_bc_parameters() const;

  protected:
  std::shared_ptr<const mesh_type> mesh() const;

  common::PtrHandle<const cell_dofs> geo_dofs() const;

  common::PtrHandle<const cell_dofs> sol_dofs() const;

  interpolation::VectorMeshFunction<Real>::ptr solution() const;

  interpolation::VectorMeshFunction<Real>::ptr residuals() const;

  typename f_space::ptr geo_space() const;

  typename f_space::ptr sol_space() const;

  /// Pointers to the functions holding solution, residuals and update
  /// coefficients

  private:
  /// Name of this bc
  std::string m_name;

  /// Name of the domain on which this bc is applied
  std::string m_domain_name;
  std::shared_ptr<const mesh_type> m_mesh;
  common::PtrHandle<const cell_dofs> m_geo_dofs;
  common::PtrHandle<const cell_dofs> m_sol_dofs;
  interpolation::VectorMeshFunction<Real>::ptr m_solution;
  interpolation::VectorMeshFunction<Real>::ptr m_residuals;

  /// Geometry space on the boundary
  typename f_space::ptr m_geo_space;

  /// Solution space on the boundary
  typename f_space::ptr m_sol_space;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
RDMBCBase<MeshConfig, Physics, BcDIM>::RDMBCBase() : m_name(""), m_domain_name("")
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
RDMBCBase<MeshConfig, Physics, BcDIM>::RDMBCBase(const std::string &name)
    : m_name(name), m_domain_name("")
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
RDMBCBase<MeshConfig, Physics, BcDIM>::~RDMBCBase()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::set_name(const std::string &bc_name)
{
  m_name = bc_name;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
const std::string &RDMBCBase<MeshConfig, Physics, BcDIM>::name() const
{
  return m_name;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
const std::string RDMBCBase<MeshConfig, Physics, BcDIM>::domain_name() const
{
  return m_domain_name;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::configure_mesh_data(
    const std::shared_ptr<const mesh_type> &mesh, const common::PtrHandle<const cell_dofs> geo_dofs,
    const common::PtrHandle<const cell_dofs> sol_dofs, const std::string &domain_name)
{
  m_mesh     = mesh;
  m_geo_dofs = geo_dofs;
  m_sol_dofs = sol_dofs;

  m_domain_name = domain_name;
  return;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::configure_spaces(const SFunc sf_type,
                                                             const PointSetID quad_type,
                                                             const Uint quadrature_order)
{
  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [quadrature_order, quad_type](const ElemShape shape,
                                                      const Uint elem_order) {
    return mesh::PointSetTag(shape, quadrature_order + 1, quad_type);
  };

  mesh::MeshBoundarySet<MeshConfig> const &mesh_boundaries = m_mesh->all_boundaries();
  mesh::BoundaryFacets<MeshConfig, BcDIM> const &geo_bdry =
      *(mesh_boundaries.domain(m_domain_name));

  m_geo_space = std::make_shared<f_space>();
  m_geo_space->set_reference_fe_values(geo_bdry, *m_geo_dofs, sf_generator, quad_generator);
  m_sol_space = std::make_shared<f_space>();
  m_sol_space->set_reference_fe_values(geo_bdry, *m_sol_dofs, sf_generator, quad_generator);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::set_solution(
    const interpolation::VectorMeshFunction<Real>::ptr solution)
{
  m_solution = solution;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::set_residuals(
    const interpolation::VectorMeshFunction<Real>::ptr residuals)
{
  m_residuals = residuals;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::set_expression(
    Real (*expr_ptr)(const math::DenseConstVecView<Real> &,
                     const interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint))
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::set_reference_state(const math::DenseDVec<Real> &state)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::set_parameter(const std::string &param_name,
                                                          const Real value)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
Uint RDMBCBase<MeshConfig, Physics, BcDIM>::dimension() const
{
  return BcDIM;
}

// ----------------------------------------------------------------------------

/*
template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::compute_residuals(
    RDMBCMetricData<MeshConfig, Physics, BcDIM> const &data,
    typename mesh::BoundaryFacets<MeshConfig>::const_dof_iterator_typed const
&sol_iter_begin, typename
mesh::BoundaryFacets<MeshConfig>::const_dof_iterator_typed const &sol_iter_end,
    RDTimeUpdate &time_update, interpolation::SolutionCache &res_cache)
{
}
*/

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::global_lsys_fin_diff_jacobian_entries(
    RDTimeUpdate &time_update, interpolation::SolutionCache &res_cache,
    interpolation::SolutionCache &res_cache_perturbed,
    std::vector<std::tuple<Uint, Uint, Real>> &mat_buffer,
    std::vector<std::tuple<Uint, Real>> &rhs_buffer)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::global_lsys_rhs_entries(
    RDTimeUpdate &time_update, interpolation::SolutionCache &res_cache,
    std::vector<std::tuple<Uint, Real>> &rhs_buffer)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCBase<MeshConfig, Physics, BcDIM>::print_bc_parameters() const
{
  // Default: no parameters, hence do nothing
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
std::shared_ptr<const typename RDMBCBase<MeshConfig, Physics, BcDIM>::mesh_type> RDMBCBase<
    MeshConfig, Physics, BcDIM>::mesh() const
{
  return m_mesh;
}

// ----------------------------------------------------------------------------
template <typename MeshConfig, typename Physics, Uint BcDIM>
common::PtrHandle<const typename RDMBCBase<MeshConfig, Physics, BcDIM>::cell_dofs> RDMBCBase<
    MeshConfig, Physics, BcDIM>::geo_dofs() const
{
  return m_geo_dofs;
}

// ----------------------------------------------------------------------------
template <typename MeshConfig, typename Physics, Uint BcDIM>
common::PtrHandle<const typename RDMBCBase<MeshConfig, Physics, BcDIM>::cell_dofs> RDMBCBase<
    MeshConfig, Physics, BcDIM>::sol_dofs() const
{
  return m_sol_dofs;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
interpolation::VectorMeshFunction<Real>::ptr RDMBCBase<MeshConfig, Physics, BcDIM>::solution() const
{
  return m_solution;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
interpolation::VectorMeshFunction<Real>::ptr RDMBCBase<MeshConfig, Physics, BcDIM>::residuals()
    const
{
  return m_residuals;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
typename RDMBCBase<MeshConfig, Physics, BcDIM>::f_space::ptr RDMBCBase<MeshConfig, Physics,
                                                                       BcDIM>::geo_space() const
{
  return m_geo_space;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
typename RDMBCBase<MeshConfig, Physics, BcDIM>::f_space::ptr RDMBCBase<MeshConfig, Physics,
                                                                       BcDIM>::sol_space() const
{
  return m_sol_space;
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
