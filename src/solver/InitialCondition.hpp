#ifndef PDEKIT_Solver_Initial_Condition_hpp
#define PDEKIT_Solver_Initial_Condition_hpp

#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/DenseDVec.hpp"

namespace pdekit
{

namespace interpolation
{
template <typename MeshConfig, Uint DIM>
class FunctionSpace;
}

namespace solver
{

template <typename MeshConfig>
class InitialCondition
{
  public:
  /// Constructor
  InitialCondition(const std::string &name, typename mesh::Tria<MeshConfig>::const_shared_ptr mesh);

  /// Destructor
  virtual ~InitialCondition();

  /// Set on which domain should this ic act
  void set_domain(const Uint dim, const std::string &domain_name);

  /// Set an analytical expression if needed for the bc
  void set_expression(Real (*expr_ptr)(
      const math::DenseConstVecView<Real> &,
      const interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint));

  /// Return the name of this ic
  const std::string &name() const;

  /// Return the mesh on which this ic is applied
  const mesh::Tria<MeshConfig> &mesh() const;

  /// Return the dimension of the active domain
  Uint dim() const;

  /// Return the name of the active domain
  const std::string domain_name() const;

  /// Set initial state for the whole domain
  void apply_values(const std::string &dof_handler_name, const math::DenseDVec<Real> &initial_c,
                    interpolation::VectorMeshFunction<Real> &function) const;

  void apply_function(const std::string &dof_handler_name,
                      interpolation::VectorMeshFunction<Real> &function) const;

  protected:
  /// Name of this ic
  const std::string m_name;

  /// Mesh to which this ic should be applied
  mesh::Tria<MeshConfig> const &m_mesh;

  /// Dimension of domain on which this ic is applied
  Uint m_dim;

  /// Name of the domain on which this ic is applied
  std::string m_domain_name;

  /// Pointer to function which computes the initial condition values
  /// @param math::VectorBlock<Real> ... vector of point coordinates
  /// @param
  /// inteprolation::FunctionSpace<MeshConfig>::Function_t::const_block_type
  /// ... solution values
  ///        given point in case the Dirichlet BC depends on solution as well
  /// @param Uint ... component (when we're solving a system of equations)

  Real (*m_expression_ptr)(const math::DenseConstVecView<Real> &,
                           const interpolation::VectorMeshFunction<Real>::const_entry_type &,
                           const Uint);
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
InitialCondition<MeshConfig>::InitialCondition(
    const std::string &name, typename mesh::Tria<MeshConfig>::const_shared_ptr mesh)
    : m_name(name), m_mesh(*mesh), m_dim(_0D), m_domain_name("Undefined")
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
InitialCondition<MeshConfig>::~InitialCondition()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void InitialCondition<MeshConfig>::set_domain(const Uint dim, const std::string &domain_name)
{
  m_dim         = dim;
  m_domain_name = domain_name;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void InitialCondition<MeshConfig>::set_expression(Real (*expr_ptr)(
    const math::DenseConstVecView<Real> &,
    const typename interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint))
{
  m_expression_ptr = expr_ptr;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const std::string &InitialCondition<MeshConfig>::name() const
{
  return m_name;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const mesh::Tria<MeshConfig> &InitialCondition<MeshConfig>::mesh() const
{
  return m_mesh;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint InitialCondition<MeshConfig>::dim() const
{
  return m_dim;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const std::string InitialCondition<MeshConfig>::domain_name() const
{
  return m_domain_name;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void InitialCondition<MeshConfig>::apply_values(
    const std::string &dof_handler_name, const math::DenseDVec<Real> &initial_c,
    interpolation::VectorMeshFunction<Real> &function) const
{
  typedef typename result_of::dof_map_t<MeshConfig> cell_dofs_type;

  cell_dofs_type const &cell_dofs = (*m_mesh.dof_storage(dof_handler_name));
  // typename result_of::geometry<MeshConfig>::type const &coords =
  // m_mesh.geometry();

  std::map<Uint, std::string> cell_tag_names;

  cell_dofs.all_tag_names(cell_tag_names);

  Uint reference_tag = 0;
  for (std::map<Uint, std::string>::const_iterator tag_iter = cell_tag_names.begin();
       tag_iter != cell_tag_names.end(); ++tag_iter)
  {
    if (tag_iter->second == m_domain_name)
    {
      reference_tag = tag_iter->first;
      break;
    }
  }

  for (const typename cell_dofs_type::const_dof_range_typed &cell_group :
       cell_dofs.all_active_dof_groups())
  {
    for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = cell_group.begin();
         cell_iter != cell_group.end(); ++cell_iter)
    {
      if (cell_iter->cell_tag() == reference_tag)
      {
        const mesh::MeshEntity elem_nodes = cell_iter->mesh_entity();
        for (Uint n = 0; n < elem_nodes.nb_vert(); ++n)
        {
          // typename
          // result_of::geometry<MeshConfig>::type::const_row_type
          // const point =
          //    coords.const_row(elem_nodes.vertex(n));

          interpolation::VectorMeshFunction<Real>::entry_type u =
              function.value(elem_nodes.vertex(n));

          for (Uint component = 0; component < u.size(); ++component)
          {
            u[component] = initial_c[component];
          }
        } // Loop over nodes of one boundary element

      } // If this is an element on boundary
    }   // Loop over all cell groups with given dimension
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void InitialCondition<MeshConfig>::apply_function(
    const std::string &dof_handler_name, interpolation::VectorMeshFunction<Real> &function) const
{
  typedef typename result_of::dof_map_t<MeshConfig> cell_dofs_type;

  cell_dofs_type const &connectivity = *(m_mesh.dof_storage(dof_handler_name));

  std::map<Uint, std::string> cell_tag_names;

  (*m_mesh.dof_storage(dof_handler_name)).all_tag_names(cell_tag_names);

  Uint reference_tag = 0;
  for (std::map<Uint, std::string>::const_iterator tag_iter = cell_tag_names.begin();
       tag_iter != cell_tag_names.end(); ++tag_iter)
  {
    if (tag_iter->second == m_domain_name)
    {
      reference_tag = tag_iter->first;
      break;
    }
  }

  math::DenseDMat<Real> temp_nodal_value(function.nb_fields(), 1);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (const typename cell_dofs_type::const_dof_range_typed &dof_group :
       connectivity.all_active_dof_groups())
  {
    for (typename cell_dofs_type::const_dof_iterator_typed dof_iter = dof_group.begin();
         dof_iter != dof_group.end(); ++dof_iter)
    {
      if (dof_iter->cell_tag() == reference_tag)
      {
        const mesh::CellTopologyView<MeshConfig> tcell_view = dof_iter->tcell();
        const mesh::MeshEntity elem_nodes                   = dof_iter->mesh_entity();

        const math::DenseConstMatView<Real> elem_coords = loc_interpolator.transfer_coords(
            tcell_view.pt_set_id(), elem_nodes.pt_set_id(), tcell_view.coordinates());

        for (Uint n = 0; n < elem_nodes.nb_vert(); ++n)
        {
          math::DenseConstVecView<Real> const point = elem_coords.row_transpose(n);

          interpolation::VectorMeshFunction<Real>::entry_type u =
              function.value(elem_nodes.vertex(n));

          for (Uint component = 0; component < function.nb_fields(); ++component)
          {
            temp_nodal_value(component, 0) = u[component];
          }

          for (Uint component = 0; component < function.nb_fields(); ++component)
          {
            u[component] = (*m_expression_ptr)(point, temp_nodal_value.const_col(0), component);
          }
        } // Loop over nodes of one boundary element

      } // If this is an element on boundary
    }   // Loop over all cell groups with given dimension
  }
}

// ----------------------------------------------------------------------------

} // namespace solver

} // namespace pdekit

#endif
