#ifndef PDEKIT_RDM_Strong_Dirichlet_BC_hpp
#define PDEKIT_RDM_Strong_Dirichlet_BC_hpp

#include "math/DenseDMatArray.hpp"
#include "solver/ElementalMatrixOperator.hpp"
#include "solver/rdm/bc/RDMBCBase.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, Uint BcDim = Physics::DIM - 1>
class StrongDirichletBC : public RDMBCBase<MeshConfig, Physics, BcDim>
{
  public:
  using rdm_bc_base = RDMBCBase<MeshConfig, Physics, BcDim>;
  using mesh_type   = typename rdm_bc_base::mesh_type;
  using cell_dofs   = typename rdm_bc_base::cell_dofs;
  /// Constructor
  StrongDirichletBC();

  /// Constructor
  StrongDirichletBC(const std::string &name);

  /// Destructor
  ~StrongDirichletBC() override;

  /// Return the type of the boundary condition
  RDMBCType bc_type() const override;

  /// Set the function space (reference elements)
  void configure_mesh_data(const std::shared_ptr<const mesh_type> &mesh,
                           const common::PtrHandle<const cell_dofs> geo_dofs,
                           const common::PtrHandle<const cell_dofs> sol_dofs,
                           const std::string &domain_name) override;

  /// Set an analytical expression if needed for the bc
  void set_expression(Real (*expr_ptr)(
      const math::DenseConstVecView<Real> &,
      const interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint)) override;

  /// This method cannot be const, because it's not const in the virtual base
  /// class
  void apply(RDTimeUpdate &time_update) override;

  private:
  /// TYPES
  using const_dof_iterator = typename mesh::BoundaryFacets<MeshConfig, BcDim>::const_dof_iterator;
  using const_dof_iterator_typed =
      typename mesh::BoundaryFacets<MeshConfig, BcDim>::const_dof_iterator_typed;
  using sol_elem_range_typed = common::IteratorRange<const_dof_iterator_typed>;

  /// Helper class to compute RHS for L2 projection of BCs
  ///
  class DirichletVals
  {
public:
    DirichletVals(Real (*expression_ptr)(
        const math::DenseConstVecView<Real> &,
        const interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint));

    ~DirichletVals() = default;

    const math::DenseConstMatView<Real> operator()(
        const math::DenseConstMatView<Real> &phys_coords);

    std::vector<Real> m_work;

    Real (*m_expression_ptr)(const math::DenseConstVecView<Real> &,
                             const interpolation::VectorMeshFunction<Real>::const_entry_type &,
                             const Uint);
  };

  struct PolyRule
  {
    mesh::sf::SFTag sf_tag(const mesh::PointSetTag &geo_tag, const mesh::PointSetTag &sol_tag) const
    {
      const mesh::sf::SFTag poly_space_tag(sol_tag.elem_shape(), SFunc::Lagrange,
                                           sol_tag.poly_order(), ModalBasis::Modal);

      return poly_space_tag;
    }

    // Return: a tuple consisting of
    // 1) Quadrature order to be used in order to construct mass matrix
    // operator
    Uint quadrature_order(const mesh::PointSetTag &geo_tag, const mesh::PointSetTag &sol_tag) const
    {
      const Uint quad_order = std::max(2 * sol_tag.poly_order(), geo_tag.poly_order());
      return quad_order;
    }
  };

  /// METHODS

  /// Compute and store Dirichlet values on given boundary
  void compute_Dirichlet_data();

  /// Compute the Dirichlet boundary condition on a single boundary face
  void compute_on_element_strong(const mesh::MeshEntity &bdry_facet,
                                 const mesh::CellGeometry<MeshConfig::GDIM> &facet_coords,
                                 const math::DenseMatView<Real> &curr_u_values,
                                 math::DenseMatView<Real> &new_u_values) const;

  /// DATA
  /// Pointer to function which computes the value of Dirichlet BC
  /// @param math::VectorBlock<Real> ... vector of point coordinates
  /// @param inteprolation::FunctionSpace<MeshConfig>::Function_t::column_type
  /// ... solution values
  ///        given point in case the Dirichlet BC depends on solution as well
  /// @param Uint ... component (when we're solving a system of equations)

  Real (*m_expression_ptr)(const math::DenseConstVecView<Real> &,
                           const interpolation::VectorMeshFunction<Real>::const_entry_type &,
                           const Uint);

  std::vector<sol_elem_range_typed> m_sol_dof_ranges;

  math::DenseDMatArray<Real> m_Dirichlet_data;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
StrongDirichletBC<MeshConfig, Physics, BcDim>::StrongDirichletBC() : rdm_bc_base()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
StrongDirichletBC<MeshConfig, Physics, BcDim>::StrongDirichletBC(const std::string &name)
    : rdm_bc_base(name)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
StrongDirichletBC<MeshConfig, Physics, BcDim>::~StrongDirichletBC()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
RDMBCType StrongDirichletBC<MeshConfig, Physics, BcDim>::bc_type() const
{
  return BC_TYPE_STRONG;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void StrongDirichletBC<MeshConfig, Physics, BcDim>::configure_mesh_data(
    const std::shared_ptr<const mesh_type> &mesh, const common::PtrHandle<const cell_dofs> geo_dofs,
    const common::PtrHandle<const cell_dofs> sol_dofs, const std::string &domain_name)
{
  rdm_bc_base::configure_mesh_data(mesh, geo_dofs, sol_dofs, domain_name);

  // Store ranges of solution space iterators
  mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry =
      *(rdm_bc_base::mesh()->all_boundaries().domain(this->domain_name()));
  bdry.all_bdry_dof_ranges(*sol_dofs, m_sol_dof_ranges);

  // Allocate memory for Dirichlet data
  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> mat_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>());
  const Uint nb_fields = Physics::NEQ;

  for (const auto &range : m_sol_dof_ranges)
  {
    for (const_dof_iterator_typed it = range.begin(); it != range.end(); ++it)
    {
      const mesh::MeshEntity bdry_facet = it->mesh_entity();
      mat_shapes->push_back(common::ArrayShape<_2D, SUint>(bdry_facet.nb_vert(), nb_fields));
    }
  } // Loop over boundary cells

  m_Dirichlet_data.allocate(std::move(mat_shapes));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void StrongDirichletBC<MeshConfig, Physics, BcDim>::set_expression(
    Real (*expr_ptr)(const math::DenseConstVecView<Real> &,
                     const interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint))
{
  m_expression_ptr = expr_ptr;
  compute_Dirichlet_data();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void StrongDirichletBC<MeshConfig, Physics, BcDim>::apply(RDTimeUpdate &time_update)
{
  using cell_dofs_type                              = typename result_of::dof_map_t<MeshConfig>;
  interpolation::VectorMeshFunction<Real> &solution = *(this->solution());

  const Uint nb_fields                           = solution.nb_fields();
  std::shared_ptr<const mesh_type> const mesh_in = rdm_bc_base::mesh();

  Uint idx = 0;
  for (const auto &range : m_sol_dof_ranges)
  {
    for (const_dof_iterator_typed it = range.begin(); it != range.end(); ++it)
    {
      const mesh::MeshEntity bdry_facet = it->mesh_entity();
      const Uint nb_facet_dof           = bdry_facet.nb_vert();

      const math::DenseConstMatView<Real> facet_Dirichlet_data =
          m_Dirichlet_data.const_mat_view(idx);

      for (Uint dof_id = 0; dof_id < nb_facet_dof; ++dof_id)
      {
        interpolation::VectorMeshFunction<Real>::entry_type u =
            solution.value(bdry_facet.vertex(dof_id));
        for (Uint f = 0; f < nb_fields; ++f)
        {
          u[f] = facet_Dirichlet_data(dof_id, f);
        }
      }
      idx++;
    } // Loop over active boundary cells
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void StrongDirichletBC<MeshConfig, Physics, BcDim>::compute_Dirichlet_data()
{
#if 1
  interpolation::VectorMeshFunction<Real> &solution = *(this->solution());

  const Uint nb_fields                           = solution.nb_fields();
  std::shared_ptr<const mesh_type> const mesh_in = rdm_bc_base::mesh();

  std::vector<Real> tmp_values_in;

  Uint idx = 0;

  for (const auto &range : m_sol_dof_ranges)
  {
    for (const_dof_iterator_typed it = range.begin(); it != range.end(); ++it)
    {
      const mesh::MeshEntity bdry_facet                     = it->mesh_entity();
      const mesh::CellGeometry<MeshConfig::GDIM> facet_geom = it->cell_geometry();

      const Uint nb_facet_dof = bdry_facet.nb_vert();
      tmp_values_in.resize(nb_facet_dof * nb_fields);

      math::DenseMatView<Real> view_in(tmp_values_in.data(), nb_fields, nb_facet_dof, nb_fields);
      math::DenseMatView<Real> view_out = m_Dirichlet_data.mat_view(idx);

      for (Uint dof_id = 0; dof_id < nb_facet_dof; ++dof_id)
      {
        interpolation::VectorMeshFunction<Real>::const_entry_type u =
            solution.const_value(bdry_facet.vertex(dof_id));
        for (Uint f = 0; f < nb_fields; ++f)
        {
          view_in(dof_id, f) = u[f];
        }
      }

      compute_on_element_strong(bdry_facet, facet_geom, view_in, view_out);
      idx++;

    } // Loop over active boundary cells
  }   // Loop over ranges

#endif

#if 0
  // ----------------
  // I) L2 projection
  // ----------------
  using cell_dofs_type = typename result_of::dof_map_t<MeshConfig>;

  mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry =
      *(rdm_bc_base::mesh()->all_boundaries().domain(this->domain_name()));

  cell_dofs_type const &geo_dofs = *rdm_bc_base::geo_dofs();
  // cell_dofs_type const &sol_dofs = *rdm_bc_base::sol_dofs();

  const_dof_iterator geo_dof_it = bdry.cbegin(geo_dofs);

  math::DenseDMatArray<Real> bdry_mass_mats;
  math::DenseDMatArray<Real> L2_rhs_values;
  DirichletVals rhs_accessor(m_expression_ptr);

  ElementalMatrixOperator<MeshConfig, BcDim> L2_operator;
  L2_operator.build_mass_matrix(geo_dof_it, m_sol_dof_ranges, PolyRule{}, InvertMassMat{ true },
                                bdry_mass_mats);
  L2_operator.build_L2_proj_rhs(geo_dof_it, m_sol_dof_ranges, PolyRule{}, rhs_accessor,
                                Physics::NEQ, L2_rhs_values);

  // Resize the container for Dirichlet data
  // m_Dirichlet_data.reshape(L2_rhs_values);

  for (Uint c = 0; c < bdry_mass_mats.size(); ++c)
  {
    const math::DenseConstMatView<Real> inv_M = bdry_mass_mats.const_mat_view(c);
    const math::DenseConstMatView<Real> rhs = L2_rhs_values.const_mat_view(c);
    math::DenseMatView<Real> Dir_bc_data = m_Dirichlet_data.mat_view(c);

    Dir_bc_data = inv_M * rhs;
  }

#endif

#if 0
  // ----------------
  // II) Collocation
  // ----------------

  using cell_dofs_type = typename result_of::dof_map_t<MeshConfig>;

  mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> const &bdry =
      *(rdm_bc_base::mesh()->all_boundaries().domain(this->domain_name()));

  cell_dofs_type const &geo_dofs = *rdm_bc_base::geo_dofs();
  // cell_dofs_type const &sol_dofs = *rdm_bc_base::sol_dofs();

  const_dof_iterator geo_it = bdry.cbegin(geo_dofs);
  DirichletVals rhs_accessor(m_expression_ptr);

  using V_map_t = common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>;
  using inv_V_map_t = common::DataMap<mesh::PointSetTagExt, math::DenseDMat<Real>>;
  V_map_t map_bdry_V;
  inv_V_map_t map_bdry_inv_V;

  mesh::StdRegion bdry_std_reg;
  PolyRule poly_rule;

  std::vector<Real> facet_coord_data;

  Uint idx = 0;

  for (const auto &range : m_sol_dof_ranges)
  {
    for (const_dof_iterator_typed sol_it = range.begin(); sol_it != range.end(); ++sol_it)
    {
      mesh::synchronize_dof_iterators(sol_it, geo_it);

      const mesh::MeshEntity bdry_facet_geo = geo_it->mesh_entity();
      const mesh::MeshEntity bdry_facet_sol = sol_it->mesh_entity();

      const mesh::PointSetTag facet_geo_tag = bdry_facet_geo.std_region_id();
      const mesh::PointSetTag facet_sol_tag = bdry_facet_sol.std_region_id();

      const Uint quad_order = poly_rule.quadrature_order(facet_geo_tag, facet_sol_tag);
      const mesh::sf::SFTag sf_tag = poly_rule.sf_tag(facet_geo_tag, facet_sol_tag);

      const mesh::PointSetTagExt tag_ext(facet_sol_tag, quad_order, mesh::CellTransform::NO_TRANS,
                                         0);

      common::PtrHandle<interpolation::FEValues> V_handle = map_bdry_V.std_region_data(tag_ext);
      common::PtrHandle<math::DenseDMat<Real>> inv_V_handle =
          map_bdry_inv_V.std_region_data(tag_ext);

      // If cached data for this standard region has not been computed yet, add it
      // to corresponding maps
      if (V_handle.is_null())
      {
        V_handle = map_bdry_V.create(tag_ext);
        (*V_handle).configure(facet_sol_tag, sf_tag);
        bdry_std_reg.change_type(facet_sol_tag);
        const math::DenseDMat<Real> &std_reg_pts = bdry_std_reg.get().coordinates();
        math::DenseDVec<Real> weights;
        weights.resize(std_reg_pts.rows());
        weights.fill(1.0);
        (*V_handle).fill_Vandermonde(std_reg_pts, weights);

        const math::DenseDMat<Real> &V = (*V_handle).Vandermonde();

        inv_V_handle = map_bdry_inv_V.create(tag_ext);
        (*inv_V_handle).resize(V.rows(), V.cols());
        V.inv(*inv_V_handle);
      }

      const Uint nb_facet_dofs = bdry_facet_sol.nb_vert();
      facet_coord_data.resize(nb_facet_dofs * MeshConfig::GDIM);

      math::DenseMatView<Real> coord_view(facet_coord_data.data(), MeshConfig::GDIM, nb_facet_dofs,
                                          MeshConfig::GDIM);

      const mesh::CellGeometry<MeshConfig::GDIM> sol_dof_coords = sol_it->cell_geometry();

      for (Uint n = 0; n < nb_facet_dofs; ++n)
      {
        const math::DenseConstVecView<Real> node_coord = sol_dof_coords.const_node_view(n);
        for (Uint d = 0; d < MeshConfig::GDIM; ++d)
        {
          coord_view(n, d) = node_coord[d];
        }
      }

      math::DenseMatView<Real> computed_Dir_data_view = m_Dirichlet_data.mat_view(idx);
      math::DenseConstMatView<Real> coord_view_const(facet_coord_data.data(), MeshConfig::GDIM,
                                                     nb_facet_dofs, MeshConfig::GDIM);

      const math::DenseConstMatView<Real> g = rhs_accessor(coord_view_const);
      computed_Dir_data_view = (*inv_V_handle) * g;

      idx++;
    } // Loop over one range
  }   // Loop over boundary ranges
#endif
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void StrongDirichletBC<MeshConfig, Physics, BcDim>::compute_on_element_strong(
    const mesh::MeshEntity &bdry_facet, const mesh::CellGeometry<MeshConfig::GDIM> &facet_coords,
    const math::DenseMatView<Real> &curr_u_values, math::DenseMatView<Real> &new_u_values) const
{
  for (Uint n = 0; n < bdry_facet.nb_vert(); ++n)
  {
    std::vector<Real> buffer;

    const math::DenseVecView<const Real> const_pt = facet_coords.const_node_view(n);
    // Workaround - the function m_expression_ptr takes a
    // math::DenseConstVecView<T> as its first argument, but
    // facet_coords.const_node_view(n) returns math::DenseVecView<const T>
    // We therefore temporarily copy the data into a buffer and then
    // generate a new const view with suitable type THIS SHOULD BE FIXED!
    buffer.resize(const_pt.size());
    for (Uint i = 0; i < const_pt.size(); ++i)
    {
      buffer[i] = const_pt[i];
    }
    const math::DenseConstVecView<Real> const_pt_buffer(buffer.data(), buffer.size());
    const math::DenseConstVecView<Real> curr_dof_value = curr_u_values.row_transpose(n);

    for (Uint component = 0; component < Physics::NEQ; ++component)
    {
      new_u_values(n, component) = (*m_expression_ptr)(const_pt_buffer, curr_dof_value, component);
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
StrongDirichletBC<MeshConfig, Physics, BcDim>::DirichletVals::DirichletVals(Real (*expression_ptr)(
    const math::DenseConstVecView<Real> &,
    const interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint))
    : m_expression_ptr(expression_ptr)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
const math::DenseConstMatView<Real> StrongDirichletBC<MeshConfig, Physics, BcDim>::
    DirichletVals::operator()(const math::DenseConstMatView<Real> &phys_coords)
{
  const Uint nb_pts    = phys_coords.rows();
  const Uint nb_fields = Physics::NEQ;
  m_work.resize(nb_pts, nb_fields);

  // Empty view of values u at quadrature points
  // This should be filled properly ...
  const math::DenseConstVecView<Real> uq;

  math::DenseMatView<Real> tmp(m_work.data(), nb_fields, nb_pts, nb_fields);

  for (Uint q = 0; q < nb_pts; ++q)
  {
    const math::DenseConstVecView<Real> const_pt = phys_coords.row_transpose(q);
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      tmp(q, eq) = (*m_expression_ptr)(const_pt, uq, eq);
    }
  }
  math::DenseConstMatView<Real> result(m_work.data(), nb_fields, nb_pts, nb_fields);
  return result;
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
