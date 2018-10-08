#ifndef RDM_Boundary_Condition_Projection_Data_hpp
#define RDM_Boundary_Condition_Projection_Data_hpp

#include "interpolation/SolutionSpaceMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "mesh/shape_function/SFTag.hpp"
#include "mesh/std_region/StdRegion.hpp"
#include "solver/ElementalMatrixOperator.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, Uint BcDim>
class RDMBCProjectionData
{
  public:
  /// Default constructor
  RDMBCProjectionData() = default;

  /// Do not allow copy construction
  RDMBCProjectionData(const RDMBCProjectionData &other) = delete;

  /// Default destructor
  ~RDMBCProjectionData() = default;

  /// Disable assignment operator
  RDMBCProjectionData &operator=(const RDMBCProjectionData &rhs) = delete;

  /// Prepare internal data
  template <typename GeoDofIterType, typename SolDofIterType>
  void setup_projection_lhs(const GeoDofIterType &geo_dofs_begin,
                            const std::vector<common::IteratorRange<SolDofIterType>> &sol_dofs);

  /// Prepare internal data
  template <typename GeoDofIterType, typename SolDofIterType>
  void setup_projection_rhs(
      const GeoDofIterType &geo_dofs_begin,
      const std::vector<common::IteratorRange<SolDofIterType>> &sol_dofs,
      interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, BcDim> const &geo_metric,
      interpolation::SolutionSpaceMetric<MeshConfig, BcDim> const &sol_metric,
      Real (*expression_ptr)(const math::DenseConstVecView<Real> &,
                             const interpolation::VectorMeshFunction<Real>::const_entry_type &,
                             const Uint));

  void project();

  const math::DenseConstMatView<Real> projected_values(const Uint idx) const;

  private:
  /// TYPES

  struct ProjectionPolyRule
  {
    // Return: a tuple consisting of
    // 1) Shape function tag of source shape function (which we're
    // projecting FROM) 2) Shape function tag of target polynomial (i.e
    // space we're projecting TO)
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

  /// DATA

  math::DenseDMatArray<Real> m_projection_ops;

  math::DenseDMatArray<Real> m_rhs;

  math::DenseDMatArray<Real> m_u;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
template <typename GeoDofIterType, typename SolDofIterType>
void RDMBCProjectionData<MeshConfig, BcDim>::setup_projection_lhs(
    const GeoDofIterType &geo_dofs_begin,
    const std::vector<common::IteratorRange<SolDofIterType>> &sol_dofs)
{
  ElementalMatrixOperator<MeshConfig, BcDim>::build_mass_matrix(
      geo_dofs_begin, sol_dofs, ProjectionPolyRule{}, InvertMassMat{true}, m_projection_ops);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
template <typename GeoDofIterType, typename SolDofIterType>
void RDMBCProjectionData<MeshConfig, BcDim>::setup_projection_rhs(
    const GeoDofIterType &geo_dofs_begin,
    const std::vector<common::IteratorRange<SolDofIterType>> &sol_dofs,
    interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, BcDim> const &geo_metric,
    interpolation::SolutionSpaceMetric<MeshConfig, BcDim> const &sol_metric,
    Real (*expression_ptr)(const math::DenseConstVecView<Real> &,
                           const interpolation::VectorMeshFunction<Real>::const_entry_type &,
                           const Uint))
{
  const Uint nb_fields = sol_metric.nb_fields();

  GeoDofIterType geo_dof_it = geo_dofs_begin;

  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> rhs_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>());

  for (const auto &sol_dof_range : sol_dofs)
  {
    for (SolDofIterType sol_dof_it = sol_dof_range.begin(); sol_dof_it != sol_dof_range.end();
         ++sol_dof_it)
    {
      const Uint nb_vert = sol_dof_it->mesh_entity().nb_vert();
      rhs_shapes->push_back(common::ArrayShape<_2D, SUint>(nb_vert, nb_fields));
    }
  }

  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> u_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>());

  u_shapes->resize(rhs_shapes->size());
  (*u_shapes) = (*rhs_shapes);

  m_rhs.allocate(std::move(rhs_shapes));
  m_u.allocate(std::move(u_shapes));

  using cell_geo_metric_t =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             BcDim>::cellwise_metric;
  using cell_sol_metric_t =
      typename interpolation::SolutionSpaceMetric<MeshConfig, BcDim>::cellwise_metric;

  Uint idx_in_metric = 0;

  for (const auto &sol_dof_range : sol_dofs)
  {
    for (SolDofIterType sol_dof_it = sol_dof_range.begin(); sol_dof_it != sol_dof_range.end();
         ++sol_dof_it)
    {
      const cell_geo_metric_t cell_geo_met = geo_metric.cellwise_values(idx_in_metric);
      const cell_sol_metric_t cell_sol_met = sol_metric.cellwise_values(idx_in_metric);

      math::DenseMatView<Real> local_rhs = m_rhs.mat_view(idx_in_metric);

      const Uint nb_qd_pts = cell_geo_met.nb_qd_pts();
      const Uint nb_nodes  = cell_sol_met.nb_dof_in_cell();

      const math::DenseConstVecView<Real> jac               = cell_geo_met.jdet();
      const math::DenseDVec<Real> &wgt                      = cell_geo_met.pt_weights();
      const math::DenseConstMatView<Real> qd_pts_phys_coord = cell_geo_met.interpolated_coords();

      const math::DenseDMat<Real> &V                    = cell_sol_met.reference_sf_values();
      const math::DenseConstMatView<Real> sol_at_qd_pts = cell_sol_met.field_values();

      local_rhs.fill(0.0);
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        const math::DenseConstVecView<Real> coord_at_qd_pt = qd_pts_phys_coord.row_transpose(q);
        const math::DenseConstVecView<Real> state_at_qd_pt = sol_at_qd_pts.row_transpose(q);

        const Real wj = jac[q] * wgt[q];

        for (Uint n = 0; n < nb_nodes; ++n)
        {
          for (Uint f = 0; f < nb_fields; ++f)
          {
            const Real fval = (*expression_ptr)(coord_at_qd_pt, state_at_qd_pt, f);
            local_rhs(n, f) += wj * fval * V(q, n);
          }
        }
      } // loop over quadrature points

      idx_in_metric++;

    } // iterator loop over one range
  }   // loop over all ranges

  // std::cout << "Number of entries to project on boundary = " <<
  // nb_entries_to_project << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
void RDMBCProjectionData<MeshConfig, BcDim>::project()
{
  for (Uint op = 0; op < m_projection_ops.size(); ++op)
  {
    math::DenseConstMatView<Real> A = m_projection_ops.const_mat_view(op);
    math::DenseConstMatView<Real> b = m_rhs.const_mat_view(op);
    math::DenseMatView<Real> u      = m_u.mat_view(op);

    u = A * b;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
const math::DenseConstMatView<Real> RDMBCProjectionData<MeshConfig, BcDim>::projected_values(
    const Uint idx) const
{
  return m_u.const_mat_view(idx);
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
