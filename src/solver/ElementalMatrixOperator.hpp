#ifndef PDEKIT_Solver_Elemental_Matrix_Operator_hpp
#define PDEKIT_Solver_Elemental_Matrix_Operator_hpp

#include <unordered_map>

#include "common/IteratorRange.hpp"
#include "common/TaggedBool.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "math/DenseDMatArray.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "solver/SolverSetupAlgorithm.hpp"

namespace pdekit
{

namespace solver
{

// ----------------------------------------------------------------------------

struct _InvertMassMatTag
{
};
using InvertMassMat = common::TaggedBool<_InvertMassMatTag>;

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM = MeshConfig::TDIM>
class ElementalMatrixOperator
{
  public:
  ElementalMatrixOperator() = default;

  ElementalMatrixOperator(const ElementalMatrixOperator &other) = delete;

  ~ElementalMatrixOperator() = default;

  ElementalMatrixOperator &operator=(const ElementalMatrixOperator &other) = delete;

  template <typename GeoDofIterType, typename SolDofIterType, typename PolySelectRule>
  static void build_mass_matrix(const GeoDofIterType &geo_dofs_begin,
                                const std::vector<common::IteratorRange<SolDofIterType>> &sol_dofs,
                                const PolySelectRule &rule, const InvertMassMat inv_mat_tag,
                                math::DenseDMatArray<Real> &ops);

  template <typename GeoDofIterType, typename SolDofIterType, typename PolySelectRule,
            typename RhsAccessor>
  static void build_L2_proj_rhs(const GeoDofIterType &geo_dofs_begin,
                                const std::vector<common::IteratorRange<SolDofIterType>> &sol_dofs,
                                const PolySelectRule &rule, RhsAccessor &rhs_accessor,
                                const Uint nb_rhs_fields,
                                math::DenseDMatArray<Real> &projected_values);
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
template <typename GeoDofIterType, typename SolDofIterType, typename PolySelectRule>
void ElementalMatrixOperator<MeshConfig, DIM>::build_mass_matrix(
    const GeoDofIterType &geo_dofs_begin,
    const std::vector<common::IteratorRange<SolDofIterType>> &sol_dofs, const PolySelectRule &rule,
    const InvertMassMat inv_mat_flag, math::DenseDMatArray<Real> &ops)
{
  std::vector<mesh::DiscreteElemKey> elem_types_geo;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> elem_types_sol;
  using mass_mat_map_t = std::unordered_map<Uint, math::DenseDMat<Real>>;

  mass_mat_map_t tmp_mat_map;
  mass_mat_map_t inv_mass_mat_map;
  mesh::StdPointSet quad;

  auto start = std::chrono::system_clock::now();

  // --------------------------------------------------------------------------
  // Loop over all cells in the solution mesh
  // Prepare temporary data structures
  // --------------------------------------------------------------------------

  SolverSetupAlgorithm::generate_fe_pairs(geo_dofs_begin, sol_dofs, rule, elem_types_geo,
                                          elem_types_sol);

  math::DenseDMat<Real> M;

  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> op_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>());

  GeoDofIterType geo_iterator = geo_dofs_begin;

  for (const auto &range : sol_dofs)
  {
    for (auto sol_iterator = range.begin(); sol_iterator != range.end(); ++sol_iterator)
    {
      mesh::synchronize_dof_iterators(sol_iterator, geo_iterator);

      const mesh::MeshEntity geo_cell = geo_iterator->mesh_entity();
      const mesh::MeshEntity sol_cell = sol_iterator->mesh_entity();

      op_shapes->push_back(common::ArrayShape<_2D, SUint>(sol_cell.nb_vert(), sol_cell.nb_vert()));

      const mesh::PointSetTag sol_cell_tag = sol_cell.pt_set_id();

      if (inv_mat_flag)
      {
        const Uint key = sol_cell_tag.store_value();

        mass_mat_map_t::iterator inv_map_it = inv_mass_mat_map.find(key);
        if (inv_map_it == inv_mass_mat_map.end())
        {
          const Uint nb_elem_vert = sol_cell.nb_vert();

          inv_mass_mat_map.insert(std::make_pair(key, M));
          inv_map_it = inv_mass_mat_map.find(key);
          inv_map_it->second.resize(nb_elem_vert, nb_elem_vert);

          tmp_mat_map.insert(std::make_pair(key, M));
          mass_mat_map_t::iterator tmp_map_it = tmp_mat_map.find(key);
          tmp_map_it->second.resize(nb_elem_vert, nb_elem_vert);
        }
      }
    } // Loop over cells of one range
  }   // Loop over ranges

  // --------------

  ops.allocate(std::move(op_shapes));

  using geo_cache_t = interpolation::GeometryCache<GeoDofIterType::geo_dim()>;
  using geo_metric_t =
      interpolation::GeometryMetric<GeoDofIterType::geo_dim(), GeoDofIterType::traversal_dim()>;

  // interpolation::GeometryCache<MeshConfig::GDIM> geo_cache;
  geo_cache_t geo_cache;
  geo_metric_t geo_metric;

  // --------------

  geo_cache.allocate(elem_types_geo.cbegin(), elem_types_geo.cend(), 1u);
  geo_metric.allocate_buffer(elem_types_geo.cbegin(), elem_types_geo.cend(), 1u);

  // --------------

  Uint mat_idx = 0;

  for (const auto &range : sol_dofs)
  {
    for (auto sol_iterator = range.begin(); sol_iterator != range.end(); ++sol_iterator)
    {
      mesh::synchronize_dof_iterators(sol_iterator, geo_iterator);

      const mesh::MeshEntity geo_cell = geo_iterator->mesh_entity();
      const mesh::MeshEntity sol_cell = sol_iterator->mesh_entity();

      const mesh::PointSetTag geo_cell_tag = geo_cell.pt_set_id();
      const mesh::PointSetTag sol_cell_tag = sol_cell.pt_set_id();

      const Uint quad_order = rule.quadrature_order(geo_cell_tag, sol_cell_tag);

      // --------------------------------

      geo_cache.flush();
      const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = geo_iterator->cell_geometry();

      const mesh::DiscreteElemKey geo_key =
          SolverSetupAlgorithm::geo_key(geo_cell_tag, sol_cell_tag, rule);
      geo_cache.push_back_to_buffer(geo_cell_coords, geo_key);

      geo_metric.empty_buffer();
      geo_metric.evaluate(geo_cache, interpolation::RebuildMetricIndex{true});

      math::DenseMatView<Real> result_mat = ops.mat_view(mat_idx);

      const common::PtrHandle<interpolation::FEValues> V_sol_ptr =
          elem_types_sol.std_region_data(mesh::PointSetTagExt(sol_cell.pt_set_id(), quad_order));
      const math::DenseDMat<Real> &V_sol = (*V_sol_ptr).Vandermonde();

      const Uint nb_vert_sol = sol_cell.nb_vert();

      const typename geo_metric_t::cellwise_metric cell_metric_geo = geo_metric.cellwise_values(0);

      const Uint n_quad_pts                     = cell_metric_geo.nb_qd_pts();
      const math::DenseConstVecView<Real> j_det = cell_metric_geo.jdet();
      const math::DenseDVec<Real> &w            = cell_metric_geo.pt_weights();

      // Mass matrix M_qq

      if (inv_mat_flag)
      {
        const Uint key = sol_cell_tag.store_value();

        mass_mat_map_t::iterator inv_mass_map_it = inv_mass_mat_map.find(key);
        mass_mat_map_t::iterator tmp_mat_map_it  = tmp_mat_map.find(key);

        math::DenseDMat<Real> &inv_mat = inv_mass_map_it->second;
        math::DenseDMat<Real> &tmp_mat = tmp_mat_map_it->second;

        tmp_mat.fill(0.0);

        for (Uint q = 0; q < n_quad_pts; ++q)
        {
          for (Uint i = 0; i < nb_vert_sol; ++i)
          {
            for (Uint j = 0; j < nb_vert_sol; ++j)
            {
              tmp_mat(i, j) += w[q] * j_det[q] * V_sol(q, i) * V_sol(q, j);
            }
          }
        }
        tmp_mat.inv(inv_mat);

        for (Uint i = 0; i < nb_vert_sol; ++i)
        {
          for (Uint j = 0; j < nb_vert_sol; ++j)
          {
            result_mat(i, j) = inv_mat(i, j);
          }
        }
      }
      else
      {
        result_mat.fill(0.0);

        for (Uint q = 0; q < n_quad_pts; ++q)
        {
          for (Uint i = 0; i < nb_vert_sol; ++i)
          {
            for (Uint j = 0; j < nb_vert_sol; ++j)
            {
              result_mat(i, j) += w[q] * j_det[q] * V_sol(q, i) * V_sol(q, j);
            }
          }
        }
      } // else

      mat_idx++;
    } // Loop over elements of one range

  } // Loop over ranges

  auto end = std::chrono::system_clock::now();

  const int elapsed_miliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  std::cout.precision(10);
  std::cout.setf(std::ios::fixed);
  std::cout << "Elemental matrix builder: mass matrices built in: " << elapsed_miliseconds << " ms"
            << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
template <typename GeoDofIterType, typename SolDofIterType, typename PolySelectRule,
          typename RhsAccessor>
void ElementalMatrixOperator<MeshConfig, DIM>::build_L2_proj_rhs(
    const GeoDofIterType &geo_dofs_begin,
    const std::vector<common::IteratorRange<SolDofIterType>> &sol_dofs, const PolySelectRule &rule,
    RhsAccessor &rhs_accessor, const Uint nb_rhs_fields, math::DenseDMatArray<Real> &L2_rhs_values)
{
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> elem_types_geo;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> elem_types_sol;

  auto start = std::chrono::system_clock::now();

  // --------------------------------------------------------------------------
  // Loop over all cells in the solution mesh
  // Prepare temporary data structures
  // --------------------------------------------------------------------------

  SolverSetupAlgorithm::generate_fe_pairs(geo_dofs_begin, sol_dofs, rule, elem_types_geo,
                                          elem_types_sol);

  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> rhs_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>());

  GeoDofIterType geo_iterator = geo_dofs_begin;

  for (const auto &range : sol_dofs)
  {
    for (auto sol_iterator = range.begin(); sol_iterator != range.end(); ++sol_iterator)
    {
      mesh::synchronize_dof_iterators(sol_iterator, geo_iterator);

      const mesh::MeshEntity geo_cell = geo_iterator->mesh_entity();
      const mesh::MeshEntity sol_cell = sol_iterator->mesh_entity();

      rhs_shapes->push_back(common::ArrayShape<_2D, SUint>(sol_cell.nb_vert(), nb_rhs_fields));
    } // Loop over cells of one range
  }   // Loop over ranges

  // --------------

  L2_rhs_values.allocate(std::move(rhs_shapes));

  using geo_cache_t = interpolation::GeometryCache<GeoDofIterType::geo_dim()>;
  using geo_metric_t =
      interpolation::GeometryMetric<GeoDofIterType::geo_dim(), GeoDofIterType::traversal_dim()>;

  geo_cache_t geo_cache;
  geo_metric_t geo_metric;

  // --------------

  geo_cache.allocate(elem_types_geo.cbegin(), elem_types_geo.cend(), 1u);
  geo_metric.allocate_buffer(elem_types_geo.cbegin(), elem_types_geo.cend(), 1u);

  // --------------

  Uint pos_idx = 0;

  for (const auto &range : sol_dofs)
  {
    for (auto sol_iterator = range.begin(); sol_iterator != range.end(); ++sol_iterator)
    {
      mesh::synchronize_dof_iterators(sol_iterator, geo_iterator);

      const mesh::MeshEntity geo_cell = geo_iterator->mesh_entity();
      const mesh::MeshEntity sol_cell = sol_iterator->mesh_entity();

      const mesh::PointSetTag geo_cell_tag = geo_cell.pt_set_id();
      const mesh::PointSetTag sol_cell_tag = sol_cell.pt_set_id();

      const Uint quad_order = rule.quadrature_order(geo_cell_tag, sol_cell_tag);

      // --------------------------------

      geo_cache.flush();

      const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = geo_iterator->cell_geometry();

      // std::cout << geo_cell_coords << std::endl;

      const mesh::PointSetTagExt tmp_key(geo_cell_tag, quad_order);
      geo_cache.push_back_to_buffer(geo_cell_coords, tmp_key);

      geo_metric.empty_buffer();
      geo_metric.evaluate(geo_cache, interpolation::RebuildMetricIndex{true});

      const common::PtrHandle<interpolation::FEValues> V_sol_ptr =
          elem_types_sol.std_region_data(mesh::PointSetTagExt(sol_cell.pt_set_id(), quad_order));
      const math::DenseDMat<Real> &V_sol = (*V_sol_ptr).Vandermonde();

      const typename geo_metric_t::cellwise_metric cell_metric_geo = geo_metric.cellwise_values(0);

      const Uint n_quad_pts                     = cell_metric_geo.nb_qd_pts();
      const math::DenseConstVecView<Real> j_det = cell_metric_geo.jdet();
      const math::DenseDVec<Real> &w            = cell_metric_geo.pt_weights();

      const math::DenseConstMatView<Real> f_rhs_at_qd_pts =
          rhs_accessor(cell_metric_geo.interpolated_coords());

      math::DenseMatView<Real> result_mat = L2_rhs_values.mat_view(pos_idx);
      result_mat.fill(0.0);

      for (Uint q = 0; q < n_quad_pts; ++q)
      {
        const Real wj = j_det[q] * w[q];

        for (Uint dof_id = 0; dof_id < sol_cell.nb_vert(); ++dof_id)
        {
          for (Uint f = 0; f < nb_rhs_fields; ++f)
          {
            result_mat(dof_id, f) += wj * V_sol(q, dof_id) * f_rhs_at_qd_pts(q, f);
          }
        }
      } // loop over quadrature points

      pos_idx++;
    } // Loop over elements of one range

  } // Loop over ranges

  auto end = std::chrono::system_clock::now();

  const int elapsed_miliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  std::cout.precision(10);
  std::cout.setf(std::ios::fixed);
  std::cout << "Elemental matrix builder: L2 projection RHS matrices built in: "
            << elapsed_miliseconds << " ms" << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace solver

} // namespace pdekit

#endif
