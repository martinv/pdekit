#ifndef PDEKIT_Solver_L2_Projection_Global_hpp
#define PDEKIT_Solver_L2_Projection_Global_hpp

#include "PDEKit_Config.hpp"

#include <ctime>
#include <memory>

#include "interpolation/GeometryMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "linear_system/LSTpetra.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseMatView.hpp"
#include "mesh/Tria.hpp"
#include "mesh/point_set/StdPointSet.hpp"

namespace pdekit
{

namespace interpolation
{

class L2ProjectionGlobal
{

  public:
  /// Constructor
  L2ProjectionGlobal();

  /// Destructor
  ~L2ProjectionGlobal();

  /// Compute the mass matrix for L2 projection
  /// This is a sparse matrix of dimension N x N, where
  /// N is the number of dof in the target mesh/space
  template <typename MeshConfig, typename PolySelectRule>
  void build_projection_operator(const mesh::Tria<MeshConfig> &mesh_in,
                                 const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
                                 const typename result_of::dof_map_t<MeshConfig> &cell_target_dofs,
                                 const PolySelectRule &rule);

  /// Project data from one mesh onto another
  /// using previously assembled projection matrix
  /// We assume that the corresponding cells in the source and target mesh
  /// have the same ids
  template <typename MeshConfig, typename PolySelectRule>
  void project(const mesh::Tria<MeshConfig> &mesh_in,
               const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
               const typename result_of::dof_map_t<MeshConfig> &cell_target_dofs,
               const PolySelectRule &rule, VectorMeshFunction<Real> const &f_src,
               VectorMeshFunction<Real> &f_tgt, const bool verbose = false);

  /// Print the structure of the projection matrix to file
  void save_projection_matrix_sparsity_pattern(const std::string &filename) const;

  private:
  /// DATA
  std::shared_ptr<ls::TpetraCrsMatrix<Real>> m_projection_matrix;

  std::shared_ptr<ls::TpetraMultiVector<Real>> m_rhs_vector;

  std::shared_ptr<ls::TpetraMultiVector<Real>> m_solution;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename PolySelectRule>
void L2ProjectionGlobal::build_projection_operator(
    const mesh::Tria<MeshConfig> &mesh_in,
    const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
    const typename result_of::dof_map_t<MeshConfig> &cell_target_dofs, const PolySelectRule &rule)
{
  if (m_projection_matrix)
  {
    m_projection_matrix.reset();
  }

  clock_t start, end;
  Real elapsed;

  using MeshType = mesh::Tria<MeshConfig>;

  using cell_dofs_type = typename MeshType::dof_storage_type;

  const Uint nb_nodes_in_sol_mesh = cell_target_dofs.nb_nodes();

  std::vector<mesh::DiscreteElemKey> elem_types_geo;
  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_sol;

  mesh::StdPointSet quad;

  // --------------------------------------------------------------------------
  // Loop over all cells and prepare the data for geometry metric
  // --------------------------------------------------------------------------

  for (Uint c = 0; c < mesh_in.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = mesh_in.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity src_cell = cell_source_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity tgt_cell = cell_target_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::PointSetTag geo_cell_tag = tcell.std_region().get().pt_set_id();
    const ElemShape eshape_geo           = geo_cell_tag.elem_shape();
    const Uint poly_order_geo            = geo_cell_tag.poly_order();

    const mesh::PointSetTag src_cell_tag = src_cell.pt_set_id();
    const mesh::PointSetTag tgt_cell_tag = tgt_cell.pt_set_id();

    const std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> space_tags =
        rule.sf_tags(geo_cell_tag, src_cell_tag, tgt_cell_tag);

    const std::tuple<Uint, Uint> quad_order_tags =
        rule.quadrature_orders(geo_cell_tag, src_cell_tag, tgt_cell_tag);

    // 1) Create FEValues object for each element type
    //    The quadrature must be such that we can integrate exactly the
    //    shape function products that define the entries in the mass
    //    matrix, but it must also be sufficient to integrate exactly on the
    //    geometrical support

    const Uint quad_order = std::get<1>(quad_order_tags);

    const mesh::PointSetTagExt geo_cell_tag_ext(geo_cell_tag, quad_order,
                                                mesh::CellTransform::NO_TRANS, 0);
    const mesh::sf::SFTag geo_sf(eshape_geo, SFunc::Lagrange, poly_order_geo, ModalBasis::Modal);
    const mesh::PointSetTag quad_tag(eshape_geo, quad_order, PointSetID::Gauss);
    const mesh::PointSetTagExt quad_tag_ext(quad_tag, P0, mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey geo_key(geo_cell_tag_ext, geo_sf, quad_tag_ext);
    mesh::add_unique_discr_elem_key(elem_types_geo, geo_key);

    quad.change_type(quad_tag);

    common::PtrHandle<FEValues> fe_values_sol_ptr =
        elem_types_sol.std_region_data(mesh::PointSetTagExt(tgt_cell.pt_set_id(), quad_order));

    if (fe_values_sol_ptr.is_null())
    {
      fe_values_sol_ptr =
          elem_types_sol.create(mesh::PointSetTagExt(tgt_cell.pt_set_id(), quad_order));
      (*fe_values_sol_ptr).configure(tgt_cell_tag, std::get<1>(space_tags));
      (*fe_values_sol_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

  } // Loop over all cells

  GeometryCache<MeshConfig::GDIM> geo_cache;
  geo_cache.allocate(elem_types_geo.cbegin(), elem_types_geo.cend(), 1);
  GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric;
  geo_metric.allocate_buffer(elem_types_geo.cbegin(), elem_types_geo.cend(), 1);

  // --------------------------------------------------------------------------
  // PREPARE THE STRUCTURE OF THE TARGET SPACE MASS MATRIX
  // --------------------------------------------------------------------------

  m_projection_matrix = std::make_shared<ls::TpetraCrsMatrix<Real>>(nb_nodes_in_sol_mesh);

  std::vector<Real> values, values_in_row;
  std::vector<Int> indices;

  // Loop over all cells going type by type and initialize the entries
  // in the stiffness matrix

  for (const typename cell_dofs_type::const_dof_range_typed &dof_group :
       cell_target_dofs.all_active_dof_groups())
  {
    const mesh::MeshEntity first_cell = (*dof_group.begin()).mesh_entity();
    // std::cout << "First cell in group = " << first_cell << std::endl;

    values_in_row.resize(first_cell.nb_vert());
    for (Uint v = 0; v < first_cell.nb_vert(); ++v)
    {
      values_in_row[v] = 0.0;
    }

    indices.resize(first_cell.nb_vert());

    for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = dof_group.begin();
         cell_iter != dof_group.end(); ++cell_iter)
    {
      const mesh::MeshEntity cell = cell_iter->mesh_entity();
      // std::cout << cell << std::endl;
      for (Uint vi = 0; vi < cell.nb_vert(); ++vi)
      {
        for (Uint vj = 0; vj < cell.nb_vert(); ++vj)
        {
          indices[vj] = cell.vertex(vj);
        }
        m_projection_matrix->insert_values_in_row(cell.vertex(vi), values_in_row, indices);
      }
    }
  } // Loop over all cell groups

  // m_projection_matrix->lock_structure();
  m_projection_matrix->fill(0.0);
  // m_projection_matrix->print_structure_to_file("sparsity_L2_projection.ps");

  // Start the actual assembly

  // This matrix block will serve as proxy to the 'values' vector
  math::DenseMatView<Real> local_stiffness;

  // --------------------------------------------------------------------------
  // SYSTEM ASSEMBLY - MAIN LOOP OVER ELEMENTS
  // --------------------------------------------------------------------------

  start = clock();

  // Loop over all cells and prepare the data for geometry metric

  for (Uint c = 0; c < mesh_in.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = mesh_in.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity src_cell = cell_source_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity tgt_cell = cell_target_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::PointSetTag geo_cell_tag = tcell.std_region().get().pt_set_id();
    const ElemShape eshape_geo           = geo_cell_tag.elem_shape();
    const Uint poly_order_geo            = geo_cell_tag.poly_order();

    const mesh::PointSetTag src_cell_tag = src_cell.pt_set_id();
    const mesh::PointSetTag tgt_cell_tag = tgt_cell.pt_set_id();

    const std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> space_tags =
        rule.sf_tags(geo_cell_tag, src_cell_tag, tgt_cell_tag);

    const std::tuple<Uint, Uint> quad_order_tags =
        rule.quadrature_orders(geo_cell_tag, src_cell_tag, tgt_cell_tag);

    const Uint quad_order = std::get<1>(quad_order_tags);

    values.resize(tgt_cell.nb_vert() * tgt_cell.nb_vert());
    values.assign(values.size(), 0.0);

    local_stiffness = math::DenseMatView<Real>(values.data(), tgt_cell.nb_vert(),
                                               tgt_cell.nb_vert(), tgt_cell.nb_vert());

    geo_cache.flush();
    const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell.coordinates();

    const mesh::PointSetTagExt geo_cell_tag_ext(geo_cell_tag, quad_order,
                                                mesh::CellTransform::NO_TRANS, 0);
    const mesh::sf::SFTag geo_sf(eshape_geo, SFunc::Lagrange, poly_order_geo, ModalBasis::Modal);
    const mesh::PointSetTag quad_tag(eshape_geo, quad_order, PointSetID::Gauss);
    const mesh::PointSetTagExt quad_tag_ext(quad_tag, P0, mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey geo_key(geo_cell_tag_ext, geo_sf, quad_tag_ext);

    geo_cache.push_back_to_buffer(geo_cell_coords, geo_key);

    geo_metric.empty_buffer();
    geo_metric.evaluate(geo_cache, interpolation::RebuildMetricIndex{true});

    const typename GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>::cellwise_metric cell_metric =
        geo_metric.cellwise_values(0);

    const common::PtrHandle<FEValues> V_sol_ptr =
        elem_types_sol.std_region_data(mesh::PointSetTagExt(tgt_cell.pt_set_id(), quad_order));
    const math::DenseDMat<Real> &V_sol = (*V_sol_ptr).Vandermonde();

    // Start assembly of the local stiffness matrix
    const Uint nb_qd_pts = cell_metric.nb_qd_pts();

    const math::DenseConstVecView<Real> j_det = cell_metric.jdet();
    const math::DenseDVec<Real> &w            = cell_metric.pt_weights();

    for (Uint q = 0; q < nb_qd_pts; ++q)
    {
      const Real wj = w[q] * j_det[q];
      for (Uint vi = 0; vi < tgt_cell.nb_vert(); ++vi)
      {
        for (Uint vj = 0; vj < tgt_cell.nb_vert(); ++vj)
        {
          local_stiffness(vi, vj) += wj * V_sol(q, vi) * V_sol(q, vj);
        } // Loop over vj
      }   // Loop over vi
    }     // Loop over quadrature points

    values_in_row.resize(tgt_cell.nb_vert());

    // Distribute to the global system

    for (Uint vi = 0; vi < tgt_cell.nb_vert(); ++vi)
    {
      for (Uint vj = 0; vj < tgt_cell.nb_vert(); ++vj)
      {
        indices[vj]       = tgt_cell.vertex(vj);
        values_in_row[vj] = local_stiffness(vi, vj);
      }
      m_projection_matrix->add_values_to_row(tgt_cell.vertex(vi), values_in_row, indices);
    }

  } // Loop over all cells

  m_projection_matrix->lock();

  end = clock();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "L2 projection matrix assembly: " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename PolySelectRule>
void L2ProjectionGlobal::project(const mesh::Tria<MeshConfig> &mesh_in,
                                 const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
                                 const typename result_of::dof_map_t<MeshConfig> &cell_target_dofs,
                                 const PolySelectRule &rule, VectorMeshFunction<Real> const &f_src,
                                 VectorMeshFunction<Real> &f_tgt, const bool verbose)
{
  // "source_mesh" is mesh FROM which we are projecting
  // Msrc is data matrix which is being projected ('source data')
  // Mtgt is the rhs of the linear system whose solution is the projected data

  clock_t start, end;
  Real elapsed;

  const Uint nb_nodes_in_tgt_mesh = cell_target_dofs.nb_nodes();

  if (m_rhs_vector == nullptr)
  {
    m_rhs_vector = std::make_shared<ls::TpetraMultiVector<Real>>(m_projection_matrix->domain_map(),
                                                                 f_src.nb_fields());
  }

  if (m_solution == nullptr)
  {
    m_solution = std::make_shared<ls::TpetraMultiVector<Real>>(m_projection_matrix->domain_map(),
                                                               f_src.nb_fields());
  }

  // The resize should maybe be done outside of the function
  f_tgt.resize(f_src.nb_fields(), nb_nodes_in_tgt_mesh);

  mesh::StdPointSet quad;

  std::vector<mesh::DiscreteElemKey> elem_types_geo;
  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_src;
  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_tgt;

  // --------------------------------------------------------------------------
  // Loop over all cells and prepare the data for geometry metric
  // --------------------------------------------------------------------------

  for (Uint c = 0; c < mesh_in.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = mesh_in.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity src_cell = cell_source_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity tgt_cell = cell_target_dofs.active_cell(mesh::ActiveIdx(c));

    // 1) Create FEValues object for each element type
    //    The quadrature must be such that we can integrate exactly the
    //    shape function products that define the entries in the mass
    //    matrix, but it must also be sufficient to integrate exactly on the
    //    geometrical support

    const mesh::PointSetTag geo_cell_tag = tcell.std_region().get().pt_set_id();
    const ElemShape eshape_geo           = geo_cell_tag.elem_shape();
    const Uint poly_order_geo            = geo_cell_tag.poly_order();

    const mesh::PointSetTag src_cell_tag = src_cell.pt_set_id();
    const mesh::PointSetTag tgt_cell_tag = tgt_cell.pt_set_id();

    const std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> space_tags =
        rule.sf_tags(geo_cell_tag, src_cell_tag, tgt_cell_tag);

    const std::tuple<Uint, Uint> quad_order_tags =
        rule.quadrature_orders(geo_cell_tag, src_cell_tag, tgt_cell_tag);

    const Uint quad_order = std::get<1>(quad_order_tags);

    const mesh::PointSetTagExt geo_pt_set_ext(geo_cell_tag, quad_order,
                                              mesh::CellTransform::NO_TRANS, 0);
    const mesh::sf::SFTag geo_sf(eshape_geo, SFunc::Lagrange, poly_order_geo, ModalBasis::Modal);
    const mesh::PointSetTag geo_quad_pt_set(eshape_geo, quad_order, PointSetID::Gauss);
    const mesh::PointSetTagExt geo_quad_pt_set_ext(geo_quad_pt_set, P0,
                                                   mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey geo_key(geo_pt_set_ext, geo_sf, geo_quad_pt_set_ext);
    mesh::add_unique_discr_elem_key(elem_types_geo, geo_key);

    quad.change_type(geo_quad_pt_set);

    common::PtrHandle<FEValues> fe_values_src_ptr =
        elem_types_src.std_region_data(mesh::PointSetTagExt(src_cell.pt_set_id(), quad_order));

    if (fe_values_src_ptr.is_null())
    {
      fe_values_src_ptr =
          elem_types_src.create(mesh::PointSetTagExt(src_cell.pt_set_id(), quad_order));
      (*fe_values_src_ptr).configure(src_cell_tag, std::get<0>(space_tags));
      (*fe_values_src_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    common::PtrHandle<FEValues> fe_values_tgt_ptr =
        elem_types_tgt.std_region_data(mesh::PointSetTagExt(tgt_cell.pt_set_id(), quad_order));

    if (fe_values_tgt_ptr.is_null())
    {
      fe_values_tgt_ptr =
          elem_types_tgt.create(mesh::PointSetTagExt(tgt_cell.pt_set_id(), quad_order));
      (*fe_values_tgt_ptr).configure(tgt_cell_tag, std::get<1>(space_tags));
      (*fe_values_tgt_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

  } // Loop over all cells

  GeometryCache<MeshConfig::GDIM> geo_cache;
  geo_cache.allocate(elem_types_geo.cbegin(), elem_types_geo.cend(), 1);
  GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric;
  geo_metric.allocate_buffer(elem_types_geo.cbegin(), elem_types_geo.cend(), 1);

  using nodal_value_type       = typename VectorMeshFunction<Real>::entry_type;
  using const_nodal_value_type = typename VectorMeshFunction<Real>::const_entry_type;

  // Reset the target data to zero
  f_tgt.fill(0.0);

  // Loop over all cells and and assemble the right-hand side
  (*m_solution).fill(0.0);
  (*m_rhs_vector).fill(0.0);

  std::vector<Real> local_rhs_raw_data;
  math::DenseMatView<Real> local_rhs;

  const Uint nb_fields = f_src.nb_fields();

  // -----------------------------------------
  // RHS ASSEMBLY - MAIN LOOP OVER ELEMENTS
  // -----------------------------------------

  start = clock();

  for (Uint c = 0; c < mesh_in.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = mesh_in.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity src_cell = cell_source_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity tgt_cell = cell_target_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::PointSetTag geo_cell_tag = tcell.std_region().get().pt_set_id();
    const mesh::PointSetTag src_cell_tag = src_cell.pt_set_id();
    const mesh::PointSetTag tgt_cell_tag = tgt_cell.pt_set_id();

    const std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> space_tags =
        rule.sf_tags(geo_cell_tag, src_cell_tag, tgt_cell_tag);

    const std::tuple<Uint, Uint> quad_order_tags =
        rule.quadrature_orders(geo_cell_tag, src_cell_tag, tgt_cell_tag);

    const Uint quad_order = std::get<1>(quad_order_tags);

    geo_cache.flush();
    const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell.coordinates();

    const ElemShape eshape_geo = geo_cell_tag.elem_shape();
    const Uint poly_order_geo  = geo_cell_tag.poly_order();

    const mesh::PointSetTagExt geo_pt_set_ext(geo_cell_tag, quad_order,
                                              mesh::CellTransform::NO_TRANS, 0);
    const mesh::sf::SFTag geo_sf(eshape_geo, SFunc::Lagrange, poly_order_geo, ModalBasis::Modal);
    const mesh::PointSetTag geo_quad_pt_set(eshape_geo, quad_order, PointSetID::Gauss);
    const mesh::PointSetTagExt geo_quad_pt_set_ext(geo_quad_pt_set, P0,
                                                   mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey geo_key(geo_pt_set_ext, geo_sf, geo_quad_pt_set_ext);

    geo_cache.push_back_to_buffer(geo_cell_coords, geo_key);

    geo_metric.empty_buffer();
    geo_metric.evaluate(geo_cache, interpolation::RebuildMetricIndex{true});

    const typename GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>::cellwise_metric cell_metric =
        geo_metric.cellwise_values(0);

    const common::PtrHandle<FEValues> V_src_ptr =
        elem_types_src.std_region_data(mesh::PointSetTagExt(src_cell.pt_set_id(), quad_order));
    const math::DenseDMat<Real> &V_src = (*V_src_ptr).Vandermonde();

    const common::PtrHandle<FEValues> V_tgt_ptr =
        elem_types_tgt.std_region_data(mesh::PointSetTagExt(tgt_cell.pt_set_id(), quad_order));
    const math::DenseDMat<Real> &V_tgt = (*V_tgt_ptr).Vandermonde();

    // Start assembly of the local stiffness matrix
    local_rhs_raw_data.resize(tgt_cell.nb_vert() * nb_fields);
    local_rhs_raw_data.assign(local_rhs_raw_data.size(), 0.0);
    local_rhs = math::DenseMatView<Real>(local_rhs_raw_data.data(), nb_fields, tgt_cell.nb_vert(),
                                         nb_fields);

    const Uint nb_qd_pts = cell_metric.nb_qd_pts();

    const math::DenseConstVecView<Real> j_det = cell_metric.jdet();
    const math::DenseDVec<Real> &w            = cell_metric.pt_weights();

    for (Uint q = 0; q < nb_qd_pts; ++q)
    {
      const Real wj = j_det[q] * w[q];
      for (Uint vi = 0; vi < tgt_cell.nb_vert(); ++vi)
      {
        for (Uint vj = 0; vj < src_cell.nb_vert(); ++vj)
        {
          const_nodal_value_type const nodal_value = f_src.const_value(src_cell.vertex(vj));
          const Real mass_matrix_factor            = wj * V_tgt(q, vi) * V_src(q, vj);

          // Accumulate multiple local RHS
          for (Uint f = 0; f < nb_fields; ++f)
          {
            local_rhs(vi, f) += mass_matrix_factor * nodal_value[f];
          }
        } // Loop over vj
      }   // Loop over vi
    }     // Loop over quadrature points

    for (Uint vi = 0; vi < tgt_cell.nb_vert(); ++vi)
    {
      for (Uint f = 0; f < nb_fields; ++f)
      {
        (*m_rhs_vector).add_value(tgt_cell.vertex(vi), local_rhs(vi, f), f);
      }
    }
  }

  end = clock();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  if (verbose)
  {
    std::cout << "L2 global projection operator build: " << elapsed << " s" << std::endl;
  }

  ls::LSTpetra<Real> lin_system;
  // lin_system.configure(m_projection_matrix, m_rhs_vector, m_solution, true,
  // false);
  lin_system.initialize_solver(m_projection_matrix, m_rhs_vector, m_solution, false);

  std::shared_ptr<ls::TrilinosPC<Real>> preconditioner = std::make_shared<ls::IfpackPC<Real>>();
  preconditioner->create("ILUT", m_projection_matrix);

  lin_system.connect_preconditioner(preconditioner);
  lin_system.update_after_mat_values_change(m_projection_matrix, m_rhs_vector, m_solution, false);
  lin_system.solve(550, 1.e-12, verbose);

  for (Uint i = 0; i < (*m_solution).size(); ++i)
  {
    nodal_value_type nodal_value = f_tgt.value(i);
    for (Uint f = 0; f < f_src.nb_fields(); ++f)
    {
      nodal_value[f] = (*m_solution).value(i, f);
    }
  }
}

// ----------------------------------------------------------------------------

void L2ProjectionGlobal::save_projection_matrix_sparsity_pattern(const std::string &filename) const
{
  m_projection_matrix->print_structure_to_file(filename);
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
