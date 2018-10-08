#ifndef PDEKIT_Solver_L2_Projection_Local_hpp
#define PDEKIT_Solver_L2_Projection_Local_hpp

#include "PDEKit_Config.hpp"

#include <ctime>
#include <memory>
#include <unordered_map>

#include "interpolation/GeometryMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDMatArray.hpp"
#include "math/DenseDVec.hpp"
#include "mesh/Tria.hpp"
#include "mesh/point_set/StdPointSet.hpp"

namespace pdekit
{

namespace interpolation
{

class L2ProjectionLocal
{

  public:
  /// Constructor
  L2ProjectionLocal();

  /// Destructor
  ~L2ProjectionLocal();

  /// Clear all data
  void clear();

  /// Compute the elementwise projection operators
  /// This is a vector of dense matrices, one matrix per cell
  template <typename MeshConfig, typename PolySelectRule>
  void build_projection_operator(const mesh::Tria<MeshConfig> &mesh_in,
                                 const result_of::dof_map_t<MeshConfig> &cell_source_dofs,
                                 const result_of::dof_map_t<MeshConfig> &cell_target_dofs,
                                 const PolySelectRule &rule);

  /// Project data from one mesh onto another
  /// using previously assembled projection matrix
  /// We assume that the corresponding cells in the source and target mesh
  /// have the same ids
  template <typename MeshConfig, typename PolySelectRule>
  void project(const mesh::Tria<MeshConfig> &mesh_in,
               const result_of::dof_map_t<MeshConfig> &cell_source_dofs,
               const result_of::dof_map_t<MeshConfig> &cell_target_dofs, const PolySelectRule &rule,
               VectorMeshFunction<Real> const &f_src, VectorMeshFunction<Real> &f_tgt,
               const bool verbose = false);

  /// Print the structure of the projection matrix to file
  void save_projection_matrix_sparsity_pattern(const std::string &filename) const;

  private:
  /// TYPES
  struct UnsignedIntPairHasher
  {
    inline std::size_t operator()(const std::pair<Uint, Uint> &key) const
    {
      return key.first ^ key.second;
    }
  };

  /// Vector of all interpolation data
  math::DenseDMatArray<Real> m_local_ops;

  std::vector<Real> m_tmp_buffer_in;
  std::vector<Real> m_tmp_buffer_out;
};

// ----------------------------------------------------------------------------

void L2ProjectionLocal::clear()
{
  m_tmp_buffer_in.resize(0);
  m_tmp_buffer_out.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename PolySelectRule>
void L2ProjectionLocal::build_projection_operator(
    const mesh::Tria<MeshConfig> &mesh_in, const result_of::dof_map_t<MeshConfig> &cell_source_dofs,
    const result_of::dof_map_t<MeshConfig> &cell_target_dofs, const PolySelectRule &rule)
{
  using target_mass_map_type = std::unordered_map<Uint, math::DenseDMat<Real>>;
  using source_mass_map_type =
      std::unordered_map<std::pair<Uint, Uint>, math::DenseDMat<Real>, UnsignedIntPairHasher>;

  std::vector<mesh::DiscreteElemKey> elem_types_geo_qq;
  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_qq;

  std::vector<mesh::DiscreteElemKey> elem_types_geo_pq;
  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_q;
  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_p;

  // These maps store Vandermonde matrices for the computation of
  // the mass matrix Mqq in the 'target' space
  // The quadrature has to be exact for polynomials of order 2 * q
  target_mass_map_type target_mass_mat_map;
  target_mass_map_type target_mass_mat_inv_map;
  source_mass_map_type src_mass_mat_map;

  mesh::StdPointSet quad;

  clock_t start, end;
  Real elapsed;

  start = clock();

  // --------------------------------------------------------------------------
  // Loop over all cells going type by type in the TARGET MESH and compute
  // Lagrange shape function values for each type and given quadrature
  // --------------------------------------------------------------------------

  math::DenseDMat<Real> M;

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
    quad.change_type(eshape_geo, quad_order, PointSetID::Gauss);

    const mesh::PointSetTagExt std_reg_tag_ext(geo_cell_tag, quad_order,
                                               mesh::CellTransform::NO_TRANS, 0);
    const mesh::sf::SFTag geo_sf_tag(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                     ModalBasis::Modal);
    const mesh::PointSetTag quad_tag(eshape_geo, quad_order, PointSetID::Gauss);
    const mesh::PointSetTagExt quad_tag_ext(quad_tag, P0, mesh::CellTransform::NO_TRANS, 0);
    const mesh::DiscreteElemKey geo_key(std_reg_tag_ext, geo_sf_tag, quad_tag_ext);

    mesh::add_unique_discr_elem_key(elem_types_geo_qq, geo_key);

    common::PtrHandle<FEValues> fe_values_qq_ptr =
        elem_types_qq.std_region_data(mesh::PointSetTagExt(tgt_cell_tag, quad_order));

    if (fe_values_qq_ptr.is_null())
    {
      fe_values_qq_ptr = elem_types_qq.create(mesh::PointSetTagExt(tgt_cell_tag, quad_order));
      (*fe_values_qq_ptr).configure(tgt_cell_tag, std::get<1>(space_tags));
      (*fe_values_qq_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    const Uint key = tgt_cell.pt_set_id().store_value();

    target_mass_map_type::iterator map_it = target_mass_mat_map.find(key);
    if (map_it == target_mass_mat_map.end())
    {
      target_mass_mat_map.insert(std::make_pair(key, M));
      map_it = target_mass_mat_map.find(key);

      target_mass_mat_inv_map.insert(std::make_pair(key, M));
      target_mass_map_type::iterator map_inv_it = target_mass_mat_inv_map.find(key);

      const Uint nb_elem_vert = tgt_cell.nb_vert();
      map_it->second.resize(nb_elem_vert, nb_elem_vert);
      map_inv_it->second.resize(nb_elem_vert, nb_elem_vert);
    }
  } // Loop over cells

  // --------------------------------------------------------------------------
  // Prepare the internal data - resize transfer matrices from p to q, i.e.
  // rectangular 'mass matrices'
  // --------------------------------------------------------------------------
  const Uint nb_cells = cell_source_dofs.nb_active_cells();

  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> op_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>());
  op_shapes->resize(nb_cells);

  for (Uint c = 0; c < nb_cells; ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = mesh_in.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell_src_p = cell_source_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell_tgt_q = cell_target_dofs.active_cell(mesh::ActiveIdx(c));

    // Store the corresponding element types. Note that we might need to
    // integrate different polynomial order combinations in source and
    // target mesh (if the mesh was adapted for example), which needs that
    // for ONE ELEMENT SHAPE, WE HAVE TO PREPARE MULTIPLE VANDERMONDE
    // MATRICES WITH DIFFERENT QUADRATURE ORDERS!

    const mesh::PointSetTag geo_cell_tag = tcell.std_region().get().pt_set_id();
    const ElemShape eshape_geo           = geo_cell_tag.elem_shape();
    const Uint poly_order_geo            = geo_cell_tag.poly_order();

    const mesh::PointSetTag src_cell_p_tag = cell_src_p.pt_set_id();
    const ElemShape eshape_src             = src_cell_p_tag.elem_shape();

    const mesh::PointSetTag tgt_cell_q_tag = cell_tgt_q.pt_set_id();
    const ElemShape eshape_tgt             = tgt_cell_q_tag.elem_shape();

    // The quadrature has to have such order that phi_q * phi_p can be
    // integrated exactly
    const std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> space_tags =
        rule.sf_tags(geo_cell_tag, src_cell_p_tag, tgt_cell_q_tag);

    const std::tuple<Uint, Uint> quad_order_tags =
        rule.quadrature_orders(geo_cell_tag, src_cell_p_tag, tgt_cell_q_tag);

    const Uint quad_order_pq = std::get<0>(quad_order_tags);

    const mesh::PointSetTagExt std_reg_tag_ext(geo_cell_tag, quad_order_pq,
                                               mesh::CellTransform::NO_TRANS, 0);
    const mesh::sf::SFTag geo_sf_tag(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                     ModalBasis::Modal);

    const mesh::PointSetTag quad_tag(eshape_geo, quad_order_pq, PointSetID::Gauss);
    const mesh::PointSetTagExt quad_tag_ext(quad_tag, P0, mesh::CellTransform::NO_TRANS, 0);
    const mesh::DiscreteElemKey geo_key(std_reg_tag_ext, geo_sf_tag, quad_tag_ext);

    mesh::add_unique_discr_elem_key(elem_types_geo_pq, geo_key);

    // a) Make sure there is Vandermonde matrix for the 'source/p-type'
    // element
    //    Note that the second entry in the key is the COMBINED QUADRATURE
    //    ORDER!
    const mesh::PointSetTagExt key_src(src_cell_p_tag, quad_order_pq);

    common::PtrHandle<FEValues> fe_values_ptr_p = elem_types_p.std_region_data(key_src);
    if (fe_values_ptr_p.is_null())
    {
      quad.change_type(eshape_src, quad_order_pq, PointSetID::Gauss);

      fe_values_ptr_p = elem_types_p.create(key_src);
      (*fe_values_ptr_p).configure(src_cell_p_tag, std::get<0>(space_tags));
      (*fe_values_ptr_p).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // b) Make sure there is Vandermonde matrix for the 'target/q-type'
    // element
    //    Note that the second entry in the key is the COMBINED QUADRATURE
    //    ORDER!
    const mesh::PointSetTagExt key_tgt(cell_tgt_q.pt_set_id(), quad_order_pq);

    common::PtrHandle<FEValues> fe_values_ptr_q = elem_types_q.std_region_data(key_tgt);
    if (fe_values_ptr_q.is_null())
    {
      quad.change_type(eshape_tgt, quad_order_pq, PointSetID::Gauss);

      fe_values_ptr_q = elem_types_q.create(key_tgt);
      (*fe_values_ptr_q).configure(tgt_cell_q_tag, std::get<1>(space_tags));
      (*fe_values_ptr_q).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // c) Count how many entries there will be to hold transfer operators
    // for all cells
    (*op_shapes)[c] = common::ArrayShape<_2D, SUint>(cell_tgt_q.nb_vert(), cell_src_p.nb_vert());

    // d) Prepare storage for the rectangular 'mass' matrices M_pq
    const std::pair<Uint, Uint> key_combined =
        std::make_pair(cell_tgt_q.pt_set_id().store_value(), src_cell_p_tag.store_value());

    source_mass_map_type::iterator map_it = src_mass_mat_map.find(key_combined);

    if (map_it == src_mass_mat_map.end())
    {
      math::DenseDMat<Real> M;
      src_mass_mat_map.insert(std::make_pair(key_combined, M));
      map_it = src_mass_mat_map.find(key_combined);
      map_it->second.resize(cell_tgt_q.nb_vert(), cell_src_p.nb_vert());
    }
  }

  // --------------

  m_local_ops.allocate(std::move(op_shapes));

  // --------------

  GeometryCache<MeshConfig::GDIM> geo_cache_qq;
  GeometryCache<MeshConfig::GDIM> geo_cache_pq;
  GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric_qq;
  GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric_pq;

  geo_cache_qq.allocate(elem_types_geo_qq.cbegin(), elem_types_geo_qq.cend(), 1u);
  geo_metric_qq.allocate_buffer(elem_types_geo_qq.cbegin(), elem_types_geo_qq.cend(), 1u);

  geo_cache_pq.allocate(elem_types_geo_pq.cbegin(), elem_types_geo_pq.cend(), 1u);
  geo_metric_pq.allocate_buffer(elem_types_geo_pq.cbegin(), elem_types_geo_pq.cend(), 1u);

  // --------------

  for (Uint c = 0; c < nb_cells; ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = mesh_in.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell_src = cell_source_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell_tgt = cell_target_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::PointSetTag std_reg_tag_geo = tcell.std_region().get().pt_set_id();

    const Uint poly_order_geo = std_reg_tag_geo.poly_order();
    const Uint poly_order_p   = cell_src.pt_set_id().poly_order();
    const Uint poly_order_q   = cell_tgt.pt_set_id().poly_order();

    const Uint quad_order_qq = std::max(poly_order_geo, 2 * poly_order_q);
    const Uint quad_order_pq = std::max(poly_order_geo, poly_order_p + poly_order_q);

    // --------------------------------

    geo_cache_qq.flush();
    const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell.coordinates();

    const ElemShape eshape_geo = std_reg_tag_geo.elem_shape();
    const mesh::PointSetTagExt std_reg_tag_qq_ext(std_reg_tag_geo, quad_order_qq,
                                                  mesh::CellTransform::NO_TRANS, 0);
    const mesh::sf::SFTag geo_sf_tag_qq(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                        ModalBasis::Modal);
    const mesh::PointSetTag quad_tag_qq(eshape_geo, quad_order_qq, PointSetID::Gauss);
    const mesh::PointSetTagExt quad_tag_qq_ext(quad_tag_qq, P0, mesh::CellTransform::NO_TRANS, 0);
    const mesh::DiscreteElemKey geo_key_qq(std_reg_tag_qq_ext, geo_sf_tag_qq, quad_tag_qq_ext);

    geo_cache_qq.push_back_to_buffer(geo_cell_coords, geo_key_qq);

    geo_metric_qq.empty_buffer();
    geo_metric_qq.evaluate(geo_cache_qq, RebuildMetricIndex{true});

    // --------------------------------

    geo_cache_pq.flush();

    const mesh::PointSetTagExt std_reg_tag_pq_ext(std_reg_tag_geo, quad_order_pq,
                                                  mesh::CellTransform::NO_TRANS, 0);
    const mesh::sf::SFTag geo_sf_tag_pq(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                        ModalBasis::Modal);
    const mesh::PointSetTag quad_tag_pq(eshape_geo, quad_order_pq, PointSetID::Gauss);
    const mesh::PointSetTagExt quad_tag_pq_ext(quad_tag_pq, P0, mesh::CellTransform::NO_TRANS, 0);
    const mesh::DiscreteElemKey geo_key_pq(std_reg_tag_pq_ext, geo_sf_tag_pq, quad_tag_pq_ext);

    geo_cache_pq.push_back_to_buffer(geo_cell_coords, geo_key_pq);

    geo_metric_pq.empty_buffer();
    geo_metric_pq.evaluate(geo_cache_pq, RebuildMetricIndex{true});

    // --------------------------------

    const Uint key_src = cell_src.pt_set_id().store_value();
    const Uint key_tgt = cell_tgt.pt_set_id().store_value();

    const std::pair<Uint, Uint> key_combined = std::make_pair(key_tgt, key_src);

    source_mass_map_type::iterator src_map_it = src_mass_mat_map.find(key_combined);

    target_mass_map_type::iterator tgt_map_it     = target_mass_mat_map.find(key_tgt);
    target_mass_map_type::iterator tgt_map_inv_it = target_mass_mat_inv_map.find(key_tgt);

    math::DenseDMat<Real> &M_qp     = src_map_it->second;
    math::DenseDMat<Real> &M_qq     = tgt_map_it->second;
    math::DenseDMat<Real> &M_qq_inv = tgt_map_inv_it->second;

    const common::PtrHandle<FEValues> V_p_ptr =
        elem_types_p.std_region_data(mesh::PointSetTagExt(cell_src.pt_set_id(), quad_order_pq));
    const math::DenseDMat<Real> &V_p = (*V_p_ptr).Vandermonde();

    const common::PtrHandle<FEValues> V_q_ptr =
        elem_types_q.std_region_data(mesh::PointSetTagExt(cell_tgt.pt_set_id(), quad_order_pq));
    const math::DenseDMat<Real> &V_q = (*V_q_ptr).Vandermonde();

    const common::PtrHandle<FEValues> V_qq_ptr =
        elem_types_qq.std_region_data(mesh::PointSetTagExt(cell_tgt.pt_set_id(), quad_order_qq));
    const math::DenseDMat<Real> &V_qq = (*V_qq_ptr).Vandermonde();

    const Uint nb_vert_src = cell_src.nb_vert();
    const Uint nb_vert_tgt = cell_tgt.nb_vert();

    const typename GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>::cellwise_metric
        cell_metric_pq = geo_metric_pq.cellwise_values(0);

    // Build the rectangular matrix M_qp
    M_qp.fill(0.0);

    const Uint n_quad_pq                         = cell_metric_pq.nb_qd_pts();
    const math::DenseConstVecView<Real> j_det_pq = cell_metric_pq.jdet();

    const math::DenseDVec<Real> &w_pq = cell_metric_pq.pt_weights();

    // Build the matrix M_qq (the mass matrix of the element in 'target'
    // space

    M_qp.fill(0.0);

    for (Uint q = 0; q < n_quad_pq; ++q)
    {
      for (Uint i = 0; i < nb_vert_tgt; ++i)
      {
        for (Uint j = 0; j < nb_vert_src; ++j)
        {
          M_qp(i, j) += w_pq[q] * j_det_pq[q] * V_q(q, i) * V_p(q, j);
        }
      }
    }

    const typename GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>::cellwise_metric
        cell_metric_qq = geo_metric_qq.cellwise_values(0);

    const Uint n_quad_qq                         = cell_metric_qq.nb_qd_pts();
    const math::DenseConstVecView<Real> j_det_qq = cell_metric_qq.jdet();
    const math::DenseDVec<Real> &w_qq            = cell_metric_qq.pt_weights();

    // Mass matrix M_qq

    M_qq.fill(0.0);

    for (Uint q = 0; q < n_quad_qq; ++q)
    {
      for (Uint i = 0; i < nb_vert_tgt; ++i)
      {
        for (Uint j = 0; j < nb_vert_tgt; ++j)
        {
          M_qq(i, j) += w_qq[q] * j_det_qq[q] * V_qq(q, i) * V_qq(q, j);
        }
      }
    }

    M_qq.inv(M_qq_inv);
    math::DenseMatView<Real> M_transfer = m_local_ops.mat_view(c);
    M_transfer                          = M_qq_inv * M_qp;
  }

  end = clock();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "L2 local projection operator build: " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename PolySelectRule>
void L2ProjectionLocal::project(const mesh::Tria<MeshConfig> &mesh_in,
                                const result_of::dof_map_t<MeshConfig> &cell_source_dofs,
                                const result_of::dof_map_t<MeshConfig> &cell_target_dofs,
                                const PolySelectRule &rule, VectorMeshFunction<Real> const &f_src,
                                VectorMeshFunction<Real> &f_tgt, const bool verbose)
{
  if (f_src.nb_fields() != f_tgt.nb_fields())
  {
    std::cerr << "L2ProjectionLocal::project: number of fields in input and "
                 "output functions "
                 "differ. Aborting."
              << std::endl;
    return;
  }

  const Uint nb_fields = f_src.nb_fields();

  math::DenseConstMatView<Real> M_in;
  math::DenseMatView<Real> M_out;

  const Uint nb_cells = cell_source_dofs.nb_active_cells();

  for (Uint c = 0; c < nb_cells; ++c)
  {
    const mesh::MeshEntity cell_src = cell_source_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell_tgt = cell_target_dofs.active_cell(mesh::ActiveIdx(c));

    const Uint nb_vert_src = cell_src.nb_vert();
    const Uint nb_vert_tgt = cell_tgt.nb_vert();

    // Size of the input and output matrix should be [ nb_fields x
    // nb_vertices ]
    m_tmp_buffer_in.resize(nb_fields * nb_vert_src);
    m_tmp_buffer_out.resize(nb_fields * nb_vert_tgt);

    // Copy the data from input vector mesh function into input buffer
    for (Uint n = 0; n < nb_vert_src; ++n)
    {
      const VectorMeshFunction<Real>::const_entry_type node_value =
          f_src.const_value(cell_src.vertex(n));

      for (Uint f = 0; f < nb_fields; ++f)
        m_tmp_buffer_in[nb_fields * n + f] = node_value[f];
    }

    M_in = math::DenseConstMatView<Real>(m_tmp_buffer_in.data(), nb_fields, nb_vert_src, nb_fields);
    M_out = math::DenseMatView<Real>(m_tmp_buffer_out.data(), nb_fields, nb_vert_tgt, nb_fields);

    const math::DenseConstMatView<Real> M_transfer = m_local_ops.const_mat_view(c);
    M_out                                          = M_transfer * M_in;

    // Copy data from output buffer into output vector mesh function
    for (Uint n = 0; n < nb_vert_tgt; ++n)
    {
      VectorMeshFunction<Real>::entry_type node_value = f_tgt.value(cell_tgt.vertex(n));

      for (Uint f = 0; f < nb_fields; ++f)
      {
        node_value[f] = M_out(n, f);
      }
    }

  } // Loop over cells
}

// ----------------------------------------------------------------------------

void L2ProjectionLocal::save_projection_matrix_sparsity_pattern(const std::string &filename) const
{
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
