#ifndef PDEKIT_Solver_Residual_Restriction_Local_hpp
#define PDEKIT_Solver_Residual_Restriction_Local_hpp

#include "PDEKit_Config.hpp"

#include <ctime>
#include <memory>
#include <unordered_map>

#include "interpolation/GeometryMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseMatView.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"
#include "mesh/Tria.hpp"
#include "mesh/point_set/StdPointSet.hpp"

namespace pdekit
{

namespace interpolation
{

class ResidualRestrictionLocal
{

  public:
  /// Constructor
  ResidualRestrictionLocal();

  /// Destructor
  ~ResidualRestrictionLocal();

  /// Compute the mass matrix for L2 projection
  /// This is a sparse matrix of dimension N x N, where
  /// N is the number of dof in the target mesh/space
  template <typename MeshConfig>
  void build_projection_operator(const mesh::Tria<MeshConfig> &mesh_in,
                                 const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
                                 const typename result_of::dof_map_t<MeshConfig> &cell_target_dofs);

  /// Project data from one mesh onto another
  /// using previously assembled projection matrix
  /// We assume that the corresponding cells in the source and target mesh
  /// have the same ids
  template <typename MeshConfig>
  void project(const mesh::Tria<MeshConfig> &mesh_in,
               const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
               const typename result_of::dof_map_t<MeshConfig> &cell_target_dofs,
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

  /// Vector of all interpolation data held in one large continuous array
  std::vector<Real> m_raw_data;
  /// Vector of proxies to the data that enables to interpret them
  /// as matrices
  std::vector<math::DenseConstMatView<Real>> m_local_ops;

  std::vector<Real> m_tmp_buffer_in;
  std::vector<Real> m_tmp_buffer_out;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void ResidualRestrictionLocal::build_projection_operator(
    const mesh::Tria<MeshConfig> &mesh_in,
    const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
    const typename result_of::dof_map_t<MeshConfig> &cell_target_dofs)
{
  typedef std::unordered_map<Uint, math::DenseDMat<Real>> source_mass_map_type;
  typedef std::unordered_map<std::pair<Uint, Uint>, math::DenseDMat<Real>, UnsignedIntPairHasher>
      target_mass_map_type;

  std::vector<mesh::DiscreteElemKey> elem_types_geo_pp;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> elem_types_pp;

  std::vector<mesh::DiscreteElemKey> elem_types_geo_pq;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> elem_types_q;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> elem_types_p;

  // These map stores Vandermonde matrices for the computation of
  // the mass matrix Mpp in the 'source' space
  // The quadrature has to be exact for polynomials of order 2 * p
  source_mass_map_type src_mass_mat_map;
  source_mass_map_type src_mass_mat_inv_map;
  target_mass_map_type target_mass_mat_map;

  ElemShape eshape_geo, eshape_src, eshape_tgt;
  Uint poly_order_geo, poly_order_src, poly_order_tgt;
  PointSetID ref_topo_geo, ref_topo_src, ref_topo_tgt;

  mesh::StdPointSet quad;

  clock_t start, end;
  Real elapsed;

  start = clock();

  // --------------------------------------------------------------------------
  // Loop over all cells going type by type in the SOURCE MESH and compute
  // Lagrange shape function values for each type and given quadrature
  // --------------------------------------------------------------------------

  math::DenseDMat<Real> M;

  for (Uint c = 0; c < mesh_in.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = mesh_in.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity src_cell = cell_source_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::PointSetTag geo_std_reg_tag = tcell.std_region().get().pt_set_id();
    mesh::PointSetTag::decompose_into_fields(geo_std_reg_tag, eshape_geo, poly_order_geo,
                                             ref_topo_geo);
    mesh::PointSetTag::decompose_into_fields(src_cell.pt_set_id(), eshape_src, poly_order_src,
                                             ref_topo_src);

    // 1) Create FEValues object for each element type
    const Uint quad_order = std::max(2 * poly_order_src, poly_order_geo);
    quad.change_type(eshape_src, quad_order, PointSetID::Gauss);

    const mesh::PointSetTagExt geo_std_reg_tag_pp_ext(geo_std_reg_tag, quad_order);
    const mesh::sf::SFTag geo_basis_pp(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                       ModalBasis::Modal);
    const mesh::PointSetTag quad_pp_tag(eshape_src, quad_order, PointSetID::Gauss);
    const mesh::PointSetTagExt quad_pp_tag_ext(quad_pp_tag, P0, mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey key_pp(geo_std_reg_tag_pp_ext, geo_basis_pp, quad_pp_tag_ext);

    add_unique_discr_elem_key(elem_types_geo_pp, key_pp);

    // ---

    common::PtrHandle<FEValues> fe_values_pp_ptr =
        elem_types_pp.std_region_data(mesh::PointSetTagExt(src_cell.pt_set_id(), quad_order));

    if (fe_values_pp_ptr.is_null())
    {
      fe_values_pp_ptr =
          elem_types_pp.create(mesh::PointSetTagExt(src_cell.pt_set_id(), quad_order));
      (*fe_values_pp_ptr)
          .configure(src_cell.pt_set_id(), mesh::sf::SFTag(eshape_src, SFunc::Lagrange,
                                                           poly_order_src, ModalBasis::Modal));
      (*fe_values_pp_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    const Uint key = src_cell.pt_set_id().store_value();

    source_mass_map_type::iterator map_it = src_mass_mat_map.find(key);
    if (map_it == src_mass_mat_map.end())
    {
      src_mass_mat_map.insert(std::make_pair(key, M));
      map_it = src_mass_mat_map.find(key);

      src_mass_mat_inv_map.insert(std::make_pair(key, M));
      source_mass_map_type::iterator map_inv_it = src_mass_mat_inv_map.find(key);

      const Uint nb_elem_vert = src_cell.nb_vert();
      map_it->second.resize(nb_elem_vert, nb_elem_vert);
      map_inv_it->second.resize(nb_elem_vert, nb_elem_vert);
    }

  } // Loop over all cell groups

  // --------------------------------------------------------------------------
  // Prepare the internal data - resize transfer matrices from p to q, i.e.
  // rectangular 'mass matrices'
  // --------------------------------------------------------------------------
  const Uint nb_cells = mesh_in.nb_active_cells();

  Uint raw_data_size = 0;

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

    const mesh::PointSetTag geo_std_reg_tag = tcell.std_region().get().pt_set_id();

    mesh::PointSetTag::decompose_into_fields(geo_std_reg_tag, eshape_geo, poly_order_geo,
                                             ref_topo_geo);

    mesh::PointSetTag::decompose_into_fields(cell_src_p.pt_set_id(), eshape_src, poly_order_src,
                                             ref_topo_src);

    mesh::PointSetTag::decompose_into_fields(cell_tgt_q.pt_set_id(), eshape_tgt, poly_order_tgt,
                                             ref_topo_tgt);

    // The quadrature has to have such order that phi_q * phi_p can be
    // integrated exactly
    const Uint quad_order_pq = std::max(poly_order_src + poly_order_tgt, poly_order_geo);

    const mesh::PointSetTagExt geo_std_reg_tag_pq_ext(geo_std_reg_tag, quad_order_pq);
    const mesh::sf::SFTag geo_basis_pq(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                       ModalBasis::Modal);
    const mesh::PointSetTag quad_pq_tag(eshape_geo, quad_order_pq, PointSetID::Gauss);
    const mesh::PointSetTagExt quad_pq_tag_ext(quad_pq_tag, P0, mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey key_pq(geo_std_reg_tag_pq_ext, geo_basis_pq, quad_pq_tag_ext);

    add_unique_discr_elem_key(elem_types_geo_pq, key_pq);

    // a) Make sure there is Vandermonde matrix for the 'source/p-type'
    // element
    //    Note that the second entry in the key is the COMBINED QUADRATURE
    //    ORDER!
    const mesh::PointSetTagExt key_src(cell_src_p.pt_set_id(), quad_order_pq);

    common::PtrHandle<FEValues> fe_values_ptr_p = elem_types_p.std_region_data(key_src);
    if (fe_values_ptr_p.is_null())
    {
      quad.change_type(eshape_src, quad_order_pq, PointSetID::Gauss);

      fe_values_ptr_p = elem_types_p.create(key_src);
      (*fe_values_ptr_p)
          .configure(cell_src_p.pt_set_id(), mesh::sf::SFTag(eshape_src, SFunc::Lagrange,
                                                             poly_order_src, ModalBasis::Modal));
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
      (*fe_values_ptr_q)
          .configure(cell_tgt_q.pt_set_id(), mesh::sf::SFTag(eshape_tgt, SFunc::Lagrange,
                                                             poly_order_tgt, ModalBasis::Modal));
      (*fe_values_ptr_q).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // c) Count how many entries there will be to hold transfer operators
    // for all cells
    raw_data_size += cell_tgt_q.nb_vert() * cell_src_p.nb_vert();

    // d) Prepare storage for the rectangular 'mass' matrices M_pq
    const std::pair<Uint, Uint> key_combined =
        std::make_pair(cell_tgt_q.pt_set_id().store_value(), cell_src_p.pt_set_id().store_value());

    target_mass_map_type::iterator map_it = target_mass_mat_map.find(key_combined);

    if (map_it == target_mass_mat_map.end())
    {
      math::DenseDMat<Real> M;
      target_mass_mat_map.insert(std::make_pair(key_combined, M));
      map_it = target_mass_mat_map.find(key_combined);
      map_it->second.resize(cell_tgt_q.nb_vert(), cell_src_p.nb_vert());
    }
  }

  // --------------

  m_raw_data.resize(raw_data_size);
  m_local_ops.resize(nb_cells);

  // --------------

  GeometryCache<MeshConfig::GDIM> geo_cache_pp;
  GeometryCache<MeshConfig::GDIM> geo_cache_pq;
  GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric_pp;
  GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric_pq;

  geo_cache_pp.allocate(elem_types_geo_pp.cbegin(), elem_types_geo_pp.cend(), 1u);
  geo_metric_pp.allocate_buffer(elem_types_geo_pp.cbegin(), elem_types_geo_pp.cend(), 1u);

  geo_cache_pq.allocate(elem_types_geo_pq.cbegin(), elem_types_geo_pq.cend(), 1u);
  geo_metric_pq.allocate_buffer(elem_types_geo_pq.cbegin(), elem_types_geo_pq.cend(), 1u);

  // --------------

  raw_data_size = 0;

  for (Uint c = 0; c < nb_cells; ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = mesh_in.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell_src = cell_source_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell_tgt = cell_target_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::PointSetTag geo_std_reg_tag = tcell.std_region().get().pt_set_id();

    const Uint poly_order_geo = geo_std_reg_tag.poly_order();
    const Uint poly_order_p   = cell_src.pt_set_id().poly_order();
    const Uint poly_order_q   = cell_tgt.pt_set_id().poly_order();

    const Uint quad_order_pp = std::max(poly_order_geo, 2 * poly_order_p);
    const Uint quad_order_pq = std::max(poly_order_geo, poly_order_p + poly_order_q);

    // --------------------------------

    const ElemShape elem_shape = tcell.pt_set_id().elem_shape();

    geo_cache_pp.flush();
    const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell.coordinates();

    const mesh::PointSetTagExt geo_std_reg_tag_pp_ext(geo_std_reg_tag, quad_order_pp,
                                                      mesh::CellTransform::NO_TRANS, 0);
    const mesh::sf::SFTag geo_basis_pp(elem_shape, SFunc::Lagrange, poly_order_geo,
                                       ModalBasis::Modal);
    const mesh::PointSetTag quad_pp_tag(elem_shape, quad_order_pp, PointSetID::Gauss);
    const mesh::PointSetTagExt quad_pp_tag_ext(quad_pp_tag, P0, mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey key_pp(geo_std_reg_tag_pp_ext, geo_basis_pp, quad_pp_tag_ext);

    geo_cache_pp.push_back_to_buffer(geo_cell_coords, key_pp);
    geo_metric_pp.empty_buffer();
    geo_metric_pp.evaluate(geo_cache_pp, RebuildMetricIndex{true});

    // --------------------------------

    geo_cache_pq.flush();
    const mesh::PointSetTagExt tmp_key_pq(geo_std_reg_tag, quad_order_pq,
                                          mesh::CellTransform::NO_TRANS, 0);

    const mesh::PointSetTagExt geo_std_reg_tag_pq_ext(geo_std_reg_tag, quad_order_pq);
    const mesh::sf::SFTag geo_basis_pq(elem_shape, SFunc::Lagrange, poly_order_geo,
                                       ModalBasis::Modal);
    const mesh::PointSetTag quad_pq_tag(elem_shape, quad_order_pq, PointSetID::Gauss);
    const mesh::PointSetTagExt quad_pq_tag_ext(quad_pq_tag, P0, mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey key_pq(geo_std_reg_tag_pq_ext, geo_basis_pq, quad_pq_tag_ext);

    geo_cache_pq.push_back_to_buffer(geo_cell_coords, key_pq);

    geo_metric_pq.empty_buffer();
    geo_metric_pq.evaluate(geo_cache_pq, RebuildMetricIndex{true});

    // --------------------------------

    const Uint key_src = cell_src.pt_set_id().store_value();
    const Uint key_tgt = cell_tgt.pt_set_id().store_value();

    source_mass_map_type::iterator src_map_it     = src_mass_mat_map.find(key_src);
    source_mass_map_type::iterator src_map_inv_it = src_mass_mat_inv_map.find(key_src);

    const std::pair<Uint, Uint> key_combined = std::make_pair(key_tgt, key_src);

    target_mass_map_type::iterator tgt_map_it = target_mass_mat_map.find(key_combined);

    math::DenseDMat<Real> &M_qp     = tgt_map_it->second;
    math::DenseDMat<Real> &M_pp     = src_map_it->second;
    math::DenseDMat<Real> &M_pp_inv = src_map_inv_it->second;

    const common::PtrHandle<FEValues> V_p_ptr =
        elem_types_p.std_region_data(mesh::PointSetTagExt(cell_src.pt_set_id(), quad_order_pq));
    const math::DenseDMat<Real> &V_p = (*V_p_ptr).Vandermonde();

    const common::PtrHandle<FEValues> V_q_ptr =
        elem_types_q.std_region_data(mesh::PointSetTagExt(cell_tgt.pt_set_id(), quad_order_pq));
    const math::DenseDMat<Real> &V_q = (*V_q_ptr).Vandermonde();

    const common::PtrHandle<FEValues> V_pp_ptr =
        elem_types_pp.std_region_data(mesh::PointSetTagExt(cell_src.pt_set_id(), quad_order_pp));
    const math::DenseDMat<Real> &V_pp = (*V_pp_ptr).Vandermonde();

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
        cell_metric_pp = geo_metric_pp.cellwise_values(0);

    const Uint n_quad_pp                         = cell_metric_pp.nb_qd_pts();
    const math::DenseConstVecView<Real> j_det_pp = cell_metric_pp.jdet();
    const math::DenseDVec<Real> &w_pp            = cell_metric_pp.pt_weights();

    // Mass matrix M_pp

    M_pp.fill(0.0);

    for (Uint q = 0; q < n_quad_pp; ++q)
    {
      for (Uint i = 0; i < nb_vert_src; ++i)
      {
        for (Uint j = 0; j < nb_vert_src; ++j)
        {
          M_pp(i, j) += w_pp[q] * j_det_pp[q] * V_pp(q, i) * V_pp(q, j);
        }
      }
    }

    M_pp.inv(M_pp_inv);

    math::DenseMatView<Real> M_transfer = math::DenseMatView<Real>(
        m_raw_data.data() + raw_data_size, nb_vert_src, nb_vert_tgt, nb_vert_src);
    M_transfer = M_qp * M_pp_inv;

    m_local_ops[c] = math::DenseConstMatView<Real>(m_raw_data.data() + raw_data_size, nb_vert_src,
                                                   nb_vert_tgt, nb_vert_src);
  }

  end = clock();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "L2 local projection operator build: " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void ResidualRestrictionLocal::project(
    const mesh::Tria<MeshConfig> &mesh_in,
    const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
    const typename result_of::dof_map_t<MeshConfig> &cell_target_dofs,
    VectorMeshFunction<Real> const &f_src, VectorMeshFunction<Real> &f_tgt, const bool verbose)
{
  const Uint nb_tgt_nodes = cell_target_dofs.nb_nodes();

  std::vector<Uint> tgt_node_multiplicity(nb_tgt_nodes, 0u);

  if (f_src.nb_fields() != f_tgt.nb_fields())
  {
    std::cerr << "ResidualRestrictionLocal::project: number of fields in input "
                 "and output functions "
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

    M_out = m_local_ops[c] * M_in;

    // Copy data from output buffer into output vector mesh function
    for (Uint n = 0; n < nb_vert_tgt; ++n)
    {
      VectorMeshFunction<Real>::entry_type node_value = f_tgt.value(cell_tgt.vertex(n));

      for (Uint f = 0; f < nb_fields; ++f)
      {
        // Here we ACCUMULATE
        node_value[f] += M_out(n, f);
      }

      tgt_node_multiplicity[cell_tgt.vertex(n)]++;
    }

  } // Loop over cells

  for (Uint n = 0; n < f_tgt.nb_entries(); ++n)
  {
    VectorMeshFunction<Real>::entry_type node_value = f_tgt.value(n);
    for (Uint f = 0; f < nb_fields; ++f)
    {
      node_value[f] *= 1. / tgt_node_multiplicity[n];
    }
  }
}

// ----------------------------------------------------------------------------

void ResidualRestrictionLocal::save_projection_matrix_sparsity_pattern(
    const std::string &filename) const
{
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
