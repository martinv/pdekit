#ifndef PDEKIT_Interpolation_Coarse_Scale_Correction_Op_Builder_hpp
#define PDEKIT_Interpolation_Coarse_Scale_Correction_Op_Builder_hpp

#include <ctime>
#include <memory>
#include <unordered_map>

#include "interpolation/GeometryMetric.hpp"
#include "linear_system/LocalTransferOps.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "mesh/Tria.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"

namespace pdekit
{

namespace interpolation
{

class CoarseScaleCorrectionOpBuilder
{

  public:
  /// Constructor
  CoarseScaleCorrectionOpBuilder();

  /// Destructor
  ~CoarseScaleCorrectionOpBuilder();

  /// Build an array of matrices for uniform restriction to
  /// lower-order polynomial space
  template <typename MeshConfig>
  static void build_restriction_ops(
      const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
      const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
      const std::vector<Uint> &cell_reordering, const Uint q,
      ls::LocalTransferOps<Real> &local_restriction_ops);

  /// Build an array of matrices for uniform restriction to
  /// lower-order polynomial space
  template <typename MeshConfig>
  static void build_prolongation_ops(
      const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
      const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
      const std::vector<Uint> &cell_reordering, const Uint q,
      ls::LocalTransferOps<Real> &local_prolongation_ops);

  private:
  /// TYPES
  struct UnsignedIntPairHasher
  {
    inline std::size_t operator()(const std::pair<Uint, Uint> &key) const
    {
      return key.first ^ key.second;
    }
  };

  /// Build an array of matrices for uniform restriction to
  /// lower-order polynomial space. The transfer operators are
  /// built using L2 projection.
  template <typename MeshConfig>
  static void build_restriction_ops_L2(
      const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
      const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
      const std::vector<Uint> &cell_reordering, const Uint q,
      ls::LocalTransferOps<Real> &local_restriction_ops);

  /// Build an array of matrices for uniform restriction to
  /// lower-order polynomial space. The transfer operators are
  /// built using L2 projection.
  template <typename MeshConfig>
  static void build_prolongation_ops_L2(
      const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
      const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
      const std::vector<Uint> &cell_reordering, const Uint q,
      ls::LocalTransferOps<Real> &local_prolongation_ops);

  /// Build an array of matrices for uniform restriction to
  /// lower-order polynomial space. The transfer operators are
  /// built using modal basis
  template <typename MeshConfig>
  static void build_restriction_ops_modal(
      const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
      const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
      const std::vector<Uint> &cell_reordering, const Uint q,
      ls::LocalTransferOps<Real> &local_restriction_ops);

  /// Build an array of matrices for uniform prolongation to
  /// high-order polynomial space. The transfer operators are
  /// built using modal basis
  template <typename MeshConfig>
  static void build_prolongation_ops_modal(
      const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
      const typename result_of::dof_map_t<MeshConfig> &cell_source_dofs,
      const std::vector<Uint> &cell_reordering, const Uint q,
      ls::LocalTransferOps<Real> &local_restriction_ops);
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void CoarseScaleCorrectionOpBuilder::build_restriction_ops(
    const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &cell_sol_dofs,
    const std::vector<Uint> &cell_reordering, const Uint q,
    ls::LocalTransferOps<Real> &local_restriction_ops)
{
  /*
  build_restriction_ops_L2<MeshConfig>(NEQ, cell_geo_dofs, cell_sol_dofs,
  cell_reordering, q, local_restriction_ops);
  */

  build_restriction_ops_modal<MeshConfig>(NEQ, cell_geo_dofs, cell_sol_dofs, cell_reordering, q,
                                          local_restriction_ops);
}

template <typename MeshConfig>
void CoarseScaleCorrectionOpBuilder::build_prolongation_ops(
    const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &cell_sol_dofs,
    const std::vector<Uint> &cell_reordering, const Uint p,
    ls::LocalTransferOps<Real> &local_prolongation_ops)
{
  /*
  build_prolongation_ops_L2<MeshConfig>(NEQ, cell_geo_dofs, cell_sol_dofs,
  cell_reordering, p, local_prolongation_ops);
  */

  build_prolongation_ops_modal<MeshConfig>(NEQ, cell_geo_dofs, cell_sol_dofs, cell_reordering, p,
                                           local_prolongation_ops);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void CoarseScaleCorrectionOpBuilder::build_restriction_ops_L2(
    const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &cell_sol_dofs,
    const std::vector<Uint> &cell_reordering, const Uint q,
    ls::LocalTransferOps<Real> &local_restriction_ops)
{
  using target_mass_map_type = std::unordered_map<Uint, math::DenseDMat<Real>>;
  using source_mass_map_type =
      std::unordered_map<std::pair<Uint, Uint>, math::DenseDMat<Real>, UnsignedIntPairHasher>;

  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_geo_qq;
  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_qq;

  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_geo_pq;
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
  mesh::StdRegion tgt_cell_q_std_reg;

  const Uint nb_cells = cell_sol_dofs.nb_active_cells();

  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> op_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>());
  op_shapes->resize(nb_cells);

  for (Uint c = 0; c < nb_cells; ++c)
  {
    const mesh::MeshEntity cell_geo   = cell_geo_dofs.active_cell(c);
    const mesh::MeshEntity cell_src_p = cell_sol_dofs.active_cell(c);

    const mesh::PointSetTag geo_cell_tag = cell_geo.pt_set_id();
    const ElemShape eshape_geo           = geo_cell_tag.elem_shape();
    const Uint poly_order_geo            = geo_cell_tag.poly_order();

    const mesh::PointSetTag src_cell_p_tag = cell_src_p.pt_set_id();
    const ElemShape eshape_sol             = src_cell_p_tag.elem_shape();
    const Uint poly_order_src              = src_cell_p_tag.poly_order();
    // For the moment we assume the reference topology in 'source'
    // and 'target' element to be the same
    const PointSetID ref_topo_tgt = src_cell_p_tag.ref_topology();
    const mesh::PointSetTag tgt_cell_q_tag(eshape_sol, q, ref_topo_tgt);
    tgt_cell_q_std_reg.change_type(tgt_cell_q_tag);
    const Uint tgt_cell_nb_vert = tgt_cell_q_std_reg.get().nb_nodes();

    // ----------------------------------------------------
    // I) Process the 'target' cell data - shape function
    // values and mass matrix with quadrature order
    // 2*q (aka 'qq')
    // ----------------------------------------------------

    const Uint quad_order_qq = std::max(2 * q, poly_order_geo);
    quad.change_type(eshape_geo, quad_order_qq, PointSetID::Gauss);

    common::PtrHandle<FEValues> fe_values_geo_qq_ptr =
        elem_types_geo_qq.std_region_data(mesh::PointSetTagExt(geo_cell_tag, quad_order_qq));

    // Geometry expansion for target cell type
    if (fe_values_geo_qq_ptr.is_null())
    {
      fe_values_geo_qq_ptr =
          elem_types_geo_qq.create(mesh::PointSetTagExt(geo_cell_tag, quad_order_qq));
      (*fe_values_geo_qq_ptr)
          .configure(geo_cell_tag, mesh::sf::SFTag(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                                   ModalBasis::Modal));
      (*fe_values_geo_qq_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // Solution expansion for target cell type
    common::PtrHandle<FEValues> fe_values_qq_ptr =
        elem_types_qq.std_region_data(mesh::PointSetTagExt(tgt_cell_q_tag, quad_order_qq));

    if (fe_values_qq_ptr.is_null())
    {
      fe_values_qq_ptr = elem_types_qq.create(mesh::PointSetTagExt(tgt_cell_q_tag, quad_order_qq));
      (*fe_values_qq_ptr)
          .configure(tgt_cell_q_tag,
                     mesh::sf::SFTag(eshape_sol, SFunc::Lagrange, q, ModalBasis::Modal));
      (*fe_values_qq_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    const Uint key = tgt_cell_q_tag.store_value();

    target_mass_map_type::iterator target_mass_map_it = target_mass_mat_map.find(key);
    if (target_mass_map_it == target_mass_mat_map.end())
    {
      target_mass_mat_map.insert(std::make_pair(key, M));
      target_mass_map_it = target_mass_mat_map.find(key);

      target_mass_mat_inv_map.insert(std::make_pair(key, M));
      target_mass_map_type::iterator map_inv_it = target_mass_mat_inv_map.find(key);

      target_mass_map_it->second.resize(tgt_cell_nb_vert, tgt_cell_nb_vert);
      map_inv_it->second.resize(tgt_cell_nb_vert, tgt_cell_nb_vert);
    }

    // ----------------------------------------------------
    // II) Process the 'source' cell data - shape function
    // values and mass matrix with quadrature order
    // p+q (aka 'pq')
    // ----------------------------------------------------

    // The quadrature has to have such order that phi_q * phi_p can be
    // integrated exactly
    const Uint quad_order_pq = std::max(poly_order_src + q, poly_order_geo);

    // Geometry expansion for source cell type
    common::PtrHandle<FEValues> fe_values_geo_pq_ptr =
        elem_types_geo_pq.std_region_data(mesh::PointSetTagExt(geo_cell_tag, quad_order_pq));

    if (fe_values_geo_pq_ptr.is_null())
    {
      quad.change_type(eshape_geo, quad_order_pq, PointSetID::Gauss);

      fe_values_geo_pq_ptr =
          elem_types_geo_pq.create(mesh::PointSetTagExt(geo_cell_tag, quad_order_pq));
      (*fe_values_geo_pq_ptr)
          .configure(geo_cell_tag, mesh::sf::SFTag(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                                   ModalBasis::Modal));
      (*fe_values_geo_pq_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // a) Make sure there is Vandermonde matrix for the 'source/p-type'
    // element
    //    Note that the second entry in the key is the COMBINED QUADRATURE
    //    ORDER!
    const mesh::PointSetTagExt key_src(src_cell_p_tag, quad_order_pq);

    common::PtrHandle<FEValues> fe_values_p_ptr = elem_types_p.std_region_data(key_src);
    if (fe_values_p_ptr.is_null())
    {
      quad.change_type(eshape_sol, quad_order_pq, PointSetID::Gauss);

      fe_values_p_ptr = elem_types_p.create(key_src);
      (*fe_values_p_ptr)
          .configure(src_cell_p_tag, mesh::sf::SFTag(eshape_sol, SFunc::Lagrange, poly_order_src,
                                                     ModalBasis::Modal));
      (*fe_values_p_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // b) Make sure there is Vandermonde matrix for the 'target/q-type'
    // element
    //    Note that the second entry in the key is the COMBINED QUADRATURE
    //    ORDER!
    const mesh::PointSetTagExt key_tgt(tgt_cell_q_tag, quad_order_pq);

    common::PtrHandle<FEValues> fe_values_ptr_q = elem_types_q.std_region_data(key_tgt);
    if (fe_values_ptr_q.is_null())
    {
      quad.change_type(eshape_sol, quad_order_pq, PointSetID::Gauss);

      fe_values_ptr_q = elem_types_q.create(key_tgt);
      (*fe_values_ptr_q)
          .configure(tgt_cell_q_tag,
                     mesh::sf::SFTag(eshape_sol, SFunc::Lagrange, q, ModalBasis::Modal));
      (*fe_values_ptr_q).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // c) Count how many entries there will be to hold transfer operators
    // for all cells
    (*op_shapes)[cell_reordering[c]] =
        common::ArrayShape<_2D, SUint>(NEQ * tgt_cell_nb_vert, NEQ * cell_src_p.nb_vert());

    // d) Prepare storage for the rectangular 'mass' matrices M_pq
    const std::pair<Uint, Uint> key_combined =
        std::make_pair(tgt_cell_q_tag.store_value(), src_cell_p_tag.store_value());

    source_mass_map_type::iterator map_it = src_mass_mat_map.find(key_combined);

    if (map_it == src_mass_mat_map.end())
    {
      math::DenseDMat<Real> M;
      src_mass_mat_map.insert(std::make_pair(key_combined, M));
      map_it = src_mass_mat_map.find(key_combined);
      map_it->second.resize(tgt_cell_nb_vert, cell_src_p.nb_vert());
    }

  } // Loop over cells

  // --------------

  std::unique_ptr<std::vector<Real>> mat_values(new std::vector<Real>());
  std::unique_ptr<std::vector<Uint>> mat_offsets(new std::vector<Uint>());

  mat_offsets->reserve(nb_cells + 1);
  mat_offsets->push_back(0);

  Uint tot_mat_storage = 0;
  for (const auto shape : (*op_shapes))
  {
    const Uint mat_block_size = shape.size(0) * shape.size(1);
    mat_offsets->push_back(mat_block_size);
    tot_mat_storage += mat_block_size;
  }

  for (Uint i = 1; i < (*mat_offsets).size(); ++i)
  {
    (*mat_offsets)[i] += (*mat_offsets)[i - 1];
  }

  mat_values->resize(tot_mat_storage);

  // --------------

  GeometryCache<MeshConfig::GDIM> geo_cache_qq;
  GeometryCache<MeshConfig::GDIM> geo_cache_pq;
  GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric_qq;
  GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric_pq;

  geo_cache_qq.allocate(elem_types_geo_qq, 1u);
  geo_metric_qq.allocate_buffer(elem_types_geo_qq, 1u);

  geo_cache_pq.allocate(elem_types_geo_pq, 1u);
  geo_metric_pq.allocate_buffer(elem_types_geo_pq, 1u);

  // --------------

  std::vector<Real> workspace;

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < nb_cells; ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_sol_dofs.tcell(c);
    const mesh::MeshEntity cell_geo                     = cell_geo_dofs.active_cell(c);
    const mesh::MeshEntity cell_src                     = cell_sol_dofs.active_cell(c);

    const ElemShape eshape_tgt    = cell_src.pt_set_id().elem_shape();
    const PointSetID ref_topo_tgt = cell_src.pt_set_id().ref_topology();
    const mesh::PointSetTag tgt_cell_q_tag(eshape_tgt, q, ref_topo_tgt);
    tgt_cell_q_std_reg.change_type(tgt_cell_q_tag);
    const Uint nb_vert_tgt = tgt_cell_q_std_reg.get().nb_nodes();

    const Uint poly_order_geo = cell_geo.pt_set_id().poly_order();
    const Uint poly_order_p   = cell_src.pt_set_id().poly_order();
    const Uint poly_order_q   = q;

    const Uint quad_order_qq = std::max(poly_order_geo, 2 * poly_order_q);
    const Uint quad_order_pq = std::max(poly_order_geo, poly_order_p + poly_order_q);

    // --------------------------------

    geo_cache_qq.flush();

    const math::DenseConstMatView<Real> geo_cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell_geo.pt_set_id(), tcell_view.coordinates());
    const mesh::PointSetTagExt tmp_key_qq(cell_geo.pt_set_id(), quad_order_qq,
                                          mesh::CellTransform::NO_TRANS, 0);
    geo_cache_qq.push_back_to_buffer(geo_cell_coords, tmp_key_qq);

    geo_metric_qq.empty_buffer();
    geo_metric_qq.evaluate(geo_cache_qq, true);

    // --------------------------------

    geo_cache_pq.flush();

    const mesh::PointSetTagExt tmp_key_pq(cell_geo.pt_set_id(), quad_order_pq,
                                          mesh::CellTransform::NO_TRANS, 0);
    geo_cache_pq.push_back_to_buffer(geo_cell_coords, tmp_key_pq);

    geo_metric_pq.empty_buffer();
    geo_metric_pq.evaluate(geo_cache_pq, true);

    // --------------------------------

    const Uint key_src = cell_src.pt_set_id().store_value();
    const Uint key_tgt = tgt_cell_q_tag.store_value();

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
        elem_types_q.std_region_data(mesh::PointSetTagExt(tgt_cell_q_tag, quad_order_pq));
    const math::DenseDMat<Real> &V_q = (*V_q_ptr).Vandermonde();

    const common::PtrHandle<FEValues> V_qq_ptr =
        elem_types_qq.std_region_data(mesh::PointSetTagExt(tgt_cell_q_tag, quad_order_qq));
    const math::DenseDMat<Real> &V_qq = (*V_qq_ptr).Vandermonde();

    const Uint nb_vert_src = cell_src.nb_vert();

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

    const Uint nb_rows_one_field = M_qq_inv.rows();
    const Uint nb_cols_one_field = M_qp.cols();
    workspace.resize(nb_rows_one_field * nb_cols_one_field);
    math::DenseMatView<Real> M_work(workspace.data(), nb_cols_one_field, nb_rows_one_field,
                                    nb_cols_one_field);
    M_work = M_qq_inv * M_qp;

    math::DenseMatView<Real> M_transfer(mat_values->data() + (*mat_offsets)[cell_reordering[c]],
                                        nb_cols_one_field * NEQ, nb_rows_one_field * NEQ,
                                        nb_cols_one_field * NEQ);

    for (Uint r = 0; r < nb_rows_one_field; ++r)
    {
      for (Uint c = 0; c < nb_cols_one_field; ++c)
      {
        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          M_transfer(NEQ * r + eq, NEQ * c + eq) = M_work(r, c);
        }
      }
    } // Loop over rows of the transfer matrix corresponding to one field

    /*
    std::cout << "*****************************************" << std::endl;

    std::cout << "Shape one field: " << M_work.rows() << " x " <<
    M_work.cols()
    << std::endl; std::cout << "Shape " << NEQ << " fields: " <<
    M_transfer.rows() << " x " << M_transfer.cols()
              << std::endl;

    std::cout << M_work << std::endl;
    std::cout << M_transfer << std::endl;
    */

  } // Loop over cells

  std::unique_ptr<std::vector<Uint>> idx_to_op_map(new std::vector<Uint>());
  idx_to_op_map->resize(nb_cells);
  std::iota(idx_to_op_map->begin(), idx_to_op_map->end(), 0);

  std::unique_ptr<common::BlockArray<Real, Uint>> mat_storage(new common::BlockArray<Real, Uint>());
  mat_storage->build_from_offsets(std::move(mat_values), std::move(mat_offsets));

  local_restriction_ops.allocate(std::move(op_shapes), std::move(mat_storage),
                                 std::move(idx_to_op_map));

  end = clock();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Coarse scale correction operator builder: computing L2 "
               "restriction matrices took : "
            << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void CoarseScaleCorrectionOpBuilder::build_prolongation_ops_L2(
    const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &cell_sol_dofs,
    const std::vector<Uint> &cell_reordering, const Uint p,
    ls::LocalTransferOps<Real> &local_prolongation_ops)
{
  using target_mass_map_type = std::unordered_map<Uint, math::DenseDMat<Real>>;
  using source_mass_map_type =
      std::unordered_map<std::pair<Uint, Uint>, math::DenseDMat<Real>, UnsignedIntPairHasher>;

  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_geo_qq;
  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_qq;

  common::DataMap<mesh::PointSetTagExt, FEValues> elem_types_geo_pq;
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
  mesh::StdRegion src_cell_p_std_reg;

  const Uint nb_cells = cell_sol_dofs.nb_active_cells();

  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> op_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>());
  op_shapes->resize(nb_cells);

  for (Uint c = 0; c < nb_cells; ++c)
  {
    const mesh::MeshEntity cell_geo   = cell_geo_dofs.active_cell(c);
    const mesh::MeshEntity cell_src_q = cell_sol_dofs.active_cell(c);

    const mesh::PointSetTag geo_cell_tag = cell_geo.pt_set_id();
    const ElemShape eshape_geo           = geo_cell_tag.elem_shape();
    const Uint poly_order_geo            = geo_cell_tag.poly_order();

    const mesh::PointSetTag sol_cell_q_tag = cell_src_q.pt_set_id();
    const ElemShape eshape_sol             = sol_cell_q_tag.elem_shape();
    const Uint poly_order_tgt              = sol_cell_q_tag.poly_order();
    const PointSetID ref_topo              = sol_cell_q_tag.ref_topology();

    const mesh::PointSetTag sol_cell_p_tag(eshape_sol, p, ref_topo);
    src_cell_p_std_reg.change_type(sol_cell_p_tag);

    const Uint src_cell_nb_vert = src_cell_p_std_reg.get().nb_nodes();
    const Uint tgt_cell_nb_vert = cell_src_q.nb_vert();

    // ----------------------------------------------------
    // I) Process the 'target' cell data - shape function
    // values and mass matrix with quadrature order
    // 2*q (aka 'qq')
    // ----------------------------------------------------

    const Uint quad_order_qq = std::max(2 * poly_order_tgt, poly_order_geo);
    quad.change_type(eshape_geo, quad_order_qq, PointSetID::Gauss);

    common::PtrHandle<FEValues> fe_values_geo_qq_ptr =
        elem_types_geo_qq.std_region_data(mesh::PointSetTagExt(geo_cell_tag, quad_order_qq));

    // Geometry expansion for target cell type
    if (fe_values_geo_qq_ptr.is_null())
    {
      fe_values_geo_qq_ptr =
          elem_types_geo_qq.create(mesh::PointSetTagExt(geo_cell_tag, quad_order_qq));
      (*fe_values_geo_qq_ptr)
          .configure(geo_cell_tag, mesh::sf::SFTag(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                                   ModalBasis::Modal));
      (*fe_values_geo_qq_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // Solution expansion for target cell type
    common::PtrHandle<FEValues> fe_values_qq_ptr =
        elem_types_qq.std_region_data(mesh::PointSetTagExt(sol_cell_q_tag, quad_order_qq));

    if (fe_values_qq_ptr.is_null())
    {
      fe_values_qq_ptr = elem_types_qq.create(mesh::PointSetTagExt(sol_cell_q_tag, quad_order_qq));
      (*fe_values_qq_ptr)
          .configure(sol_cell_q_tag, mesh::sf::SFTag(eshape_sol, SFunc::Lagrange, poly_order_tgt,
                                                     ModalBasis::Modal));
      (*fe_values_qq_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    const Uint key = sol_cell_q_tag.store_value();

    target_mass_map_type::iterator target_mass_map_it = target_mass_mat_map.find(key);
    if (target_mass_map_it == target_mass_mat_map.end())
    {
      target_mass_mat_map.insert(std::make_pair(key, M));
      target_mass_map_it = target_mass_mat_map.find(key);

      target_mass_mat_inv_map.insert(std::make_pair(key, M));
      target_mass_map_type::iterator map_inv_it = target_mass_mat_inv_map.find(key);

      target_mass_map_it->second.resize(tgt_cell_nb_vert, tgt_cell_nb_vert);
      map_inv_it->second.resize(tgt_cell_nb_vert, tgt_cell_nb_vert);
    }

    // ----------------------------------------------------
    // II) Process the 'source' cell data - shape function
    // values and mass matrix with quadrature order
    // p+q (aka 'pq')
    // ----------------------------------------------------

    // The quadrature has to have such order that phi_q * phi_p can be
    // integrated exactly
    const Uint quad_order_pq = std::max(poly_order_tgt + p, poly_order_geo);

    // Geometry expansion for source cell type
    common::PtrHandle<FEValues> fe_values_geo_pq_ptr =
        elem_types_geo_pq.std_region_data(mesh::PointSetTagExt(geo_cell_tag, quad_order_pq));

    if (fe_values_geo_pq_ptr.is_null())
    {
      quad.change_type(eshape_geo, quad_order_pq, PointSetID::Gauss);

      fe_values_geo_pq_ptr =
          elem_types_geo_pq.create(mesh::PointSetTagExt(geo_cell_tag, quad_order_pq));
      (*fe_values_geo_pq_ptr)
          .configure(geo_cell_tag, mesh::sf::SFTag(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                                   ModalBasis::Modal));
      (*fe_values_geo_pq_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // a) Make sure there is Vandermonde matrix for the 'source/p-type'
    // element
    //    Note that the second entry in the key is the COMBINED QUADRATURE
    //    ORDER!
    const mesh::PointSetTagExt key_src(sol_cell_p_tag, quad_order_pq);

    common::PtrHandle<FEValues> fe_values_p_ptr = elem_types_p.std_region_data(key_src);
    if (fe_values_p_ptr.is_null())
    {
      quad.change_type(eshape_sol, quad_order_pq, PointSetID::Gauss);

      fe_values_p_ptr = elem_types_p.create(key_src);
      (*fe_values_p_ptr)
          .configure(sol_cell_p_tag,
                     mesh::sf::SFTag(eshape_sol, SFunc::Lagrange, p, ModalBasis::Modal));
      (*fe_values_p_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // b) Make sure there is Vandermonde matrix for the 'target/q-type'
    // element
    //    Note that the second entry in the key is the COMBINED QUADRATURE
    //    ORDER!
    const mesh::PointSetTagExt key_tgt(sol_cell_q_tag, quad_order_pq);

    common::PtrHandle<FEValues> fe_values_ptr_q = elem_types_q.std_region_data(key_tgt);
    if (fe_values_ptr_q.is_null())
    {
      quad.change_type(eshape_sol, quad_order_pq, PointSetID::Gauss);

      fe_values_ptr_q = elem_types_q.create(key_tgt);
      (*fe_values_ptr_q)
          .configure(sol_cell_q_tag, mesh::sf::SFTag(eshape_sol, SFunc::Lagrange, poly_order_tgt,
                                                     ModalBasis::Modal));
      (*fe_values_ptr_q).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
    }

    // c) Count how many entries there will be to hold transfer operators
    // for all cells
    (*op_shapes)[cell_reordering[c]] =
        common::ArrayShape<_2D, SUint>(NEQ * tgt_cell_nb_vert, NEQ * src_cell_nb_vert);

    // d) Prepare storage for the rectangular 'mass' matrices M_pq
    const std::pair<Uint, Uint> key_combined =
        std::make_pair(sol_cell_q_tag.store_value(), sol_cell_p_tag.store_value());

    source_mass_map_type::iterator map_it = src_mass_mat_map.find(key_combined);

    if (map_it == src_mass_mat_map.end())
    {
      math::DenseDMat<Real> M;
      src_mass_mat_map.insert(std::make_pair(key_combined, M));
      map_it = src_mass_mat_map.find(key_combined);
      map_it->second.resize(tgt_cell_nb_vert, src_cell_nb_vert);
    }

  } // Loop over cells

  // --------------

  std::unique_ptr<std::vector<Real>> mat_values(new std::vector<Real>());
  std::unique_ptr<std::vector<Uint>> mat_offsets(new std::vector<Uint>());

  mat_offsets->reserve(nb_cells + 1);
  mat_offsets->push_back(0);

  Uint tot_mat_storage = 0;
  for (const auto shape : (*op_shapes))
  {
    const Uint mat_block_size = shape.size(0) * shape.size(1);
    mat_offsets->push_back(mat_block_size);
    tot_mat_storage += mat_block_size;
  }

  for (Uint i = 1; i < (*mat_offsets).size(); ++i)
  {
    (*mat_offsets)[i] += (*mat_offsets)[i - 1];
  }

  mat_values->resize(tot_mat_storage);

  // --------------

  GeometryCache<MeshConfig::GDIM> geo_cache_qq;
  GeometryCache<MeshConfig::GDIM> geo_cache_pq;
  GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric_qq;
  GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric_pq;

  geo_cache_qq.allocate(elem_types_geo_qq, 1u);
  geo_metric_qq.allocate_buffer(elem_types_geo_qq, 1u);

  geo_cache_pq.allocate(elem_types_geo_pq, 1u);
  geo_metric_pq.allocate_buffer(elem_types_geo_pq, 1u);

  // --------------

  std::vector<Real> workspace;
  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < nb_cells; ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_sol_dofs.tcell(c);
    const mesh::MeshEntity cell_geo                     = cell_geo_dofs.active_cell(c);
    const mesh::MeshEntity cell_tgt                     = cell_sol_dofs.active_cell(c);

    const mesh::PointSetTag tgt_cell_q_tag = cell_tgt.pt_set_id();
    const ElemShape eshape                 = tgt_cell_q_tag.elem_shape();
    const PointSetID ref_topo_tgt          = tgt_cell_q_tag.ref_topology();
    const mesh::PointSetTag src_cell_p_tag(eshape, p, ref_topo_tgt);
    src_cell_p_std_reg.change_type(src_cell_p_tag);

    const Uint nb_vert_src = src_cell_p_std_reg.get().nb_nodes();
    const Uint nb_vert_tgt = cell_tgt.nb_vert();

    const Uint poly_order_geo = cell_geo.pt_set_id().poly_order();
    const Uint poly_order_p   = p;
    const Uint poly_order_q   = cell_tgt.pt_set_id().poly_order();

    const Uint quad_order_qq = std::max(poly_order_geo, 2 * poly_order_q);
    const Uint quad_order_pq = std::max(poly_order_geo, poly_order_p + poly_order_q);

    // --------------------------------

    geo_cache_qq.flush();

    const math::DenseConstMatView<Real> geo_cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell_geo.pt_set_id(), tcell_view.coordinates());

    const mesh::PointSetTagExt tmp_key_qq(cell_geo.pt_set_id(), quad_order_qq,
                                          mesh::CellTransform::NO_TRANS, 0);
    geo_cache_qq.push_back_to_buffer(geo_cell_coords, tmp_key_qq);

    geo_metric_qq.empty_buffer();
    geo_metric_qq.evaluate(geo_cache_qq, true);

    // --------------------------------

    geo_cache_pq.flush();

    const mesh::PointSetTagExt tmp_key_pq(cell_geo.pt_set_id(), quad_order_pq,
                                          mesh::CellTransform::NO_TRANS, 0);
    geo_cache_pq.push_back_to_buffer(geo_cell_coords, tmp_key_pq);

    geo_metric_pq.empty_buffer();
    geo_metric_pq.evaluate(geo_cache_pq, true);

    // --------------------------------

    const Uint key_src = src_cell_p_tag.store_value();
    const Uint key_tgt = tgt_cell_q_tag.store_value();

    const std::pair<Uint, Uint> key_combined = std::make_pair(key_tgt, key_src);

    source_mass_map_type::iterator src_map_it = src_mass_mat_map.find(key_combined);

    target_mass_map_type::iterator tgt_map_it     = target_mass_mat_map.find(key_tgt);
    target_mass_map_type::iterator tgt_map_inv_it = target_mass_mat_inv_map.find(key_tgt);

    math::DenseDMat<Real> &M_qp = src_map_it->second;

    math::DenseDMat<Real> &M_qq     = tgt_map_it->second;
    math::DenseDMat<Real> &M_qq_inv = tgt_map_inv_it->second;

    const common::PtrHandle<FEValues> V_p_ptr =
        elem_types_p.std_region_data(mesh::PointSetTagExt(src_cell_p_tag, quad_order_pq));
    const math::DenseDMat<Real> &V_p = (*V_p_ptr).Vandermonde();

    const common::PtrHandle<FEValues> V_q_ptr =
        elem_types_q.std_region_data(mesh::PointSetTagExt(tgt_cell_q_tag, quad_order_pq));
    const math::DenseDMat<Real> &V_q = (*V_q_ptr).Vandermonde();

    const common::PtrHandle<FEValues> V_qq_ptr =
        elem_types_qq.std_region_data(mesh::PointSetTagExt(tgt_cell_q_tag, quad_order_qq));
    const math::DenseDMat<Real> &V_qq = (*V_qq_ptr).Vandermonde();

    const typename GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>::cellwise_metric
        cell_metric_pq = geo_metric_pq.cellwise_values(0);

    // Build the rectangular matrix M_qp
    M_qp.fill(0.0);

    const Uint n_quad_pq                         = cell_metric_pq.nb_qd_pts();
    const math::DenseConstVecView<Real> j_det_pq = cell_metric_pq.jdet();

    const math::DenseDVec<Real> &w_pq = cell_metric_pq.pt_weights();

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

    // Build the matrix M_qq (the mass matrix of the element in 'target'
    // space

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

    const Uint nb_rows_one_field = M_qq_inv.rows();
    const Uint nb_cols_one_field = M_qp.cols();
    workspace.resize(nb_rows_one_field * nb_cols_one_field);
    math::DenseMatView<Real> M_work(workspace.data(), nb_cols_one_field, nb_rows_one_field,
                                    nb_cols_one_field);
    M_work = M_qq_inv * M_qp;

    math::DenseMatView<Real> M_transfer(mat_values->data() + (*mat_offsets)[cell_reordering[c]],
                                        nb_cols_one_field * NEQ, nb_rows_one_field * NEQ,
                                        nb_cols_one_field * NEQ);

    for (Uint r = 0; r < nb_rows_one_field; ++r)
    {
      for (Uint c = 0; c < nb_cols_one_field; ++c)
      {
        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          M_transfer(NEQ * r + eq, NEQ * c + eq) = M_work(r, c);
        }
      }
    } // Loop over rows of the transfer matrix corresponding to one field
  }   // Loop over cells

  std::unique_ptr<std::vector<Uint>> idx_to_op_map(new std::vector<Uint>());
  idx_to_op_map->resize(nb_cells);
  std::iota(idx_to_op_map->begin(), idx_to_op_map->end(), 0);

  std::unique_ptr<common::BlockArray<Real, Uint>> mat_storage(new common::BlockArray<Real, Uint>());
  mat_storage->build_from_offsets(std::move(mat_values), std::move(mat_offsets));

  local_prolongation_ops.allocate(std::move(op_shapes), std::move(mat_storage),
                                  std::move(idx_to_op_map));

  end = clock();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Coarse scale correction operator builder: computing L2 "
               "prolongation matrices took : "
            << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void CoarseScaleCorrectionOpBuilder::build_restriction_ops_modal(
    const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &cell_sol_dofs,
    const std::vector<Uint> &cell_reordering, const Uint q,
    ls::LocalTransferOps<Real> &local_restriction_ops)
{
  const Uint nb_cells = cell_sol_dofs.nb_active_cells();

  std::unique_ptr<std::vector<Uint>> idx_to_op_map(new std::vector<Uint>());
  idx_to_op_map->resize(nb_cells);

  std::vector<mesh::PointSetTag> source_cell_tags;

  clock_t start, end;
  Real elapsed;

  start = clock();

  for (Uint ac = 0; ac < cell_sol_dofs.nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity sol_cell      = cell_sol_dofs.active_cell(mesh::ActiveIdx(ac));
    const mesh::PointSetTag sol_cell_tag = sol_cell.pt_set_id();
    Uint search_idx                      = source_cell_tags.size();
    for (Uint j = 0; j < source_cell_tags.size(); ++j)
    {
      if (source_cell_tags[j] == sol_cell_tag)
      {
        search_idx = j;
      }
    }

    if (search_idx == source_cell_tags.size())
    {
      source_cell_tags.push_back(sol_cell_tag);
    }
    (*idx_to_op_map)[cell_reordering[ac]] = search_idx;
  }

  mesh::StdRegion region_fine;
  mesh::StdRegion region_coarse;

  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> op_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>());
  op_shapes->resize(source_cell_tags.size());

  for (Uint i = 0; i < source_cell_tags.size(); ++i)
  {
    const mesh::PointSetTag fine_tag = source_cell_tags[i];
    const mesh::PointSetTag coarse_tag(fine_tag.elem_shape(), q, fine_tag.ref_topology());

    region_fine.change_type(fine_tag);
    region_coarse.change_type(coarse_tag);

    (*op_shapes)[i] = common::ArrayShape<_2D, SUint>(region_coarse.get().nb_nodes() * NEQ,
                                                     region_fine.get().nb_nodes() * NEQ);
  }

  std::unique_ptr<std::vector<Real>> mat_values(new std::vector<Real>());
  std::unique_ptr<std::vector<Uint>> mat_offsets(new std::vector<Uint>());

  mat_offsets->reserve(source_cell_tags.size() + 1);
  mat_offsets->push_back(0);

  Uint tot_mat_storage = 0;
  for (const auto shape : (*op_shapes))
  {
    const Uint mat_block_size = shape.size(0) * shape.size(1);
    mat_offsets->push_back(mat_block_size);
    tot_mat_storage += mat_block_size;
  }

  for (Uint i = 1; i < (*mat_offsets).size(); ++i)
  {
    (*mat_offsets)[i] += (*mat_offsets)[i - 1];
  }

  mat_values->resize(tot_mat_storage);
  std::fill(mat_values->begin(), mat_values->end(), 0);

  mesh::sf::ShapeFunction sf_Lagrange;
  mesh::sf::ShapeFunction sf_modal;

  // Transpose of Lagrange Vandermonde matrix
  math::DenseDMat<Real> V_Lagrange_T;
  // Transpose of modal Vandermonde matrix
  math::DenseDMat<Real> V_modal_T;

  // Matrix representing change of basis modal->Lagrange
  math::DenseDMat<Real> C_modal_to_Lagrange;

  // Prolongation matrix
  math::DenseDMat<Real> prolongation_mat;

  // Restriction matrix (transpose of prolongation matrix)
  math::DenseDMat<Real> restriction_mat;

  math::DenseDMat<Real> V_tmp;

  for (Uint i = 0; i < source_cell_tags.size(); ++i)
  {
    const mesh::PointSetTag fine_tag = source_cell_tags[i];
    const mesh::PointSetTag coarse_tag(fine_tag.elem_shape(), q, fine_tag.ref_topology());

    region_fine.change_type(fine_tag);
    region_coarse.change_type(coarse_tag);

    const mesh::sf::SFTag sf_tag_Lagrange(fine_tag.elem_shape(), SFunc::Lagrange, q,
                                          ModalBasis::Modal);
    const mesh::sf::SFTag sf_tag_modal(fine_tag.elem_shape(), SFunc::Modal, q, ModalBasis::Modal);

    sf_Lagrange.change_type(fine_tag, sf_tag_Lagrange);
    sf_modal.change_type(fine_tag, sf_tag_modal);

    // Compute Vandermonde matrix of coarse-space Lagrange SF in coarse
    // space DOFs
    sf_Lagrange.get().compute_ref_values(region_coarse.get().coordinates(), V_tmp);

    // Compute the transpose
    V_Lagrange_T.resize(V_tmp.cols(), V_tmp.rows());
    for (Uint r = 0; r < V_tmp.rows(); ++r)
    {
      for (Uint c = 0; c < V_tmp.cols(); ++c)
      {
        V_Lagrange_T(c, r) = V_tmp(r, c);
      }
    }

    // Compute Vandermonde matrix of coarse-space modal basis in coarse
    // space DOFs
    sf_modal.get().compute_ref_values(region_coarse.get().coordinates(), V_tmp);

    // Compute the transpose
    V_modal_T.resize(V_tmp.cols(), V_tmp.rows());
    for (Uint r = 0; r < V_tmp.rows(); ++r)
    {
      for (Uint c = 0; c < V_tmp.cols(); ++c)
      {
        V_modal_T(c, r) = V_tmp(r, c);
      }
    }

    // Evaluate the matrix of coefficients of basis change from Lagrange to
    // modal basis V_Lagrange_T = C_modal_to_Lagrange * V_modal_T
    // C_Lagrange_to_modal = V_modal_T * inv(V_Lagrange_T)
    V_modal_T.inv(V_tmp);
    C_modal_to_Lagrange = V_Lagrange_T * V_tmp;

    // Now transpose the values in C_modal_to_Lagrange
    V_tmp = C_modal_to_Lagrange;

    for (Uint r = 0; r < V_tmp.rows(); ++r)
    {
      for (Uint c = 0; c < V_tmp.cols(); ++c)
      {
        C_modal_to_Lagrange(c, r) = V_tmp(r, c);
      }
    }

    // Re-compute Vandermonde matrix of modal basis: take coarse-space modes
    // and evaluate them in fine-space DOFs
    sf_modal.get().compute_ref_values(region_fine.get().coordinates(), V_tmp);

    prolongation_mat.resize(V_tmp.rows(), C_modal_to_Lagrange.cols());
    prolongation_mat = V_tmp * C_modal_to_Lagrange;

    // Compute restriction matrix as transpose of prolongation matrix
    restriction_mat.resize(prolongation_mat.cols(), prolongation_mat.rows());
    for (Uint r = 0; r < prolongation_mat.rows(); ++r)
    {
      for (Uint c = 0; c < prolongation_mat.cols(); ++c)
      {
        restriction_mat(c, r) = prolongation_mat(r, c);
      }
    }

    /*
    // Normalize rows of restriction matrix
    for (Uint r = 0; r < restriction_mat.rows(); ++r)
    {
      Real col_sum = 0.0;
      for (Uint c = 0; c < restriction_mat.cols(); ++c)
      {
        col_sum += restriction_mat(r, c);
      }
      if (std::abs(col_sum) > 1.e-14)
      {
        for (Uint c = 0; c < restriction_mat.cols(); ++c)
        {
          restriction_mat(r, c) *= 1. / col_sum;
        }
      }
    }
    */

    const Uint nb_rows_one_field = restriction_mat.rows();
    const Uint nb_cols_one_field = restriction_mat.cols();

    math::DenseMatView<Real> M_transfer(mat_values->data() + (*mat_offsets)[i],
                                        nb_cols_one_field * NEQ, nb_rows_one_field * NEQ,
                                        nb_cols_one_field * NEQ);

    for (Uint r = 0; r < nb_rows_one_field; ++r)
    {
      for (Uint c = 0; c < nb_cols_one_field; ++c)
      {
        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          M_transfer(NEQ * r + eq, NEQ * c + eq) = restriction_mat(r, c);
        }
      }
    } // Loop over rows of the transfer matrix corresponding to one field

  } // Loop over source cell tags

  std::unique_ptr<common::BlockArray<Real, Uint>> mat_storage(new common::BlockArray<Real, Uint>());
  mat_storage->build_from_offsets(std::move(mat_values), std::move(mat_offsets));

  local_restriction_ops.allocate(std::move(op_shapes), std::move(mat_storage),
                                 std::move(idx_to_op_map));

  end = clock();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Coarse scale correction operator builder: computing modal "
               "restriction matrices took : "
            << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void CoarseScaleCorrectionOpBuilder::build_prolongation_ops_modal(
    const Uint NEQ, const typename result_of::dof_map_t<MeshConfig> &cell_geo_dofs,
    const typename result_of::dof_map_t<MeshConfig> &cell_sol_dofs,
    const std::vector<Uint> &cell_reordering, const Uint q,
    ls::LocalTransferOps<Real> &local_restriction_ops)
{
  const Uint nb_cells = cell_sol_dofs.nb_active_cells();

  std::unique_ptr<std::vector<Uint>> idx_to_op_map(new std::vector<Uint>());
  idx_to_op_map->resize(nb_cells);

  std::vector<mesh::PointSetTag> source_cell_tags;

  clock_t start, end;
  Real elapsed;

  start = clock();

  for (Uint ac = 0; ac < cell_sol_dofs.nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity sol_cell      = cell_sol_dofs.active_cell(mesh::ActiveIdx(ac));
    const mesh::PointSetTag sol_cell_tag = sol_cell.pt_set_id();
    Uint search_idx                      = source_cell_tags.size();
    for (Uint j = 0; j < source_cell_tags.size(); ++j)
    {
      if (source_cell_tags[j] == sol_cell_tag)
      {
        search_idx = j;
      }
    }

    if (search_idx == source_cell_tags.size())
    {
      source_cell_tags.push_back(sol_cell_tag);
    }
    (*idx_to_op_map)[cell_reordering[ac]] = search_idx;
  }

  mesh::StdRegion region_fine;
  mesh::StdRegion region_coarse;

  std::unique_ptr<std::vector<common::ArrayShape<_2D, SUint>>> op_shapes(
      new std::vector<common::ArrayShape<_2D, SUint>>());
  op_shapes->resize(source_cell_tags.size());

  for (Uint i = 0; i < source_cell_tags.size(); ++i)
  {
    const mesh::PointSetTag fine_tag = source_cell_tags[i];
    const mesh::PointSetTag coarse_tag(fine_tag.elem_shape(), q, fine_tag.ref_topology());

    region_fine.change_type(fine_tag);
    region_coarse.change_type(coarse_tag);

    (*op_shapes)[i] = common::ArrayShape<_2D, SUint>(region_fine.get().nb_nodes() * NEQ,
                                                     region_coarse.get().nb_nodes() * NEQ);
  }

  std::unique_ptr<std::vector<Real>> mat_values(new std::vector<Real>());
  std::unique_ptr<std::vector<Uint>> mat_offsets(new std::vector<Uint>());

  mat_offsets->reserve(source_cell_tags.size() + 1);
  mat_offsets->push_back(0);

  Uint tot_mat_storage = 0;
  for (const auto shape : (*op_shapes))
  {
    const Uint mat_block_size = shape.size(0) * shape.size(1);
    mat_offsets->push_back(mat_block_size);
    tot_mat_storage += mat_block_size;
  }

  for (Uint i = 1; i < (*mat_offsets).size(); ++i)
  {
    (*mat_offsets)[i] += (*mat_offsets)[i - 1];
  }

  mat_values->resize(tot_mat_storage);
  std::fill(mat_values->begin(), mat_values->end(), 0);

  mesh::sf::ShapeFunction sf_Lagrange;
  mesh::sf::ShapeFunction sf_modal;

  // Transpose of Lagrange Vandermonde matrix
  math::DenseDMat<Real> V_Lagrange_T;
  // Transpose of modal Vandermonde matrix
  math::DenseDMat<Real> V_modal_T;

  // Matrix representing change of basis modal->Lagrange
  math::DenseDMat<Real> C_modal_to_Lagrange;

  // Prolongation matrix
  math::DenseDMat<Real> prolongation_mat;

  math::DenseDMat<Real> V_tmp;

  for (Uint i = 0; i < source_cell_tags.size(); ++i)
  {
    const mesh::PointSetTag fine_tag = source_cell_tags[i];
    const mesh::PointSetTag coarse_tag(fine_tag.elem_shape(), q, fine_tag.ref_topology());

    region_fine.change_type(fine_tag);
    region_coarse.change_type(coarse_tag);

    const mesh::sf::SFTag sf_tag_Lagrange(fine_tag.elem_shape(), SFunc::Lagrange, q,
                                          ModalBasis::Modal);
    const mesh::sf::SFTag sf_tag_modal(fine_tag.elem_shape(), SFunc::Modal, q, ModalBasis::Modal);

    sf_Lagrange.change_type(fine_tag, sf_tag_Lagrange);
    sf_modal.change_type(fine_tag, sf_tag_modal);

    // Compute Vandermonde matrix of coarse-space Lagrange SF in coarse
    // space DOFs
    sf_Lagrange.get().compute_ref_values(region_coarse.get().coordinates(), V_tmp);

    // Compute the transpose
    V_Lagrange_T.resize(V_tmp.cols(), V_tmp.rows());
    for (Uint r = 0; r < V_tmp.rows(); ++r)
    {
      for (Uint c = 0; c < V_tmp.cols(); ++c)
      {
        V_Lagrange_T(c, r) = V_tmp(r, c);
      }
    }

    // Compute Vandermonde matrix of coarse-space modal basis in coarse
    // space DOFs
    sf_modal.get().compute_ref_values(region_coarse.get().coordinates(), V_tmp);

    // Compute the transpose
    V_modal_T.resize(V_tmp.cols(), V_tmp.rows());
    for (Uint r = 0; r < V_tmp.rows(); ++r)
    {
      for (Uint c = 0; c < V_tmp.cols(); ++c)
      {
        V_modal_T(c, r) = V_tmp(r, c);
      }
    }

    // Evaluate the matrix of coefficients of basis change from Lagrange to
    // modal basis V_Lagrange_T = C_modal_to_Lagrange * V_modal_T
    // C_Lagrange_to_modal = V_modal_T * inv(V_Lagrange_T)
    V_modal_T.inv(V_tmp);
    C_modal_to_Lagrange = V_Lagrange_T * V_tmp;

    // Now transpose the values in C_modal_to_Lagrange
    V_tmp = C_modal_to_Lagrange;

    for (Uint r = 0; r < V_tmp.rows(); ++r)
    {
      for (Uint c = 0; c < V_tmp.cols(); ++c)
      {
        C_modal_to_Lagrange(c, r) = V_tmp(r, c);
      }
    }

    // Re-compute Vandermonde matrix of modal basis: take coarse-space modes
    // and evaluate them in fine-space DOFs
    sf_modal.get().compute_ref_values(region_fine.get().coordinates(), V_tmp);

    prolongation_mat.resize(V_tmp.rows(), C_modal_to_Lagrange.cols());
    prolongation_mat = V_tmp * C_modal_to_Lagrange;

    const Uint nb_rows_one_field = prolongation_mat.rows();
    const Uint nb_cols_one_field = prolongation_mat.cols();

    math::DenseMatView<Real> M_transfer(mat_values->data() + (*mat_offsets)[i],
                                        nb_cols_one_field * NEQ, nb_rows_one_field * NEQ,
                                        nb_cols_one_field * NEQ);

    for (Uint r = 0; r < nb_rows_one_field; ++r)
    {
      for (Uint c = 0; c < nb_cols_one_field; ++c)
      {
        for (Uint eq = 0; eq < NEQ; ++eq)
        {
          M_transfer(NEQ * r + eq, NEQ * c + eq) = prolongation_mat(r, c);
        }
      }
    } // Loop over rows of the transfer matrix corresponding to one field

  } // Loop over source cell tags

  std::unique_ptr<common::BlockArray<Real, Uint>> mat_storage(new common::BlockArray<Real, Uint>());
  mat_storage->build_from_offsets(std::move(mat_values), std::move(mat_offsets));

  local_restriction_ops.allocate(std::move(op_shapes), std::move(mat_storage),
                                 std::move(idx_to_op_map));

  end = clock();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Coarse scale correction operator builder: computing modal "
               "prolongation matrices took : "
            << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
