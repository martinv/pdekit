#ifndef PDEKIT_Solver_FE_Interior_Solver_CG_HDG_hpp
#define PDEKIT_Solver_FE_Interior_Solver_CG_HDG_hpp

/// ***********************************************
/// TODO
/// ***********************************************
/// - Where possible, don't access system matrix
///   and solution by pointer, but through a reference
///   instead
/// - Fix the process of generating a triangulation from
///   cell buffer: why do material ids need to be set
///   separately

/// Standard template library headers
#include <array>
#include <cmath>
#include <ctime>
#include <forward_list>
#include <iostream>
#include <map>

/// PDEKIT headers
#include "common/Constants.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/MathConstants.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class InteriorSolverCGHDG
{
  public:
  /// TYPEDEFS
  using tria_t       = typename mesh::Tria<MeshConfig>;
  using dof_map_t    = typename pdekit::result_of::dof_map_t<MeshConfig>;
  using multi_edge_t = mesh::internal::TriaFacets<MeshConfig>;

  /// Default constructor
  InteriorSolverCGHDG();

  /// Destructor
  ~InteriorSolverCGHDG() = default;

  void setup(const tria_t &tria, const dof_map_t &global_dofs,
             const std::vector<mesh::ActiveIdx> &cell_ids);

  template <typename BdryDofIterator>
  void add_boundary(const common::IteratorRange<BdryDofIterator> &global_boundary,
                    const std::string &name);

  void add_boundary(const tria_t &tria_global, const multi_edge_t &global_boundary,
                    const LeftRightOrientation orient, const std::string &name);

  template <typename RHS>
  void assemble(const RHS &rhs);

  template <typename DirichletBC>
  void add_weak_bc(const std::string &boundary_name, const DirichletBC &dirichlet_bc);

  void solve();

  const tria_t &local_mesh() const;

  const dof_map_t &local_dofs() const;

  const interpolation::VectorMeshFunction<Real> &solution() const;

  void write_to_gmsh(const std::string &filename) const;

  private:
  struct ActiveIdxHash
  {
    std::size_t operator()(mesh::ActiveIdx const &idx) const noexcept
    {
      std::size_t h1 = std::hash<mesh::ActiveIdx::value_type>{}(idx.id());
      return h1;
    }
  };

  /*
  struct ActiveIdxCompare
  {
    bool operator()(mesh::ActiveIdx const &idx_l, mesh::ActiveIdx const &idx_r) const noexcept
    {
      return idx_l.id() < idx_r.id();
    }
  };
  */

  using cell_buffer_t = mesh::CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM>;
  // using active_cell_idx_set_t = std::set<mesh::ActiveIdx, ActiveIdxCompare>;
  using glob_to_loc_map_t = std::unordered_map<mesh::ActiveIdx, mesh::ActiveIdx, ActiveIdxHash>;
  using cell_geo_metric_type =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             MeshConfig::GDIM>::cellwise_metric;
  using facet_geo_metric_type =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             MeshConfig::GDIM - 1>::cellwise_metric;

  static std::unique_ptr<cell_buffer_t> create_sub_mesh_buffer(
      const tria_t &tria, const dof_map_t &dofs, const std::vector<mesh::ActiveIdx> &cell_ids,
      const Uint material_id);

  template <typename RHS>
  void setup_solver(const RHS &rhs);

  template <typename RHS>
  void assemble_interior(const RHS &rhs);

  template <typename BdryDofIterator>
  void prepare_boundary_metric_data(const common::IteratorRange<BdryDofIterator> &boundary);

  static void build_adj_volume_matrices(const cell_geo_metric_type &cell_geo_metric,
                                        math::DenseDMat<Real> &M_mat,
                                        std::array<math::DenseDMat<Real>, MeshConfig::GDIM> &D_mat);

  static void build_E_matrices(const cell_geo_metric_type &cell_geo_metric,
                               const facet_geo_metric_type &facet_geo_metric,
                               math::DenseDMat<Real> &E_mat,
                               std::array<math::DenseDMat<Real>, MeshConfig::GDIM> &E_tilde_mat);

  static void build_F_matrices(const cell_geo_metric_type &cell_geo_metric,
                               const facet_geo_metric_type &facet_geo_metric,
                               math::DenseDMat<Real> &F_mat,
                               std::array<math::DenseDMat<Real>, MeshConfig::GDIM> &F_tilde_mat);

  void solve_ls(std::shared_ptr<ls::TpetraCrsMatrix<Real>> &sys_A,
                std::shared_ptr<ls::TpetraMultiVector<Real>> &sys_RHS,
                interpolation::VectorMeshFunction<Real> &solution);

  // Triangulation corresponding to one patch in global mesh
  tria_t m_tria;

  glob_to_loc_map_t m_glob_to_loc_map;

  common::BlockArray<Uint, Uint> m_bdry_face_local_id;

  interpolation::GeometryCache<MeshConfig::GDIM> m_geo_cache_elems;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, MeshConfig::GDIM>
      m_geo_metric_elems;

  /// Cache and metric for shape function gradients on Dirichlet boundary
  interpolation::GeometryCache<MeshConfig::GDIM> m_geo_cache_elems_on_bdry;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, MeshConfig::GDIM>
      m_geo_metric_elems_on_bdry;

  /// Cache and metric for interpolation on boundary traces
  interpolation::GeometryCache<MeshConfig::GDIM> m_geo_cache_faces_on_bdry;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, MeshConfig::GDIM - 1>
      m_geo_metric_faces_on_bdry;

  Uint m_quadrature_order;

  std::shared_ptr<ls::TpetraCrsMatrix<Real>> m_sys_A;
  std::shared_ptr<ls::TpetraMultiVector<Real>> m_sys_RHS;

  interpolation::VectorMeshFunction<Real> m_solution;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
InteriorSolverCGHDG<MeshConfig>::InteriorSolverCGHDG() : m_tria("tria"), m_solution("", "solution")
{
  m_quadrature_order = 0;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void InteriorSolverCGHDG<MeshConfig>::setup(const tria_t &tria, const dof_map_t &global_dofs,
                                            const std::vector<mesh::ActiveIdx> &cell_ids)
{
  const Uint material_id = 1;

  using cell_buffer_t = mesh::CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM>;
  std::unique_ptr<cell_buffer_t> buffer =
      create_sub_mesh_buffer(tria, global_dofs, cell_ids, material_id);

  buffer->remove_dof_numbering_gaps();

  m_tria.create_from_cells(*buffer);

  common::PtrHandle<dof_map_t> local_dofs = m_tria.create_dof_storage("dofs");
  (*local_dofs).create_from_cells(*buffer);

  const std::vector<std::pair<Uint, std::string>> cell_tag_names = {
      std::make_pair(material_id, "interior")};
  std::vector<Uint> tag_values((*local_dofs).nb_nodes());
  tag_values.assign(tag_values.size(), material_id);

  (*local_dofs).tag_all_active_cells(cell_tag_names, tag_values);

  Int i = 0;

  for (auto id : cell_ids)
  {
    m_glob_to_loc_map[id] = mesh::ActiveIdx(i++);
  }

  /*
  std::cout << "Interior solver HDG: interior to global cell id mapping:" << std::endl;
  for (auto id_tuple : m_glob_to_loc_map)
  {
    std::cout << "{" << std::get<0>(id_tuple) << "," << std::get<1>(id_tuple) << "}" << std::endl;
    const mesh::CellTopologyView<MeshConfig> tcell_loc  = m_tria.active_cell(std::get<1>(id_tuple));
    const mesh::CellTopologyView<MeshConfig> tcell_glob = tria.active_cell(std::get<0>(id_tuple));

    const auto loc_coord  = tcell_loc.coordinates();
    const auto glob_coord = tcell_loc.coordinates();

    std::cout << loc_coord << std::endl;
    std::cout << glob_coord << std::endl;

    for (Uint c = 0; c < loc_coord.size(); ++c)
    {
      const auto node_loc  = loc_coord.const_node_view(c);
      const auto node_glob = glob_coord.const_node_view(c);

      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        if (std::abs(node_loc[d] - node_glob[d]) > 1.e-8)
        {
          std::cerr << "PROBLEM" << std::endl;
        }
      }
    }

    std::cout << "========================================================" << std::endl;
  }
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename BdryDofIterator>
void InteriorSolverCGHDG<MeshConfig>::add_boundary(
    const common::IteratorRange<BdryDofIterator> &global_boundary, const std::string &name)
{
  const BdryDofIterator bdry_begin = global_boundary.begin();
  const BdryDofIterator bdry_end   = global_boundary.end();

  std::unique_ptr<std::vector<mesh::IncidenceEntry>> bdry_incidences(
      new std::vector<mesh::IncidenceEntry>());

  for (BdryDofIterator it = bdry_begin; it != bdry_end; ++it)
  {
    const mesh::CellTopologyView<MeshConfig> glob_tcell_view = it->tcell();
    const mesh::ActiveIdx local_active_idx = m_glob_to_loc_map[glob_tcell_view.active_idx()];
    const mesh::CellTopologyView<MeshConfig> loc_tcell_view = m_tria.active_cell(local_active_idx);
    bdry_incidences->push_back(
        mesh::IncidenceEntry(loc_tcell_view.linear_pos_idx().id(), it->local_id()));
  }

  mesh::MeshBoundarySet<MeshConfig> &boundaries = m_tria.all_boundaries();
  const Uint num_current_boundaries             = boundaries.nb_domains();
  typename mesh::MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr bdry_facets =
      boundaries.create(MeshConfig::TDIM - 1, name);

  // The '1' below accounts for the interior domain, which has material id = 1
  const Uint material_id = num_current_boundaries + 2;
  bdry_facets->emplace_bdry_cell_ids(std::move(bdry_incidences), material_id);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void InteriorSolverCGHDG<MeshConfig>::add_boundary(const tria_t &tria_global,
                                                   const multi_edge_t &global_boundary,
                                                   const LeftRightOrientation orient,
                                                   const std::string &name)
{
  std::unique_ptr<std::vector<mesh::IncidenceEntry>> bdry_incidences(
      new std::vector<mesh::IncidenceEntry>());

  for (Uint bdry_block = 0; bdry_block < global_boundary.nb_all_facets(); ++bdry_block)
  {
    const mesh::TraceIncidences facet_incidences =
        global_boundary.facet_data(mesh::FlatIdx(bdry_block));

    const mesh::FlatIdx global_flat_id(facet_incidences.cell_id(orient));
    const Uint local_id                                      = facet_incidences.local_id(orient);
    const mesh::CellTopologyView<MeshConfig> glob_tcell_view = tria_global.cell(global_flat_id);
    const mesh::ActiveIdx local_active_id = m_glob_to_loc_map[glob_tcell_view.active_idx()];
    const mesh::CellTopologyView<MeshConfig> loc_tcell_view = m_tria.active_cell(local_active_id);
    bdry_incidences->push_back(
        mesh::IncidenceEntry(loc_tcell_view.linear_pos_idx().id(), local_id));
  }

  mesh::MeshBoundarySet<MeshConfig> &boundaries = m_tria.all_boundaries();
  const Uint num_current_boundaries             = boundaries.nb_domains();
  typename mesh::MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr bdry_facets =
      boundaries.create(MeshConfig::TDIM - 1, name);

  // The '1' below accounts for the interior domain, which has material id = 1
  const Uint material_id = num_current_boundaries + 2;
  bdry_facets->emplace_bdry_cell_ids(std::move(bdry_incidences), material_id);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS>
void InteriorSolverCGHDG<MeshConfig>::assemble(const RHS &rhs)
{
  setup_solver(rhs);
  assemble_interior(rhs);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename DirichletBC>
void InteriorSolverCGHDG<MeshConfig>::add_weak_bc(const std::string &boundary_name,
                                                  const DirichletBC &dirichlet_bc)
{
  const dof_map_t &cell_dofs = *(m_tria.dof_storage("dofs"));

  const mesh::MeshBoundarySet<MeshConfig> &boundaries = m_tria.all_boundaries();
  const typename mesh::MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr boundary_ptr =
      boundaries.domain(boundary_name);

  using bdry_dof_iterator_t =
      typename mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_dof_iterator;

  const bdry_dof_iterator_t bdry_begin = boundary_ptr->cbegin(cell_dofs);
  const bdry_dof_iterator_t bdry_end   = boundary_ptr->cend(cell_dofs);

  prepare_boundary_metric_data(common::make_iter_range(bdry_begin, bdry_end));
  // prepare_boundary_metric_data(tria, {boundary_name});

  constexpr Uint TDIM = MeshConfig::TDIM;

  clock_t start, end;
  Real elapsed;

  std::vector<Real> values;
  std::vector<Int> indices;

  m_sys_A->unlock();

  // -----------------------------------------
  // SYSTEM ASSEMBLY - LOOP OVER BOUNDARIES
  // -----------------------------------------

  start = clock();

  // First prepare volume data of cells adjacent to Dirichlet boundaries
  m_geo_cache_elems.flush();

  const auto basis_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  const auto eval_pt_set_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::PointSetTag(shape, m_quadrature_order, PointSetID::Gauss);
  };

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (bdry_dof_iterator_t bdry_it = bdry_begin; bdry_it != bdry_end; ++bdry_it)
  {
    const mesh::CellTopologyView<MeshConfig> adj_tcell_view = bdry_it->tcell();
    const mesh::MeshEntity cell = cell_dofs.active_cell(adj_tcell_view.active_idx());

    const mesh::PointSetTag std_reg_tag = cell.pt_set_id();
    const ElemShape cell_shape          = std_reg_tag.elem_shape();
    const Uint cell_order               = std_reg_tag.poly_order();
    const mesh::PointSetTagExt std_reg_tag_ext(std_reg_tag, P0, mesh::CellTransform::NO_TRANS, 0u);
    const mesh::sf::SFTag basis_tag = basis_generator(cell_shape, cell_order);

    const mesh::PointSetTag eval_pt_set_tag = eval_pt_set_generator(cell_shape, cell_order);
    const mesh::PointSetTagExt eval_pt_set_tag_ext(eval_pt_set_tag, P0,
                                                   mesh::CellTransform::NO_TRANS, 0u);

    const mesh::DiscreteElemKey cell_key(std_reg_tag_ext, basis_tag, eval_pt_set_tag_ext);

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        adj_tcell_view.pt_set_id(), cell.pt_set_id(), adj_tcell_view.coordinates());

    m_geo_cache_elems.push_back_to_buffer(cell_coords, cell_key);
  }

  m_geo_metric_elems.empty_buffer();
  m_geo_metric_elems.evaluate(m_geo_cache_elems, interpolation::RebuildMetricIndex{true});

  // Prepare face terms on Dirichlet boundaries

  m_geo_metric_elems_on_bdry.evaluate(m_geo_cache_elems_on_bdry,
                                      interpolation::RebuildMetricIndex{true});
  m_geo_metric_faces_on_bdry.evaluate(m_geo_cache_faces_on_bdry,
                                      interpolation::RebuildMetricIndex{true});

  // Prepare auxiliary data structures

  // Mass matrix and its inverse
  math::DenseDMat<Real> M_mat;
  math::DenseDMat<Real> M_mat_inv;

  // Derivative matrix
  std::array<math::DenseDMat<Real>, MeshConfig::GDIM> D_mat;
  std::array<math::DenseDMat<Real>, MeshConfig::GDIM> D_matT;
  // Edge mass matrices
  math::DenseDMat<Real> E_mat;
  math::DenseDMat<Real> E_mat_sum;
  std::array<math::DenseDMat<Real>, MeshConfig::GDIM> E_tilde_mat;
  std::array<math::DenseDMat<Real>, MeshConfig::GDIM> E_tilde_mat_sum;

  // Edge mass matrices 'F'
  math::DenseDMat<Real> F_mat;
  math::DenseDVec<Real> F_dot_lam_sum;
  std::array<math::DenseDMat<Real>, MeshConfig::GDIM> F_tilde_mat;
  std::array<math::DenseDVec<Real>, MeshConfig::GDIM> F_tilde_dot_lam_sum;

  math::DenseDMat<Real> A_weak;
  math::DenseDMat<Real> mat_tmp;
  math::DenseDVec<Real> vec_lambda;
  math::DenseDVec<Real> vec_tmp;

  math::DenseDVec<Real> b_weak;

  Uint idx_in_metric = 0;
  mesh::StdRegion std_reg;

  Uint bdry_block = 0;

  for (bdry_dof_iterator_t bdry_it = bdry_begin; bdry_it != bdry_end; ++bdry_it)
  {
    // SUB-STEP 1: get adjacent volume element and compute its mass and
    // derivative matrices
    const mesh::CellTopologyView<MeshConfig> adjacent_tcell_view = bdry_it->tcell();
    const mesh::MeshEntity adjacent_volume_cell =
        cell_dofs.active_cell(adjacent_tcell_view.active_idx());

    const Uint nb_adj_cell_dof = adjacent_volume_cell.nb_vert();

    M_mat.resize(nb_adj_cell_dof, nb_adj_cell_dof);
    M_mat.fill(0.0);
    M_mat_inv.resize(nb_adj_cell_dof, nb_adj_cell_dof);

    for (Uint d = 0; d < TDIM; ++d)
    {
      D_mat[d].resize(nb_adj_cell_dof, nb_adj_cell_dof);
      D_mat[d].fill(0.0);
      D_matT[d].resize(nb_adj_cell_dof, nb_adj_cell_dof);
      D_matT[d].fill(0.0);
    }

    cell_geo_metric_type const cell_geo_metric = m_geo_metric_elems.cellwise_values(bdry_block);

    build_adj_volume_matrices(cell_geo_metric, M_mat, D_mat);

    // Fill the transpose of derivative matrix D
    for (Uint d = 0; d < TDIM; ++d)
    {
      const math::DenseDMat<Real> &D = D_mat[d];
      math::DenseDMat<Real> &DT      = D_matT[d];

      for (Uint i = 0; i < D.rows(); ++i)
      {
        for (Uint j = 0; j < D.cols(); ++j)
        {
          DT(j, i) = D(i, j);
        }
      }
    }

    // Compute the inverse mass matrix
    M_mat.inv(M_mat_inv);

    /*
    std::cout << "M = \n" << M_mat << std::endl;
    for (Uint dim = 0; dim < TDIM; ++dim)
    {
      std::cout << "D[" << dim << "] = \n" << D_mat[dim] << std::endl;
    }

    std::cout << "\n\n" << std::endl;
    */

    // SUB-STEP 2: get all local Dirichlet faces of given element and
    // compute their local mass matrices (with support on faces)

    E_mat_sum.resize(nb_adj_cell_dof, nb_adj_cell_dof);
    E_mat.resize(nb_adj_cell_dof, nb_adj_cell_dof);
    for (Uint dim = 0; dim < TDIM; ++dim)
    {
      E_tilde_mat[dim].resize(nb_adj_cell_dof, nb_adj_cell_dof);
    }

    for (Uint dim = 0; dim < TDIM; ++dim)
    {
      E_tilde_mat_sum[dim].resize(nb_adj_cell_dof, nb_adj_cell_dof);
      E_tilde_mat_sum[dim].fill(0.0);
    }

    E_mat_sum.fill(0.0);

    F_dot_lam_sum.resize(nb_adj_cell_dof);
    F_dot_lam_sum.fill(0.0);

    for (Uint dim = 0; dim < TDIM; ++dim)
    {
      F_tilde_dot_lam_sum[dim].resize(nb_adj_cell_dof);
      F_tilde_dot_lam_sum[dim].fill(0.0);
    }

    A_weak.resize(nb_adj_cell_dof, nb_adj_cell_dof);
    A_weak.fill(0.0);
    mat_tmp.resize(nb_adj_cell_dof, nb_adj_cell_dof);

    b_weak.resize(nb_adj_cell_dof);
    b_weak.fill(0.0);
    vec_tmp.resize(nb_adj_cell_dof);

    const common::ArrayView<const Uint, _1D, Uint> local_elem_faces =
        m_bdry_face_local_id.const_block(bdry_block);

    // ------------
    // Loop over local boundary facets of one element
    // ------------

    for (Uint f = 0; f < local_elem_faces.size(); ++f)
    {
      cell_geo_metric_type const cell_geo_metric_on_trace =
          m_geo_metric_elems_on_bdry.cellwise_values(idx_in_metric);

      facet_geo_metric_type const facet_geo_metric =
          m_geo_metric_faces_on_bdry.cellwise_values(idx_in_metric);

      const Uint loc_face_id = local_elem_faces[f];

      /*
      const mesh::PointSetTag adjacent_elem_tag =
      adjacent_volume_cell.std_region_id(); const mesh::PointSetTagExt
      adjacent_elem_tag_ext(adjacent_elem_tag, P0,
                                                        mesh::CellTransform::NO_TRANS,
      loc_face_id);
      */

      const mesh::MeshEntity boundary_elem =
          adjacent_volume_cell.sub_entity(MeshConfig::TDIM - 1, loc_face_id);
      const mesh::PointSetTagExt boundary_elem_tag_ext(boundary_elem.pt_set_id(), P0,
                                                       mesh::CellTransform::NO_TRANS, 0);

      // const Uint nb_bdry_cell_dof = boundary_elem.nb_vert();

      E_mat.fill(0.0);
      for (Uint dim = 0; dim < TDIM; ++dim)
      {
        E_tilde_mat[dim].fill(0.0);
      }

      build_E_matrices(cell_geo_metric_on_trace, facet_geo_metric, E_mat, E_tilde_mat);

      E_mat_sum += E_mat;

      for (Uint dim = 0; dim < MeshConfig::GDIM; ++dim)
      {
        E_tilde_mat_sum[dim] += E_tilde_mat[dim];
      }

      const mesh::PointSetTag adjacent_elem_tag = adjacent_volume_cell.pt_set_id();
      std_reg.change_type(adjacent_elem_tag);

      std::shared_ptr<mesh::StdRegionEntity const> facet =
          std_reg.get().elem_entity(MeshConfig::GDIM - 1, loc_face_id);

      const Uint nb_facet_dof = (*facet).nb_vert();

      F_mat.resize(nb_adj_cell_dof, nb_facet_dof);
      F_mat.fill(0.0);

      for (Uint dim = 0; dim < MeshConfig::GDIM; ++dim)
      {
        F_tilde_mat[dim].resize(nb_adj_cell_dof, nb_facet_dof);
        F_tilde_mat[dim].fill(0.0);
      }

      const math::DenseConstMatView<Real> facet_coords = loc_interpolator.transfer_coords(
          adjacent_tcell_view.pt_set_id(MeshConfig::TDIM - 1, loc_face_id),
          boundary_elem.pt_set_id(),
          adjacent_tcell_view.coordinates(MeshConfig::TDIM - 1, loc_face_id));

      vec_lambda.resize(nb_facet_dof);

      for (Uint v = 0; v < boundary_elem.nb_vert(); ++v)
      {
        vec_lambda[v] = dirichlet_bc.eval(facet_coords.row_transpose(v));
      }

      build_F_matrices(cell_geo_metric_on_trace, facet_geo_metric, F_mat, F_tilde_mat);

      vec_tmp = F_mat * vec_lambda;
      F_dot_lam_sum += vec_tmp;

      for (Uint dim = 0; dim < MeshConfig::GDIM; ++dim)
      {
        vec_tmp = F_tilde_mat[dim] * vec_lambda;
        F_tilde_dot_lam_sum[dim] += vec_tmp;
      }

      idx_in_metric++;

    } // Loop over local face ids of one block

    const Real tau = 1.0;
    A_weak         = tau * E_mat_sum;

    b_weak = tau * F_dot_lam_sum;

    for (Uint dim = 0; dim < MeshConfig::GDIM; ++dim)
    {
      // -------------------------
      // Original - do not modify:
      // -------------------------
#if 0
         mat_tmp = D_mat[dim] * M_mat_inv * E_tilde_mat_sum[dim] -
                   E_tilde_mat_sum[dim] * M_mat_inv * D_mat[dim];
#endif

      // mat_tmp = E_tilde_mat_sum[dim] * M_mat_inv * (E_tilde_mat_sum[dim] -
      // D_mat[dim]);

      // Form 1 - fully symmetrized:
#if 1
      mat_tmp = E_tilde_mat_sum[dim] * M_mat_inv * E_tilde_mat_sum[dim] -
                D_matT[dim] * M_mat_inv * E_tilde_mat_sum[dim] -
                E_tilde_mat_sum[dim] * M_mat_inv * D_mat[dim];
      // + D_matT[dim] * M_mat_inv * D_mat[dim];
#endif

      // Form 2 - simplified by merging 2 terms using discrete integration by parts:
#if 0
         mat_tmp = D_mat[dim] * M_mat_inv * E_tilde_mat_sum[dim] -
                   E_tilde_mat_sum[dim] * M_mat_inv * D_mat[dim] +
                   D_matT[dim] * M_mat_inv * D_mat[dim];
#endif

      // Form 3 - original 'interior solve' without any modifications
#if 0
         mat_tmp = (E_tilde_mat_sum[dim] - D_matT[dim]) * M_mat_inv * D_matT[dim];
#endif

      // Form 4 - simplified
#if 0
         mat_tmp = D_mat[dim] * M_mat_inv * D_matT[dim];
#endif
      // std::cout << "mat_tmp = \n" << mat_tmp << std::endl;
      A_weak += mat_tmp;

      // vec_tmp = D_mat[dim] * M_mat_inv * F_tilde_dot_lam_sum[dim];
      vec_tmp = (E_tilde_mat_sum[dim] - D_matT[dim]) * M_mat_inv * F_tilde_dot_lam_sum[dim];
      b_weak += vec_tmp;
    }
    // std::cout << "BC weak = \n" << bc_weak << std::endl;

    indices.resize(nb_adj_cell_dof);
    values.resize(nb_adj_cell_dof);

    for (Uint vi = 0; vi < nb_adj_cell_dof; ++vi)
    {
      for (Uint vj = 0; vj < nb_adj_cell_dof; ++vj)
      {
        indices[vj] = adjacent_volume_cell.vertex(vj);
        values[vj]  = A_weak(vi, vj);
      }
      m_sys_A->add_values_to_row(adjacent_volume_cell.vertex(vi), values, indices);
      m_sys_RHS->add_value(adjacent_volume_cell.vertex(vi), b_weak[vi], 0);
    }

    bdry_block++;

  } // Loop over Dirichlet boundary facet blocks

  m_sys_A->lock();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Global system assembly (bc): " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void InteriorSolverCGHDG<MeshConfig>::solve()
{
  solve_ls(m_sys_A, m_sys_RHS, m_solution);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const typename InteriorSolverCGHDG<MeshConfig>::tria_t &InteriorSolverCGHDG<
    MeshConfig>::local_mesh() const
{
  return m_tria;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const typename InteriorSolverCGHDG<MeshConfig>::dof_map_t &InteriorSolverCGHDG<
    MeshConfig>::local_dofs() const
{
  common::PtrHandle<const dof_map_t> dofs = m_tria.dof_storage("dofs");
  return *dofs;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const interpolation::VectorMeshFunction<Real> &InteriorSolverCGHDG<MeshConfig>::solution() const
{
  return m_solution;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void InteriorSolverCGHDG<MeshConfig>::write_to_gmsh(const std::string &filename) const
{
  mesh::gmsh::GmshWriter gmsh_writer;
  gmsh_writer.write_mesh_to_file(m_tria, "dofs", filename);
  gmsh_writer.append_nodal_function_to_file(m_tria, filename, m_solution, "u_num");
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
std::unique_ptr<typename InteriorSolverCGHDG<MeshConfig>::cell_buffer_t> InteriorSolverCGHDG<
    MeshConfig>::create_sub_mesh_buffer(const tria_t &tria, const dof_map_t &dofs,
                                        const std::vector<mesh::ActiveIdx> &cell_ids,
                                        const Uint material_id)
{
  mesh::adapt::LocalInterpolator loc_interpolator;

  std::vector<Uint> vert_buffer;
  std::vector<Real> coord_buffer;

  Uint id_in_buffer = 0;

  std::unique_ptr<cell_buffer_t> cell_buffer(new cell_buffer_t());
  for (auto active_id : cell_ids)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = tria.active_cell(active_id);
    const mesh::MeshEntity active_cell                  = dofs.active_cell(active_id);

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), active_cell.pt_set_id(), tcell_view.coordinates());
    const Uint nb_cell_dof = active_cell.nb_vert();
    vert_buffer.resize(nb_cell_dof);

    for (Uint v = 0; v < nb_cell_dof; ++v)
    {
      vert_buffer[v] = active_cell.vertex(v);
    }

    coord_buffer.resize(cell_coords.rows() * cell_coords.cols());
    Uint insert_id = 0;

    for (Uint r = 0; r < cell_coords.rows(); ++r)
    {
      for (Uint c = 0; c < cell_coords.cols(); ++c)
      {
        coord_buffer[insert_id++] = cell_coords(r, c);
      }
    }

    cell_buffer->push_back_cell(id_in_buffer++, active_cell.pt_set_id(), vert_buffer, material_id,
                                coord_buffer);
  }

  return cell_buffer;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS>
void InteriorSolverCGHDG<MeshConfig>::setup_solver(const RHS &rhs)
{
  clock_t start, end;
  Real elapsed;

  const dof_map_t &cell_dofs = *(m_tria.dof_storage("dofs"));

  const Uint nb_nodes_in_mesh = cell_dofs.nb_nodes();

  m_geo_cache_elems.clear();
  m_geo_metric_elems.clear();

  // Setup metric containers for boundary terms
  // prepare_boundary_metric_data(tria, boundary_names);

  // Create system matrix
  std::cout << "Allocating global system matrix" << std::endl;
  m_sys_A = std::make_shared<ls::TpetraCrsMatrix<Real>>(nb_nodes_in_mesh);
  std::cout << " ... done" << std::endl;
  m_sys_RHS = std::make_shared<ls::TpetraMultiVector<Real>>(m_sys_A->map());

  std::vector<Real> values, boundary_values;
  std::vector<Int> indices, boundary_indices;

  m_quadrature_order = P1;

  start = clock();

  // Loop over all cells going type by type find the maximum quadrature order
  // The initialize the geometry metric
  for (const typename dof_map_t::const_dof_range_typed &cell_group :
       cell_dofs.all_active_dof_groups())
  {
    const mesh::MeshEntity first_cell = (*cell_group.begin()).mesh_entity();

    const mesh::PointSetTag first_cell_tag = first_cell.pt_set_id();
    const Uint poly_order                  = first_cell_tag.poly_order();

    // quadrature_order = std::max(quadrature_order, 2 * (poly_order - 1));
    m_quadrature_order = std::max(m_quadrature_order, 2 * poly_order);
  }

  const auto basis_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  const auto eval_pt_set_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::PointSetTag(shape, m_quadrature_order, PointSetID::Gauss);
  };

  typename interpolation::FunctionSpace<MeshConfig>::ptr fs_geometry =
      std::make_shared<interpolation::FunctionSpace<MeshConfig>>();
  fs_geometry->set_reference_fe_values(cell_dofs.as_range(), basis_generator,
                                       eval_pt_set_generator);

  m_geo_cache_elems.allocate(fs_geometry->discrete_elements().cbegin(),
                             fs_geometry->discrete_elements().cend(), cell_dofs.nb_active_cells());
  m_geo_metric_elems.allocate_buffer(fs_geometry->discrete_elements().cbegin(),
                                     fs_geometry->discrete_elements().cend(),
                                     cell_dofs.nb_active_cells());

  // fs_geometry->print_element_types();

  // ------------------------------------------------------
  // PART I:
  // Loop over all cells going type by type and:
  // 1) Push back the cell coordinates into geometry metric
  // 2) Initialize the entries in the stiffness matrix
  // ------------------------------------------------------

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (const typename dof_map_t::const_dof_range_typed &dof_group :
       cell_dofs.all_active_dof_groups())
  {
    const mesh::MeshEntity first_cell = (*dof_group.begin()).mesh_entity();
    // std::cout << "First cell in group = " << first_cell << std::endl;

    values.resize(first_cell.nb_vert());
    for (Uint v = 0; v < first_cell.nb_vert(); ++v)
    {
      values[v] = 0.0;
    }

    indices.resize(first_cell.nb_vert());

    for (typename dof_map_t::const_dof_iterator_typed cell_iter = dof_group.begin();
         cell_iter != dof_group.end(); ++cell_iter)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_view = cell_iter->tcell();
      const mesh::MeshEntity cell                         = cell_iter->mesh_entity();

      const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
          tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

      const mesh::PointSetTag std_reg_tag = cell.pt_set_id();
      const ElemShape cell_shape          = std_reg_tag.elem_shape();
      const Uint cell_order               = std_reg_tag.poly_order();
      const mesh::PointSetTagExt std_reg_tag_ext(std_reg_tag, P0, mesh::CellTransform::NO_TRANS,
                                                 0u);
      const mesh::sf::SFTag basis_tag = basis_generator(cell_shape, cell_order);

      const mesh::PointSetTag eval_pt_set_tag = eval_pt_set_generator(cell_shape, cell_order);
      const mesh::PointSetTagExt eval_pt_set_tag_ext(eval_pt_set_tag, P0,
                                                     mesh::CellTransform::NO_TRANS, 0u);

      const mesh::DiscreteElemKey cell_key(std_reg_tag_ext, basis_tag, eval_pt_set_tag_ext);

      m_geo_cache_elems.push_back_to_buffer(cell_coords, cell_key);

      // std::cout << cell << std::endl;
      for (Uint vi = 0; vi < cell.nb_vert(); ++vi)
      {

        for (Uint vj = 0; vj < cell.nb_vert(); ++vj)
        {
          indices[vj] = cell.vertex(vj);
        }
        m_sys_A->insert_values_in_row(cell.vertex(vi), values, indices);
        // (*sys_RHS).insert_value(cell.vertex(vi), 0.0, 0);
      }
    }
  } // Loop over all cell groups

  // Prepare the right-hand side
  m_sys_RHS->fill(0.0);

  //  m_sys_A->lock_structure();
  //  m_sys_A->print_structure_to_file("sparsity_Poisson.ps");

  m_sys_A->lock();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Sparse matrix setup (weak): " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS>
void InteriorSolverCGHDG<MeshConfig>::assemble_interior(const RHS &rhs)
{
  constexpr Uint TDIM = MeshConfig::TDIM;

  clock_t start, end;
  Real elapsed;

  math::DenseDMat<Real> local_stiffness;
  math::DenseDVec<Real> local_rhs;

  std::vector<Real> values;
  std::vector<Int> indices;

  std::array<math::DenseConstMatView<Real>, MeshConfig::GDIM> sf_deriv_phys;

  m_sys_A->unlock();
  // m_sys_A.fill(0.0);

  const dof_map_t &cell_dofs = *(m_tria.dof_storage("dofs"));

  for (const typename dof_map_t::const_dof_range_typed &cell_group :
       cell_dofs.all_active_dof_groups())
  {
    const mesh::MeshEntity first_cell = (*cell_group.begin()).mesh_entity();

    values.resize(first_cell.nb_vert());
    indices.resize(first_cell.nb_vert());

    // Resize local stiffness matrix for given element type
    const Uint n_loc_dof = first_cell.nb_vert();
    local_stiffness.resize(n_loc_dof, n_loc_dof);

    // Resize local rhs for given element type
    local_rhs.resize(first_cell.nb_vert());

    // -----------------------------------------
    // SYSTEM ASSEMBLY - MAIN LOOP OVER ELEMENTS
    // -----------------------------------------

    start = clock();

    m_geo_metric_elems.evaluate(m_geo_cache_elems, interpolation::RebuildMetricIndex{true});

    Uint idx_in_metric = 0;

    for (typename dof_map_t::const_dof_iterator_typed cell_iter = cell_group.begin();
         cell_iter != cell_group.end(); ++cell_iter)
    {
      const mesh::MeshEntity cell = cell_iter->mesh_entity();

      cell_geo_metric_type const cell_geo_metric =
          m_geo_metric_elems.cellwise_values(idx_in_metric);
      idx_in_metric++;

      math::DenseConstMatView<Real> const qd_pts_phys       = cell_geo_metric.interpolated_coords();
      math::DenseConstVecView<Real> const jacobian_in_qd_pt = cell_geo_metric.jdet();
      math::DenseDVec<Real> const &weight_in_qd_pt          = cell_geo_metric.pt_weights();

      for (Uint d = 0; d < MeshConfig::TDIM; ++d)
      {
        sf_deriv_phys[d] = cell_geo_metric.sf_derivatives(d);
      }

      math::DenseDMat<Real> const &V = cell_geo_metric.reference_sf_values();

      const Uint nb_qd_pts = jacobian_in_qd_pt.size();

      // Start assembly of the local stiffness matrix

      local_stiffness.fill(0.0);
      local_rhs.fill(0.0);

      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        const Real wq = jacobian_in_qd_pt[q] * weight_in_qd_pt[q];
        for (Uint vi = 0; vi < cell.nb_vert(); ++vi)
        {
          for (Uint vj = 0; vj < cell.nb_vert(); ++vj)
          {
            for (Uint dim = 0; dim < TDIM; ++dim)
            {
              local_stiffness(vi, vj) += wq * sf_deriv_phys[dim](q, vi) * sf_deriv_phys[dim](q, vj);
            }
          } // Loop over vj

          local_rhs[vi] += wq * V(q, vi) * rhs.eval(qd_pts_phys.row_transpose(q));

        } // Loop over vi
      }   // Loop over quadrature points

      // ------------------------------

      // std::cout << local_stiffness << std::endl;

      // Distribute to the global system

      for (Uint vi = 0; vi < cell.nb_vert(); ++vi)
      {
        for (Uint vj = 0; vj < cell.nb_vert(); ++vj)
        {
          indices[vj] = cell.vertex(vj);
          values[vj]  = local_stiffness(vi, vj);
        }
        m_sys_A->add_values_to_row(cell.vertex(vi), values, indices);
        m_sys_RHS->add_value(cell.vertex(vi), local_rhs[vi], 0);
      }
    }
  } // Loop over all cell groups

  m_sys_A->lock();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Global system assembly (interior): " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename BdryDofIterator>
void InteriorSolverCGHDG<MeshConfig>::prepare_boundary_metric_data(
    const common::IteratorRange<BdryDofIterator> &boundary)
{
  const dof_map_t &cell_dofs = *(m_tria.dof_storage("dofs"));

  Uint nb_active_cells_on_Dir_bdry = 0;
  Uint nb_Dir_bdry_faces           = 0;

  // Stores number of weak boundary faces of ANY ACTIVE elements
  // We will later on build a filtered vector which stores the weak
  // boundary faces only for ACTIVE elements that actually are incident
  // to weak boundaries
  std::vector<std::forward_list<Uint>> elem_to_weak_facets_map;
  elem_to_weak_facets_map.resize(cell_dofs.nb_active_cells());

  const BdryDofIterator bdry_begin = boundary.begin();
  const BdryDofIterator bdry_end   = boundary.end();

  for (BdryDofIterator bdry_it = bdry_begin; bdry_it != bdry_end; ++bdry_it)
  {
    nb_Dir_bdry_faces++;

    const mesh::CellTopologyView<MeshConfig> adj_tcell_view = bdry_it->tcell();
    const mesh::MeshEntity cell = cell_dofs.active_cell(adj_tcell_view.active_idx());
    const Uint local_id         = bdry_it->local_id();
    elem_to_weak_facets_map[cell.idx()].push_front(local_id);
  }

  // Count how many active cells touch weak boundaries
  for (const auto &list : elem_to_weak_facets_map)
  {
    if (!list.empty())
    {
      nb_active_cells_on_Dir_bdry++;
    }
  }

  std::cout << "Number of Dirichlet faces on boundary = " << nb_Dir_bdry_faces << std::endl;
  std::cout << "Number of cells touching Dirichlet boundary = " << nb_active_cells_on_Dir_bdry
            << std::endl;

  std::unique_ptr<std::vector<Uint>> bdry_face_local_id(new std::vector<Uint>());
  bdry_face_local_id->reserve(nb_Dir_bdry_faces);
  bdry_face_local_id->resize(0);

  // This vector nb_bdry_faces_per_elem_offsets holds number of weak boundary
  // faces of each element THAT IS INCIDENT TO WEAK BOUNDARY Corresponding
  // active element ids are stored in m_bdry_elem_id
  std::unique_ptr<std::vector<Uint>> bdry_faces_per_elem_offsets(new std::vector<Uint>());
  bdry_faces_per_elem_offsets->reserve(nb_active_cells_on_Dir_bdry + 1);
  bdry_faces_per_elem_offsets->resize(0);
  bdry_faces_per_elem_offsets->push_back(0);

  for (BdryDofIterator bdry_it = bdry_begin; bdry_it != bdry_end; ++bdry_it)
  {
    const mesh::CellTopologyView<MeshConfig> adj_tcell_view = bdry_it->tcell();
    const mesh::MeshEntity cell = cell_dofs.active_cell(adj_tcell_view.active_idx());

    if (!elem_to_weak_facets_map[cell.idx()].empty())
    {
      Uint nb_local_faces = 0;
      for (auto item : elem_to_weak_facets_map[cell.idx()])
      {
        bdry_face_local_id->push_back(item);
        nb_local_faces++;
      }
      bdry_faces_per_elem_offsets->push_back(nb_local_faces);
      // m_bdry_elem_id.push_back(cell.idx());
    }
  }

  for (Uint i = 1; i < bdry_faces_per_elem_offsets->size(); ++i)
  {
    (*bdry_faces_per_elem_offsets)[i] += (*bdry_faces_per_elem_offsets)[i - 1];
  }

  m_bdry_face_local_id.build_from_offsets(std::move(bdry_face_local_id),
                                          std::move(bdry_faces_per_elem_offsets));

  // ---

  std::vector<mesh::DiscreteElemKey> boundary_fe_type_map;
  std::vector<mesh::DiscreteElemKey> boundary_adjacent_vol_fe_type_map;

  /*
  const std::vector<std::shared_ptr<mesh::BoundaryFacets<MeshConfig,
  MeshConfig::TDIM-1>>> &boundary_domains = boundaries.all_domains();
  */

  Uint bdry_block = 0;
  for (BdryDofIterator bdry_it = bdry_begin; bdry_it != bdry_end; ++bdry_it)
  {
    const mesh::CellTopologyView<MeshConfig> adj_tcell_view = bdry_it->tcell();
    const mesh::MeshEntity adjacent_volume_cell =
        cell_dofs.active_cell(adj_tcell_view.active_idx());

    const common::ArrayView<const Uint, _1D, Uint> local_elem_faces =
        m_bdry_face_local_id.const_block(bdry_block);

    bdry_block++;

    for (Uint f = 0; f < local_elem_faces.size(); ++f)
    {
      const Uint loc_face_id = local_elem_faces[f];

      const mesh::PointSetTag adjacent_elem_tag = adjacent_volume_cell.pt_set_id();
      const mesh::PointSetTagExt adjacent_elem_tag_ext(adjacent_elem_tag, P0,
                                                       mesh::CellTransform::NO_TRANS, loc_face_id);
      const mesh::sf::SFTag vol_cell_sf_tag(adjacent_elem_tag.elem_shape(), SFunc::Lagrange,
                                            adjacent_elem_tag.poly_order(), ModalBasis::Modal);
      const mesh::PointSetTag adjacent_elem_quad_tag(adjacent_elem_tag.elem_shape(),
                                                     2. * adjacent_elem_tag.poly_order(),
                                                     PointSetID::FaceGauss);
      const mesh::PointSetTagExt adjacent_elem_quad_tag_ext(
          adjacent_elem_quad_tag, P0, mesh::CellTransform::RESTRICT_TO_CODIM_1, loc_face_id);

      const mesh::DiscreteElemKey adjacent_elem_key(adjacent_elem_tag_ext, vol_cell_sf_tag,
                                                    adjacent_elem_quad_tag_ext);

      bool key_found = false;
      for (const auto &elem : boundary_adjacent_vol_fe_type_map)
      {
        if (elem == adjacent_elem_key)
        {
          key_found = true;
          break;
        }
      }
      if (!key_found)
      {
        boundary_adjacent_vol_fe_type_map.push_back(adjacent_elem_key);
      }

      // ---

      const mesh::MeshEntity boundary_elem =
          adjacent_volume_cell.sub_entity(MeshConfig::TDIM - 1, loc_face_id);

      const mesh::PointSetTag boundary_elem_tag = boundary_elem.pt_set_id();
      const mesh::PointSetTagExt boundary_elem_tag_ext(boundary_elem_tag, P0,
                                                       mesh::CellTransform::NO_TRANS, 0);

      const mesh::sf::SFTag boundary_elem_sf_tag(boundary_elem_tag.elem_shape(), SFunc::Lagrange,
                                                 boundary_elem_tag.poly_order(), ModalBasis::Modal);
      const mesh::PointSetTag boundary_elem_quad_tag(
          boundary_elem_tag.elem_shape(), 2. * boundary_elem_tag.poly_order(), PointSetID::Gauss);

      const mesh::PointSetTagExt boundary_elem_quad_tag_ext(boundary_elem_quad_tag, P0,
                                                            mesh::CellTransform::NO_TRANS, 0);

      const mesh::DiscreteElemKey boundary_elem_key(boundary_elem_tag_ext, boundary_elem_sf_tag,
                                                    boundary_elem_quad_tag_ext);

      key_found = false;
      for (const auto &elem : boundary_fe_type_map)
      {
        if (elem == boundary_elem_key)
        {
          key_found = true;
          break;
        }
      }
      if (!key_found)
      {
        boundary_fe_type_map.push_back(boundary_elem_key);
      }
    }
  }

  std::cout << "Number of values in map: " << boundary_fe_type_map.size() << std::endl;
  for (const auto &elem : boundary_fe_type_map)
  {
    std::cout << elem << std::endl;
  }

  const Uint nb_bdry_cells = m_bdry_face_local_id.size();

  m_geo_cache_elems_on_bdry.clear();
  m_geo_metric_elems_on_bdry.clear();

  m_geo_cache_elems_on_bdry.allocate(boundary_adjacent_vol_fe_type_map.cbegin(),
                                     boundary_adjacent_vol_fe_type_map.cend(), nb_bdry_cells);
  m_geo_metric_elems_on_bdry.allocate_buffer(boundary_adjacent_vol_fe_type_map.cbegin(),
                                             boundary_adjacent_vol_fe_type_map.cend(),
                                             nb_bdry_cells);

  m_geo_cache_faces_on_bdry.clear();
  m_geo_metric_faces_on_bdry.clear();

  m_geo_cache_faces_on_bdry.allocate(boundary_fe_type_map.cbegin(), boundary_fe_type_map.cend(),
                                     nb_bdry_cells);
  m_geo_metric_faces_on_bdry.allocate_buffer(boundary_fe_type_map.cbegin(),
                                             boundary_fe_type_map.cend(), nb_bdry_cells);

  mesh::adapt::LocalInterpolator loc_interpolator;

  bdry_block = 0;

  for (BdryDofIterator bdry_it = bdry_begin; bdry_it != bdry_end; ++bdry_it)
  {
    const mesh::CellTopologyView<MeshConfig> adj_tcell_view = bdry_it->tcell();
    const mesh::MeshEntity adjacent_volume_cell =
        cell_dofs.active_cell(adj_tcell_view.active_idx());

    const common::ArrayView<const Uint, _1D, Uint> local_elem_faces =
        m_bdry_face_local_id.const_block(bdry_block);

    bdry_block++;

    for (Uint f = 0; f < local_elem_faces.size(); ++f)
    {
      const Uint loc_face_id = local_elem_faces[f];

      const mesh::PointSetTag adjacent_elem_tag = adjacent_volume_cell.pt_set_id();
      const mesh::PointSetTagExt adjacent_elem_tag_ext(adjacent_elem_tag, P0,
                                                       mesh::CellTransform::NO_TRANS, loc_face_id);

      const mesh::sf::SFTag vol_cell_sf_tag(adjacent_elem_tag.elem_shape(), SFunc::Lagrange,
                                            adjacent_elem_tag.poly_order(), ModalBasis::Modal);
      const mesh::PointSetTag adjacent_elem_quad_tag(adjacent_elem_tag.elem_shape(),
                                                     2. * adjacent_elem_tag.poly_order(),
                                                     PointSetID::FaceGauss);
      const mesh::PointSetTagExt adjacent_elem_quad_tag_ext(
          adjacent_elem_quad_tag, P0, mesh::CellTransform::RESTRICT_TO_CODIM_1, loc_face_id);

      const mesh::DiscreteElemKey adjacent_elem_key(adjacent_elem_tag_ext, vol_cell_sf_tag,
                                                    adjacent_elem_quad_tag_ext);

      // ---

      const mesh::MeshEntity boundary_elem =
          adjacent_volume_cell.sub_entity(MeshConfig::TDIM - 1, loc_face_id);
      const mesh::PointSetTag boundary_elem_tag = boundary_elem.pt_set_id();
      const mesh::PointSetTagExt boundary_elem_tag_ext(boundary_elem_tag, P0,
                                                       mesh::CellTransform::NO_TRANS, 0);

      const mesh::sf::SFTag boundary_elem_sf_tag(boundary_elem_tag.elem_shape(), SFunc::Lagrange,
                                                 boundary_elem_tag.poly_order(), ModalBasis::Modal);
      const mesh::PointSetTag boundary_elem_quad_tag(
          boundary_elem_tag.elem_shape(), 2. * boundary_elem_tag.poly_order(), PointSetID::Gauss);

      const mesh::PointSetTagExt boundary_elem_quad_tag_ext(boundary_elem_quad_tag, P0,
                                                            mesh::CellTransform::NO_TRANS, 0);

      const mesh::DiscreteElemKey boundary_elem_key(boundary_elem_tag_ext, boundary_elem_sf_tag,
                                                    boundary_elem_quad_tag_ext);

      const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
          adj_tcell_view.pt_set_id(), adjacent_volume_cell.pt_set_id(),
          adj_tcell_view.coordinates());

      m_geo_cache_elems_on_bdry.push_back_to_buffer(cell_coords, adjacent_elem_key);

      const math::DenseConstMatView<Real> facet_coords = loc_interpolator.transfer_coords(
          adj_tcell_view.pt_set_id(MeshConfig::TDIM - 1, loc_face_id), boundary_elem.pt_set_id(),
          adj_tcell_view.coordinates(MeshConfig::TDIM - 1, loc_face_id));

      m_geo_cache_faces_on_bdry.push_back_to_buffer(facet_coords, boundary_elem_key);
    }
  }

  /*
  std::cout << "PRINTING CONTENTS OF CACHE FOR ADJACENT BOUNDARY ELEMENTS:" <<
  std::endl; m_geo_cache_elems_on_bdry.print_types(); std::cout << "PRINTING
  CONTENTS OF METRIC FOR ADJACENT BOUNDARY ELEMENTS:" << std::endl;
  m_geo_metric_elems_on_bdry.print_types();
  std::cout << "PRINTING CONTENTS OF CACHE FOR BOUNDARY FACES:" << std::endl;
  m_geo_cache_faces_on_bdry.print_types();
  std::cout << "PRINTING CONTENTS OF METRIC FOR BOUNDARY FACES:" << std::endl;
  m_geo_metric_faces_on_bdry.print_types();
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline void InteriorSolverCGHDG<MeshConfig>::build_adj_volume_matrices(
    const cell_geo_metric_type &cell_geo_metric, math::DenseDMat<Real> &M_mat,
    std::array<math::DenseDMat<Real>, MeshConfig::GDIM> &D_mat)
{
  // math::DenseConstMatView<Real> const qd_pts_phys =
  // cell_geo_metric.interpolated_coords();
  math::DenseConstVecView<Real> const jacobian_in_qd_pt = cell_geo_metric.jdet();
  math::DenseDVec<Real> const &weight_in_qd_pt          = cell_geo_metric.pt_weights();

  // Vandermonde matrix for shape function values
  math::DenseDMat<Real> const &V = cell_geo_metric.reference_sf_values();

  std::array<math::DenseConstMatView<Real>, MeshConfig::GDIM> dV_phys;

  // Vandermonde matrices for shape function derivatives
  for (Uint d = 0; d < MeshConfig::TDIM; ++d)
  {
    dV_phys[d] = cell_geo_metric.sf_derivatives(d);
  }

  const Uint nb_local_dof = M_mat.rows();
  const Uint nb_qd_pts    = jacobian_in_qd_pt.size();

  for (Uint q = 0; q < nb_qd_pts; ++q)
  {
    const Real wq = jacobian_in_qd_pt[q] * weight_in_qd_pt[q];
    for (Uint vi = 0; vi < nb_local_dof; ++vi)
    {
      for (Uint vj = 0; vj < nb_local_dof; ++vj)
      {
        M_mat(vi, vj) += wq * V(q, vi) * V(q, vj);
        for (Uint dim = 0; dim < MeshConfig::GDIM; ++dim)
        {
          D_mat[dim](vi, vj) += wq * V(q, vi) * dV_phys[dim](q, vj);
        }
      } // Loop over vj

    } // Loop over vi
  }   // Loop over quadrature points
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline void InteriorSolverCGHDG<MeshConfig>::build_E_matrices(
    const cell_geo_metric_type &cell_geo_metric_on_trace,
    const facet_geo_metric_type &facet_geo_metric, math::DenseDMat<Real> &E_mat,
    std::array<math::DenseDMat<Real>, MeshConfig::GDIM> &E_tilde_mat)
{
  math::DenseConstMatView<Real> const qd_pts_phys_cell =
      cell_geo_metric_on_trace.interpolated_coords();
  math::DenseConstMatView<Real> const qd_pts_phys_facet = facet_geo_metric.interpolated_coords();

  math::DenseConstVecView<Real> const jacobian_in_qd_pt = facet_geo_metric.jdet();
  math::DenseDVec<Real> const &weight_in_qd_pt          = facet_geo_metric.pt_weights();
  math::DenseConstMatView<Real> const facet_normals     = facet_geo_metric.normals();

  // Vandermonde matrix for shape function values
  math::DenseDMat<Real> const &V = cell_geo_metric_on_trace.reference_sf_values();

  const Uint nb_cell_dof = E_mat.rows();

  const Uint nb_qd_pts = jacobian_in_qd_pt.size();

  for (Uint q = 0; q < nb_qd_pts; ++q)
  {
    const Real wq = jacobian_in_qd_pt[q] * weight_in_qd_pt[q];
    for (Uint vi = 0; vi < nb_cell_dof; ++vi)
    {
      for (Uint vj = 0; vj < nb_cell_dof; ++vj)
      {
        const Real term = wq * V(q, vi) * V(q, vj);
        E_mat(vi, vj) += term;
        for (Uint dim = 0; dim < MeshConfig::GDIM; ++dim)
        {
          E_tilde_mat[dim](vi, vj) += facet_normals(q, dim) * term;
        }
      } // Loop over vj

    } // Loop over vi
  }   // Loop over quadrature points

  /*
  std::cout <<
  "---------------------------------------------------------------"
  << std::endl;

  std::cout << "Quadrature points in cell = \n" << qd_pts_phys_cell <<
  std::endl; std::cout << "Quadrature points in facet = \n" <<
  qd_pts_phys_facet
  << std::endl;

  std::cout << "E = \n" << E_mat << std::endl;
  for (Uint dim = 0; dim < MeshConfig::GDIM; ++dim)
  {
    std::cout << "tilde(E)[" << dim << "] = \n" << E_tilde_mat[dim] <<
  std::endl;
  }

  std::cout <<
  "---------------------------------------------------------------"
  << std::endl;
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline void InteriorSolverCGHDG<MeshConfig>::build_F_matrices(
    const cell_geo_metric_type &cell_geo_metric_on_trace,
    const facet_geo_metric_type &facet_geo_metric, math::DenseDMat<Real> &F_mat,
    std::array<math::DenseDMat<Real>, MeshConfig::GDIM> &F_tilde_mat)
{
  math::DenseConstMatView<Real> const qd_pts_phys_cell =
      cell_geo_metric_on_trace.interpolated_coords();
  math::DenseConstMatView<Real> const qd_pts_phys_facet = facet_geo_metric.interpolated_coords();

  math::DenseConstVecView<Real> const jacobian_in_qd_pt = facet_geo_metric.jdet();
  math::DenseDVec<Real> const &weight_in_qd_pt          = facet_geo_metric.pt_weights();
  math::DenseConstMatView<Real> const facet_normals     = facet_geo_metric.normals();

  // Vandermonde matrix for shape function values
  math::DenseDMat<Real> const &V_cell  = cell_geo_metric_on_trace.reference_sf_values();
  math::DenseDMat<Real> const &V_facet = facet_geo_metric.reference_sf_values();

  const Uint nb_cell_dof  = cell_geo_metric_on_trace.nb_dof_in_cell();
  const Uint nb_facet_dof = facet_geo_metric.nb_dof_in_cell();

  const Uint nb_qd_pts = jacobian_in_qd_pt.size();

  for (Uint q = 0; q < nb_qd_pts; ++q)
  {
    const Real wq = jacobian_in_qd_pt[q] * weight_in_qd_pt[q];
    for (Uint vi = 0; vi < nb_cell_dof; ++vi)
    {
      for (Uint vj = 0; vj < nb_facet_dof; ++vj)
      {
        const Real term = wq * V_cell(q, vi) * V_facet(q, vj);
        F_mat(vi, vj) += term;
        for (Uint dim = 0; dim < MeshConfig::GDIM; ++dim)
        {
          F_tilde_mat[dim](vi, vj) += facet_normals(q, dim) * term;
        }
      } // Loop over vj

    } // Loop over vi
  }   // Loop over quadrature points

  /*
  std::cout <<
  "---------------------------------------------------------------"
  << std::endl;

  std::cout << "Quadrature points in cell = \n" << qd_pts_phys_cell <<
  std::endl; std::cout << "Quadrature points in facet = \n" <<
  qd_pts_phys_facet
  << std::endl;

  std::cout << "F = \n" << F_mat << std::endl;
  for (Uint dim = 0; dim < MeshConfig::GDIM; ++dim)
  {
    std::cout << "tilde(F)[" << dim << "] = \n" << F_tilde_mat[dim] <<
  std::endl;
  }

  std::cout <<
  "---------------------------------------------------------------"
  << std::endl;
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void InteriorSolverCGHDG<MeshConfig>::solve_ls(
    std::shared_ptr<ls::TpetraCrsMatrix<Real>> &sys_A,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &sys_RHS,
    interpolation::VectorMeshFunction<Real> &solution)
{
  std::shared_ptr<ls::TpetraMultiVector<Real>> x(new ls::TpetraMultiVector<Real>(sys_A->map()));

  ls::LSTpetra<Real> lin_system(ls::SparseSolverType::eIterative);
  // lin_system.configure(sys_A, sys_RHS, x, true,
  // false);
  lin_system.initialize_solver(sys_A, sys_RHS, x, false);

  std::shared_ptr<ls::TrilinosPC<Real>> preconditioner = std::make_shared<ls::IfpackPC<Real>>();
  preconditioner->create("ILUT", sys_A);

  lin_system.connect_preconditioner(preconditioner);
  lin_system.update_after_mat_values_change(sys_A, sys_RHS, x, false);
  lin_system.solve(550, 1.e-12);

  // Create a mesh function to copy the result

  typedef typename pdekit::result_of::dof_map_t<MeshConfig> cell_dofs_type;
  const Uint nb_dofs = sys_RHS->size();

  solution.resize(1, nb_dofs);
  solution.fill(0.0);

  for (Uint n = 0; n < nb_dofs; ++n)
  {
    typename interpolation::VectorMeshFunction<Real>::entry_type nodal_value = solution.value(n);
    nodal_value[0]                                                           = (*x).value(n, 0);
  }
}

// ----------------------------------------------------------------------------

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
