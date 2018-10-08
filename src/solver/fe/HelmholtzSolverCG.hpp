#ifndef PDEKIT_Solver_FE_Helmholtz_Solver_CG_hpp
#define PDEKIT_Solver_FE_Helmholtz_Solver_CG_hpp

/// Standard template library headers
#include <array>
#include <cmath>
#include <ctime>
#include <forward_list>
#include <iostream>
#include <map>

/// PDEKIT headers
#include "common/MPI/MPIEnv.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/MathConstants.hpp"
#include "mesh/MeshConfig.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class HelmholtzSolverCG
{
  public:
  /// TYPEDEFS
  using mesh_type      = typename mesh::Tria<MeshConfig>;
  using cell_dofs_type = typename pdekit::result_of::dof_map_t<MeshConfig>;

  /// Default constructor
  HelmholtzSolverCG();

  /// Destructor
  ~HelmholtzSolverCG();

  template <typename RHS, typename DirichletBC>
  void assemble(const std::string &method, const mesh_type &in_mesh, const RHS &rhs,
                const std::vector<std::string> &boundary_names, const DirichletBC &dirichlet_bc,
                std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
                std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs);

  private:
  using cell_geo_metric_type =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             MeshConfig::GDIM>::cellwise_metric;
  using facet_geo_metric_type =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             MeshConfig::GDIM - 1>::cellwise_metric;

  template <typename RHS, typename DirichletBC>
  void setup_solver_standard(const mesh_type &in_mesh, const RHS &rhs,
                             const std::vector<std::string> &boundary_names,
                             const DirichletBC &dirichlet_bc,
                             std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
                             std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs);

  template <typename RHS, typename DirichletBC>
  void setup_solver_Nitsche(const mesh_type &in_mesh, const RHS &rhs,
                            const std::vector<std::string> &boundary_names,
                            const DirichletBC &dirichlet_bc,
                            std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
                            std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs);

  template <typename RHS, typename DirichletBC>
  void setup_solver_weak(const mesh_type &in_mesh, const RHS &rhs,
                         const std::vector<std::string> &boundary_names,
                         const DirichletBC &dirichlet_bc,
                         std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
                         std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs);

  template <typename RHS, typename DirichletBC>
  void assemble_system_standard(const mesh_type &in_mesh, const RHS &rhs,
                                const std::vector<std::string> &boundary_names,
                                const DirichletBC &dirichlet_bc,
                                std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
                                std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs);

  template <typename RHS, typename DirichletBC>
  void assemble_system_Nitsche(const mesh_type &in_mesh, const RHS &rhs,
                               const std::vector<std::string> &boundary_names,
                               const DirichletBC &dirichlet_bc,
                               std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
                               std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs);

  template <typename RHS, typename DirichletBC>
  void assemble_system_weak(const mesh_type &in_mesh, const RHS &rhs,
                            const std::vector<std::string> &boundary_names,
                            const DirichletBC &dirichlet_bc,
                            std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
                            std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs);

  void prepare_boundary_metric_data(const mesh_type &in_mesh,
                                    const std::vector<std::string> &boundary_names);

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

  std::vector<bool> m_is_on_boundary;
  common::BlockArray<Uint, Uint> m_bdry_face_local_id;

  // Stores active index of all elements incident to weak boundary
  std::vector<Uint> m_bdry_elem_id;

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
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
HelmholtzSolverCG<MeshConfig>::HelmholtzSolverCG()
{
  m_quadrature_order = 0;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
HelmholtzSolverCG<MeshConfig>::~HelmholtzSolverCG()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS, typename DirichletBC>
void HelmholtzSolverCG<MeshConfig>::assemble(
    const std::string &method, const mesh_type &in_mesh, const RHS &rhs,
    const std::vector<std::string> &boundary_names, const DirichletBC &dirichlet_bc,
    std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs)
{
  if (method == "standard")

  {
    setup_solver_standard(in_mesh, rhs, boundary_names, dirichlet_bc, global_stiffness_matrix,
                          global_rhs);
    assemble_system_standard(in_mesh, rhs, boundary_names, dirichlet_bc, global_stiffness_matrix,
                             global_rhs);
  }
  else if (method == "Nitsche")
  {
    setup_solver_Nitsche(in_mesh, rhs, boundary_names, dirichlet_bc, global_stiffness_matrix,
                         global_rhs);
    assemble_system_Nitsche(in_mesh, rhs, boundary_names, dirichlet_bc, global_stiffness_matrix,
                            global_rhs);
  }
  else if (method == "weak")
  {
    setup_solver_weak(in_mesh, rhs, boundary_names, dirichlet_bc, global_stiffness_matrix,
                      global_rhs);

    assemble_system_weak(in_mesh, rhs, boundary_names, dirichlet_bc, global_stiffness_matrix,
                         global_rhs);
  }
  else
  {
    std::cerr << "Unknown method for system assembly specified. Aborting." << std::endl;
    return;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS, typename DirichletBC>
void HelmholtzSolverCG<MeshConfig>::setup_solver_standard(
    const mesh_type &in_mesh, const RHS &rhs, const std::vector<std::string> &boundary_names,
    const DirichletBC &dirichlet_bc,
    std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs)
{
  clock_t start, end;
  Real elapsed;

  typedef typename pdekit::result_of::dof_map_t<MeshConfig> cell_dofs_type;
  const cell_dofs_type &cell_dofs = *(in_mesh.dof_storage("geo_dofs"));

  const Uint nb_nodes_in_mesh = cell_dofs.nb_nodes();

  m_geo_cache_elems.clear();
  m_geo_metric_elems.clear();

  // Create a global vector which indicates which nodes are on boundary
  // Set the boundary condition
  m_is_on_boundary.resize(nb_nodes_in_mesh);
  m_is_on_boundary.assign(nb_nodes_in_mesh, false);

  const mesh::MeshBoundarySet<MeshConfig> &boundaries = in_mesh.all_boundaries();

  for (auto &name : boundary_names)
  {
    const typename mesh::MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr bdomain =
        boundaries.domain(name);
    for (Uint c = 0; c < (*bdomain).nb_active_cells(); ++c)
    {
      const mesh::MeshEntity bcell = (*bdomain).active_cell(cell_dofs, mesh::ActiveIdx(c));

      for (Uint v = 0; v < bcell.nb_vert(); ++v)
      {
        m_is_on_boundary[bcell.vertex(v)] = true;
      }
    }
  }

  /*
  const std::vector<std::shared_ptr<mesh::BoundaryFacets<MeshConfig,
  MeshConfig::TDIM-1>>> &boundary_domains = boundaries.all_domains();

  for (Uint dom = 0; dom < boundary_domains.size(); ++dom)
  {
    const mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM-1> &bdomain =
  (*boundary_domains[dom]);

    for (Uint c = 0; c < bdomain.nb_active_cells(); ++c)
    {
      const mesh::MeshEntity bcell = bdomain.active_cell(cell_dofs, c);

      for (Uint v = 0; v < bcell.nb_vert(); ++v)
      {
        m_is_on_boundary[bcell.vertex(v)] = true;
      }
    }
  }
  */

  // Create system matrix
  std::cout << "Allocating global system matrix" << std::endl;
  global_stiffness_matrix = std::make_shared<ls::TpetraCrsMatrix<Real>>(nb_nodes_in_mesh);
  std::cout << " ... done" << std::endl;
  global_rhs = std::make_shared<ls::TpetraMultiVector<Real>>(global_stiffness_matrix->map());

  std::vector<Real> values, boundary_values;
  std::vector<Int> indices, boundary_indices;

  m_quadrature_order = P1;

  start = clock();

  // Loop over all cells going type by type find the maximum quadrature order
  // The initialize the geometry metric
  for (const typename cell_dofs_type::const_dof_range_typed &cell_group :
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

  for (const typename cell_dofs_type::const_dof_range_typed &dof_group :
       cell_dofs.all_active_dof_groups())
  {
    const mesh::MeshEntity first_cell = (*dof_group.begin()).mesh_entity();
    // std::cout << "First cell in group = " << first_cell << std::endl;

    values.resize(first_cell.nb_vert());
    for (Uint v = 0; v < first_cell.nb_vert(); ++v)
    {
      values[v] = 0.0;
    }
    boundary_values.resize(1);
    boundary_values[0] = 0.0;

    indices.resize(first_cell.nb_vert());
    boundary_indices.resize(1);

    for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = dof_group.begin();
         cell_iter != dof_group.end(); ++cell_iter)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_view = cell_iter->tcell();
      const mesh::MeshEntity cell                         = cell_iter->mesh_entity();

      const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
          tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

      const ElemShape cell_shape = cell.pt_set_id().elem_shape();
      const Uint cell_order      = cell.pt_set_id().poly_order();

      const mesh::PointSetTagExt cell_tag_ext(cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS,
                                              0u);
      const mesh::sf::SFTag basis_tag         = basis_generator(cell_shape, cell_order);
      const mesh::PointSetTag eval_pt_set_tag = eval_pt_set_generator(cell_shape, cell_order);
      const mesh::PointSetTagExt eval_pt_set_tag_ext(eval_pt_set_tag, P0,
                                                     mesh::CellTransform::NO_TRANS, 0u);

      const mesh::DiscreteElemKey cell_key(cell_tag_ext, basis_tag, eval_pt_set_tag_ext);

      m_geo_cache_elems.push_back_to_buffer(cell_coords, cell_key);

      // std::cout << cell << std::endl;
      for (Uint vi = 0; vi < cell.nb_vert(); ++vi)
      {
        if (m_is_on_boundary[cell.vertex(vi)] == false)
        {
          for (Uint vj = 0; vj < cell.nb_vert(); ++vj)
          {
            indices[vj] = cell.vertex(vj);
          }
          global_stiffness_matrix->insert_values_in_row(cell.vertex(vi), values, indices);
          // (*global_rhs).insert_value(cell.vertex(vi), 0.0, 0);
        } // If this is not a boundary vertex
#if 0
          else
          {
            boundary_indices[0] = cell.vertex(vi);
            boundary_values[0] = 1.0;
            global_stiffness_matrix->insert_values_in_row(cell.vertex(vi), boundary_values,
                                                          boundary_indices);
            (*global_rhs)(cell.vertex(vi), 0) = 0.0;
          } // If this IS a boundary vertex
#endif
      }
    }
  } // Loop over all cell groups

  for (Uint v = 0; v < m_is_on_boundary.size(); ++v)
  {
    if (m_is_on_boundary[v])
    {
      boundary_indices[0] = v;
      boundary_values[0]  = 1.0;
      global_stiffness_matrix->insert_values_in_row(v, boundary_values, boundary_indices);
      //(*global_rhs).insert_value(v, 0.0, 0);
    }
  }

  // Prepare the right-hand side
  global_rhs->fill(0.0);

  //  global_stiffness_matrix->lock_structure();
  //  global_stiffness_matrix->print_structure_to_file("sparsity_Poisson.ps");

  global_stiffness_matrix->lock();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Sparse matrix setup: " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS, typename DirichletBC>
void HelmholtzSolverCG<MeshConfig>::setup_solver_Nitsche(
    const mesh_type &in_mesh, const RHS &rhs, const std::vector<std::string> &boundary_names,
    const DirichletBC &dirichlet_bc,
    std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs)
{
  clock_t start, end;
  Real elapsed;

  typedef typename pdekit::result_of::dof_map_t<MeshConfig> cell_dofs_type;
  const cell_dofs_type &cell_dofs = *(in_mesh.dof_storage("geo_dofs"));

  const Uint nb_nodes_in_mesh = cell_dofs.nb_nodes();

  // Create system matrix
  std::cout << "Allocating global system matrix" << std::endl;
  global_stiffness_matrix = std::make_shared<ls::TpetraCrsMatrix<Real>>(nb_nodes_in_mesh);
  std::cout << " ... done" << std::endl;
  global_rhs = std::make_shared<ls::TpetraMultiVector<Real>>(global_stiffness_matrix->map());

  m_geo_cache_elems.clear();
  m_geo_metric_elems.clear();

  std::vector<Real> values;
  std::vector<Int> indices;

  m_quadrature_order = P1;

  start = clock();

  // Loop over all cells going type by type find the maximum quadrature order
  // The initialize the geometry metric
  for (const typename cell_dofs_type::const_dof_range_typed &cell_group :
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

  for (const typename cell_dofs_type::const_dof_range_typed &dof_group :
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

    for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = dof_group.begin();
         cell_iter != dof_group.end(); ++cell_iter)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_view = cell_iter->tcell();
      const mesh::MeshEntity cell                         = cell_iter->mesh_entity();

      const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
          tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

      const ElemShape cell_shape = cell.pt_set_id().elem_shape();
      const Uint cell_order      = cell.pt_set_id().poly_order();
      const mesh::PointSetTagExt std_reg_tag_ext(cell.pt_set_id(), P0,
                                                 mesh::CellTransform::NO_TRANS, 0u);
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
        global_stiffness_matrix->insert_values_in_row(cell.vertex(vi), values, indices);
        (*global_rhs).insert_value(cell.vertex(vi), 0.0, 0);
      }
    }
  } // Loop over all cell groups

  // Prepare the right-hand side
  global_rhs->fill(0.0);

  //  global_stiffness_matrix->lock_structure();
  //  global_stiffness_matrix->print_structure_to_file("sparsity_Poisson.ps");

  global_stiffness_matrix->lock();

  // SETUP DATA FOR COMPUTATION OF BOUNDARY MATRICES
  std::vector<mesh::DiscreteElemKey> boundary_fe_type_map;
  std::vector<mesh::DiscreteElemKey> boundary_adjacent_vol_fe_type_map;

  const mesh::MeshBoundarySet<MeshConfig> &boundaries = in_mesh.all_boundaries();

  const std::vector<std::shared_ptr<mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>>>
      &boundary_domains = boundaries.all_domains();

  Uint nb_bdry_cells = 0;

  for (Uint dom = 0; dom < boundary_domains.size(); ++dom)
  {
    const mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> &bdomain =
        (*boundary_domains[dom]);

    for (Uint c = 0; c < bdomain.nb_active_cells(); ++c)
    {
      nb_bdry_cells++;

      const mesh::MeshEntity adjacent_volume_cell =
          bdomain.adjacent_active_vol_cell(cell_dofs, mesh::ActiveIdx(c));
      const mesh::PointSetTag adjacent_elem_tag = adjacent_volume_cell.pt_set_id();
      const mesh::PointSetTagExt adjacent_elem_tag_ext(adjacent_elem_tag, P0,
                                                       mesh::CellTransform::NO_TRANS,
                                                       bdomain.local_id(mesh::ActiveIdx(c)));

      const mesh::sf::SFTag adjacent_elem_sf_tag(adjacent_elem_tag.elem_shape(), SFunc::Lagrange,
                                                 adjacent_elem_tag.poly_order(), ModalBasis::Modal);
      const mesh::PointSetTag adjacent_elem_quad_tag(adjacent_elem_tag.elem_shape(),
                                                     2. * adjacent_elem_tag.poly_order(),
                                                     PointSetID::FaceGauss);
      const Uint local_id = bdomain.local_id(mesh::ActiveIdx(c));
      const mesh::PointSetTagExt adjacent_elem_quad_tag_ext(
          adjacent_elem_quad_tag, P0, mesh::CellTransform::RESTRICT_TO_CODIM_1, local_id);
      const mesh::DiscreteElemKey adjacent_elem_key(adjacent_elem_tag_ext, adjacent_elem_sf_tag,
                                                    adjacent_elem_quad_tag_ext);

      bool key_found = false;

      for (const auto &cell : boundary_adjacent_vol_fe_type_map)
      {
        if (cell == adjacent_elem_key)
        {
          key_found = true;
          break;
        }
      }

      if (!key_found)
      {
        boundary_adjacent_vol_fe_type_map.push_back(adjacent_elem_key);
      }

      const mesh::MeshEntity boundary_elem = bdomain.active_cell(cell_dofs, mesh::ActiveIdx(c));
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

      for (const auto &cell : boundary_fe_type_map)
      {
        if (boundary_elem_key == cell)
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
  for (const auto &key : boundary_fe_type_map)
  {
    std::cout << key << std::endl;
  }

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

  Uint item_nb = 1;

  for (Uint dom = 0; dom < boundary_domains.size(); ++dom)
  {
    const mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> &bdomain =
        (*boundary_domains[dom]);

    for (Uint c = 0; c < bdomain.nb_active_cells(); ++c)
    {
      const mesh::CellTopologyView<MeshConfig> adj_tcell_view =
          bdomain.adjacent_tcell(mesh::ActiveIdx(c));
      const mesh::MeshEntity adjacent_volume_cell =
          bdomain.adjacent_active_vol_cell(cell_dofs, mesh::ActiveIdx(c));
      const mesh::PointSetTag adjacent_elem_tag = adjacent_volume_cell.pt_set_id();
      const mesh::PointSetTagExt adjacent_elem_tag_ext(adjacent_elem_tag, P0,
                                                       mesh::CellTransform::NO_TRANS,
                                                       bdomain.local_id(mesh::ActiveIdx(c)));

      const mesh::sf::SFTag adjacent_elem_sf_tag(adjacent_elem_tag.elem_shape(), SFunc::Lagrange,
                                                 adjacent_elem_tag.poly_order(), ModalBasis::Modal);
      const mesh::PointSetTag adjacent_elem_quad_tag(adjacent_elem_tag.elem_shape(),
                                                     2. * adjacent_elem_tag.poly_order(),
                                                     PointSetID::FaceGauss);
      const Uint local_id = bdomain.local_id(mesh::ActiveIdx(c));
      const mesh::PointSetTagExt adjacent_elem_quad_tag_ext(
          adjacent_elem_quad_tag, P0, mesh::CellTransform::RESTRICT_TO_CODIM_1, local_id);
      const mesh::DiscreteElemKey adjacent_elem_key(adjacent_elem_tag_ext, adjacent_elem_sf_tag,
                                                    adjacent_elem_quad_tag_ext);

      const mesh::MeshEntity boundary_elem = bdomain.active_cell(cell_dofs, mesh::ActiveIdx(c));

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
          adj_tcell_view.pt_set_id(MeshConfig::TDIM - 1, local_id), boundary_elem.pt_set_id(),
          adj_tcell_view.coordinates(MeshConfig::TDIM - 1, local_id));

      m_geo_cache_faces_on_bdry.push_back_to_buffer(facet_coords, boundary_elem_key);

      item_nb++;
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

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Sparse matrix setup (Nitsche): " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS, typename DirichletBC>
void HelmholtzSolverCG<MeshConfig>::setup_solver_weak(
    const mesh_type &in_mesh, const RHS &rhs, const std::vector<std::string> &boundary_names,
    const DirichletBC &dirichlet_bc,
    std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs)
{
  clock_t start, end;
  Real elapsed;

  typedef typename pdekit::result_of::dof_map_t<MeshConfig> cell_dofs_type;
  const cell_dofs_type &cell_dofs = *(in_mesh.dof_storage("geo_dofs"));

  const Uint nb_nodes_in_mesh = cell_dofs.nb_nodes();

  m_geo_cache_elems.clear();
  m_geo_metric_elems.clear();

  // Setup metric containers for boundary terms
  prepare_boundary_metric_data(in_mesh, boundary_names);

  // Create system matrix
  std::cout << "Allocating global system matrix" << std::endl;
  global_stiffness_matrix = std::make_shared<ls::TpetraCrsMatrix<Real>>(nb_nodes_in_mesh);
  std::cout << " ... done" << std::endl;
  global_rhs = std::make_shared<ls::TpetraMultiVector<Real>>(global_stiffness_matrix->map());

  std::vector<Real> values, boundary_values;
  std::vector<Int> indices, boundary_indices;

  m_quadrature_order = P1;

  start = clock();

  // Loop over all cells going type by type find the maximum quadrature order
  // The initialize the geometry metric
  for (const typename cell_dofs_type::const_dof_range_typed &cell_group :
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

  for (const typename cell_dofs_type::const_dof_range_typed &dof_group :
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

    for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = dof_group.begin();
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
        global_stiffness_matrix->insert_values_in_row(cell.vertex(vi), values, indices);
        // (*global_rhs).insert_value(cell.vertex(vi), 0.0, 0);
      }
    }
  } // Loop over all cell groups

  // Prepare the right-hand side
  global_rhs->fill(0.0);

  //  global_stiffness_matrix->lock_structure();
  //  global_stiffness_matrix->print_structure_to_file("sparsity_Poisson.ps");

  global_stiffness_matrix->lock();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Sparse matrix setup (weak): " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS, typename DirichletBC>
void HelmholtzSolverCG<MeshConfig>::assemble_system_standard(
    const mesh_type &in_mesh, const RHS &rhs, const std::vector<std::string> &boundary_names,
    const DirichletBC &dirichlet_bc,
    std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs)
{
  // ------------------------------------------------------
  // PART II:
  // Perform the actual assembly
  // ------------------------------------------------------

  const Uint TDIM = MeshConfig::TDIM;

  clock_t start, end;
  Real elapsed;

  math::DenseDMat<Real> local_stiffness;
  math::DenseDVec<Real> local_rhs;

  math::DenseDMat<Real> local_stiffness_wD;
  math::DenseDMat<Real> M, invM;
  std::array<math::DenseDMat<Real>, MeshConfig::TDIM> D;
  std::array<math::DenseDMat<Real>, MeshConfig::TDIM> DT;

  std::vector<Real> values;
  std::vector<Int> indices;

  typedef typename pdekit::result_of::dof_map_t<MeshConfig> cell_dofs_type;
  typedef typename interpolation::GeometryMetric<
      MeshConfig::GDIM, MeshConfig::TDIM, MeshConfig::GDIM>::cellwise_metric cell_geo_metric_type;

  std::array<math::DenseConstMatView<Real>, TDIM> sf_deriv_phys;

  global_stiffness_matrix->unlock();
  // m_global_stiffness_matrix->fill(0.0);

  const cell_dofs_type &cell_dofs = *(in_mesh.dof_storage("geo_dofs"));

  // Initialize the clock in case there are no cell groups to loop over
  start = clock();

  for (const typename cell_dofs_type::const_dof_range_typed &cell_group :
       cell_dofs.all_active_dof_groups())
  {
    const mesh::MeshEntity first_cell = (*cell_group.begin()).mesh_entity();

    values.resize(first_cell.nb_vert());
    indices.resize(first_cell.nb_vert());

    // Resize local stiffness matrix for given element type
    const Uint n_loc_dof = first_cell.nb_vert();
    local_stiffness.resize(n_loc_dof, n_loc_dof);

    local_stiffness_wD.resize(n_loc_dof, n_loc_dof);
    M.resize(n_loc_dof, n_loc_dof);
    invM.resize(n_loc_dof, n_loc_dof);

    for (Uint d = 0; d < MeshConfig::TDIM; ++d)
    {
      D[d].resize(n_loc_dof, n_loc_dof);
      DT[d].resize(n_loc_dof, n_loc_dof);
    }

    // Resize local rhs for given element type
    local_rhs.resize(first_cell.nb_vert());

    // -----------------------------------------
    // SYSTEM ASSEMBLY - MAIN LOOP OVER ELEMENTS
    // -----------------------------------------

    start = clock();

    m_geo_metric_elems.evaluate(m_geo_cache_elems, interpolation::RebuildMetricIndex{true});

    Uint idx_in_metric = 0;

    for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = cell_group.begin();
         cell_iter != cell_group.end(); ++cell_iter)
    {
      const mesh::MeshEntity cell = cell_iter->mesh_entity();

      cell_geo_metric_type const cell_geo_metric =
          m_geo_metric_elems.cellwise_values(idx_in_metric);
      idx_in_metric++;

      math::DenseConstMatView<Real> qd_pts_phys       = cell_geo_metric.interpolated_coords();
      math::DenseConstVecView<Real> jacobian_in_qd_pt = cell_geo_metric.jdet();
      math::DenseDVec<Real> const &weight_in_qd_pt    = cell_geo_metric.pt_weights();

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

      local_stiffness_wD.fill(0.0);
      M.fill(0.0);
      invM.fill(0.0);

      for (Uint d = 0; d < MeshConfig::TDIM; ++d)
      {
        D[d].fill(0.0);
        DT[d].fill(0.0);
      }

      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        const Real wq = jacobian_in_qd_pt[q] * weight_in_qd_pt[q];
        for (Uint vi = 0; vi < cell.nb_vert(); ++vi)
        {
          for (Uint vj = 0; vj < cell.nb_vert(); ++vj)
          {
            M(vi, vj) += wq * V(q, vi) * V(q, vj);

            for (Uint d = 0; d < MeshConfig::TDIM; ++d)
            {
              D[d](vi, vj) = wq * V(q, vi) * sf_deriv_phys[d](q, vj);
            }
          } // Loop over vj
        }   // Loop over vi

      } // Loop over quadrature points

      for (Uint vi = 0; vi < cell.nb_vert(); ++vi)
      {
        for (Uint vj = 0; vj < cell.nb_vert(); ++vj)
        {
          for (Uint d = 0; d < MeshConfig::TDIM; ++d)
          {
            DT[d](vi, vj) = D[d](vj, vi);
          }
        }
      }

      M.inv(invM);

      for (Uint d = 0; d < MeshConfig::TDIM; ++d)
      {
        local_stiffness_wD += (DT[d] * invM) * D[d];
      }

      /*
      for (Uint r = 0; r < local_stiffness_wD.rows(); ++r)
      {
        for (Uint c = 0; c < local_stiffness_wD.cols(); ++c)
        {
          local_stiffness_wD(r, c) *= 3.0;
        }
      }
      */

      /*
      std::cout << "K:" << std::endl << local_stiffness << std::endl;
      std::cout << "KwD:" << std::endl << local_stiffness_wD << std::endl;
      std::cout << "V:" << std::endl << V << std::endl;
      std::cout << "M:" << std::endl << M << std::endl;
      std::cout << "invM:" << std::endl << invM << std::endl;
      // std::cout << "D0:" << std::endl << D[0] << std::endl;
      // std::cout << "D1:" << std::endl << D[1] << std::endl;
      std::cout <<
      "--------------------------------------------------------" <<
      std::endl;
      */

      // ------------------------------

      // std::cout << local_stiffness << std::endl;

      // Distribute to the global system

      for (Uint vi = 0; vi < cell.nb_vert(); ++vi)
      {
        if (m_is_on_boundary[cell.vertex(vi)] == false)
        {
          for (Uint vj = 0; vj < cell.nb_vert(); ++vj)
          {
            indices[vj] = cell.vertex(vj);
            values[vj]  = local_stiffness(vi, vj);
          }
          global_stiffness_matrix->add_values_to_row(cell.vertex(vi), values, indices);
          (*global_rhs).add_value(cell.vertex(vi), local_rhs[vi], 0);
        }
      }
    }

  } // Loop over all cell groups

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell                         = cell_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {

      if (m_is_on_boundary[cell.vertex(n)])
      {
        const math::DenseConstVecView<Real> node_coords = cell_coords.row_transpose(n);
        /*
        Real bc_val = 1.0;
        for (Uint d = 0; d < node_coords.size(); ++d)
        {
          bc_val += std::sin(math::pi * node_coords[d]);
        }
        */
        const Real bc_val = dirichlet_bc.eval(node_coords);
        (*global_rhs).insert_value(cell.vertex(n), bc_val, 0);
      }
    }
  }

  global_stiffness_matrix->lock();
  global_stiffness_matrix->print_structure_to_file("sparsity_Poisson.ps");

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Global system assembly: " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS, typename DirichletBC>
void HelmholtzSolverCG<MeshConfig>::assemble_system_Nitsche(
    const mesh_type &in_mesh, const RHS &rhs, const std::vector<std::string> &boundary_names,
    const DirichletBC &dirichlet_bc,
    std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs)
{
  // ------------------------------------------------------
  // PART II:
  // Perform the actual assembly
  // ------------------------------------------------------

  const Uint TDIM = MeshConfig::TDIM;

  clock_t start, end;
  Real elapsed;

  math::DenseDMat<Real> local_stiffness;
  math::DenseDVec<Real> local_rhs;

  std::vector<Real> values;
  std::vector<Int> indices;

  typedef typename pdekit::result_of::dof_map_t<MeshConfig> cell_dofs_type;
  typedef typename interpolation::GeometryMetric<
      MeshConfig::GDIM, MeshConfig::TDIM, MeshConfig::GDIM>::cellwise_metric cell_geo_metric_type;

  std::array<math::DenseConstMatView<Real>, TDIM> sf_deriv_phys;

  global_stiffness_matrix->unlock();
  // global_stiffness_matrix->fill(0.0);

  const cell_dofs_type &cell_dofs = *(in_mesh.dof_storage("geo_dofs"));

  for (const typename cell_dofs_type::const_dof_range_typed &cell_group :
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

    for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = cell_group.begin();
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
        global_stiffness_matrix->add_values_to_row(cell.vertex(vi), values, indices);
        (*global_rhs).add_value(cell.vertex(vi), local_rhs[vi], 0);
      }
    }
  } // Loop over all cell groups

  // -----------------------------------------
  // SYSTEM ASSEMBLY - LOOP OVER BOUNDARIES
  // -----------------------------------------

  m_geo_metric_elems_on_bdry.evaluate(m_geo_cache_elems_on_bdry,
                                      interpolation::RebuildMetricIndex{true});
  m_geo_metric_faces_on_bdry.evaluate(m_geo_cache_faces_on_bdry,
                                      interpolation::RebuildMetricIndex{true});

  const mesh::MeshBoundarySet<MeshConfig> &boundaries = in_mesh.all_boundaries();

  const std::vector<std::shared_ptr<mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>>>
      &boundary_domains = boundaries.all_domains();

  math::DenseSVec<Real, MeshConfig::GDIM> grad_phi_i_ref;
  math::DenseSVec<Real, MeshConfig::GDIM> grad_phi_i_phys;
  math::DenseSVec<Real, MeshConfig::GDIM> grad_phi_j_ref;
  math::DenseSVec<Real, MeshConfig::GDIM> grad_phi_j_phys;

  std::array<math::DenseConstMatView<Real>, MeshConfig::GDIM> dV_ref;

  Uint idx_in_bdry = 0;

  for (Uint dom = 0; dom < boundary_domains.size(); ++dom)
  {
    const mesh::BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1> &bdomain =
        (*boundary_domains[dom]);

    for (Uint c = 0; c < bdomain.nb_active_cells(); ++c)
    {
      const mesh::MeshEntity adjacent_volume_cell =
          bdomain.adjacent_active_vol_cell(cell_dofs, mesh::ActiveIdx(c));
      const mesh::MeshEntity boundary_elem = bdomain.active_cell(cell_dofs, mesh::ActiveIdx(c));

      local_stiffness.resize(adjacent_volume_cell.nb_vert(), adjacent_volume_cell.nb_vert());
      local_stiffness.fill(0.0);
      local_rhs.resize(adjacent_volume_cell.nb_vert());
      local_rhs.fill(0.0);

      const typename interpolation::GeometryMetric<
          MeshConfig::GDIM, MeshConfig::TDIM, MeshConfig::GDIM>::cellwise_metric adj_cell_geo_met =
          m_geo_metric_elems_on_bdry.cellwise_values(idx_in_bdry);

      const typename interpolation::GeometryMetric<
          MeshConfig::GDIM, MeshConfig::TDIM, MeshConfig::GDIM - 1>::cellwise_metric facet_geo_met =
          m_geo_metric_faces_on_bdry.cellwise_values(idx_in_bdry);

      const math::DenseConstMatView<Real> qd_pts_phys       = facet_geo_met.interpolated_coords();
      math::DenseConstVecView<Real> const jacobian_in_qd_pt = facet_geo_met.jdet();
      math::DenseDVec<Real> const &weight_in_qd_pt          = facet_geo_met.pt_weights();
      math::DenseConstMatView<Real> const normals_in_qd_pt  = facet_geo_met.normals();

      math::DenseDMat<Real> const &V_ref = adj_cell_geo_met.reference_sf_values();
      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        dV_ref[d] = adj_cell_geo_met.sf_derivatives(d);
      }

      const Real eta = 1.e6;

      for (Uint q = 0; q < facet_geo_met.nb_qd_pts(); ++q)
      {
        const Real wq = jacobian_in_qd_pt[q] * weight_in_qd_pt[q];

        const math::DenseConstMatView<Real> inv_jacobi_at_qd_pt = adj_cell_geo_met.inv_jacobi(q);

        for (Uint i = 0; i < adjacent_volume_cell.nb_vert(); ++i)
        {
          for (Uint d = 0; d < MeshConfig::GDIM; ++d)
          {
            grad_phi_i_ref[d] = dV_ref[d](q, i);
          }
          grad_phi_i_phys = inv_jacobi_at_qd_pt * grad_phi_i_ref;

          Real grad_phi_i_norm = 0.0;

          for (Uint d = 0; d < MeshConfig::GDIM; ++d)
          {
            grad_phi_i_norm += grad_phi_i_phys[d] * normals_in_qd_pt(q, d);
          }

          for (Uint j = 0; j < adjacent_volume_cell.nb_vert(); ++j)
          {
            for (Uint d = 0; d < MeshConfig::GDIM; ++d)
            {
              grad_phi_j_ref[d] = dV_ref[d](q, j);
            }
            grad_phi_j_phys = inv_jacobi_at_qd_pt * grad_phi_j_ref;

            Real grad_phi_j_norm = 0.0;

            for (Uint d = 0; d < MeshConfig::GDIM; ++d)
            {
              grad_phi_j_norm += grad_phi_j_phys[d] * normals_in_qd_pt(q, d);
            }

            local_stiffness(i, j) +=
                wq * (-V_ref(q, i) * grad_phi_j_norm - grad_phi_i_norm * V_ref(q, j) +
                      eta * V_ref(q, i) * V_ref(q, j));
          }

          const Real g = dirichlet_bc.eval(qd_pts_phys.row_transpose(q));
          local_rhs[i] += wq * (-g * grad_phi_i_norm + eta * g * V_ref(q, i));
        }
      } // Loop over quadrature points

      indices.resize(adjacent_volume_cell.nb_vert());
      values.resize(adjacent_volume_cell.nb_vert());

      for (Uint vi = 0; vi < adjacent_volume_cell.nb_vert(); ++vi)
      {
        for (Uint vj = 0; vj < adjacent_volume_cell.nb_vert(); ++vj)
        {
          indices[vj] = adjacent_volume_cell.vertex(vj);
          values[vj]  = local_stiffness(vi, vj);
        }
        global_stiffness_matrix->add_values_to_row(adjacent_volume_cell.vertex(vi), values,
                                                   indices);
        (*global_rhs).add_value(adjacent_volume_cell.vertex(vi), local_rhs[vi], 0);
      }

      idx_in_bdry++;
    }
  }

  global_stiffness_matrix->lock();
  // global_stiffness_matrix->print_structure_to_file("sparsity_Poisson.ps");

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Global system assembly (Nitsche): " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS, typename DirichletBC>
void HelmholtzSolverCG<MeshConfig>::assemble_system_weak(
    const mesh_type &in_mesh, const RHS &rhs, const std::vector<std::string> &boundary_names,
    const DirichletBC &dirichlet_bc,
    std::shared_ptr<ls::TpetraCrsMatrix<Real>> &global_stiffness_matrix,
    std::shared_ptr<ls::TpetraMultiVector<Real>> &global_rhs)
{
  const Uint TDIM = MeshConfig::TDIM;

  clock_t start, end;
  Real elapsed;

  math::DenseDMat<Real> local_stiffness;
  math::DenseDVec<Real> local_rhs;

  std::vector<Real> values;
  std::vector<Int> indices;

  using cell_dofs_type = typename pdekit::result_of::dof_map_t<MeshConfig>;

  std::array<math::DenseConstMatView<Real>, MeshConfig::GDIM> sf_deriv_phys;

  global_stiffness_matrix->unlock();
  // global_stiffness_matrix->fill(0.0);

  const cell_dofs_type &cell_dofs = *(in_mesh.dof_storage("geo_dofs"));

  for (const typename cell_dofs_type::const_dof_range_typed &cell_group :
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

    for (typename cell_dofs_type::const_dof_iterator_typed cell_iter = cell_group.begin();
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
        global_stiffness_matrix->add_values_to_row(cell.vertex(vi), values, indices);
        (*global_rhs).add_value(cell.vertex(vi), local_rhs[vi], 0);
      }
    }
  } // Loop over all cell groups

  // -----------------------------------------
  // SYSTEM ASSEMBLY - LOOP OVER BOUNDARIES
  // -----------------------------------------

  // First prepare volume data of cells adjacent to Dirichlet boundaries
  m_geo_cache_elems.flush();

  const auto basis_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  const auto eval_pt_set_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::PointSetTag(shape, m_quadrature_order, PointSetID::Gauss);
  };

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint ac : m_bdry_elem_id)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(ac));
    const mesh::MeshEntity cell = cell_dofs.active_cell(mesh::ActiveIdx(ac));

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
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

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

  for (Uint bdry_block = 0; bdry_block < m_bdry_face_local_id.nb_blocks(); ++bdry_block)
  {
    // SUB-STEP 1: get adjacent volume element and compute its mass and
    // derivative matrices
    const Uint elem_id = m_bdry_elem_id[bdry_block];
    const mesh::CellTopologyView<MeshConfig> adjacent_tcell_view =
        cell_dofs.tcell(mesh::ActiveIdx(elem_id));
    const mesh::MeshEntity adjacent_volume_cell = cell_dofs.active_cell(mesh::ActiveIdx(elem_id));

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
      global_stiffness_matrix->add_values_to_row(adjacent_volume_cell.vertex(vi), values, indices);
      (*global_rhs).add_value(adjacent_volume_cell.vertex(vi), b_weak[vi], 0);
    }
  } // Loop over Dirichlet boundary facet blocks

  global_stiffness_matrix->lock();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "Global system assembly (weak): " << elapsed << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void HelmholtzSolverCG<MeshConfig>::prepare_boundary_metric_data(
    const mesh_type &in_mesh, const std::vector<std::string> &boundary_names)
{
  typedef typename pdekit::result_of::dof_map_t<MeshConfig> cell_dofs_type;
  const cell_dofs_type &cell_dofs = *(in_mesh.dof_storage("geo_dofs"));

  Uint nb_active_cells_on_Dir_bdry = 0;
  Uint nb_Dir_bdry_faces           = 0;

  // Stores number of weak boundary faces of ANY ACTIVE elements
  // We will later on build a filtered vector which stores the weak
  // boundary faces only for ACTIVE elements that actually are incident
  // to weak boundaries
  std::vector<std::forward_list<Uint>> elem_to_weak_facets_map;
  elem_to_weak_facets_map.resize(cell_dofs.nb_active_cells());

  const mesh::MeshBoundarySet<MeshConfig> &boundaries = in_mesh.all_boundaries();

  for (auto &name : boundary_names)
  {
    const typename mesh::MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr bdomain =
        boundaries.domain(name);

    nb_Dir_bdry_faces += (*bdomain).nb_active_cells();

    for (Uint c = 0; c < (*bdomain).nb_active_cells(); ++c)
    {
      const mesh::MeshEntity cell =
          (*bdomain).adjacent_active_vol_cell(cell_dofs, mesh::ActiveIdx(c));
      const Uint local_id = (*bdomain).local_id(mesh::ActiveIdx(c));

      elem_to_weak_facets_map[cell.idx()].push_front(local_id);
    }
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

  m_bdry_elem_id.reserve(nb_active_cells_on_Dir_bdry);
  m_bdry_elem_id.resize(0);

  for (Uint i = 0; i < elem_to_weak_facets_map.size(); ++i)
  {
    if (!elem_to_weak_facets_map[i].empty())
    {
      Uint nb_local_faces = 0;
      for (auto item : elem_to_weak_facets_map[i])
      {
        bdry_face_local_id->push_back(item);
        nb_local_faces++;
      }
      bdry_faces_per_elem_offsets->push_back(nb_local_faces);
      m_bdry_elem_id.push_back(i);
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

  for (Uint bdry_block = 0; bdry_block < m_bdry_face_local_id.nb_blocks(); ++bdry_block)
  {
    const Uint elem_id                          = m_bdry_elem_id[bdry_block];
    const mesh::MeshEntity adjacent_volume_cell = cell_dofs.active_cell(mesh::ActiveIdx(elem_id));

    const common::ArrayView<const Uint, _1D, Uint> local_elem_faces =
        m_bdry_face_local_id.const_block(bdry_block);

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

  for (Uint bdry_block = 0; bdry_block < m_bdry_face_local_id.nb_blocks(); ++bdry_block)
  {
    const Uint elem_id = m_bdry_elem_id[bdry_block];

    const mesh::CellTopologyView<MeshConfig> adj_tcell_view =
        cell_dofs.tcell(mesh::ActiveIdx(elem_id));
    const mesh::MeshEntity adjacent_volume_cell = cell_dofs.active_cell(mesh::ActiveIdx(elem_id));

    const common::ArrayView<const Uint, _1D, Uint> local_elem_faces =
        m_bdry_face_local_id.const_block(bdry_block);

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
inline void HelmholtzSolverCG<MeshConfig>::build_adj_volume_matrices(
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
inline void HelmholtzSolverCG<MeshConfig>::build_E_matrices(
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
inline void HelmholtzSolverCG<MeshConfig>::build_F_matrices(
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

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
