#ifndef PDEKIT_DG_Elliptic_Solver_hpp
#define PDEKIT_DG_Elliptic_Solver_hpp

#include <array>

#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/SolutionSpaceMetric.hpp"
#include "linear_system/LSTpetra.hpp"
#include "mesh/Tria.hpp"
#include "mesh/containers/DofMap.hpp"
#include "mesh/point_set/StdPointSet.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class DGEllipticSolver
{
  public:
  /// TYPEDEFS

  using f_space = interpolation::FunctionSpace<MeshConfig>;

  /// Constructor
  DGEllipticSolver();

  /// Copy constructor: deleted
  DGEllipticSolver(const DGEllipticSolver &rhs) = delete;

  /// Assignement operator: deleted
  DGEllipticSolver &operator=(const DGEllipticSolver &rhs) = delete;

  /// Destructor
  ~DGEllipticSolver();

  /// Configure the solver
  void setup(const mesh::Tria<MeshConfig> &input_mesh, const mesh::DofMap<MeshConfig> &geo_dofs,
             const mesh::DofMap<MeshConfig> &sol_dofs);

  /// Fill the stiffness matrix and solve it
  void solve(const mesh::Tria<MeshConfig> &input_mesh, const mesh::DofMap<MeshConfig> &geo_dofs,
             const mesh::DofMap<MeshConfig> &sol_dofs);

  private:
  enum
  {
    ELEM_DIM  = MeshConfig::TDIM,
    FACET_DIM = MeshConfig::TDIM - 1,
    EDGE_DIM  = _1D,
  };

  /// TYPES AND TYPEDEFS

  /// Geometry cache and metric
  using geo_cache_type = interpolation::GeometryCache<MeshConfig::GDIM>;
  using geo_metric_type =
      interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, MeshConfig::TDIM>;

  /// Solution cache and metric
  using sol_cache_type  = interpolation::SolutionCache;
  using sol_metric_type = interpolation::SolutionSpaceMetric<MeshConfig, MeshConfig::TDIM>;

  using Vandermonde_mat_array = std::array<math::DenseDMat<Real>, MeshConfig::TDIM>;

  using fe_map_iterator = common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>::iterator;
  using sf_deriv_map_iterator =
      typename common::DataMap<mesh::PointSetTagExt, Vandermonde_mat_array>::iterator;

  /// METHODS
  /// Configure the solver for looping over cells
  void setup_cells(const mesh::Tria<MeshConfig> &input_mesh,
                   const mesh::DofMap<MeshConfig> &geo_dofs,
                   const mesh::DofMap<MeshConfig> &sol_dofs);

  /// Configure the solver for looping over edges/faces
  void setup_facets(const mesh::Tria<MeshConfig> &input_mesh,
                    const mesh::DofMap<MeshConfig> &geo_dofs,
                    const mesh::DofMap<MeshConfig> &sol_dofs);

  /// Given metric data for one cell and SOLUTION shape functions in reference
  /// space, compute SOLUTION shape function derivatives in physical space
  void compute_sf_deriv_phys(const typename geo_metric_type::cellwise_metric &geo_metric,
                             const interpolation::FEValues &sol_fe_values,
                             Vandermonde_mat_array &sol_sf_der_phys);

  Uint compute_quad_order(const mesh::PointSetTag geo_tag, const mesh::PointSetTag sol_tag) const;

  /// PRIVATE DATA
  std::vector<bool> m_is_on_boundary;

  std::shared_ptr<ls::TpetraMultiVector<Real>> m_rhs_b;

  std::shared_ptr<ls::TpetraMultiVector<Real>> m_solution_X;

  std::shared_ptr<ls::TpetraCrsMatrix<Real>> m_matrix_A;

  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> m_geo_cell_fe_map;
  // mesh::StdRegionDataMap<interpolation::FEValues> m_geo_facet_fe_map;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> m_geo_facet_L_fe_map;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> m_geo_facet_R_fe_map;

  // For each element type present in the mesh, this map holds
  // one 'FEValues', which has the Vandermonde matrix for shape functions
  // and their derivatives in REFERENCE space
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> m_sol_cell_fe_map;
  // mesh::StdRegionDataMap<interpolation::FEValues> m_sol_facet_fe_map;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> m_sol_facet_L_fe_map;
  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> m_sol_facet_R_fe_map;

  // This map has 2 (in 2D) or 3 (in 3D) little matrices holding
  // shape function derivatives in physical space
  common::DataMap<mesh::PointSetTagExt, Vandermonde_mat_array> m_sol_cell_sf_deriv_map;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
DGEllipticSolver<MeshConfig>::DGEllipticSolver()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
DGEllipticSolver<MeshConfig>::~DGEllipticSolver()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DGEllipticSolver<MeshConfig>::setup(const mesh::Tria<MeshConfig> &input_mesh,
                                         const mesh::DofMap<MeshConfig> &geo_dofs,
                                         const mesh::DofMap<MeshConfig> &sol_dofs)
{
  setup_cells(input_mesh, geo_dofs, sol_dofs);
  setup_facets(input_mesh, geo_dofs, sol_dofs);

  const Uint nb_nodes       = sol_dofs.nb_nodes();
  const Uint nb_rhs_columns = 1;

  m_matrix_A   = std::make_shared<ls::TpetraCrsMatrix<Real>>(nb_nodes);
  m_rhs_b      = std::make_shared<ls::TpetraMultiVector<Real>>(m_matrix_A->map(), nb_rhs_columns);
  m_solution_X = std::make_shared<ls::TpetraMultiVector<Real>>(m_matrix_A->map(), nb_rhs_columns);

  // Here we are "preparing/initializing" the matrix A:

  std::vector<Real> values;
  values.resize(1);
  values[0] = 0.0;
  std::vector<Int> indices;
  indices.resize(1);

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = sol_dofs.active_cell(mesh::ActiveIdx(c));
    // std::cout << "Cell [" << cell.idx() << "] = " << cell << std::endl;

    for (Uint i = 0; i < cell.nb_vert(); i++)
    {
      const Int row_index = static_cast<Int>(cell.vertex(i));
      for (Uint j = 0; j < cell.nb_vert(); j++)
      {
        indices[0] = static_cast<Int>(cell.vertex(j));
        m_matrix_A->insert_values_in_row(row_index, values, indices);
      }
    }
  }

  m_matrix_A->lock();

  /// Prepare Trilinos system matrix
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DGEllipticSolver<MeshConfig>::solve(const mesh::Tria<MeshConfig> &input_mesh,
                                         const mesh::DofMap<MeshConfig> &geo_dofs,
                                         const mesh::DofMap<MeshConfig> &sol_dofs)
{
  // ----------------------------------------
  // Phase 1: compute element jacobians
  // using geometry cache and geometry metric
  // ----------------------------------------
  geo_cache_type geo_cache;
  geo_metric_type geo_metric;

  geo_cache.allocate(m_geo_cell_fe_map, geo_dofs.nb_active_cells());
  geo_metric.allocate_buffer(m_geo_cell_fe_map, geo_dofs.nb_active_cells());

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint ac = 0; ac < geo_dofs.nb_active_cells(); ++ac)
  {
    // Get information about geometry element ...
    const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dofs.tcell(mesh::ActiveIdx(ac));
    const mesh::MeshEntity geo_cell                     = geo_dofs.active_cell(ac);
    const mesh::PointSetTag geo_cell_tag                = geo_cell.pt_set_id();

    const math::DenseConstMatView<Real> geo_cell_coord = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), geo_cell_tag, tcell_view.coordinates());

    const mesh::PointSetTagExt geo_cell_key(geo_cell_tag, geo_cell_tag.poly_order(),
                                            mesh::CellTransform::NO_TRANS, 0u);

    geo_cache.push_back_to_buffer(geo_cell_coord, geo_cell_key);
  }

  geo_metric.evaluate(geo_cache);

  // ----------------------------------------
  // Phase 2: loop over elements and for each
  // element, compute the sf. derivatives in
  // physical space
  // ----------------------------------------

  typedef typename geo_metric_type::cellwise_metric cellwise_metric;
  math::DenseDMat<Real> local_stiff_matrix;

  for (Uint ac = 0; ac < geo_dofs.nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity geo_cell      = geo_dofs.active_cell(ac);
    const mesh::MeshEntity sol_cell      = sol_dofs.active_cell(ac);
    const mesh::PointSetTag sol_cell_tag = sol_cell.pt_set_id();

    const Uint required_quad_order = compute_quad_order(geo_cell.pt_set_id(), sol_cell_tag);

    const mesh::PointSetTagExt sol_cell_key(sol_cell_tag, required_quad_order,
                                            mesh::CellTransform::NO_TRANS, 0u);

    cellwise_metric cell_geo_metric = geo_metric.cellwise_values(ac);
    local_stiff_matrix.resize(geo_cell.nb_vert(), geo_cell.nb_vert());

    const math::DenseDVec<Real> qd_weights          = cell_geo_metric.pt_weights();
    const math::DenseConstVecView<Real> qd_jacobian = cell_geo_metric.jdet();
    local_stiff_matrix.fill(0.0);

    // Get the finite element values in reference space and compute
    // shape function derivatives in physical space
    common::PtrHandle<interpolation::FEValues> sol_fe_values =
        m_sol_cell_fe_map.std_region_data(sol_cell_key);

    common::PtrHandle<Vandermonde_mat_array> fe_deriv_phys_ptr =
        m_sol_cell_sf_deriv_map.std_region_data(sol_cell_key);

    // Get a reference instead of using pointer - this is just for
    // convenience
    // ...
    Vandermonde_mat_array &sol_fe_deriv_phys = (*fe_deriv_phys_ptr);

    compute_sf_deriv_phys(cell_geo_metric, sol_fe_values, sol_fe_deriv_phys);

    // Use the SOLUTION shape function derivatives in physical space
    // to assemble the first part of the weak form (grad phi_i * grad phi_j)
    for (Uint q = 0; q < cell_geo_metric.nb_qd_pts(); q++)
    {
      for (Uint sf_i = 0; sf_i < cell_geo_metric.nb_dof_in_cell(); sf_i++)
      {
        for (Uint sf_j = 0; sf_j < cell_geo_metric.nb_dof_in_cell(); sf_j++)
        {
          Real grad_sum = 0.0;
          for (Uint dim = 0; dim < MeshConfig::TDIM; dim++)
          {
            // fe_deriv_phys[dim] is Vandermonde matrix of size (nb.
            // qd. pts, nb. shape functions) that contains
            // derivatives of shape functions with respect to dim
            // (dim = x,y,z)
            grad_sum += sol_fe_deriv_phys[dim](q, sf_i) * sol_fe_deriv_phys[dim](q, sf_j);
          }
          local_stiff_matrix(sf_i, sf_j) += grad_sum * qd_jacobian[q] * qd_weights[q];
        }
      }
    }

    std::vector<Real> values;
    values.resize(1);
    values[0] = 0.0;
    std::vector<Int> indices;
    indices.resize(1);

    for (Uint i = 0; i < geo_cell.nb_vert(); i++)
    {
      const Int row_index = static_cast<Int>(geo_cell.vertex(i));
      if (m_is_on_boundary[geo_cell.vertex(i)] == true)
      {
        indices.resize(1);
        values.resize(1);

        indices[0] = row_index;
        values[0]  = 1.0;
        m_matrix_A->insert_values_in_row(row_index, values, indices);
      }
      else
      {
        indices.resize(geo_cell.nb_vert());
        values.resize(geo_cell.nb_vert());
        for (Uint j = 0; j < geo_cell.nb_vert(); j++)
        {
          indices[j] = static_cast<Int>(geo_cell.vertex(j));
          values[j]  = local_stiff_matrix(i, j);
        }
        m_matrix_A->add_values_to_row(row_index, values, indices);
      }
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DGEllipticSolver<MeshConfig>::setup_cells(const mesh::Tria<MeshConfig> &input_mesh,
                                               const mesh::DofMap<MeshConfig> &geo_dofs,
                                               const mesh::DofMap<MeshConfig> &sol_dofs)
{
  // Clear map for geometry space reference elements
  m_geo_cell_fe_map.clear();

  // Clear map for solution space reference elements
  m_sol_cell_fe_map.clear();

  mesh::StdPointSet elem_quad;

  for (Uint ac = 0; ac < geo_dofs.nb_active_cells(); ++ac)
  {
    // Get information about geometry element ...
    const mesh::MeshEntity geo_cell      = geo_dofs.active_cell(mesh::ActiveIdx(ac));
    const mesh::PointSetTag geo_cell_tag = geo_cell.pt_set_id();

    const mesh::PointSetTagExt geo_cell_key(geo_cell_tag, geo_cell_tag.poly_order(),
                                            mesh::CellTransform::NO_TRANS, 0u);

    // ... and corresponding solution element

    const mesh::MeshEntity sol_cell      = sol_dofs.active_cell(mesh::ActiveIdx(ac));
    const mesh::PointSetTag sol_cell_tag = sol_cell.pt_set_id();

    // Compute minimum order for quadrature order
    // 2*(sol_order -1) because of the grad*grad operator
    const Uint required_quad_order = compute_quad_order(geo_cell_tag, sol_cell_tag);

    const mesh::PointSetTagExt sol_cell_key(sol_cell_tag, required_quad_order,
                                            mesh::CellTransform::NO_TRANS, 0u);

    // Create geometry and solution element values
    common::PtrHandle<interpolation::FEValues> geo_cell_fe =
        m_geo_cell_fe_map.std_region_data(geo_cell_key);

    if (geo_cell_fe.is_null())
    {
      geo_cell_fe = m_geo_cell_fe_map.create(geo_cell_key);
      const mesh::sf::SFTag geo_cell_sf_tag(geo_cell_tag.elem_shape(), SFunc::Lagrange,
                                            geo_cell_tag.poly_order(), ModalBasis::Modal);

      (*geo_cell_fe).configure(geo_cell_tag, geo_cell_sf_tag);
      const mesh::PointSetTag quad_tag =
          mesh::PointSetTag(geo_cell_tag.elem_shape(), required_quad_order, PointSetID::Gauss);
      elem_quad.change_type(quad_tag);
      (*geo_cell_fe).fill_Vandermonde(elem_quad.get().coordinates(), elem_quad.get().weights());
      (*geo_cell_fe).print();
    }

    common::PtrHandle<interpolation::FEValues> sol_cell_fe =
        m_sol_cell_fe_map.std_region_data(sol_cell_key);

    if (sol_cell_fe.is_null())
    {
      // Prepare shape function values in reference space
      sol_cell_fe = m_sol_cell_fe_map.create(sol_cell_key);
      const mesh::sf::SFTag sol_cell_sf_tag(sol_cell_tag.elem_shape(), SFunc::Lagrange,
                                            sol_cell_tag.poly_order(), ModalBasis::Modal);

      (*sol_cell_fe).configure(sol_cell_tag, sol_cell_sf_tag);
      const mesh::PointSetTag quad_tag =
          mesh::PointSetTag(sol_cell_tag.elem_shape(), required_quad_order, PointSetID::Gauss);
      elem_quad.change_type(quad_tag);
      (*sol_cell_fe).fill_Vandermonde(elem_quad.get().coordinates(), elem_quad.get().weights());
      (*sol_cell_fe).print();

      // Resize matrices that will hold shape function DERIVATIVES in
      // PHYSICAL SPACE
      common::PtrHandle<Vandermonde_mat_array> sol_cell_sf_deriv =
          m_sol_cell_sf_deriv_map.create(sol_cell_key);

      const Uint nb_qd_pts        = elem_quad.get().size();
      const Uint nb_nodes_in_elem = (*sol_cell_fe).nb_nodes();
      for (Uint dim = 0; dim < MeshConfig::TDIM; ++dim)
      {
        (*sol_cell_sf_deriv)[dim].resize(nb_qd_pts, nb_nodes_in_elem);
      }
    }
  }

  std::cout << "Print content of geo_cell_fe_map: " << std::endl;
  std::cout << "**********************************" << std::endl;
  for (fe_map_iterator iter = m_geo_cell_fe_map.begin(); iter != m_geo_cell_fe_map.end(); ++iter)
  {
    std::cout << "Type = " << iter.key_value() << std::endl;
    std::cout << "FE values: " << std::endl;
    (*iter.data_ptr()).print();
  }
  std::cout << "**********************************" << std::endl;
  std::cout << "Print content of sol_cell_fe_map: " << std::endl;
  std::cout << "**********************************" << std::endl;
  for (fe_map_iterator iter = m_sol_cell_fe_map.begin(); iter != m_sol_cell_fe_map.end(); ++iter)
  {
    std::cout << "Type = " << iter.key_value() << std::endl;
    std::cout << "FE values: " << std::endl;
    (*iter.data_ptr()).print();
  }
  std::cout << "**********************************" << std::endl;
  std::cout << "Print content of m_sol_cell_sf_deriv_map: " << std::endl;
  std::cout << "**********************************" << std::endl;
  for (sf_deriv_map_iterator iter = m_sol_cell_sf_deriv_map.begin();
       iter != m_sol_cell_sf_deriv_map.end(); ++iter)
  {
    std::cout << "Type = " << iter.key_value() << std::endl;
    std::cout << "FE values: " << std::endl;
    const Vandermonde_mat_array &deriv_array = (*iter.data_ptr());
    for (Uint d = 0; d < deriv_array.size(); ++d)
    {
      std::cout << deriv_array[d] << std::endl;
    }
  }
  std::cout << "**********************************" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DGEllipticSolver<MeshConfig>::setup_facets(const mesh::Tria<MeshConfig> &input_mesh,
                                                const mesh::DofMap<MeshConfig> &geo_dofs,
                                                const mesh::DofMap<MeshConfig> &sol_dofs)
{
  // Clear map for geometry space reference elements
  m_geo_facet_L_fe_map.clear();
  m_geo_facet_R_fe_map.clear();

  // Clear map for solution space reference elements
  m_sol_facet_L_fe_map.clear();
  m_sol_facet_R_fe_map.clear();

  mesh::StdPointSet facet_quad;

  for (Uint f = 0; f < input_mesh.active_skeleton_size(FACET_DIM); ++f)
  {
    const mesh::TraceIncidences facet =
        input_mesh.active_skeleton_entry(FACET_DIM, mesh::ActiveIdx(f));

    const Uint local_idx_L = 0;
    const Uint local_idx_R = (facet.size() == 2) ? 1 : 0;

    // Get the cell and facet on the left- and right-hand side of interface
    // in GEOMETRY space

    ///*** LEFT****
    const mesh::CellTopologyView<MeshConfig> tcell_L =
        input_mesh.cell(mesh::FlatIdx(facet.cell_id(local_idx_L)));
    const mesh::ActiveIdx active_cell_id_L = tcell_L.active_idx();
    const mesh::MeshEntity geo_cell_L      = geo_dofs.active_cell(active_cell_id_L);
    const mesh::PointSetTag geo_cell_L_tag = geo_cell_L.pt_set_id();
    const mesh::MeshEntity geo_facet_L =
        geo_cell_L.sub_entity(FACET_DIM, facet.local_id(local_idx_L));

    const mesh::EntityRealignCode pcode_L = facet.permutation(local_idx_L).get().code();
    const mesh::PointSetTagExt geo_facet_L_key(geo_facet_L.pt_set_id(), P0, pcode_L.adapt_op_id(),
                                               pcode_L.local_pos_in_parent());

    ///*** RIGHT ****
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        input_mesh.cell(mesh::FlatIdx(facet.cell_id(local_idx_R)));
    const mesh::ActiveIdx active_cell_id_R = tcell_R.active_idx();
    const mesh::MeshEntity geo_cell_R      = geo_dofs.active_cell(active_cell_id_R);
    const mesh::PointSetTag geo_cell_R_tag = geo_cell_R.pt_set_id();
    const mesh::MeshEntity geo_facet_R =
        geo_cell_R.sub_entity(FACET_DIM, facet.local_id(local_idx_R));

    const mesh::EntityRealignCode pcode_R = facet.permutation(local_idx_R).get().code();
    const mesh::PointSetTagExt geo_facet_R_key(geo_facet_R.pt_set_id(), P0, pcode_R.adapt_op_id(),
                                               pcode_R.local_pos_in_parent());

    // Get the cell and facet on the left- and right-hand side of interface
    // in SOLUTION space

    ///*** LEFT****
    const mesh::MeshEntity sol_cell_L      = sol_dofs.active_cell(active_cell_id_L);
    const mesh::PointSetTag sol_cell_L_tag = sol_cell_L.pt_set_id();
    const mesh::MeshEntity sol_facet_L =
        sol_cell_L.sub_entity(FACET_DIM, facet.local_id(local_idx_L));
    const mesh::PointSetTagExt sol_facet_L_key(sol_facet_L.pt_set_id(), P0, pcode_L.adapt_op_id(),
                                               pcode_L.local_pos_in_parent());

    ///*** RIGHT ****
    const mesh::MeshEntity sol_cell_R      = sol_dofs.active_cell(active_cell_id_R);
    const mesh::PointSetTag sol_cell_R_tag = sol_cell_R.pt_set_id();
    const mesh::MeshEntity sol_facet_R =
        sol_cell_R.sub_entity(FACET_DIM, facet.local_id(local_idx_R));
    const mesh::PointSetTagExt sol_facet_R_key(sol_facet_R.pt_set_id(), P0, pcode_R.adapt_op_id(),
                                               pcode_R.local_pos_in_parent());

    // Compute minimum order for quadrature order
    // 2*(sol_order -1) because of the grad*grad operator

    const Uint required_quad_order_L = compute_quad_order(geo_cell_L_tag, sol_cell_L_tag);
    const Uint required_quad_order_R = compute_quad_order(geo_cell_R_tag, sol_cell_R_tag);
    const Uint required_quad_order   = std::max(required_quad_order_L, required_quad_order_R);

    // Create geometry and solution element values for facet
    // Attach the Vandermonde matrix to the pointer geo_facet_L_fe
    // *** LEFT facet *****
    common::PtrHandle<interpolation::FEValues> geo_facet_L_fe =
        m_geo_facet_L_fe_map.std_region_data(geo_facet_L_key);
    if (geo_facet_L_fe.is_null())
    {
      geo_facet_L_fe = m_geo_facet_L_fe_map.create(geo_facet_L_key);
      const mesh::sf::SFTag geo_facet_sf_tag(geo_cell_L_tag.elem_shape(), SFunc::Lagrange,
                                             geo_cell_L_tag.poly_order(), ModalBasis::Modal);

      (*geo_facet_L_fe).configure(geo_cell_L_tag, geo_facet_sf_tag);
      const mesh::PointSetTag quad_tag =
          mesh::PointSetTag(geo_cell_L_tag.elem_shape(), required_quad_order, PointSetID::Gauss);
      facet_quad.change_type(quad_tag);
      (*geo_facet_L_fe)
          .fill_Vandermonde(facet_quad.get().coordinates(), facet_quad.get().weights());
      (*geo_facet_L_fe).print();
    }

    // *** RIGHT facet *****
    common::PtrHandle<interpolation::FEValues> geo_facet_R_fe =
        m_geo_facet_R_fe_map.std_region_data(geo_facet_R_key);
    if (geo_facet_R_fe.is_null())
    {
      geo_facet_R_fe = m_geo_facet_R_fe_map.create(geo_facet_R_key);
      const mesh::sf::SFTag geo_facet_sf_tag(geo_cell_R_tag.elem_shape(), SFunc::Lagrange,
                                             geo_cell_R_tag.poly_order(), ModalBasis::Modal);

      (*geo_facet_R_fe).configure(geo_cell_R_tag, geo_facet_sf_tag);
      const mesh::PointSetTag quad_tag =
          mesh::PointSetTag(geo_cell_R_tag.elem_shape(), required_quad_order, PointSetID::Gauss);
      facet_quad.change_type(quad_tag);
      (*geo_facet_R_fe)
          .fill_Vandermonde(facet_quad.get().coordinates(), facet_quad.get().weights());
      (*geo_facet_R_fe).print();
    }
    // Create geometry and solution element values for facet
    // Attach the Vandermonde matrix to the pointer sol_facet_L_fe
    // *** LEFT facet *****
    common::PtrHandle<interpolation::FEValues> sol_facet_L_fe =
        m_sol_facet_L_fe_map.std_region_data(sol_facet_L_key);
    if (sol_facet_L_fe.is_null())
    {
      sol_facet_L_fe = m_sol_facet_L_fe_map.create(sol_facet_L_key);
      const mesh::sf::SFTag sol_facet_sf_tag(sol_cell_L_tag.elem_shape(), SFunc::Lagrange,
                                             sol_cell_L_tag.poly_order(), ModalBasis::Modal);

      (*sol_facet_L_fe).configure(sol_cell_L_tag, sol_facet_sf_tag);
      const mesh::PointSetTag quad_tag =
          mesh::PointSetTag(sol_cell_L_tag.elem_shape(), required_quad_order, PointSetID::Gauss);
      facet_quad.change_type(quad_tag);
      (*sol_facet_L_fe)
          .fill_Vandermonde(facet_quad.get().coordinates(), facet_quad.get().weights());
      (*sol_facet_L_fe).print();
    }

    // Create geometry and solution element values for facet
    // Attach the Vandermonde matrix to the pointer sol_facet_R_fe
    // *** RIGHT facet *****
    common::PtrHandle<interpolation::FEValues> sol_facet_R_fe =
        m_sol_facet_R_fe_map.std_region_data(sol_facet_R_key);
    if (sol_facet_R_fe.is_null())
    {
      sol_facet_R_fe = m_sol_facet_R_fe_map.create(sol_facet_R_key);
      const mesh::sf::SFTag sol_facet_sf_tag(sol_cell_R_tag.elem_shape(), SFunc::Lagrange,
                                             sol_cell_R_tag.poly_order(), ModalBasis::Modal);

      (*sol_facet_R_fe).configure(sol_cell_R_tag, sol_facet_sf_tag);
      const mesh::PointSetTag quad_tag =
          mesh::PointSetTag(sol_cell_R_tag.elem_shape(), required_quad_order, PointSetID::Gauss);
      facet_quad.change_type(quad_tag);
      (*sol_facet_R_fe)
          .fill_Vandermonde(facet_quad.get().coordinates(), facet_quad.get().weights());
      (*sol_facet_R_fe).print();
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DGEllipticSolver<MeshConfig>::compute_sf_deriv_phys(
    const typename geo_metric_type::cellwise_metric &geo_metric,
    const interpolation::FEValues &sol_fe_values, Vandermonde_mat_array &sol_sf_der_phys)
{
  math::DenseSVec<Real, ELEM_DIM> der_sf_Phys;
  math::DenseSVec<Real, ELEM_DIM> der_sf_Ref;

  const Uint dof_in_elem = sol_fe_values.nb_nodes();

  for (Uint q = 0; q < geo_metric.nb_qd_pts(); q++)
  {
    for (Uint n = 0; n < dof_in_elem; ++n)
    {
      for (Uint dim = 0; dim < ELEM_DIM; ++dim)
      {
        const math::DenseDMat<Real> dV = sol_fe_values.deriv_Vandermonde(dim);
        der_sf_Ref[dim]                = dV(q, n);
      }
      const math::DenseConstMatView<Real> j_inv = geo_metric.inv_jacobi(q);
      der_sf_Phys                               = j_inv * der_sf_Ref;

      for (Uint dim = 0; dim < ELEM_DIM; ++dim)
      {
        sol_sf_der_phys[dim](q, n) = der_sf_Phys[dim];
      }
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline Uint DGEllipticSolver<MeshConfig>::compute_quad_order(const mesh::PointSetTag geo_tag,
                                                             const mesh::PointSetTag sol_tag) const
{
  // Compute minimum order for quadrature order
  // 2*(sol_order -1) because of the grad*grad operator
  return std::max(geo_tag.poly_order(), 2 * (sol_tag.poly_order() - 1));
}

// ----------------------------------------------------------------------------

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
