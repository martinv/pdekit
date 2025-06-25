#ifndef PDEKIT_HDG_Helmholtz_Solver_HDG_hpp
#define PDEKIT_HDG_Helmholtz_Solver_HDG_hpp

#include <array>

#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/SolutionSpaceMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "linear_system/LSTpetra.hpp"
#include "mesh/CellBuffer.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/Tria.hpp"
#include "mesh/containers/DofMap.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "solver/fe/assembly/HDGCellWorker.hpp"
#include "solver/fe/assembly/HDGTraceWorker.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class HelmholtzSolverHDG
{
  public:
  /// TYPEDEFS
  using tria_t    = typename mesh::Tria<MeshConfig>;
  using dof_map_t = typename pdekit::result_of::dof_map_t<MeshConfig>;
  using f_space   = interpolation::FunctionSpace<MeshConfig>;

  /// Constructor
  HelmholtzSolverHDG();

  /// Copy constructor: deleted
  HelmholtzSolverHDG(const HelmholtzSolverHDG &rhs) = delete;

  /// Assignement operator: deleted
  HelmholtzSolverHDG &operator=(const HelmholtzSolverHDG &rhs) = delete;

  /// Destructor
  ~HelmholtzSolverHDG();

  /// Configure the solver
  void setup(const tria_t &tria, const mesh::DofMap<MeshConfig> &dofs);

  /// Set the solution function
  void set_solution(const interpolation::VectorMeshFunction<Real>::ptr &solution);

  /// Return the solution function
  interpolation::VectorMeshFunction<Real>::ptr solution() const;

  /// Fill the stiffness matrix and solve it
  void solve(const tria_t &input_mesh, const mesh::DofMap<MeshConfig> &sol_dofs);

  private:
  enum
  {
    ELEM_DIM  = MeshConfig::TDIM,
    FACET_DIM = MeshConfig::TDIM - 1,
    EDGE_DIM  = _1D,
  };

  /// TYPES AND TYPEDEFS

  struct ComputeQuadOrder
  {
    Uint operator()(const mesh::PointSetTag geo_tag, const mesh::PointSetTag sol_tag) const;
  };

  using fe_values_t = interpolation::FEValues;

  /// Geometry cache and metric
  using geo_cache_type  = interpolation::GeometryCache<MeshConfig::GDIM>;
  using geo_metric_type = interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>;

  /// Solution cache and metric
  using sol_cache_type  = interpolation::SolutionCache;
  using sol_metric_type = interpolation::SolutionSpaceMetric<MeshConfig, MeshConfig::TDIM>;

  using Vandermonde_mat_array = std::array<math::DenseDMat<Real>, MeshConfig::TDIM>;

  using fe_map_iterator = common::DataMap<mesh::PointSetTagExt, fe_values_t>::iterator;
  using sf_deriv_map_iterator =
      typename common::DataMap<mesh::PointSetTagExt, Vandermonde_mat_array>::iterator;

  /// METHODS

  /// Configure the solver for looping over edges/faces
  void setup_facets(const tria_t &input_mesh, const mesh::DofMap<MeshConfig> &geo_dofs,
                    const mesh::DofMap<MeshConfig> &sol_dofs);

  /// Given metric data for one cell and SOLUTION shape functions in reference
  /// space, compute SOLUTION shape function derivatives in physical space
  void compute_sf_deriv_phys(const typename geo_metric_type::cellwise_metric &geo_metric,
                             const fe_values_t &sol_fe_values,
                             Vandermonde_mat_array &sol_sf_der_phys);

  void interior_solve(const tria_t &input_mesh, const mesh::DofMap<MeshConfig> &geo_dofs,
                      const mesh::DofMap<MeshConfig> &sol_dofs,
                      const common::BlockArray<Real, Uint> &lambda_trace,
                      interpolation::VectorMeshFunction<Real> &interior_solution) const;

  /// PRIVATE DATA

  ComputeQuadOrder quad_order_ftor;

  std::vector<bool> m_is_on_boundary;

  std::shared_ptr<ls::TpetraMultiVector<Real>> m_rhs_b;

  std::shared_ptr<ls::TpetraMultiVector<Real>> m_solution_X;

  std::shared_ptr<ls::TpetraCrsMatrix<Real>> m_matrix_A;

  interpolation::VectorMeshFunction<Real> m_trace_solution;
  interpolation::VectorMeshFunction<Real>::ptr m_interior_solution;

  detail::HDGCellWorker<MeshConfig> m_cell_worker;
  detail::HDGTraceWorker<MeshConfig> m_facet_worker;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
HelmholtzSolverHDG<MeshConfig>::HelmholtzSolverHDG() : m_trace_solution("", "trace_solution")
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
HelmholtzSolverHDG<MeshConfig>::~HelmholtzSolverHDG()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void HelmholtzSolverHDG<MeshConfig>::setup(const tria_t &tria, const mesh::DofMap<MeshConfig> &dofs)
{
  m_cell_worker.setup_cells(tria, dofs, quad_order_ftor);
  m_facet_worker.setup_traces(tria, dofs, quad_order_ftor);

  const Uint nb_nodes       = dofs.nb_nodes();
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

  for (Uint c = 0; c < dofs.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = dofs.active_cell(mesh::ActiveIdx(c));
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
void HelmholtzSolverHDG<MeshConfig>::set_solution(
    const interpolation::VectorMeshFunction<Real>::ptr &solution)
{
  m_interior_solution = solution;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
interpolation::VectorMeshFunction<Real>::ptr HelmholtzSolverHDG<MeshConfig>::solution() const
{
  return m_interior_solution;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void HelmholtzSolverHDG<MeshConfig>::solve(const tria_t &input_mesh,
                                           const mesh::DofMap<MeshConfig> &sol_dofs)
{

#if 0
  // ----------------------------------------
  // Phase 1: compute element jacobians
  // using geometry cache and geometry metric
  // ----------------------------------------
  geo_cache_type geo_cache;
  geo_metric_type geo_metric;


  geo_cache.allocate(m_geo_cell_fe_map, geo_dofs.nb_active_cells());
  geo_metric.allocate_buffer(m_geo_cell_fe_map, geo_dofs.nb_active_cells());

  for (Uint ac = 0; ac < geo_dofs.nb_active_cells(); ++ac)
  {
    // Get information about geometry element ...
    const mesh::MeshEntity geo_cell = geo_dofs.active_cell(mesh::ActiveIdx(ac));
    const mesh::PointSetTag geo_cell_tag = geo_cell.std_region_id();
    const mesh::DofCoordinates<MeshConfig::GDIM> geo_cell_coord =
        geo_dofs.active_cell_coords(mesh::ActiveIdx(ac));
    const mesh::PointSetTagExt geo_cell_key(geo_cell_tag, geo_cell_tag.poly_order(),
                                             mesh::CellTransform::NO_TRANS, 0u);

    geo_cache.push_back_to_buffer(geo_cell_coord, geo_cell_key);
  }

  geo_metric.evaluate(geo_cache, interpolation::RebuildMetricIndex::True);

  // ----------------------------------------
  // Phase 2: loop over elements and for each
  // element, compute the sf. derivatives in
  // physical space
  // ----------------------------------------

  typedef typename geo_metric_type::cellwise_metric cellwise_metric;
  math::DenseDMat<Real> local_stiff_matrix;

  for (Uint ac = 0; ac < geo_dofs.nb_active_cells(); ++ac)
  {
    const mesh::MeshEntity geo_cell = geo_dofs.active_cell(mesh::ActiveIdx(ac));
    const mesh::MeshEntity sol_cell = sol_dofs.active_cell(mesh::ActiveIdx(ac));
    const mesh::PointSetTag sol_cell_tag = sol_cell.std_region_id();

    const Uint required_quad_order = compute_quad_order(geo_cell.std_region_id(), sol_cell_tag);

    const mesh::PointSetTagExt sol_cell_key(sol_cell_tag, required_quad_order,
                                             mesh::CellTransform::NO_TRANS, 0u);

    const cellwise_metric cell_geo_metric = geo_metric.cellwise_values(ac);
    local_stiff_matrix.resize(geo_cell.nb_vert(), geo_cell.nb_vert());

    const math::DenseDVec<Real> qd_weights = cell_geo_metric.pt_weights();
    const math::DenseConstVecView<Real> qd_jacobian = cell_geo_metric.jdet();
    local_stiff_matrix.fill(0.0);

    // Get the finite element values in reference space and compute
    // shape function derivatives in physical space
    common::PtrHandle<fe_values_t> sol_fe_values = m_sol_cell_fe_map.std_region_data(sol_cell_key);


    common::PtrHandle<Vandermonde_mat_array> fe_deriv_phys_ptr =
        m_sol_cell_sf_deriv_map.std_region_data(sol_cell_key);

    // Get a reference instead of using pointer - this is just for convenience ...
    Vandermonde_mat_array &sol_fe_deriv_phys = (*fe_deriv_phys_ptr);

    compute_sf_deriv_phys(cell_geo_metric, *sol_fe_values, sol_fe_deriv_phys);


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
            // fe_deriv_phys[dim] is Vandermonde matrix of size (nb. qd. pts, nb. shape functions)
            // that contains derivatives of shape functions with respect to dim (dim = x,y,z)
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
        values[0] = 1.0;
        m_matrix_A->insert_values_in_row(row_index, values, indices);
      }
      else
      {
        indices.resize(geo_cell.nb_vert());
        values.resize(geo_cell.nb_vert());
        for (Uint j = 0; j < geo_cell.nb_vert(); j++)
        {
          indices[j] = static_cast<Int>(geo_cell.vertex(j));
          values[j] = local_stiff_matrix(i, j);
        }
        m_matrix_A->add_values_to_row(row_index, values, indices);
      }
    }
  }

#endif

#if 0
  // HDG interior solve
  const Uint nb_trace_blocks = input_mesh.skeleton_size(MeshConfig::TDIM - 1);

  // Count the number of DOFs on skeleton
  Uint nb_node_entries = 0;

  for (Uint b = 0; b < nb_trace_blocks; ++b)
  {
    const mesh::TraceIncidences trace =
        input_mesh.active_skeleton_entry(FACET_DIM, mesh::ActiveIdx(b));
    const mesh::CellTopologyView<MeshConfig> tcell =
        input_mesh.cell(mesh::FlatIdx(trace.cell_id(0)));

    const mesh::ActiveIdx active_cell_id = tcell.active_idx();
    const mesh::MeshEntity sol_cell      = sol_dofs.active_cell(mesh::ActiveIdx(active_cell_id));
    const mesh::MeshEntity sol_facet     = sol_cell.sub_entity(FACET_DIM, trace.local_id(0));

    nb_node_entries += sol_facet.nb_vert();
  }

  // Create trace DOFs
  // mesh::CellBuffer<MeshConfig::GDIM, FACET_DIM> facet_buffer;
  mesh::CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> facet_buffer;

  std::cout << "Number of skeleton entries = " << nb_trace_blocks << std::endl;
  std::cout << "Number of DOFs on skeleton = " << nb_node_entries << std::endl;
  facet_buffer.reserve(nb_node_entries, nb_trace_blocks);

  std::unique_ptr<std::vector<Real>> lambda_values(new std::vector<Real>());
  std::unique_ptr<std::vector<Uint>> lambda_block_sizes(new std::vector<Uint>());

  lambda_block_sizes->resize(nb_trace_blocks);

  // FIXME - THIS IS JUST FOR TESTING
  // Fill the blocks with some values of lambda
  mesh::EntityDofRealign sol_facet_permutation;

  std::vector<Uint> facet_dofs;
  std::vector<Real> facet_coords_vec;

  std::vector<Real> all_facet_coords;
  all_facet_coords.reserve(MeshConfig::GDIM * nb_node_entries);
  all_facet_coords.resize(0);

  Uint curr_dof_id = 0;

  for (Uint b = 0; b < nb_trace_blocks; ++b)
  {
    const mesh::TraceIncidences trace =
        input_mesh.active_skeleton_entry(FACET_DIM, mesh::ActiveIdx(b));
    const mesh::EntityRealignCode pcode = trace.permutation(0).get().code();

    const mesh::CellTopologyView<MeshConfig> tcell =
        input_mesh.cell(mesh::FlatIdx(trace.cell_id(0)));

    const mesh::ActiveIdx active_cell_id = tcell.active_idx();
    const mesh::MeshEntity sol_cell      = sol_dofs.active_cell(mesh::ActiveIdx(active_cell_id));
    const mesh::MeshEntity sol_facet     = sol_cell.sub_entity(FACET_DIM, trace.local_id(0));

    sol_facet_permutation.change_type(sol_facet.pt_set_id(), trace.permutation(0).get().code());

    const mesh::DofCoordinates<MeshConfig::GDIM> facet_coords =
        sol_dofs.active_cell_coords(sol_facet);

    const Uint nb_dofs_in_facet = sol_facet.nb_vert();

    facet_dofs.resize(nb_dofs_in_facet);
    facet_coords_vec.resize(MeshConfig::GDIM * nb_dofs_in_facet);

    for (Uint i = 0; i < nb_dofs_in_facet; ++i)
    {
      // facet_dofs[i] = sol_facet.vertex(i);
      facet_dofs[i] = curr_dof_id++;

      math::DenseConstVecView<Real> node_coords = facet_coords.c(i);
      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        facet_coords_vec[MeshConfig::GDIM * i + d] = node_coords[d];
        all_facet_coords.push_back(node_coords[d]);
      }
    }
    facet_buffer.push_back_cell(b, sol_facet.pt_set_id(), facet_dofs, 0, facet_coords_vec);
  }

  mesh::DofMap<MeshConfig> trace_dofs(nullptr, "lambda_dofs");
  trace_dofs.init(FACET_DIM);
  trace_dofs.create_from_cells(facet_buffer, all_facet_coords);

  // Fill trace dof values - this is just for debugging
  m_trace_solution.resize(1, nb_node_entries);

  for (Uint t = 0; t < trace_dofs.nb_active_cells(); ++t)
  {
    const mesh::MeshEntity trace_ent = trace_dofs.active_cell(mesh::ActiveIdx(t));
    const mesh::DofCoordinates<MeshConfig::GDIM> trace_ent_coords =
        trace_dofs.active_cell_coords(mesh::ActiveIdx(t));
    for (Uint n = 0; n < trace_ent.nb_vert(); ++n)
    {
      const math::DenseConstVecView<Real> node_coords = trace_ent_coords.c(n);
      interpolation::VectorMeshFunction<Real>::entry_type trace_val =
          m_trace_solution.value(trace_ent.vertex(n));

      trace_val[0] = 1.0;

      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        trace_val[0] *= std::sin(math::pi * node_coords[d]);
      }
    }
  }

  m_interior_solution->resize(1, sol_dofs.nb_nodes());
  // m_interior_solution->fill(0.0);

  std::vector<Real> mass_mat_data;
  std::vector<Real> deriv_mat_data;
  std::vector<Real> E_mat_data;
  std::vector<Real> E_tilde_mat_data;
  std::vector<Real> F_mat_data;
  std::vector<Real> F_tilde_mat_data;

  // Derivative matrices
  math::DenseConstMatView<Real> M_mat;
  std::vector<math::DenseConstMatView<Real>> DT_mat(MeshConfig::GDIM);
  std::vector<math::DenseConstMatView<Real>> E_mat;
  std::vector<math::DenseConstMatView<Real>> E_tilde_mat;
  std::vector<math::DenseConstMatView<Real>> F_mat;
  std::vector<math::DenseConstMatView<Real>> F_tilde_mat;

  for (Uint ac = 0; ac < sol_dofs.nb_active_cells(); ++ac)
  {
    const mesh::CellTopologyView<MeshConfig> tcell = input_mesh.active_cell(mesh::ActiveIdx(ac));

    const mesh::ActiveIdx active_cell_id = tcell.active_idx();
    const mesh::MeshEntity sol_cell      = sol_dofs.active_cell(mesh::ActiveIdx(active_cell_id));
  }

#endif

  // mesh::send_to_output_stream(std::cout, trace_dofs) << std::endl;

  /*
  std::ofstream outfile;
  outfile.open("traces.dat");

  for (Uint t = 0; t < trace_dofs.nb_active_cells(); ++t)
  {
    const mesh::MeshEntity trace_ent =
  trace_dofs.active_cell(mesh::ActiveIdx(t)); const
  mesh::DofCoordinates<MeshConfig::GDIM> trace_ent_coords =
        trace_dofs.active_cell_coords(mesh::ActiveIdx(t));
    for (Uint n = 0; n < trace_ent.nb_vert(); ++n)
    {
      const math::DenseConstVecView<Real> node_coords = trace_ent_coords.c(n);
      outfile << trace_ent.vertex(n);
      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        outfile << " " << node_coords[d];
      }
      outfile << std::endl;
    }
  }
  outfile.close();
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void HelmholtzSolverHDG<MeshConfig>::compute_sf_deriv_phys(
    const typename geo_metric_type::cellwise_metric &geo_metric, const fe_values_t &sol_fe_values,
    Vandermonde_mat_array &sol_sf_der_phys)
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
inline Uint HelmholtzSolverHDG<MeshConfig>::ComputeQuadOrder::operator()(
    const mesh::PointSetTag geo_tag, const mesh::PointSetTag sol_tag) const
{
  // Compute minimum order for quadrature order
  // 2*(sol_order -1) because of the grad*grad operator
  return std::max(geo_tag.poly_order(), 2 * (sol_tag.poly_order() - 1));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void HelmholtzSolverHDG<MeshConfig>::interior_solve(
    const tria_t &input_mesh, const mesh::DofMap<MeshConfig> &geo_dofs,
    const mesh::DofMap<MeshConfig> &sol_dofs, const common::BlockArray<Real, Uint> &lambda_trace,
    interpolation::VectorMeshFunction<Real> &interior_solution) const
{
}

// ----------------------------------------------------------------------------

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
