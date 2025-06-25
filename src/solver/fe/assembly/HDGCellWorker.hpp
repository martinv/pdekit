#ifndef PDEKIT_Solver_HDG_Cell_Worker_hpp
#define PDEKIT_Solver_HDG_Cell_Worker_hpp

#include <memory>

#include "common/DataMap.hpp"
#include "interpolation/FEValues.hpp"
#include "interpolation/SolutionSpaceMetric.hpp"
#include "mesh/Tria.hpp"
#include "mesh/containers/DofMap.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "mesh/std_region/PointSetTagExt.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

namespace detail
{

template <typename MeshConfig>
class HDGCellWorker
{
  public:
  /// TYPEDEFS
  using tria_t    = typename mesh::Tria<MeshConfig>;
  using dof_map_t = typename pdekit::result_of::dof_map_t<MeshConfig>;

  /// Default constructor
  HDGCellWorker();

  /// Default destructor
  ~HDGCellWorker();

  /// Configure the solver for looping over cells
  template <typename QuadOrderRule>
  void setup_cells(const tria_t &tria, const mesh::DofMap<MeshConfig> &dofs,
                   const QuadOrderRule &quad_order_rule);

  template <typename QuadOrderRule>
  void compute_metric_data(const tria_t &tria, const mesh::DofMap<MeshConfig> &dofs,
                           const QuadOrderRule &quad_order_rule);

  template <typename QuadOrderRule>
  void interior_solve(const tria_t &tria, const mesh::DofMap<MeshConfig> &cell_dofs,
                      const mesh::DofMap<MeshConfig> &trace_dofs,
                      const interpolation::VectorMeshFunction<Real> &lambda_trace,
                      const QuadOrderRule &quad_order_rule);

  private:
  /// TYPES AND TYPEDEFS

  /// Type of function space
  using f_space_t = typename interpolation::FunctionSpace<MeshConfig>;

  using fe_values_t = interpolation::FEValues;

  using Vandermonde_mat_array = std::array<math::DenseDMat<Real>, MeshConfig::TDIM>;

  using fe_map_iterator =
      common::DataMap<mesh::PointSetTagExtPair, std::tuple<fe_values_t, fe_values_t>>::iterator;
  using sf_deriv_map_iterator =
      typename common::DataMap<mesh::PointSetTagExt, Vandermonde_mat_array>::iterator;

  /// Geometry cache and metric
  using geo_cache_type  = interpolation::GeometryCache<MeshConfig::GDIM>;
  using geo_metric_type = interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM>;
  /// Solution cache and metric
  using sol_cache_type  = interpolation::SolutionCache;
  using sol_metric_type = interpolation::SolutionSpaceMetric<MeshConfig, MeshConfig::GDIM>;

  enum Field
  {
    GEO = 0,
    SOL = 1
  };

  /// DATA
  // For each element type present in the mesh, this map holds
  // one 'FEValues', which has the Vandermonde matrix for shape functions
  // and their derivatives in REFERENCE space
  common::DataMap<mesh::PointSetTagExtPair, std::tuple<fe_values_t, fe_values_t>> m_cell_fe_map;

  // This map has 2 (in 2D) or 3 (in 3D) little matrices holding
  // shape function derivatives in physical space
  common::DataMap<mesh::PointSetTagExt, Vandermonde_mat_array> m_sol_cell_sf_deriv_map;

  typename f_space_t::ptr m_geo_cell_space;
  typename f_space_t::ptr m_sol_cell_space;

  /// Geometry cache
  geo_cache_type m_geo_cache;

  /// Geometry metric
  geo_metric_type m_geo_metric;

  /// Solution cache
  sol_cache_type m_sol_cache;

  /// Interpolated values and derivatives of the solution u_h
  sol_metric_type m_sol_metric;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
HDGCellWorker<MeshConfig>::HDGCellWorker()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
HDGCellWorker<MeshConfig>::~HDGCellWorker()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename QuadOrderRule>
void HDGCellWorker<MeshConfig>::setup_cells(const tria_t &tria,
                                            const mesh::DofMap<MeshConfig> &dofs,
                                            const QuadOrderRule &quad_order_rule)
{
  // m_geo_cell_space = std::make_shared<f_space_t>();
  // m_sol_cell_space = std::make_shared<f_space_t>();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename QuadOrderRule>
void HDGCellWorker<MeshConfig>::compute_metric_data(const tria_t &tria,
                                                    const mesh::DofMap<MeshConfig> &dofs,
                                                    const QuadOrderRule &quad_order_rule)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename QuadOrderRule>
void HDGCellWorker<MeshConfig>::interior_solve(
    const tria_t &tria, const mesh::DofMap<MeshConfig> &cell_dofs,
    const mesh::DofMap<MeshConfig> &trace_dofs,
    const interpolation::VectorMeshFunction<Real> &lambda_trace,
    const QuadOrderRule &quad_order_rule)
{
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
