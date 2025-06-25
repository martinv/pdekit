#ifndef PDEKIT_Mesh_Adaptation_Mesh_Adapt_Step_hpp
#define PDEKIT_Mesh_Adaptation_Mesh_Adapt_Step_hpp

#include <algorithm>
#include <queue>
#include <vector>

#include "mesh/CellTransform.hpp"
#include "mesh/containers/CellPath.hpp"
#include "mesh/view/CellTopologyView.hpp"

namespace pdekit
{

namespace mesh
{

template <typename MeshConfig>
class Tria;

enum class AdaptationType : short unsigned
{
  h  = 0,
  p  = 1,
  hp = 2
};

enum class h_AdaptStrategy : short unsigned
{
  w_hanging_nodes = 0,
  coarsen         = 1,
  red_green       = 2
};

class MeshAdaptStep
{

  public:
  MeshAdaptStep();

  ~MeshAdaptStep();

  /// Return adaptation type
  AdaptationType adapt_type() const;

  /// Return h_adaptation strategy
  h_AdaptStrategy h_adapt_strategy() const;

  /// Define what h-type adaptation
  /// operation should be applied to each cell in mesh
  template <typename MeshConfig>
  void define_h_adapt_ops(Tria<MeshConfig> const &mesh, std::vector<CellTransform> const &adapt_op,
                          const h_AdaptStrategy strategy, const bool store_cell_paths = false);

  /// Define what h-type adaptation
  /// operation should be applied to each cell in mesh
  template <typename MeshConfig>
  void define_h_adapt_ops(Tria<MeshConfig> const &mesh, std::vector<CellTransform> const &adapt_op,
                          std::vector<Real> const &indicator, const Real ratio,
                          const h_AdaptStrategy strategy, const bool store_cell_paths = false);

  /// Define what p-type adaptation
  /// operation should be applied to each cell in mesh
  void define_p_adapt_ops(const std::vector<Uint> &p_order);

  /// Define what p-type adaptation
  /// operation should be applied to each cell in mesh
  void define_p_adapt_ops(const std::vector<Uint> &old_p_order,
                          const std::vector<Uint> &proposed_p_order,
                          const std::vector<Real> &quality_measure, const Real refinement_ratio);

  /// Return a vector of adaptation operations
  const std::vector<CellTransform> &adapt_ops() const;

  /// Return a vector of polynomial orders
  const std::vector<Uint> &cell_poly_orders() const;

  /// Return a value indicating whether this is a pure coarsening operation
  bool is_pure_coarsening() const;

  /// Return the number of cells to which coarsening will be applied
  Uint nb_coarsening_ops() const;

  /// Print operations
  void print_ops() const;

  private:
  /// Type of adaptation procedure:
  /// 1) h-adaptation
  /// 2) p-adaptation
  /// 3) hp-adaptation
  AdaptationType m_adapt_type;

  /// Type of h-adaptation strategy:
  h_AdaptStrategy m_h_strategy;

  /// Adaptation operation applied to each cell
  std::vector<CellTransform> m_cell_transform;

  /// Adaptation operation applied to each cell. Stores the path to the cell
  std::vector<std::tuple<CellTransform, CellPath>> m_cell_transform_w_path;

  /// Polynomial order of each cell (after adaptation)
  std::vector<Uint> m_cell_p_order;

  /// Additional flags with information about the character
  /// of the refinement operation
  Uint m_nb_coarsening_ops;
  Uint m_nb_do_nothing_ops;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshAdaptStep::define_h_adapt_ops(Tria<MeshConfig> const &mesh,
                                       std::vector<CellTransform> const &adapt_op,
                                       const h_AdaptStrategy strategy, const bool store_cell_paths)
{
  std::vector<Real> indicator;
  define_h_adapt_ops(mesh, adapt_op, indicator, 1.0, strategy, store_cell_paths);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshAdaptStep::define_h_adapt_ops(Tria<MeshConfig> const &mesh,
                                       std::vector<CellTransform> const &adapt_op,
                                       std::vector<Real> const &indicator, const Real ratio,
                                       const h_AdaptStrategy strategy, const bool store_cell_paths)
{
  m_adapt_type = AdaptationType::h;
  m_h_strategy = strategy;
  m_cell_transform.resize(mesh.nb_active_cells());

  if (ratio < 0.0 || ratio > 1.0)
  {
    m_cell_transform.assign(mesh.nb_active_cells(), CellTransform::NO_TRANS);
    std::cerr << "MeshAdaptSchedule:: Refinement ratio must have a value "
                 "between 0 and 1"
              << std::endl;
    return;
  }

  if (ratio < 1.0)
  {
    std::vector<std::pair<Real, Uint>> adapt_candidates;
    adapt_candidates.resize(indicator.size());

    for (Uint i = 0; i < adapt_candidates.size(); ++i)
    {
      adapt_candidates[i].first  = indicator[i];
      adapt_candidates[i].second = i;
    }

    std::sort(adapt_candidates.begin(), adapt_candidates.end());

    const Uint threshold = std::round(ratio * adapt_op.size());

    for (Uint i = 0; i < threshold; ++i)
    {
      m_cell_transform[i] = CellTransform::NO_TRANS;
    }
    for (Uint i = threshold; i < m_cell_transform.size(); ++i)
    {
      m_cell_transform[adapt_candidates[i].second] = adapt_op[adapt_candidates[i].second];
    }
  }
  else
  {
    std::copy(adapt_op.begin(), adapt_op.end(), m_cell_transform.begin());
  }

  // Let the mesh make sure that the adaptation field will generate a valid
  // mesh This eventually modifies some entries in 'm_adapt_op'
  mesh.verify_h_adapt_ops(m_h_strategy, m_cell_transform);

  m_nb_coarsening_ops = 0;
  m_nb_do_nothing_ops = 0;

  for (Uint i = 0; i < m_cell_transform.size(); ++i)
  {
    if (m_cell_transform[i] == CellTransform::COARSEN)
    {
      m_nb_coarsening_ops++;
    }
    else if (m_cell_transform[i] == CellTransform::NO_TRANS)
    {
      m_nb_do_nothing_ops++;
    }
  }

  if (store_cell_paths)
  {
    CellTopologyView<MeshConfig> zero_level_cell;
    std::vector<Uint> path_entries;

    m_cell_transform_w_path.reserve(m_cell_transform.size());
    m_cell_transform_w_path.resize(0);

    for (Uint ac = 0; ac < m_cell_transform.size(); ++ac)
    {
      const CellTopologyView<MeshConfig> curr_cell = mesh.active_cell(ActiveIdx(ac));
      mesh.path(curr_cell, zero_level_cell, path_entries);

      const CellPath cp(FlatIdx(zero_level_cell.linear_pos_idx()), path_entries);
      m_cell_transform_w_path.push_back(std::make_tuple(m_cell_transform[ac], cp));
    }
  }
  else
  {
    m_cell_transform_w_path.resize(0);
  }

  // Set the p-adaptation vector size to 0 - only h adaptation is considered
  m_cell_p_order.resize(0);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
