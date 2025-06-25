#ifndef PDEKIT_Mesh_Tria_hpp
#define PDEKIT_Mesh_Tria_hpp

#include <memory>
#include <string>

#include "common/Component.hpp"
#include "common/Meta.hpp"
#include "graph/Graph.hpp"
#include "math/DenseSVec.hpp"
#include "mesh/adaptation/LocalInterpolator.hpp"
#include "mesh/containers/DofMap.hpp"
#include "mesh/containers/MeshStatistics.hpp"
#include "mesh/containers/MeshTopologyAdaptAlgorithm.hpp"
#include "mesh/containers/TopologyAlgorithms.hpp"
#include "mesh/iterators/BdryDofIterator.hpp"
#include "mesh/iterators/CellTopologyIterator.hpp"
#include "mesh/view/TraceDofView.hpp"

namespace pdekit
{

namespace mesh
{

template <typename MeshConfig>
class Tria : public common::Component
{
  public:
  static const Uint GDIM = MeshConfig::GDIM;
  static const Uint TDIM = MeshConfig::TDIM;

  enum
  {
    CELL   = TDIM,
    FACET  = TDIM - 1,
    EDGE   = 1,
    VERTEX = 0
  };

  /// TYPEDEFS

  using shared_ptr       = std::shared_ptr<Tria>;
  using const_shared_ptr = std::shared_ptr<Tria const>;

  using dof_storage_type = typename result_of::dof_map_t<MeshConfig>;

  /// TYPEDEFS
  using mesh_config = MeshConfig;

  using const_cell_iterator =
      CellTopologyIterator<CellTopologyView<MeshConfig>, CellTopologyIterFilterDefault>;
  using const_active_cell_iterator =
      CellTopologyIterator<CellTopologyView<MeshConfig>, CellTopologyIterFilterActive>;

  using const_skeleton_iterator = typename internal::TriaFacets<MeshConfig>::const_iterator;
  // using const_bdry_topo_iterator = typename
  // MeshBoundarySet::const_bdry_topo_iterator;

  using const_trace_dof_iterator = BdryDofIterator<
      TraceDofView<dof_storage_type, typename internal::TriaFacets<MeshConfig>, ViewIsConst>>;

  /// METHODS

  /// Construct a mesh with given name
  explicit Tria(const std::string &name);

  /// Disabled copy constructor
  Tria(const Tria &other) = delete;

  /// Destroy the mesh object
  ~Tria() override;

  /// Return the topology dimension
  constexpr Uint topo_dim() const;

  /// Return the geometrical dimension
  constexpr Uint geo_dim() const;

  /// Disabled assignent operator
  Tria &operator=(const Tria &rhs) = delete;

  /// Return the type name of this class
  std::string derived_type_name() const override;

  /// Clear all internal data
  void clear();

  /// Return the mesh domains
  MeshBoundarySet<MeshConfig> &all_boundaries();

  /// Return the mesh domains, const version
  const MeshBoundarySet<MeshConfig> &all_boundaries() const;

  /// Return all boundaries wrapped in a pointer handle
  common::PtrHandle<MeshBoundarySet<MeshConfig>> all_boundaries_ptr();

  /// Return all boundaries wrapped in a pointer handle, const version
  common::PtrHandle<MeshBoundarySet<MeshConfig> const> all_boundaries_ptr() const;

  /// Return the first incidence dimension of this connectivity
  Uint dim() const;

  /// Return the number of cells
  Uint nb_active_cells() const;

  /// Get the number of entries in skeleton
  Uint skeleton_size(const Uint skeleton_dim = MeshConfig::TDIM - 1) const;

  /// Get the number of entries in skeleton
  Uint active_skeleton_size(const Uint skeleton_dim = MeshConfig::TDIM - 1) const;

  /// Return one topological cell (which is a proxy class)
  const CellTopologyView<MeshConfig> cell(const FlatIdx idx) const;

  /// Return one (proxy) cell from the container
  const CellTopologyView<MeshConfig> active_cell(const ActiveIdx idx) const;

  /// Return all active neighbours of one active cell
  /// The first entry in each tuple <TopologyCell, Uint, Uint> is an adjacent
  /// cell The second entry is a local id of a facet in that cell which
  /// touches 'tcell' The third entry is a local id of a facet in tcell which
  /// touches the neighbour
  const std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> active_neighbours(
      const CellTopologyView<MeshConfig> tcell) const;

  /// Generate the path to a cell
  /// @param start_tcell     ... topology cell whose path we're trying to
  /// generate
  /// @param zero_level_cell ... root cell of the generated path
  /// @param path_entries    ... child indices of the path, going from
  /// zero_level_cell
  ///                            to path_entries. The first entry of
  ///                            path_entries corresponds to local child id of
  ///                            cell on level 1
  void path(const CellTopologyView<MeshConfig> &start_tcell,
            CellTopologyView<MeshConfig> &zero_level_cell, std::vector<Uint> &path_entries) const;

  /// Get the leaf element of given path
  /// @param  path  ... path to the leaf
  /// @return topology cell which is the last element in the path
  const CellTopologyView<MeshConfig> path_leaf(const CellPath &path) const;

  /// Return one skeleton entry
  /// @param: skeleton_dim ... which skeleton to pick: facets or edges (in
  /// case
  ///                          MeshConfig::TDIM == 3D)
  /// @note:  In 2 dimensions, there is only one skeleton: edges
  /// @param: entry_idx ... index of skeleton entry to choose
  const TraceIncidences skeleton_entry(const Uint skeleton_dim, const FlatIdx entry_idx) const;

  /// Return one skeleton entry
  /// @param: skeleton_dim ... which skeleton to pick: facets or edges (in
  /// case
  ///                          MeshConfig::TDIM == 3D)
  /// @note:  In 2 dimensions, there is only one skeleton: edges
  /// @param: entry_idx ... index of skeleton entry to choose
  const TraceIncidences active_skeleton_entry(const Uint skeleton_dim,
                                              const ActiveIdx entry_idx) const;

  /// Fill the cell connectivity from cells stored in buffer
  /// This resets all previously allocated data
  void create_from_cells(const CellBuffer<GDIM, TDIM> &cell_buffer,
                         const bool build_edge_skeleton = false);

  /// Rebuild a skeleton. This is meant primarily to build edge skeleton in 3D
  /// mesh, when facet skeleton is already present by default
  void rebuild_skeleton(const Uint skeleton_dim,
                        std::unique_ptr<std::vector<IncidenceEntry>> &&incidences,
                        std::unique_ptr<std::vector<EntityDofRealign>> &&facet_permutations,
                        std::unique_ptr<std::vector<Uint>> &&facet_data_offsets);

  /// Return the number of all cells (active or not)
  Uint nb_all_cells_in_all_levels() const;

  /// Return the number of levels in the cells
  Uint nb_levels() const;

  /// Static method for copying of AdaptiveCells
  /// @todo: MAKE THIS METHOD PRIVATE AND ACCESSIBLE ONLY BY FRIENDS
  static void clone(Tria const &cells_in, Tria &cells_out);

  /// Change the type of each standard region
  /// @param rule ... a function that takes a PointSetTag as an argument
  ///                 and returns another PointSetTag
  ///                 The input is the existing standard region id
  ///                 The output is the new (desired) standard region tag
  template <typename StdRegChangeRule>
  void change_std_region_types(const StdRegChangeRule &rule);

  /// Perform a geometrical transformation of the mesh
  /// @param rule .. a void function which takes two parameters:
  ///                1) an input vector view (DenseConstVecView) representing
  ///                   the coordinates of one node
  ///                2) an output vector view (DenseConstVecView) representing
  ///                   the new (transformed) coordinates of the same node
  template <typename GeoTransformRule>
  void geo_transform(const GeoTransformRule &rule);

  /// Perform a geometrical transformation of one cell in the mesh
  /// @param idx  .. index of cell whose coordinates should be transformed
  /// @param rule .. a void function which takes two parameters:
  ///                1) an input vector view (DenseConstVecView) representing
  ///                   the coordinates of one node
  ///                2) an output vector view (DenseConstVecView) representing
  ///                   the new (transformed) coordinates of the same node
  template <typename GeoTransformRule>
  void geo_transform(const ActiveIdx idx, const GeoTransformRule &rule);

  /// Adapt the topology
  /// @todo: MAKE THIS METHOD PRIVATE AND ACCESSIBLE ONLY BY FRIENDS
  /// @param cell_ops ... vector of cell operations indicating what to do with
  /// each cell:
  ///                     refine, unrefine, do nothing, ...
  void adapt(MeshAdaptSequence<MeshConfig> &adapter);

  /// Verify that the prescribed h-adaptation operations will produce a valid
  /// mesh
  void verify_h_adapt_ops(const h_AdaptStrategy strategy,
                          std::vector<CellTransform> &adapt_ops) const;

  /// Find all elements incident to each (P1) node in the mesh
  void compute_node_to_cell_connectivity(
      common::BlockArray<std::tuple<Uint, Uint>, Uint> &connectivity,
      common::BlockArray<Uint, Uint> &cell_to_p1_node);

  /// Build mesh dual graph
  template <typename VertType>
  void build_dual_graph_undirected(common::BlockArray<VertType, Uint> &crs_graph) const;

  /// Create an iterator pointing to the beginning of the cell array, const
  /// version
  const_cell_iterator begin_cells() const;

  /// Return const iterator regardless of whether Cells is a constant object
  /// or not
  const_cell_iterator cbegin_cells() const;

  /// Create an iterator pointing 1 position after the end of the cell array
  const_cell_iterator end_cells() const;

  /// Return const iterator regardless of whether Cells is a constant object
  /// or not
  const_cell_iterator cend_cells() const;

  /// Create an iterator pointing to the first active cell, const
  /// version
  const_active_cell_iterator begin_cells_active() const;

  /// Return const iterator regardless of whether Cells is a constant object
  /// or not The iterator will point to the first active cell
  const_active_cell_iterator cbegin_cells_active() const;

  /// Create an iterator pointing 1 position after the end of the cell array
  const_active_cell_iterator end_cells_active() const;

  /// Return const iterator regardless of whether Cells is a constant object
  /// or not
  const_active_cell_iterator cend_cells_active() const;

  /// Return a range of all active cells
  common::IteratorRange<const_active_cell_iterator> as_active_cell_range() const;

  /// Return const iterator regardless of whether Cells is a constant
  /// object or not. Iterator points to the first skeleton entity, i.e.
  /// to the first edge or facet (depending on dim)
  const_skeleton_iterator cbegin_skeleton(const Uint dim) const;

  /// Return const iterator regardless of whether Cells is a constant
  /// object or not. Iterator points to the first position after
  /// the last skeleton entry (edge or facet, depending on dim)
  const_skeleton_iterator cend_skeleton(const Uint dim) const;

  /// Return const iterator over all degrees of freedom on skeleton of dimension dim
  const_trace_dof_iterator cbegin_skeleton_dofs(const Uint dim, const dof_storage_type &dofs) const;

  /// Return const iterator over all degrees of freedom on skeleton of dimension dim
  const_trace_dof_iterator cend_skeleton_dofs(const Uint dim, const dof_storage_type &dofs) const;

  /// Get the statistics of the mesh
  const MeshStatistics &statistics() const;

  /// Print cell storage information (debugging)
  void print_cell_data() const;

  /// Print all skeleton data
  void print_complete_skeleton(const Uint dim);

  /// Write the dual graph to a gmsh file using given DOF handler
  void write_dual_graph_to_gmsh(const std::string &dof_handler_name,
                                const std::string &outfilename) const;

  /// Dump cell coordinates to file (for debugging)
  void print_cell_coords_to_file(const std::string &filename) const;

  void color_cells_red_green(std::vector<CellTransform> &adapt_op, std::vector<Uint> &colors) const;

  void rebuild_active_skeleton();

  /// Get a pointer to dof handler with given name
  common::PtrHandle<dof_storage_type> dof_storage(const std::string &name);

  /// Get a pointer to dof handler with given name, const version
  common::PtrHandle<dof_storage_type const> dof_storage(const std::string &name) const;

  /// Create a new dof handler
  common::PtrHandle<dof_storage_type> create_dof_storage(const std::string &name);

  /// Print the names of all stored handlers (for debugging)
  void print_dof_handler_names() const;

  private:
  /// DATA

  /// TYPEDEFS
  using skeleton_data = internal::TriaFacets<MeshConfig>;

  // Find active cells adjacent to center cell. If necessary, descend in the
  // refinement hierarchy to reach the adjacent cells
  void add_active_adjacent_cells(const CellTopologyView<MeshConfig> center_cell,
                                 const TraceIncidences &facet_block,
                                 std::vector<CellTopologyView<MeshConfig>> &adjacent_cells) const;

  // Update the statistics describing this mesh. This should be called after
  // every mesh-modifying operation
  void update_statistics();

  /// Name of this container
  std::string m_name;

  /// Array holding the skeleton data - this will be edges in 2D
  /// (and the size of m_skeleton array should be 1), or edges and faces
  /// in 3D (in that case, the size of m_skeleton will be 2), because
  /// the first slot in m_skeleton will be used for edges, while the
  /// second slot will be used for faces
  std::array<internal::TriaFacets<MeshConfig>, common::StaticMax<MeshConfig::TDIM - 1, 1>::value>
      m_skeleton;

  /// Facets on level zero - these are stored upon loading the mesh
  /// and don't change during adaptation
  internal::TriaFacets<MeshConfig> m_zero_level_facets;

  /// Levels which contain layers of refined cells
  internal::TriaCells<MeshConfig> m_cells;

  /// Boundaries of the mesh
  MeshBoundarySet<MeshConfig> m_boundary_set;

  /// Next available cell index that can be used
  Uint m_next_avail_cell_id;

  /// Mesh statistics
  MeshStatistics m_stats;

  /// Dof handlers
  std::vector<std::unique_ptr<dof_storage_type>> m_dofs;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Tria<MeshConfig>::Tria(const std::string &name) : Component(name)
{
  for (Uint i = 0; i < m_skeleton.size(); ++i)
  {
    m_skeleton[i].set_dim(i + 1);
  }

  m_zero_level_facets.set_dim(MeshConfig::FACET_DIM);

  m_boundary_set.init(m_cells);

  m_next_avail_cell_id = 0u;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Tria<MeshConfig>::~Tria()
{
  clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
constexpr Uint Tria<MeshConfig>::topo_dim() const
{
  return MeshConfig::TDIM;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
constexpr Uint Tria<MeshConfig>::geo_dim() const
{
  return MeshConfig::GDIM;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
std::string Tria<MeshConfig>::derived_type_name() const
{
  return "Mesh";
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::clear()
{
  for (Uint i = 0; i < m_skeleton.size(); ++i)
  {
    m_skeleton[i].clear();
  }
  m_zero_level_facets.clear();

  m_cells.clear();

  m_next_avail_cell_id = 0;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
MeshBoundarySet<MeshConfig> &Tria<MeshConfig>::all_boundaries()
{
  return m_boundary_set;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const MeshBoundarySet<MeshConfig> &Tria<MeshConfig>::all_boundaries() const
{
  return m_boundary_set;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
common::PtrHandle<MeshBoundarySet<MeshConfig>> Tria<MeshConfig>::all_boundaries_ptr()
{
  return common::PtrHandle<MeshBoundarySet<MeshConfig>>(&m_boundary_set);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
common::PtrHandle<MeshBoundarySet<MeshConfig> const> Tria<MeshConfig>::all_boundaries_ptr() const
{
  return common::PtrHandle<MeshBoundarySet<MeshConfig> const>(&m_boundary_set);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline Uint Tria<MeshConfig>::dim() const
{
  return MeshConfig::TDIM;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint Tria<MeshConfig>::nb_active_cells() const
{
  return m_cells.nb_active_cells();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint Tria<MeshConfig>::skeleton_size(const Uint skeleton_dim) const
{
  return m_skeleton[skeleton_dim - 1].nb_all_facets();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint Tria<MeshConfig>::active_skeleton_size(const Uint skeleton_dim) const
{
  return m_skeleton[skeleton_dim - 1].nb_active_facets();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const CellTopologyView<MeshConfig> Tria<MeshConfig>::cell(const FlatIdx idx) const
{
  return m_cells.cell(idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const CellTopologyView<MeshConfig> Tria<MeshConfig>::active_cell(const ActiveIdx idx) const
{
  return m_cells.active_cell(idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const std::vector<std::tuple<CellTopologyView<MeshConfig>, Uint, Uint>> Tria<
    MeshConfig>::active_neighbours(const CellTopologyView<MeshConfig> tcell) const
{
  const skeleton_data &skeleton_facets = m_skeleton[MeshConfig::TDIM - 2];
  return internal::MeshTopologyAdaptAlgorithm<MeshConfig>::active_neighbours(
      m_cells, skeleton_facets, tcell);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::path(const CellTopologyView<MeshConfig> &start_tcell,
                            CellTopologyView<MeshConfig> &zero_level_cell,
                            std::vector<Uint> &path_entries) const
{
  m_cells.path(start_tcell, zero_level_cell, path_entries);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const CellTopologyView<MeshConfig> Tria<MeshConfig>::path_leaf(const CellPath &path) const
{
  return m_cells.path_leaf(path);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const TraceIncidences Tria<MeshConfig>::skeleton_entry(const Uint skeleton_dim,
                                                       const FlatIdx entry_idx) const
{
  const TraceIncidences block(m_skeleton[skeleton_dim - 1].facet_data(entry_idx));
  return block;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const TraceIncidences Tria<MeshConfig>::active_skeleton_entry(const Uint skeleton_dim,
                                                              const ActiveIdx entry_idx) const
{
  TraceIncidences block(m_skeleton[skeleton_dim - 1].active_facet_data(entry_idx));
  return block;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::create_from_cells(const CellBuffer<GDIM, TDIM> &cell_buffer,
                                         const bool build_edge_skeleton)
{
  clear();

  // Data arrays to hold raw facet information
  internal::TriaFacets<MeshConfig> &facets = m_skeleton.back();

  // Compute the facet skeleton
  facets.create_from_cell_buffer(cell_buffer);

  // Store facets on level 0
  internal::TriaFacets<MeshConfig>::clone(facets, m_zero_level_facets);

  if ((MeshConfig::TDIM == 3) && build_edge_skeleton)
  {
    internal::TriaFacets<MeshConfig> &edges = m_skeleton.front();

    // Compute the edge skeleton
    edges.create_from_cell_buffer(cell_buffer);
  }

  std::vector<StdRegion> cell_types(cell_buffer.nb_active_cells());
  std::vector<Uint> cell_face(facets.nb_stored_entries(), 0u);
  std::vector<Uint> cell_face_offsets(cell_buffer.nb_active_cells() + 1, 0u);

  // Build the cell index
  for (Uint c = 0; c < cell_buffer.nb_active_cells(); ++c)
  {
    const MeshEntity cell = cell_buffer.active_cell(ActiveIdx(c));

    cell_types[c].change_type(cell.pt_set_id());
    // Here, cell_face_offsets[c+1] stores the number of faces in cell [c]
    cell_face_offsets[c + 1] = cell.nb_sub_elements(MeshConfig::TDIM - 1);
  }

  for (Uint c = 1; c < cell_face_offsets.size(); ++c)
  {
    cell_face_offsets[c] += cell_face_offsets[c - 1];
  }

  for (Uint f = 0; f < facets.nb_active_facets(); ++f)
  {
    const TraceIncidences facet = facets.active_facet_data(ActiveIdx(f));
    for (Uint idx_in_facet = 0; idx_in_facet < facet.size(); ++idx_in_facet)
    {
      const Uint cell_idx                               = facet.cell_id(idx_in_facet);
      const SUint local_id                              = facet.local_id(idx_in_facet);
      cell_face[cell_face_offsets[cell_idx] + local_id] = f;
    }
  }

  m_cells.emplace_cells(cell_types, cell_face, cell_buffer.cell_coordinates());

  m_next_avail_cell_id = cell_buffer.nb_active_cells();

  update_statistics();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::rebuild_skeleton(
    const Uint skeleton_dim, std::unique_ptr<std::vector<IncidenceEntry>> &&incidences,
    std::unique_ptr<std::vector<EntityDofRealign>> &&facet_permutations,
    std::unique_ptr<std::vector<Uint>> &&facet_data_offsets)
{
  if ((skeleton_dim < _1D) || (skeleton_dim > _2D))
  {
    std::cerr << "Mesh::rebuild_skeleton: can only build skeleton with "
                 "skeleton_dim = 1 or 2"
              << std::endl;
    return;
  }

  skeleton_data &skeleton = m_skeleton[skeleton_dim - 1];
  skeleton.create_from_raw_data(std::move(incidences), std::move(facet_permutations),
                                std::move(facet_data_offsets));

  if (skeleton_dim == MeshConfig::FACET_DIM)
  {
    m_zero_level_facets.clear();
    internal::TriaFacets<MeshConfig>::clone(skeleton, m_zero_level_facets);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint Tria<MeshConfig>::nb_all_cells_in_all_levels() const
{
  return m_cells.nb_all_cells();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint Tria<MeshConfig>::nb_levels() const
{
  return m_cells.nb_all_levels();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::clone(Tria const &mesh_in, Tria &mesh_out)
{
  // 1) Copy/clone the skeleton: facets and edges
  for (Uint i = 0; i < mesh_in.m_skeleton.size(); ++i)
  {
    internal::TriaFacets<MeshConfig>::clone(mesh_in.m_skeleton[i], mesh_out.m_skeleton[i]);
  }

  internal::TriaFacets<MeshConfig>::clone(mesh_in.m_zero_level_facets,
                                          mesh_out.m_zero_level_facets);

  // 2) Copy/clone all cells
  mesh_out.m_cells.clear();
  internal::TriaCells<MeshConfig>::clone(mesh_in.m_cells, mesh_out.m_cells);

  // 3) Update the next available cell id
  mesh_out.m_next_avail_cell_id = mesh_in.m_next_avail_cell_id;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename StdRegChangeRule>
void Tria<MeshConfig>::change_std_region_types(const StdRegChangeRule &rule)
{
  m_cells.change_std_region_types(rule);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename GeoTransformRule>
void Tria<MeshConfig>::geo_transform(const GeoTransformRule &rule)
{
  m_cells.geo_transform(rule);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename GeoTransformRule>
void Tria<MeshConfig>::geo_transform(const ActiveIdx idx, const GeoTransformRule &rule)
{
  m_cells.geo_transform(idx, rule);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::adapt(MeshAdaptSequence<MeshConfig> &adapter)
{
  if (adapter.adapt_type() != AdaptationType::h)
  {
    std::cerr << "Mesh::adapt: only h-adaptation is allowed for mesh topology" << std::endl;
  }

  std::vector<CellTransform> const &prescribed_adapt_ops = adapter.adapt_ops();

  internal::TriaFacets<MeshConfig> &skeleton_facets = m_skeleton[MeshConfig::TDIM - 2];

  // Check whether only removal operations were scheduled
  if (adapter.is_pure_coarsening())
  {
    internal::MeshTopologyAdaptAlgorithm<MeshConfig>::remove_cells(
        m_cells, m_zero_level_facets, skeleton_facets, m_boundary_set, prescribed_adapt_ops);
    return;
  }

  if ((adapter.nb_coarsening_ops() > 0) && !adapter.is_pure_coarsening())
  {
    std::cerr << "Mesh: adapt: can't coarsen and refine at the same time.\n"
              << "Please specify only refinement or only coarsening operations "
                 "at a time."
              << std::endl;
    return;
  }

  if (adapter.h_adapt_strategy() == h_AdaptStrategy::w_hanging_nodes)
  {
    internal::MeshTopologyAdaptAlgorithm<MeshConfig>::adapt_w_hanging_nodes(
        m_cells, skeleton_facets, m_boundary_set, prescribed_adapt_ops);
  }
  else if (adapter.h_adapt_strategy() == h_AdaptStrategy::red_green)
  {
    internal::MeshTopologyAdaptAlgorithm<MeshConfig>::adapt_red_green(
        m_cells, m_zero_level_facets, skeleton_facets, m_boundary_set, prescribed_adapt_ops);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::verify_h_adapt_ops(const h_AdaptStrategy strategy,
                                          std::vector<CellTransform> &adapt_ops) const
{
  internal::TriaFacets<MeshConfig> const &skeleton_facets = m_skeleton[MeshConfig::TDIM - 2];
  if (strategy == h_AdaptStrategy::coarsen)
  {
    internal::MeshTopologyAdaptAlgorithm<MeshConfig>::prepare_for_cell_removal(
        m_cells, skeleton_facets, adapt_ops);
  }
  else if (strategy == h_AdaptStrategy::w_hanging_nodes)
  {
    internal::MeshTopologyAdaptAlgorithm<MeshConfig>::prepare_cells_for_hanging_node_adapt(
        m_cells, skeleton_facets, adapt_ops);
  }
  else if (strategy == h_AdaptStrategy::red_green)
  {
    /*
    std::vector<Uint> colors;
    internal::MeshTopologyAdaptAlgorithm<MeshConfig>::prepare_cells_for_red_green_adapt(
        m_cells, skeleton_facets, adapt_ops, colors);
    */

    // internal::MeshTopologyAdaptAlgorithm<MeshConfig>::prepare_cells_for_hanging_node_adapt(
    //    m_cells, skeleton_facets, adapt_ops);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::compute_node_to_cell_connectivity(
    common::BlockArray<std::tuple<Uint, Uint>, Uint> &connectivity,
    common::BlockArray<Uint, Uint> &cell_to_p1_node)
{
  internal::TriaFacets<MeshConfig> &skeleton_facets = m_skeleton[MeshConfig::TDIM - 2];
  TopologyAlgorithms::compute_node_to_cell_connectivity(m_cells, skeleton_facets, connectivity,
                                                        cell_to_p1_node);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename VertType>
void Tria<MeshConfig>::build_dual_graph_undirected(
    common::BlockArray<VertType, Uint> &crs_graph) const
{
  internal::TriaFacets<MeshConfig> const &skeleton_facets = m_skeleton[MeshConfig::TDIM - 2];
  graph::Graph<VertType> graph(m_cells.nb_active_cells());

  for (Uint f = 0; f < skeleton_facets.nb_active_facets(); ++f)
  {
    const mesh::TraceIncidences facet_block = skeleton_facets.active_facet_data(ActiveIdx(f));
    if (facet_block.size() == 2)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_L =
          m_cells.cell(FlatIdx(facet_block.cell_id(0)));
      const mesh::CellTopologyView<MeshConfig> tcell_R =
          m_cells.cell(FlatIdx(facet_block.cell_id(1)));

      /// @TODO: for active facets, this is not necessary?
      if ((tcell_L.status() == EntityStatus::Active) && (tcell_R.status() == EntityStatus::Active))
      {
        graph.insert_edge(tcell_L.active_idx().id(), tcell_R.active_idx().id());
        graph.insert_edge(tcell_R.active_idx().id(), tcell_L.active_idx().id());
      }
    } // if
  }

  graph.compress_to_crs(crs_graph);
  crs_graph.sort_blocks();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_cell_iterator Tria<MeshConfig>::begin_cells() const
{
  return const_cell_iterator(m_cells, 0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_cell_iterator Tria<MeshConfig>::cbegin_cells() const
{
  return const_cell_iterator(m_cells, 0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_cell_iterator Tria<MeshConfig>::end_cells() const
{
  return const_cell_iterator(m_cells, m_cells.nb_all_cells());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_cell_iterator Tria<MeshConfig>::cend_cells() const
{
  return const_cell_iterator(m_cells, m_cells.nb_all_cells());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_active_cell_iterator Tria<MeshConfig>::begin_cells_active() const
{
  return const_active_cell_iterator(m_cells, m_cells.first_active_cell_pos());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_active_cell_iterator Tria<MeshConfig>::cbegin_cells_active() const
{
  return const_active_cell_iterator(m_cells, m_cells.first_active_cell_pos());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_active_cell_iterator Tria<MeshConfig>::end_cells_active() const
{
  return const_active_cell_iterator(m_cells, m_cells.nb_all_cells());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_active_cell_iterator Tria<MeshConfig>::cend_cells_active() const
{
  return const_active_cell_iterator(m_cells, m_cells.nb_all_cells());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
common::IteratorRange<typename Tria<MeshConfig>::const_active_cell_iterator> Tria<
    MeshConfig>::as_active_cell_range() const
{
  return common::IteratorRange<const_active_cell_iterator>(cbegin_cells_active(),
                                                           cend_cells_active());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_skeleton_iterator Tria<MeshConfig>::cbegin_skeleton(
    const Uint dim) const
{
  const_skeleton_iterator it(m_skeleton[dim - 1], 0);
  return it;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_skeleton_iterator Tria<MeshConfig>::cend_skeleton(
    const Uint dim) const
{
  const_skeleton_iterator it(m_skeleton[dim - 1], m_skeleton[dim - 1].nb_all_facets());
  return it;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_trace_dof_iterator Tria<MeshConfig>::cbegin_skeleton_dofs(
    const Uint dim, const dof_storage_type &dofs) const
{
  const_trace_dof_iterator it(dofs, m_skeleton[dim - 1], ActiveIdx(0));
  return it;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename Tria<MeshConfig>::const_trace_dof_iterator Tria<MeshConfig>::cend_skeleton_dofs(
    const Uint dim, const dof_storage_type &dofs) const
{

  const_trace_dof_iterator it(dofs, m_skeleton[dim - 1],
                              ActiveIdx(m_skeleton[dim - 1].nb_all_facets()));
  return it;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const MeshStatistics &Tria<MeshConfig>::statistics() const
{
  return m_stats;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::print_cell_data() const
{
  std::cout << m_cells << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::print_complete_skeleton(const Uint dim)
{
  std::cout << m_skeleton[dim - 1];
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::write_dual_graph_to_gmsh(const std::string &dof_handler_name,
                                                const std::string &outfilename) const
{
  // Gravity centers of all active cells
  std::vector<math::DenseSVec<Real, MeshConfig::GDIM>> gravity_centers;
  gravity_centers.resize(m_cells.nb_active_cells());

  common::PtrHandle<dof_storage_type const> dofs = this->dof_storage(dof_handler_name);

  adapt::LocalInterpolator loc_interpolator;

  for (Uint ac = 0; ac < m_cells.nb_active_cells(); ++ac)
  {
    const CellTopologyView<MeshConfig> tcell_view = (*dofs).tcell(ActiveIdx(ac));
    const MeshEntity active_cell                  = (*dofs).active_cell(ActiveIdx(ac));

    const math::DenseConstMatView<Real> cell_coord = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), active_cell.pt_set_id(), tcell_view.coordinates());

    math::DenseSVec<Real, MeshConfig::GDIM> &gc = gravity_centers[ac];
    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      gc[d] = 0.0;
    }

    for (Uint v = 0; v < cell_coord.rows(); ++v)
    {
      gc += cell_coord.row_transpose(v);
    }

    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      gc[d] /= cell_coord.rows();
    }
  }

  std::ofstream out_stream;
  out_stream.open(outfilename.c_str());

  out_stream << "$MeshFormat" << std::endl;
  out_stream << "2.2 0 8" << std::endl;
  out_stream << "$EndMeshFormat" << std::endl;
  out_stream << "$PhysicalNames" << std::endl;
  out_stream << "1\n1000 1000 \"DualGraph\"" << std::endl;
  out_stream << "$EndPhysicalNames" << std::endl;
  out_stream << "$Nodes" << std::endl;
  out_stream << gravity_centers.size() << std::endl;

  out_stream.precision(15);

  for (Uint i = 0; i < gravity_centers.size(); ++i)
  {
    out_stream << i + 1;
    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      out_stream << " " << gravity_centers[i][d];
    }

    for (Uint d = MeshConfig::GDIM; d < 3; ++d)
    {
      out_stream << " 0.0";
    }

    out_stream << std::endl;
  }

  out_stream << "$EndNodes\n$Elements" << std::endl;

  internal::TriaFacets<MeshConfig> const &skeleton_facets = m_skeleton[MeshConfig::TDIM - 2];

  Uint skeleton_entry_nr = 0;

  for (Uint f = 0; f < skeleton_facets.nb_active_facets(); ++f)
  {
    const mesh::TraceIncidences facet_block = skeleton_facets.active_facet_data(ActiveIdx(f));
    if (facet_block.size() == 2)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_L =
          m_cells.cell(FlatIdx(facet_block.cell_id(0)));
      const mesh::CellTopologyView<MeshConfig> tcell_R =
          m_cells.cell(FlatIdx(facet_block.cell_id(1)));

      if ((tcell_L.status() == EntityStatus::Active) && (tcell_R.status() == EntityStatus::Active))
      {
        skeleton_entry_nr++;
      }
    }
  }

  out_stream << skeleton_entry_nr << std::endl;

  skeleton_entry_nr = 0;

  for (Uint f = 0; f < skeleton_facets.nb_active_facets(); ++f)
  {
    const mesh::TraceIncidences facet_block = skeleton_facets.active_facet_data(ActiveIdx(f));
    if (facet_block.size() == 2)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_L =
          m_cells.cell(FlatIdx(facet_block.cell_id(0)));
      const mesh::CellTopologyView<MeshConfig> tcell_R =
          m_cells.cell(FlatIdx(facet_block.cell_id(1)));

      if ((tcell_L.status() == EntityStatus::Active) && (tcell_R.status() == EntityStatus::Active))
      {
        out_stream << skeleton_entry_nr + 1 << " 1 2 1000 1000 " << tcell_L.active_idx().id() + 1
                   << " " << tcell_R.active_idx().id() + 1 << std::endl;
        skeleton_entry_nr++;
      }
    }
  }

  out_stream << "$EndElements" << std::endl;

  out_stream.close();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::print_cell_coords_to_file(const std::string &filename) const
{
  m_cells.print_cell_coords_to_file(filename);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::color_cells_red_green(std::vector<CellTransform> &adapt_op,
                                             std::vector<Uint> &colors) const
{
  internal::TriaFacets<MeshConfig> const &skeleton_facets = m_skeleton[MeshConfig::TDIM - 2];

  internal::MeshTopologyAdaptAlgorithm<MeshConfig>::generate_red_adapt_ops(m_cells, skeleton_facets,
                                                                           adapt_op, colors);

  for (Uint i = 0; i < adapt_op.size(); ++i)
  {
    if (adapt_op[i] == CellTransform::UNIFORM_REFINE)
    {
      colors[i] = 1;
    }
    else if (CellTransformTraits::is_aniso_refinement(adapt_op[i]))
    {
      colors[i] = 100 + CellTransformTraits::aniso_refinement_face_id(adapt_op[i]);
    }
    else
    {
      colors[i] = 300;
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::rebuild_active_skeleton()
{
  internal::TriaFacets<MeshConfig> &skeleton_facets = m_skeleton[MeshConfig::TDIM - 2];
  std::cout << "Size of original facets = " << skeleton_facets.nb_all_facets() << std::endl;

  // internal::TriaFacets<MeshConfig> new_facets;
  // new_facets.set_dim(MeshConfig::TDIM - 1);
  internal::MeshTopologyAdaptAlgorithm<MeshConfig>::rebuild_active_skeleton(
      m_cells, m_zero_level_facets, skeleton_facets);

  std::cout << "Size of new facets = " << skeleton_facets.nb_all_facets() << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
common::PtrHandle<typename Tria<MeshConfig>::dof_storage_type> Tria<MeshConfig>::dof_storage(
    const std::string &name)
{
  for (Uint i = 0; i < m_dofs.size(); ++i)
  {
    if (m_dofs[i]->name() == name)
    {
      return common::PtrHandle<dof_storage_type>(m_dofs[i].get());
    }
  }
  return common::PtrHandle<dof_storage_type>();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
common::PtrHandle<typename Tria<MeshConfig>::dof_storage_type const> Tria<MeshConfig>::dof_storage(
    const std::string &name) const
{
  for (Uint i = 0; i < m_dofs.size(); ++i)
  {
    if (m_dofs[i]->name() == name)
    {
      return common::PtrHandle<dof_storage_type const>(m_dofs[i].get());
    }
  }
  return common::PtrHandle<dof_storage_type const>();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
common::PtrHandle<typename Tria<MeshConfig>::dof_storage_type> Tria<MeshConfig>::create_dof_storage(
    const std::string &name)
{
  for (Uint i = 0; i < m_dofs.size(); ++i)
  {
    if (m_dofs[i]->name() == name)
    {
      std::cerr << "Mesh::create_dof_storage: dof handler with name " << name << " already exists!"
                << std::endl;
      return common::PtrHandle<dof_storage_type>(m_dofs[i].get());
    }
  }

  std::unique_ptr<dof_storage_type> new_dofs(new dof_storage_type(this, name));
  new_dofs->init(MeshConfig::TDIM);
  m_dofs.push_back(std::move(new_dofs));

  common::PtrHandle<dof_storage_type> dofs_ptr_handle(m_dofs.back().get());
  return dofs_ptr_handle;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::print_dof_handler_names() const
{
  std::cout << "Mesh [" << name() << "] has " << m_dofs.size() << " dof handlers: " << std::endl;
  for (Uint i = 0; i < m_dofs.size(); ++i)
  {
    std::cout << "  " << m_dofs[i]->name() << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::add_active_adjacent_cells(
    const CellTopologyView<MeshConfig> center_cell, const TraceIncidences &facet_block,
    std::vector<CellTopologyView<MeshConfig>> &adjacent_cells) const
{
  const skeleton_data &skeleton_facets = m_skeleton[MeshConfig::TDIM - 2];
  internal::MeshTopologyAdaptAlgorithm<MeshConfig>::add_active_adjacent_cells(
      m_cells, skeleton_facets, center_cell, facet_block, adjacent_cells);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void Tria<MeshConfig>::update_statistics()
{
  m_stats.tot_nb_cells  = m_cells.nb_all_cells();
  m_stats.tot_nb_facets = 0;

  for (Uint d = 0; d < m_skeleton.size(); ++d)
  {
    m_stats.tot_nb_facets += m_skeleton[d].nb_all_facets();
  }

  m_stats.tot_nb_cells_in_level.resize(m_cells.nb_all_levels(), 0);
  for (Uint c = 0; c < m_cells.nb_all_cells(); ++c)
  {
    const CellTopologyView<MeshConfig> tcell = m_cells.cell(FlatIdx(c));
    m_stats.tot_nb_cells_in_level[tcell.refinement_level()]++;
  }
}

// ----------------------------------------------------------------------------
} // namespace mesh

} // namespace pdekit

#endif
