#ifndef PDEKIT_Mesh_View_Cell_Topology_View_hpp
#define PDEKIT_Mesh_View_Cell_Topology_View_hpp

#include "common/ArrayView.hpp"
#include "math/DenseConstMatView.hpp"
#include "mesh/CellGeometry.hpp"
#include "mesh/EntityStatus.hpp"
#include "mesh/MeshIndex.hpp"
#include "mesh/adaptation/CellAdaptOp.hpp"
#include "mesh/std_region/StdRegion.hpp"

// This class acts as a proxy collecting data corresponding
// to one cell in CellConnectivity

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

namespace internal
{
template <typename MeshConfig>
class TriaCells;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class CellTopologyView
{
  public:
  // Topology cell as a 'view into container' is always const - treats the
  // container as immutable and never changes its contents
  static const bool is_immutable_proxy = true;

  using container_type     = internal::TriaCells<MeshConfig>;
  using container_ptr_type = const container_type *;
  using container_ref_type = const container_type &;

  enum
  {
    GDIM = MeshConfig::GDIM
  };

  enum
  {
    TDIM = MeshConfig::TDIM
  };

  /// Default constructor
  CellTopologyView();

  /// Construct from cell container and given index
  CellTopologyView(container_ptr_type cells, const FlatIdx linear_pos_idx);

  /// Copy constructor
  CellTopologyView(const CellTopologyView &other_cell);

  /// Assignment operator
  CellTopologyView &operator=(const CellTopologyView &other_cell);

  /// Default destructor
  ~CellTopologyView();

  CellTopologyView &increment();

  CellTopologyView &decrement();

  bool equal(const CellTopologyView &other_cell) const;

  /// Return the index of this cell
  const FlatIdx linear_pos_idx() const;

  /// Return the index among active cells
  const ActiveIdx active_idx() const;

  /// Return the indexes of faces of this cell
  const common::ArrayView<const Uint, _1D, Uint> incident_facets() const;

  /// Return the shape of this cell
  ElemShape cell_shape() const;

  /// Return the shape of this cell
  PointSetTag pt_set_id() const;

  /// Return the shape of this cell's sub-entity
  PointSetTag pt_set_id(const Uint dim, const Uint idx) const;

  /// Return the cell type of this cell
  StdRegion std_region() const;

  /// Return the type of a sub-entity
  std::shared_ptr<const StdRegionEntity> sub_entity(const Uint dim, const Uint idx) const;

  /// Get the type of adaptation operation that was applied to this cell
  const adapt::CellAdaptOp cell_adapt_op() const;

  /// Return on which level this cell lies
  Uint refinement_level() const;

  /// Get the status of this cell
  EntityStatus status() const;

  /// Return the index of parent of this cell
  const CellTopologyView<MeshConfig> parent() const;

  /// Return the indices of children of this cell
  const std::vector<CellTopologyView<MeshConfig>> children() const;

  /// Return the number of children of this cell
  Uint nb_children() const;

  /// Return the number of vertices
  Uint nb_vert() const;

  /// Return the coordinates of one cell
  const mesh::CellGeometry<MeshConfig::GDIM> coordinates() const;

  /// Return the coordinates of one cell
  const mesh::CellGeometry<MeshConfig::GDIM> coordinates(const Uint dim, const Uint idx) const;

  /// Return the number of all cells (on all levels, active or not)
  /// in TriaCells
  Uint nb_all_cells() const;

  int position() const;

  private:
  /// DATA
  /// Pointer to underlying connectivity container
  container_ptr_type m_mesh_cells;

  /// Index of this cell
  FlatIdx m_linear_idx;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
CellTopologyView<MeshConfig>::CellTopologyView()
    : m_mesh_cells(nullptr), m_linear_idx(INVALID_CELL_ID)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
CellTopologyView<MeshConfig>::CellTopologyView(container_ptr_type cells,
                                               const FlatIdx linear_pos_idx)
    : m_mesh_cells(cells), m_linear_idx(linear_pos_idx)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
CellTopologyView<MeshConfig>::CellTopologyView(const CellTopologyView &other_cell)
    : m_mesh_cells(other_cell.m_mesh_cells), m_linear_idx(other_cell.m_linear_idx)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
CellTopologyView<MeshConfig> &CellTopologyView<MeshConfig>::operator=(
    const CellTopologyView<MeshConfig> &other_cell)
{
  m_mesh_cells = other_cell.m_mesh_cells;
  m_linear_idx = other_cell.m_linear_idx;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
CellTopologyView<MeshConfig>::~CellTopologyView()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
CellTopologyView<MeshConfig> &CellTopologyView<MeshConfig>::increment()
{
  m_linear_idx++;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
CellTopologyView<MeshConfig> &CellTopologyView<MeshConfig>::decrement()
{
  m_linear_idx--;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
bool CellTopologyView<MeshConfig>::equal(const CellTopologyView<MeshConfig> &other_cell) const
{
  return (m_mesh_cells == other_cell.m_mesh_cells) && (m_linear_idx == other_cell.m_linear_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const FlatIdx CellTopologyView<MeshConfig>::linear_pos_idx() const
{
  return m_linear_idx;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const ActiveIdx CellTopologyView<MeshConfig>::active_idx() const
{
  return m_mesh_cells->active_idx(m_linear_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const common::ArrayView<const Uint, _1D, Uint> CellTopologyView<MeshConfig>::incident_facets() const
{
  const common::ArrayView<const Uint, _1D, Uint> faces(m_mesh_cells->cell_faces(m_linear_idx));
  return faces;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
ElemShape CellTopologyView<MeshConfig>::cell_shape() const
{
  const StdRegion std_region = m_mesh_cells->cell_type(m_linear_idx);
  return std_region.get().pt_set_id().elem_shape();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
PointSetTag CellTopologyView<MeshConfig>::pt_set_id() const
{
  const StdRegion std_region = m_mesh_cells->cell_type(m_linear_idx);
  return std_region.get().pt_set_id();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
PointSetTag CellTopologyView<MeshConfig>::pt_set_id(const Uint dim, const Uint idx) const
{
  const StdRegion std_region = m_mesh_cells->cell_type(m_linear_idx);
  return std_region.get().elem_entity(dim, idx)->pt_set_id();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
StdRegion CellTopologyView<MeshConfig>::std_region() const
{
  return m_mesh_cells->cell_type(m_linear_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
std::shared_ptr<const StdRegionEntity> CellTopologyView<MeshConfig>::sub_entity(
    const Uint dim, const Uint idx) const
{
  const StdRegion std_reg = m_mesh_cells->cell_type(m_linear_idx);
  return std_reg.get().elem_entity(dim, idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const adapt::CellAdaptOp CellTopologyView<MeshConfig>::cell_adapt_op() const
{
  return m_mesh_cells->cell_adapt_op(m_linear_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint CellTopologyView<MeshConfig>::refinement_level() const
{
  return m_mesh_cells->refinement_level(m_linear_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
EntityStatus CellTopologyView<MeshConfig>::status() const
{
  return m_mesh_cells->cell_status(m_linear_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const CellTopologyView<MeshConfig> CellTopologyView<MeshConfig>::parent() const
{
  return m_mesh_cells->cell_parent(*this);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const std::vector<CellTopologyView<MeshConfig>> CellTopologyView<MeshConfig>::children() const
{
  return m_mesh_cells->cell_children(m_linear_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint CellTopologyView<MeshConfig>::nb_children() const
{
  return m_mesh_cells->nb_children(m_linear_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint CellTopologyView<MeshConfig>::nb_vert() const
{
  const StdRegion std_reg = m_mesh_cells->cell_type(m_linear_idx);
  return std_reg.get().nb_nodes();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const mesh::CellGeometry<MeshConfig::GDIM> inline CellTopologyView<MeshConfig>::coordinates() const
{
  return m_mesh_cells->cell_geometry(m_linear_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const mesh::CellGeometry<MeshConfig::GDIM> inline CellTopologyView<MeshConfig>::coordinates(
    const Uint dim, const Uint idx) const
{
  return m_mesh_cells->cell_geometry(m_linear_idx, dim, idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint CellTopologyView<MeshConfig>::nb_all_cells() const
{
  return m_mesh_cells->nb_all_cells();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
int CellTopologyView<MeshConfig>::position() const
{
  return static_cast<int>(m_linear_idx.id());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
std::ostream &operator<<(std::ostream &out_stream, const CellTopologyView<MeshConfig> &topo_cell)
{
  out_stream << "Topology cell [" << topo_cell.linear_pos_idx() << "]\n";
  out_stream << "   Linear pos. idx = " << topo_cell.linear_pos_idx() << "\n";
  out_stream << "   Active pos. idx = " << topo_cell.active_idx() << "\n";
  out_stream << "   Cell type:  " << topo_cell.std_region().get().pt_set_id().as_string() << "\n";
  out_stream << "   Parent index: " << topo_cell.parent().linear_pos_idx() << "\n";
  out_stream << "   Children:";
  const std::vector<CellTopologyView<MeshConfig>> cell_children = topo_cell.children();
  for (Uint c = 0; c < cell_children.size(); ++c)
  {
    out_stream << " " << cell_children[c].linear_pos_idx();
  }
  out_stream << "\n";

  out_stream << "   Faces: " << topo_cell.incident_facets() << "\n";
  out_stream << "   Status: " << topo_cell.status() << "\n";
  out_stream << "   Refinement level: " << topo_cell.refinement_level() << "\n";
  return out_stream;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
