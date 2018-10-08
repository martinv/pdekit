#ifndef PDEKIT_Mesh_Cell_Buffer_hpp
#define PDEKIT_Mesh_Cell_Buffer_hpp

#include "common/BlockArray.hpp"
#include "mesh/MeshEntity.hpp"
#include "mesh/MeshIndex.hpp"
#include "mesh/iterators/CellDofIterator.hpp"
#include "mesh/view/BufferDofView.hpp"
#include <memory>
#include <vector>

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
class CellBuffer
{
  public:
  using iterator = CellDofIterator<BufferDofView<CellBuffer<GeoDim, TDim>, ViewIsNotConst>>;
  using const_dof_iterator = CellDofIterator<BufferDofView<CellBuffer<GeoDim, TDim>, ViewIsConst>>;

  enum
  {
    GDIM = GeoDim
  };

  enum
  {
    TDIM = TDim
  };

  /// Constructor
  CellBuffer();

  /// Copy constructor
  CellBuffer(const CellBuffer &buffer_rhs);

  /// Destructor
  ~CellBuffer();

  /// Assignment operator
  CellBuffer &operator=(const CellBuffer &buffer_rhs);

  /// Return geometry dimension
  Uint geo_dim() const;

  /// Return topological dimension of stored cells
  Uint topo_dim() const;

  /// Return the number of node_id entries
  Uint nb_node_entries() const;

  /// Return the number of cells
  Uint nb_active_cells() const;

  /// Return the capacity of the node_id array
  Uint capacity_node_entries() const;

  /// Resize the node storage and the mesh entity storage
  void reserve(const Uint nb_node_entries, const Uint nb_active_cells);

  /// Clear all internal data
  void clear();

  /// Insert a new cell
  /// @param cell_nb       ... cell index (number)
  /// @param cell_type     ... defines the standard region of this cell
  /// @param cell_vertices ... node indices of the cell
  /// @param material_id   ... physical tag/material tag of the cell
  /// @param cell_coords   ... nodal coordinates of the cell. The coordinates
  /// are ordered
  ///                          so that first come (x0,y0,z0) coordinates of
  ///                          node 0 of the cell, then (x1,y1, z1), i.e. the
  ///                          coordinates of node 1 etc.
  void push_back_cell(const Uint cell_nb, const PointSetTag cell_type,
                      const std::vector<Uint> &cell_vertices, const Uint material_id,
                      const std::vector<Real> &cell_coords);

  /// Insert a group of cells. This method resets all internal data
  void emplace_cell_data(std::unique_ptr<std::vector<Uint>> cell_connections,
                         std::unique_ptr<std::vector<StdRegion>> cell_types,
                         std::unique_ptr<std::vector<Uint>> material_ids,
                         std::unique_ptr<common::BlockArray<Real, Uint>> cell_coords);

  /// Make sure that the there are no gaps in numbering of DOFs
  void remove_dof_numbering_gaps();

  /// Print all the cells stored
  void print(const bool display_cell_types = false) const;

  /// Return the i-th cell
  const MeshEntity active_cell(const ActiveIdx cell_nb) const;

  /// Return the coordinate array of all nodes of one cell
  const common::ArrayView<const Real, _1D, Uint> active_cell_coords(const ActiveIdx cell_idx) const;

  /// Fill entity so that it corresponds to the cell on position i
  void fill_cell(const ActiveIdx idx, MeshEntity &entity) const;

  /// Return the material id ('tag') of the i-th cell
  Uint active_cell_tag(const ActiveIdx cell_nb) const;

  /// Return array with coordinates of all nodes in every cell
  const common::BlockArray<Real, Uint> &cell_coordinates() const;

  /// Create an iterator pointing to the beginning of the cell array, const
  /// version
  const_dof_iterator begin() const;

  /// Return const iterator regardless of whether Cells is a constant object
  /// or not
  const_dof_iterator cbegin() const;

  /// Create an iterator pointing 1 position after the end of the cell array
  const_dof_iterator end() const;

  /// Return const iterator regardless of whether Cells is a constant object
  /// or not
  const_dof_iterator cend() const;

  private:
  /// Linearized connectivity of all cells
  std::unique_ptr<common::BlockArray<Uint, Uint>> m_node_ids;

  /// Type of each cell
  std::unique_ptr<std::vector<mesh::StdRegion>> m_cell_types;

  /// Index of each cell
  std::unique_ptr<std::vector<Uint>> m_cell_ids;

  /// Material id of each cell
  std::unique_ptr<std::vector<Uint>> m_cell_tag;

  /// Coordinates of nodes of each cell
  std::unique_ptr<common::BlockArray<Real, Uint>> m_cell_coords;
};

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
CellBuffer<GeoDim, TDim>::CellBuffer()
    : m_node_ids(new common::BlockArray<Uint, Uint>()), m_cell_types(new std::vector<StdRegion>()),
      m_cell_ids(new std::vector<Uint>()), m_cell_tag(new std::vector<Uint>()),
      m_cell_coords(new common::BlockArray<Real, Uint>())
{
  m_node_ids->resize(0, 0);
  m_cell_ids->resize(0);
  m_cell_tag->resize(0);
  m_cell_coords->resize(0, 0);
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
CellBuffer<GeoDim, TDim>::CellBuffer(const CellBuffer &buffer_rhs)
{
  m_node_ids->resize(buffer_rhs.m_node_ids->size(), buffer_rhs.m_node_ids->nb_blocks());
  *m_node_ids = *(buffer_rhs.m_node_ids);

  m_cell_types->resize(buffer_rhs.m_cell_types->size());
  *m_cell_types = *(buffer_rhs.m_cell_types);

  m_cell_ids->resize(buffer_rhs.m_cell_ids->size());
  *m_cell_ids = *(buffer_rhs.m_cell_ids);

  m_cell_tag->resize(buffer_rhs.m_cell_tag->size());
  *m_cell_tag = *(buffer_rhs.m_cell_tag);

  m_cell_coords->resize(buffer_rhs.m_cell_coords->size(), buffer_rhs.m_cell_coords->nb_blocks());
  *m_cell_coords = (*buffer_rhs.m_cell_coords);
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
CellBuffer<GeoDim, TDim>::~CellBuffer()
{
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
CellBuffer<GeoDim, TDim> &CellBuffer<GeoDim, TDim>::operator=(const CellBuffer &buffer_rhs)
{
  m_node_ids->resize(buffer_rhs.m_node_ids->size(), buffer_rhs.m_node_ids->nb_blocks());
  *m_node_ids = *(buffer_rhs.m_node_ids);

  m_cell_types->resize(buffer_rhs.m_cell_types->size());
  *m_cell_types = *(buffer_rhs.m_cell_types);

  m_cell_ids->resize(buffer_rhs.m_cell_ids->size());
  *m_cell_ids = *(buffer_rhs.m_cell_ids);

  m_cell_tag->resize(buffer_rhs.m_cell_tag->size());
  *m_cell_tag = *(buffer_rhs).m_cell_tag;

  m_cell_coords->resize(buffer_rhs.m_cell_coords->size(), buffer_rhs.m_cell_coords->nb_blocks());
  *m_cell_coords = (*buffer_rhs.m_cell_coords);

  return *this;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
Uint CellBuffer<GeoDim, TDim>::geo_dim() const
{
  return GDIM;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
Uint CellBuffer<GeoDim, TDim>::topo_dim() const
{
  return TDIM;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
Uint CellBuffer<GeoDim, TDim>::nb_node_entries() const
{
  return m_node_ids->size();
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
Uint CellBuffer<GeoDim, TDim>::nb_active_cells() const
{
  return m_node_ids->nb_blocks();
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
Uint CellBuffer<GeoDim, TDim>::capacity_node_entries() const
{
  return m_node_ids->capacity();
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
void CellBuffer<GeoDim, TDim>::reserve(const Uint nb_node_entries, const Uint nb_cells)
{
  m_node_ids->reserve(nb_node_entries, nb_cells);
  m_cell_types->reserve(nb_cells);
  m_cell_ids->reserve(nb_cells);
  m_cell_tag->reserve(nb_cells);
  m_cell_coords->reserve(GDIM * nb_node_entries, nb_cells);
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
void CellBuffer<GeoDim, TDim>::clear()
{
  m_node_ids->clear();
  m_cell_types->clear();
  m_cell_ids->clear();
  m_cell_tag->clear();
  m_cell_coords->clear();
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
void CellBuffer<GeoDim, TDim>::push_back_cell(const Uint cell_nb, const PointSetTag cell_type,
                                              const std::vector<Uint> &cell_vertices,
                                              const Uint material_id,
                                              const std::vector<Real> &cell_coords)
{
  common::ArrayView<const Uint, _1D, Uint> vert_ids(cell_vertices.data(), cell_vertices.size());

  m_node_ids->create_back_block(vert_ids.size());
  m_node_ids->fill_last_block(vert_ids);

  StdRegion new_cell(cell_type);
  m_cell_types->push_back(new_cell);
  m_cell_ids->push_back(cell_nb);
  m_cell_tag->push_back(material_id);

  const common::ArrayView<const Real, _1D, Uint> block_values_view(cell_coords.data(),
                                                                   cell_coords.size());

  m_cell_coords->create_back_block(block_values_view.size());
  m_cell_coords->fill_last_block(block_values_view);
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
void CellBuffer<GeoDim, TDim>::emplace_cell_data(
    std::unique_ptr<std::vector<Uint>> cell_connections,
    std::unique_ptr<std::vector<StdRegion>> cell_types,
    std::unique_ptr<std::vector<Uint>> material_ids,
    std::unique_ptr<common::BlockArray<Real, Uint>> cell_coords)
{
  m_node_ids->clear();
  m_node_ids->resize(0, 0);
  m_node_ids->reserve(cell_connections->size(), cell_types->size());

  Uint fill_pos = 0;

  std::vector<Uint> single_cell_nodes;

  for (Uint c = 0; c < cell_types->size(); ++c)
  {
    single_cell_nodes.resize(0);
    const StdRegion cell_std_reg((*cell_types)[c]);

    for (Uint v = 0; v < cell_std_reg.get().nb_nodes(); ++v)
    {
      single_cell_nodes.push_back((*cell_connections)[fill_pos++]);
    }

    const common::ArrayView<const Uint, _1D, Uint> cell_verts(single_cell_nodes.data(),
                                                              single_cell_nodes.size());

    m_node_ids->create_back_block(cell_verts.size());
    m_node_ids->fill_last_block(cell_verts);

    // m_cell_types->push_back(cell_std_reg);
    m_cell_ids->push_back(c);
  }

  std::swap(m_cell_types, cell_types);

  std::swap(m_cell_tag, material_ids);

  std::swap(m_cell_coords, cell_coords);
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
void CellBuffer<GeoDim, TDim>::remove_dof_numbering_gaps()
{
  std::map<Uint, Uint> dof_remap;

  for (Uint b = 0; b < m_node_ids->nb_blocks(); ++b)
  {
    const common::ArrayView<Uint, _1D, Uint> dof_block = m_node_ids->block(b);
    for (Uint i = 0; i < dof_block.size(); ++i)
    {
      const Uint one_dof_old = dof_block[i];
      dof_remap[one_dof_old] = 0;
    }
  }

  Uint next_free_id = 0;

  for (auto it = dof_remap.begin(); it != dof_remap.end(); ++it)
  {
    it->second = next_free_id++;
  }

  std::vector<Uint> new_dof_ids;

  // Now change the numbering of nodes in m_node_ids;
  for (Uint b = 0; b < m_node_ids->nb_blocks(); ++b)
  {
    const common::ArrayView<Uint, _1D, Uint> dof_block = m_node_ids->block(b);
    new_dof_ids.resize(dof_block.size());
    for (Uint i = 0; i < dof_block.size(); ++i)
    {
      const auto it  = dof_remap.find(dof_block[i]);
      new_dof_ids[i] = it->second;
    }

    const common::ArrayView<const Uint, _1D, Uint> new_dof_ids_view(new_dof_ids.data(),
                                                                    new_dof_ids.size());
    m_node_ids->insert_block(b, new_dof_ids_view);
  } // Loop over blocks in m_node_ids
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
void CellBuffer<GeoDim, TDim>::print(const bool display_cell_types) const
{
  for (Uint i = 0; i < m_cell_types->size(); ++i)
  {
    const MeshEntity c = active_cell(ActiveIdx(i));
    std::cout << "[" << c.idx() << "]  " << c;
    if (display_cell_types)
    {
      std::cout << " " << c.pt_set_id().as_string();
    }
    std::cout << " tag = " << (*m_cell_tag)[i] << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
inline const MeshEntity CellBuffer<GeoDim, TDim>::active_cell(const ActiveIdx cell_nb) const
{
  const MeshEntity buffer_cell(m_node_ids->const_block(cell_nb.id()), (*m_cell_ids)[cell_nb.id()],
                               (*m_cell_types)[cell_nb.id()]);

  return buffer_cell;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
const common::ArrayView<const Real, _1D, Uint> CellBuffer<GeoDim, TDim>::active_cell_coords(
    const ActiveIdx cell_idx) const
{
  return m_cell_coords->const_block(cell_idx.id());
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
inline void CellBuffer<GeoDim, TDim>::fill_cell(const ActiveIdx idx, MeshEntity &entity) const
{
  /*
  entity.set_limits(m_node_ids->const_block(idx.id()));
  entity.set_idx(idx.id());
  entity.set_type((*m_cell_types)[idx.id()]);
  */
  entity.reinit(m_node_ids->const_block(idx.id()), idx.id(), (*m_cell_types)[idx.id()]);
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
inline Uint CellBuffer<GeoDim, TDim>::active_cell_tag(const ActiveIdx cell_nb) const
{
  return (*m_cell_tag)[cell_nb.id()];
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
const common::BlockArray<Real, Uint> &CellBuffer<GeoDim, TDim>::cell_coordinates() const
{
  return *m_cell_coords;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
typename CellBuffer<GeoDim, TDim>::const_dof_iterator CellBuffer<GeoDim, TDim>::begin() const
{
  return const_dof_iterator(*this, ActiveIdx(0));
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
typename CellBuffer<GeoDim, TDim>::const_dof_iterator CellBuffer<GeoDim, TDim>::cbegin() const
{
  return const_dof_iterator(*this, ActiveIdx(0));
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
typename CellBuffer<GeoDim, TDim>::const_dof_iterator CellBuffer<GeoDim, TDim>::end() const
{
  return const_dof_iterator(*this, ActiveIdx(m_node_ids->nb_blocks()));
}

// ----------------------------------------------------------------------------

template <Uint GeoDim, Uint TDim>
typename CellBuffer<GeoDim, TDim>::const_dof_iterator CellBuffer<GeoDim, TDim>::cend() const
{
  return const_dof_iterator(*this, ActiveIdx(m_node_ids->nb_blocks()));
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
