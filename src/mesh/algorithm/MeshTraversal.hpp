#ifndef PDEKIT_Mesh_Algorithm_Mesh_Traversal_hpp
#define PDEKIT_Mesh_Algorithm_Mesh_Traversal_hpp

#include "common/Constants.hpp"
#include "common/IteratorRange.hpp"
#include "common/Range1D.hpp"
#include "mesh/MeshEntity.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class DofMap;

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class InteriorRandAccTraversal
{
  public:
  InteriorRandAccTraversal(const mesh::DofMap<MeshConfig> &dofs,
                           const common::Range1D<Uint> &cell_range)
      : m_dofs(dofs), m_cell_range(cell_range), m_current(cell_range.lbound())
  {
  }

  InteriorRandAccTraversal(const InteriorRandAccTraversal &other) = delete;

  ~InteriorRandAccTraversal() = default;

  InteriorRandAccTraversal &operator=(const InteriorRandAccTraversal &other) = delete;

  void init()
  {
    m_current = m_cell_range.lbound();
  }

  mesh::MeshEntity make_next()
  {
    const mesh::MeshEntity cell = m_dofs.active_cell(mesh::ActiveIdx(m_current));
    m_current++;

    return cell;
  }

  bool reached_end() const
  {
    return m_current == (m_cell_range.ubound() + 1);
  }

  private:
  const mesh::DofMap<MeshConfig> &m_dofs;
  const common::Range1D<Uint> m_cell_range;
  Uint m_current;
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
