#ifndef PDEKIT_Mesh_Mesh_Entity_Iterator_hpp
#define PDEKIT_Mesh_Mesh_Entity_Iterator_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

class MeshEntity;

// ----------------------------------------------------------------------------

class MeshEntityIterator
{
  public:
  // Default constructor
  explicit MeshEntityIterator();

  // Construct from mesh entity and vertex index
  explicit MeshEntityIterator(const MeshEntity &entity, const Uint vert_id);

  void swap(MeshEntityIterator &other) noexcept;

  // Pre-increment
  MeshEntityIterator &operator++();

  // Post-increment
  MeshEntityIterator operator++(int);

  // two-way comparison: v.begin() == v.cbegin() and vice versa
  bool operator==(const MeshEntityIterator &rhs) const;

  bool operator!=(const MeshEntityIterator &rhs) const;

  Uint operator*() const;

  Uint operator->() const;

  private:
  MeshEntity const *m_entity;
  Uint m_vert_id;
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

// ----------------------------------------------------------------------------

#endif
