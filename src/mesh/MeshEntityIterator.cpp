#include "mesh/MeshEntityIterator.hpp"
#include "mesh/MeshEntity.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

MeshEntityIterator::MeshEntityIterator() : m_entity(nullptr), m_vert_id(0)
{
}

// ----------------------------------------------------------------------------

MeshEntityIterator::MeshEntityIterator(const MeshEntity &entity, const Uint vert_id)
    : m_entity(&entity), m_vert_id(vert_id)
{
}

// ----------------------------------------------------------------------------

void MeshEntityIterator::swap(MeshEntityIterator &other) noexcept
{
  std::swap(m_entity, other.m_entity);
  std::swap(m_vert_id, other.m_vert_id);
}

// ----------------------------------------------------------------------------

MeshEntityIterator &MeshEntityIterator::operator++()
{
  m_vert_id++;
  return *this;
}

// ----------------------------------------------------------------------------

MeshEntityIterator MeshEntityIterator::operator++(int)
{
  MeshEntityIterator tmp(*this);
  m_vert_id++;
  return tmp;
}

// ----------------------------------------------------------------------------

bool MeshEntityIterator::operator==(const MeshEntityIterator &rhs) const
{
  return (m_entity == rhs.m_entity) && (m_vert_id == rhs.m_vert_id);
}

// ----------------------------------------------------------------------------

bool MeshEntityIterator::operator!=(const MeshEntityIterator &rhs) const
{
  return !this->operator==(rhs);
}

// ----------------------------------------------------------------------------

Uint MeshEntityIterator::operator*() const
{
  return m_entity->vertex(m_vert_id);
}

// ----------------------------------------------------------------------------

Uint MeshEntityIterator::operator->() const
{
  return m_entity->vertex(m_vert_id);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

// ----------------------------------------------------------------------------
