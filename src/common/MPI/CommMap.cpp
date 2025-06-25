#include "common/MPI/CommMap.hpp"

namespace pdekit
{

namespace common
{

namespace mpi
{

// ----------------------------------------------------------------------------

CommMap::CommMap(const std::string &name) : m_name(name)
{
  m_local_to_global_id_map.resize(0);
  m_global_to_local_id_map.clear();
}

// ----------------------------------------------------------------------------

CommMap::~CommMap()
{
}

// ----------------------------------------------------------------------------

const std::string &CommMap::name() const
{
  return m_name;
}

// ----------------------------------------------------------------------------

Uint CommMap::size() const
{
  return m_local_to_global_id_map.size();
}

// ----------------------------------------------------------------------------

void CommMap::clear()
{
  m_local_to_global_id_map.clear();
  m_global_to_local_id_map.clear();
}

// ----------------------------------------------------------------------------

void CommMap::fill(std::vector<idx_type> &global_id)
{
  m_local_to_global_id_map.swap(global_id);
}

// ----------------------------------------------------------------------------

} // namespace mpi

} // namespace common

} // namespace pdekit
