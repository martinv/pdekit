#include "mesh/CellMarker.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

CellMarker::CellMarker() : m_marker_type("")
{
}

// ----------------------------------------------------------------------------

CellMarker::~CellMarker()
{
}

// ----------------------------------------------------------------------------

void CellMarker::copy(const CellMarker &source, CellMarker &target)
{
  target.m_marker_type = source.m_marker_type;
  target.m_names.clear();

  // Copy the maps
  target.m_names = source.m_names;

  // Copy the values
  target.m_values.resize(source.m_values.size());
  std::copy(source.m_values.begin(), source.m_values.end(), target.m_values.begin());
}

// ----------------------------------------------------------------------------

void CellMarker::set_type_name(const std::string &type_name)
{
  m_marker_type = type_name;
}

// ----------------------------------------------------------------------------

void CellMarker::fill(const std::vector<std::pair<Uint, std::string>> &marker_names,
                      std::vector<Uint> &marker_values)
{
  m_names.clear();
  for (Uint i = 0; i < marker_names.size(); ++i)
  {
    m_names[marker_names[i].first] = marker_names[i].second;
  }
  m_values.swap(marker_values);
}

// ----------------------------------------------------------------------------

void CellMarker::fill(const std::map<Uint, std::string> &marker_names,
                      std::vector<Uint> &marker_values)
{
  m_names.clear();
  for (std::map<Uint, std::string>::const_iterator it = marker_names.cbegin();
       it != marker_names.cend(); ++it)
  {
    m_names[it->first] = it->second;
  }
  m_values.swap(marker_values);
}

// ----------------------------------------------------------------------------

void CellMarker::push_back_values(std::vector<Uint> const &new_values)
{
  const Uint old_size = m_values.size();
  const Uint new_size = old_size + new_values.size();

  m_values.resize(new_size);
  for (Uint i = 0; i < new_values.size(); ++i)
  {
    m_values[old_size + i] = new_values[i];
  }
}

// ----------------------------------------------------------------------------

const std::string CellMarker::type_name() const
{
  return m_marker_type;
}

// ----------------------------------------------------------------------------

const std::string CellMarker::name(const Uint c) const
{
  std::map<Uint, std::string>::const_iterator value_to_name_iter = m_names.find(m_values[c]);
  if (value_to_name_iter != m_names.end())
  {
    return value_to_name_iter->second;
  }
  return "";
}

// ----------------------------------------------------------------------------

Uint CellMarker::nb_values() const
{
  return m_values.size();
}

// ----------------------------------------------------------------------------

Uint CellMarker::nb_types() const
{
  return m_names.size();
}

// ----------------------------------------------------------------------------

void CellMarker::all_marker_names(std::map<Uint, std::string> &marker_names) const
{
  // marker_names.resize(m_names.size());
  // std::copy(m_names.begin(),m_names.end(),marker_names.begin());

  marker_names = m_names;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
