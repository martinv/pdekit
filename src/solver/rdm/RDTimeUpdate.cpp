#include "solver/rdm/RDTimeUpdate.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

// ----------------------------------------------------------------------------

RDTimeUpdate::RDTimeUpdate()
    : m_dual_nodal_volume("", "dual_volume"), m_wave_speeds_per_node("", "nodal_wave_speeds")
{
}

// ----------------------------------------------------------------------------

RDTimeUpdate::~RDTimeUpdate()
{
}

// ----------------------------------------------------------------------------

const interpolation::ScalarMeshFunction<Real> &RDTimeUpdate::nodal_dual_volume() const
{
  return m_dual_nodal_volume;
}

// ----------------------------------------------------------------------------

const interpolation::ScalarMeshFunction<Real> &RDTimeUpdate::nodal_wave_speed() const
{
  return m_wave_speeds_per_node;
}

// ----------------------------------------------------------------------------

void RDTimeUpdate::compute_local_time_step(const Real CFL,
                                           interpolation::ScalarMeshFunction<Real> &time_step)
{
  if (time_step.nb_entries() != m_dual_nodal_volume.nb_entries())
  {
    time_step.resize(m_dual_nodal_volume.nb_entries());
  }

  for (Uint n = 0; n < time_step.nb_entries(); ++n)
  {
    time_step[n] = CFL * m_dual_nodal_volume[n] / m_wave_speeds_per_node[n];
  }
}

// ----------------------------------------------------------------------------

void RDTimeUpdate::reset_wave_speeds()
{
  m_wave_speeds_per_node.fill(0.0);
}

// ----------------------------------------------------------------------------

void RDTimeUpdate::accumulate_nodal_wave_speed(const Uint node_idx, const Real ws)
{
  m_wave_speeds_per_node[node_idx] += ws;
}

// ----------------------------------------------------------------------------

void RDTimeUpdate::accumulate_wave_speeds(const std::vector<std::tuple<Uint, Real>> &ws_buffer)
{
  std::lock_guard<std::mutex> lock(m_mutex);

  for (Uint i = 0; i < ws_buffer.size(); ++i)
  {
    m_wave_speeds_per_node[std::get<0>(ws_buffer[i])] += std::get<1>(ws_buffer[i]);
  }
}

// ----------------------------------------------------------------------------

Real RDTimeUpdate::wave_speed(const Uint node_idx) const
{
  return m_wave_speeds_per_node[node_idx];
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit
