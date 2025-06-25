#include "interpolation/ResidualRestrictionLocal.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

ResidualRestrictionLocal::ResidualRestrictionLocal()
{
  m_raw_data.resize(0);
  m_local_ops.resize(0);
  m_tmp_buffer_in.resize(0);
  m_tmp_buffer_out.resize(0);
}

// ----------------------------------------------------------------------------

ResidualRestrictionLocal::~ResidualRestrictionLocal()
{
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit
