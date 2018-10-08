#include "interpolation/L2ProjectionLocal.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

L2ProjectionLocal::L2ProjectionLocal()
{
  m_tmp_buffer_in.resize(0);
  m_tmp_buffer_out.resize(0);
}

// ----------------------------------------------------------------------------

L2ProjectionLocal::~L2ProjectionLocal()
{
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit
