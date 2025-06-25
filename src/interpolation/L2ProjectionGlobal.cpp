#include "interpolation/L2ProjectionGlobal.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

L2ProjectionGlobal::L2ProjectionGlobal()
    : m_projection_matrix(nullptr), m_rhs_vector(nullptr), m_solution(nullptr)
{
}

// ----------------------------------------------------------------------------

L2ProjectionGlobal::~L2ProjectionGlobal()
{
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit
