#include "mesh/MeshPredicates.hpp"

namespace pdekit
{

namespace mesh
{

// ============================================================================

AllCells::AllCells()
{
}

AllCells::~AllCells()
{
}

// ============================================================================

AllEntities::AllEntities()
{
}

AllEntities::~AllEntities()
{
}

// ============================================================================

CellGroup::CellGroup()
    : m_reference_etype(CellGroup::ConfigParamT1()), m_current_etype(CellGroup::ConfigParamT1())
{
}

CellGroup::~CellGroup()
{
}

// ============================================================================

} // namespace mesh

} // namespace pdekit
