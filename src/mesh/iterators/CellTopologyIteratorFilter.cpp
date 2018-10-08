
#include "mesh/iterators/CellTopologyIteratorFilter.hpp"
#include "mesh/EntityStatus.hpp"

namespace pdekit
{

namespace mesh
{

CellTopologyIterFilterTyped::CellTopologyIterFilterTyped() : m_cell_type()
{
}

CellTopologyIterFilterTyped::CellTopologyIterFilterTyped(const PointSetTag std_reg_tag)
    : m_cell_type(std_reg_tag)
{
}

} // namespace mesh

} // namespace pdekit
