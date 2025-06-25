#include "mesh/iterators/DofIteratorFilter.hpp"

namespace pdekit
{

namespace mesh
{

DofIterFilterTyped::DofIterFilterTyped() : m_cell_type()
{
}

DofIterFilterTyped::DofIterFilterTyped(const PointSetTag std_reg_tag) : m_cell_type(std_reg_tag)
{
}

} // namespace mesh

} // namespace pdekit
