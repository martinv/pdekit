#ifndef PDEKIT_Mesh_Iterators_Cell_Topology_Iterator_Filter_hpp
#define PDEKIT_Mesh_Iterators_Cell_Topology_Iterator_Filter_hpp

#include "mesh/EntityStatus.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Default filter - accepts TopologyCell as a view of any entry in TriaCells
// ----------------------------------------------------------------------------

class CellTopologyIterFilterDefault
{
  public:
  CellTopologyIterFilterDefault() = default;

  ~CellTopologyIterFilterDefault() = default;

  template <typename ViewType>
  inline constexpr bool filter_pass(const ViewType &view) const
  {
    return true;
  }
};

// ----------------------------------------------------------------------------
// Filter based on types of standard regions
// ----------------------------------------------------------------------------

class CellTopologyIterFilterActive
{
  public:
  CellTopologyIterFilterActive() = default;

  CellTopologyIterFilterActive(const CellTopologyIterFilterActive &other_filter) = default;

  ~CellTopologyIterFilterActive() = default;

  CellTopologyIterFilterActive &operator=(const CellTopologyIterFilterActive &other_filter) =
      default;

  template <typename TopoCellType>
  inline const bool filter_pass(const TopoCellType &cell) const
  {
    return (cell.status() == EntityStatus::Active);
  }
};

// ----------------------------------------------------------------------------
// Filter based on types of standard regions
// ----------------------------------------------------------------------------

class CellTopologyIterFilterTyped
{
  public:
  CellTopologyIterFilterTyped();

  CellTopologyIterFilterTyped(const PointSetTag std_reg_tag);

  CellTopologyIterFilterTyped(const CellTopologyIterFilterTyped &other_filter) = default;

  ~CellTopologyIterFilterTyped() = default;

  CellTopologyIterFilterTyped &operator=(const CellTopologyIterFilterTyped &other_filter) = default;

  template <typename ViewType>
  inline const bool filter_pass(const ViewType &view) const
  {
    return (view.cell_type().std_region_id() == m_cell_type);
  }

  private:
  PointSetTag m_cell_type;
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
