#ifndef PDEKIT_Mesh_Iterators_Dof_Iterator_Filter_hpp
#define PDEKIT_Mesh_Iterators_Dof_Iterator_Filter_hpp

#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Default filter - accepts dof view of any entry in DofMap (i.e. no
// filtering)
// ----------------------------------------------------------------------------

class DofIterFilterDefault
{
  public:
  DofIterFilterDefault() = default;

  ~DofIterFilterDefault() = default;

  template <typename ViewType>
  inline constexpr bool filter_pass(const ViewType &view) const
  {
    return true;
  }
};

// ----------------------------------------------------------------------------
// Filter based on types of standard regions
// ----------------------------------------------------------------------------

class DofIterFilterTyped
{
  public:
  DofIterFilterTyped();

  DofIterFilterTyped(const PointSetTag std_reg_tag);

  DofIterFilterTyped(const DofIterFilterTyped &other_filter) = default;

  ~DofIterFilterTyped() = default;

  DofIterFilterTyped &operator=(const DofIterFilterTyped &other_filter) = default;

  template <typename ViewType>
  inline const bool filter_pass(const ViewType &view) const
  {
    return (view.pt_set_id() == m_cell_type);
  }

  private:
  PointSetTag m_cell_type;
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
