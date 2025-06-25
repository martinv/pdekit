#ifndef PDEKIT_Mesh_View_Buffer_Dof_View_hpp
#define PDEKIT_Mesh_View_Buffer_Dof_View_hpp

#include "mesh/CellGeometry.hpp"
#include "mesh/MeshEntity.hpp"
#include "mesh/MeshIndex.hpp"
#include "mesh/view/ViewTraits.hpp"

namespace pdekit
{

namespace mesh
{

// This is a proxy class for accessing information related
// to degrees of freedom stored in CellBuffer and related to
// one element

template <typename DofContainer, typename ConstTrait>
class BufferDofView
{
  public:
  static_assert(std::is_same<ConstTrait, ViewIsConst>::value ||
                    std::is_same<ConstTrait, ViewIsNotConst>::value,
                "Use ViewIsConst or ViewIsNotConst as constness trait.");

  static const bool is_immutable_proxy = std::is_same<ConstTrait, ViewIsConst>::value;
  using container_type =
      typename std::conditional<std::is_same<ConstTrait, ViewIsConst>::value,
                                const typename std::remove_const<DofContainer>::type,
                                typename std::remove_const<DofContainer>::type>::type;

  using container_ptr_type = container_type *;
  using container_ref_type = container_type &;

  enum
  {
    GDIM = DofContainer::GDIM
  };

  enum
  {
    TDIM = DofContainer::TDIM
  };

  BufferDofView() = default;

  BufferDofView(container_ref_type dof_container, const ActiveIdx active_idx);

  BufferDofView(const BufferDofView &other_view);

  ~BufferDofView();

  BufferDofView &operator=(const BufferDofView &other_view);

  BufferDofView &increment();

  BufferDofView &decrement();

  static std::ptrdiff_t distance(const BufferDofView &first, const BufferDofView &last);

  bool idx_less(const ActiveIdx idx) const;

  bool idx_greater(const ActiveIdx idx) const;

  bool idx_less_or_equal(const ActiveIdx idx) const;

  bool idx_greater_or_equal(const ActiveIdx idx) const;

  bool equal(const BufferDofView &other_view) const;

  void synchronize_with(const BufferDofView<DofContainer, ConstTrait> &other_view);

  const MeshEntity mesh_entity() const;

  const CellGeometry<GDIM> cell_geometry() const;

  PointSetTag cell_type() const;

  /// Return the tag of underlying cell
  Uint cell_tag() const;

  void print() const;

  private:
  /// Pointer to the group of cells over which we iterate
  container_ptr_type m_dof_container;

  /// Position of this view
  ActiveIdx m_active_idx;
};

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
BufferDofView<DofContainer, ConstTrait>::BufferDofView(container_ref_type dof_container,
                                                       const ActiveIdx idx)
    : m_dof_container(&dof_container), m_active_idx(idx)
{
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
BufferDofView<DofContainer, ConstTrait>::BufferDofView(const BufferDofView &other_view)
    : m_dof_container(other_view.m_dof_container), m_active_idx(other_view.m_active_idx)
{
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
BufferDofView<DofContainer, ConstTrait>::~BufferDofView()
{
  m_dof_container = nullptr;
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
BufferDofView<DofContainer, ConstTrait> &BufferDofView<DofContainer, ConstTrait>::operator=(
    const BufferDofView &other_view)
{
  m_dof_container = other_view.m_dof_container;
  m_active_idx    = other_view.m_active_idx;
  return *this;
}
// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
BufferDofView<DofContainer, ConstTrait> &BufferDofView<DofContainer, ConstTrait>::increment()
{
  m_active_idx++;
  return *this;
}
// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
BufferDofView<DofContainer, ConstTrait> &BufferDofView<DofContainer, ConstTrait>::decrement()
{
  m_active_idx--;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline std::ptrdiff_t BufferDofView<DofContainer, ConstTrait>::distance(const BufferDofView &first,
                                                                        const BufferDofView &last)
{
  return last.m_active_idx.id() - first.m_active_idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline bool BufferDofView<DofContainer, ConstTrait>::idx_less(const ActiveIdx idx) const
{
  return m_active_idx.id() < idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline bool BufferDofView<DofContainer, ConstTrait>::idx_greater(const ActiveIdx idx) const
{
  return m_active_idx.id() > idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline bool BufferDofView<DofContainer, ConstTrait>::idx_less_or_equal(const ActiveIdx idx) const
{
  return m_active_idx.id() <= idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
inline bool BufferDofView<DofContainer, ConstTrait>::idx_greater_or_equal(const ActiveIdx idx) const
{
  return m_active_idx.id() >= idx.id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
bool BufferDofView<DofContainer, ConstTrait>::equal(const BufferDofView &other_view) const
{
  return (m_dof_container == other_view.m_dof_container) &&
         (m_active_idx == other_view.m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
void BufferDofView<DofContainer, ConstTrait>::synchronize_with(
    const BufferDofView<DofContainer, ConstTrait> &other_view)
{
  /// @todo Check that both views point to the same container instance
  m_active_idx = other_view.m_active_idx;
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
const MeshEntity BufferDofView<DofContainer, ConstTrait>::mesh_entity() const
{
  return m_dof_container->active_cell(m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
const CellGeometry<BufferDofView<DofContainer, ConstTrait>::GDIM> BufferDofView<
    DofContainer, ConstTrait>::cell_geometry() const
{
  const MeshEntity entity = m_dof_container->active_cell(m_active_idx);
  const common::ArrayView<const Real, _1D, Uint> cell_coord_array =
      m_dof_container->active_cell_coords(m_active_idx);
  return CellGeometry<GDIM>(cell_coord_array, entity);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
PointSetTag BufferDofView<DofContainer, ConstTrait>::cell_type() const
{
  const MeshEntity entity = m_dof_container->active_cell(m_active_idx);
  return entity.pt_set_id();
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
Uint BufferDofView<DofContainer, ConstTrait>::cell_tag() const
{
  return m_dof_container->active_cell_tag(m_active_idx);
}

// ----------------------------------------------------------------------------

template <typename DofContainer, typename ConstTrait>
void BufferDofView<DofContainer, ConstTrait>::print() const
{
  std::cout << "DofView: curr. pos/nb. entries in container: [" << m_active_idx.id() << "/"
            << m_dof_container->nb_active_cells() << "]" << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
