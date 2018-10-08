#ifndef PDEKIT_RD_Key_Generators_hpp
#define PDEKIT_RD_Key_Generators_hpp

#include <iostream>

#include "common/Constants.hpp"
#include "common/IteratorRange.hpp"
#include "common/Range1D.hpp"
#include "mesh/DiscreteElemKey.hpp"
#include "mesh/MeshEntity.hpp"

namespace pdekit
{

namespace mesh
{
template <typename MeshConfig>
class DofMap;
}

namespace solver
{

namespace rdm
{

// ----------------------------------------------------------------------------

template <typename DofIterator>
class DiscreteElemKeyGenerator
{
  public:
  DiscreteElemKeyGenerator(const common::IteratorRange<DofIterator> &dofs, SFunc sf_type,
                           PointSetID quad_type, Uint quad_order)
      : m_dof_range(dofs), m_current(dofs.begin()), m_sf_type(sf_type), m_quad_type(quad_type),
        m_quad_order(quad_order)
  {
  }

  DiscreteElemKeyGenerator(const DiscreteElemKeyGenerator &other) = delete;

  ~DiscreteElemKeyGenerator() = default;

  DiscreteElemKeyGenerator &operator=(const DiscreteElemKeyGenerator &other) = delete;

  void init()
  {
    m_current = m_dof_range.begin();
  }

  mesh::DiscreteElemKey make_next()
  {
    const mesh::MeshEntity cell      = m_current->mesh_entity();
    const mesh::PointSetTag cell_tag = cell.pt_set_id();

    const ElemShape cell_shape = cell_tag.elem_shape();
    const Uint cell_order      = cell_tag.poly_order();

    const mesh::PointSetTagExt cell_tag_ext(cell_tag, P0, mesh::CellTransform::NO_TRANS, 0u);

    const mesh::sf::SFTag cell_basis(cell_shape, m_sf_type, cell_order, ModalBasis::Modal);

    const mesh::PointSetTag cell_quad(cell_shape, m_quad_order, m_quad_type);
    const mesh::PointSetTagExt cell_quad_ext(cell_quad, P0, mesh::CellTransform::NO_TRANS, 0u);

    ++m_current;
    return mesh::DiscreteElemKey{cell_tag_ext, cell_basis, cell_quad_ext};
  }

  bool reached_end() const
  {
    return m_current == m_dof_range.end();
  }

  private:
  const common::IteratorRange<DofIterator> m_dof_range;
  DofIterator m_current;
  const SFunc m_sf_type;
  const PointSetID m_quad_type;
  const Uint m_quad_order;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class DiscreteElemKeyGenerator2
{
  public:
  DiscreteElemKeyGenerator2(const mesh::DofMap<MeshConfig> &dofs,
                            const common::Range1D<Uint> &cell_range, SFunc sf_type,
                            PointSetID quad_type, Uint quad_order)
      : m_dofs(dofs), m_cell_range(cell_range), m_sf_type(sf_type), m_quad_type(quad_type),
        m_quad_order(quad_order), m_current(cell_range.lbound())
  {
  }

  DiscreteElemKeyGenerator2(const DiscreteElemKeyGenerator2 &other) = delete;

  ~DiscreteElemKeyGenerator2() = default;

  DiscreteElemKeyGenerator2 &operator=(const DiscreteElemKeyGenerator2 &other) = delete;

  void init()
  {
    m_current = m_cell_range.lbound();
  }

  mesh::DiscreteElemKey make_next()
  {
    const mesh::MeshEntity cell      = m_dofs.active_cell(mesh::ActiveIdx(m_current));
    const mesh::PointSetTag cell_tag = cell.pt_set_id();

    const ElemShape cell_shape = cell_tag.elem_shape();
    const Uint cell_order      = cell_tag.poly_order();

    const mesh::PointSetTagExt cell_tag_ext(cell_tag, P0, mesh::CellTransform::NO_TRANS, 0u);

    const mesh::sf::SFTag cell_basis(cell_shape, m_sf_type, cell_order, ModalBasis::Modal);

    const mesh::PointSetTag cell_quad(cell_shape, m_quad_order, m_quad_type);
    const mesh::PointSetTagExt cell_quad_ext(cell_quad, P0, mesh::CellTransform::NO_TRANS, 0u);
    m_current++;

    return mesh::DiscreteElemKey{cell_tag_ext, cell_basis, cell_quad_ext};
  }

  bool reached_end() const
  {
    return m_current == (m_cell_range.ubound() + 1);
  }

  private:
  const mesh::DofMap<MeshConfig> &m_dofs;
  const common::Range1D<Uint> m_cell_range;
  const SFunc m_sf_type;
  const PointSetID m_quad_type;
  const Uint m_quad_order;
  Uint m_current;
};

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
