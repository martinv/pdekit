#include "interpolation/ElementAdaptInterpolator.hpp"
#include "math/DenseMatView.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

ElementAdaptInterpolator::ElementAdaptInterpolator()
{
}

// ----------------------------------------------------------------------------

ElementAdaptInterpolator::~ElementAdaptInterpolator()
{
}

// ----------------------------------------------------------------------------

void ElementAdaptInterpolator::compute_child_values(mesh::PointSetTag const parent_elem_type,
                                                    mesh::adapt::CellAdaptOpTag const adapt_op_tag,
                                                    math::DenseDMat<Real> const &parent_values)
{
  interpolator_key search_key(parent_elem_type, adapt_op_tag);

  map_iterator interp_mat_it = m_interp_matrices.find(search_key);

  if (interp_mat_it == m_interp_matrices.end())
  {
    interp_mat_it = build_interpolation_work_matrices(parent_elem_type, search_key, adapt_op_tag);
  } // If interpolation matrices not found and had to be created

  std::vector<math::DenseDMat<Real>> &interp_matrices = interp_mat_it->second;

  mesh::StdRegion child_std_region;

  const mesh::adapt::CellAdaptOp adapt_operation(adapt_op_tag);
  m_nb_children = adapt_operation.get().nb_child_elems();

  m_nb_fields = parent_values.cols();
  m_buffer_offsets.resize(m_nb_children + 1);
  m_buffer_offsets[0] = 0;

  for (Uint c = 0; c < m_nb_children; ++c)
  {
    const mesh::PointSetTag child_tag(adapt_operation.get().child_elem_shape(c),
                                      parent_elem_type.poly_order(),
                                      parent_elem_type.ref_topology());
    child_std_region.change_type(child_tag);
    m_buffer_offsets[c + 1] = child_std_region.get().nb_nodes() * m_nb_fields;
  }

  for (Uint i = 1; i < m_buffer_offsets.size(); ++i)
  {
    m_buffer_offsets[i] += m_buffer_offsets[i - 1];
  }
  m_adapted_coord_buffer.resize(m_buffer_offsets.back());

  for (Uint c = 0; c < interp_matrices.size(); ++c)
  {
    const Uint nb_child_nodes = (m_buffer_offsets[c + 1] - m_buffer_offsets[c]) / m_nb_fields;
    math::DenseMatView<Real> adapted_values(m_adapted_coord_buffer.data() + m_buffer_offsets[c],
                                            m_nb_fields, nb_child_nodes, m_nb_fields);
    adapted_values = interp_matrices[c] * parent_values;
  }
}

// ----------------------------------------------------------------------------

Uint ElementAdaptInterpolator::nb_child_value_blocks() const
{
  return m_nb_children;
}

// ----------------------------------------------------------------------------

const math::DenseConstMatView<Real> ElementAdaptInterpolator::interpolated_child_values(
    const Uint c) const
{
  const Uint nb_child_nodes = (m_buffer_offsets[c + 1] - m_buffer_offsets[c]) / m_nb_fields;
  const math::DenseConstMatView<Real> adapted_values(m_adapted_coord_buffer.data() +
                                                         m_buffer_offsets[c],
                                                     m_nb_fields, nb_child_nodes, m_nb_fields);
  return adapted_values;
}

// ----------------------------------------------------------------------------

typename ElementAdaptInterpolator::map_iterator ElementAdaptInterpolator::
    build_interpolation_work_matrices(mesh::PointSetTag const parent_elem_type,
                                      interpolator_key const search_key,
                                      mesh::adapt::CellAdaptOpTag const adapt_op_tag)
{
  // ************************************************************************
  // Create the interpolation matrices that will compute values
  // on child elements from the values of parent element
  // ************************************************************************

  std::pair<map_iterator, bool> out = m_interp_matrices.insert(
      std::pair<interpolator_key, std::vector<math::DenseDMat<Real>>>(search_key, {}));

  mesh::adapt::CellAdaptOp adapt_operation(adapt_op_tag);

  // Here the child types will be changed to have the same polynomial
  // order and reference topology type as the parent:
  const Uint nb_children = adapt_operation.get().nb_child_elems();

  std::vector<mesh::StdRegion> real_child_types(nb_children);
  for (Uint c = 0; c < nb_children; ++c)
  {
    const mesh::PointSetTag child_tag(adapt_operation.get().child_elem_shape(c),
                                      parent_elem_type.poly_order(),
                                      parent_elem_type.ref_topology());

    real_child_types[c].change_type(child_tag);
  }

  // Ask the local adaptation strategy to compute child coordinates
  std::vector<math::DenseDMat<Real>> child_coords;
  adapt_operation.get().compute_child_coords(parent_elem_type, child_coords);

  mesh::sf::ShapeFunction parent_basis;
  const mesh::sf::SFTag parent_basis_tag(parent_elem_type.elem_shape(), SFunc::Lagrange,
                                         parent_elem_type.poly_order(), ModalBasis::Modal);
  parent_basis.change_type(parent_elem_type, parent_basis_tag);

  // Update the iterator: now it does not point to the 'end' of the container
  map_iterator interp_mat_it                          = out.first;
  std::vector<math::DenseDMat<Real>> &interp_matrices = interp_mat_it->second;
  interp_matrices.resize(adapt_operation.get().nb_child_elems());

  for (Uint c = 0; c < nb_children; ++c)
  {
    // child_coord[c].rows() is the number of nodes in child element
    // Resize is done automatically when the shape function computes ref.
    // values interp_matrices[c].resize(child_coords[c].rows(),
    // child_coords[c].rows());
    parent_basis.get().compute_ref_values(child_coords[c], interp_matrices[c]);
  }
  return interp_mat_it;
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit
