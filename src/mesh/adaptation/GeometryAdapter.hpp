#ifndef PDEKIT_Mesh_Adaptation_Geometry_Adapter_hpp
#define PDEKIT_Mesh_Adaptation_Geometry_Adapter_hpp

#include "unordered_map"
#include "vector"

#include "math/DenseDMat.hpp"
#include "mesh/CellGeometry.hpp"
#include "mesh/DofCoordinates.hpp"
#include "mesh/adaptation/CellAdaptOp.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class GeometryAdapter
{
  public:
  /// Default constructor
  GeometryAdapter();

  /// Deleted copy constructor
  GeometryAdapter(const GeometryAdapter &other_adapter) = delete;

  /// Default desctructor
  ~GeometryAdapter();

  /// Deleted assignment operator
  GeometryAdapter &operator=(const GeometryAdapter &other_adapter) = delete;

  /// Compute the coordinates of child elements after
  /// splitting given parent element
  /// @param parent_elem_type ... element type that is being split. This is to
  ///                             detect polynomial order and interpolation
  ///                             point distribution
  /// @param adapt_op_tag ... determines how is the element split (type of
  /// adaptation)
  /// @param cell_coordinates ... coordinates of parent element (before
  /// adaptation)
  const common::ArrayView<const math::DenseDMat<Real>, _1D, Uint> compute_child_coords(
      PointSetTag const parent_elem_type, CellAdaptOpTag const adapt_op_tag,
      DofCoordinates<MeshConfig::GDIM> const &cell_coords);

  /// Compute the coordinates of child elements after
  /// splitting given parent element
  /// @param parent_elem_type ... element type that is being split. This is to
  ///                             detect polynomial order and interpolation
  ///                             point distribution
  /// @param adapt_op_tag ... determines how is the element split (type of
  /// adaptation)
  /// @param cell_coordinates ... coordinates of parent element (before
  /// adaptation)
  const common::ArrayView<const math::DenseDMat<Real>, _1D, Uint> compute_child_coords(
      PointSetTag const parent_elem_type, CellAdaptOpTag const adapt_op_tag,
      mesh::CellGeometry<MeshConfig::GDIM> const &cell_coords);

  private:
  /// TYPES

  typedef std::pair<PointSetTag, CellAdaptOpTag> interpolator_key;

  struct AdapterMapKeyHasher
  {
    inline std::size_t operator()(const interpolator_key &key) const
    {
      return key.first.store_value() ^ key.second.store_value();
    }
  };

  typedef std::unordered_map<interpolator_key, std::vector<math::DenseDMat<Real>>,
                             AdapterMapKeyHasher>
      interpolator_map;

  typedef typename interpolator_map::iterator map_iterator;
  typedef typename interpolator_map::const_iterator const_map_iterator;

  /// METHODS

  map_iterator build_interpolation_work_matrices(PointSetTag const parent_elem_type,
                                                 interpolator_key const search_key,
                                                 CellAdaptOpTag const adapt_op_tag,
                                                 const Uint nb_nodes);

  /// DATA
  /// Interpolation matrices are Vandermonde matrices of Lagrange shape
  /// functions of the parent element (including the correct polynomial order)
  /// evalated at the coordinates of each child element in the reference space
  interpolator_map m_interp_matrices;

  /// Convention: the first n matrices in the std::vector hold the coordinates
  /// of the children, the last matrix has the coordinates of the original
  /// (parent) element
  interpolator_map m_adapted_coord_buffer;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
GeometryAdapter<MeshConfig>::GeometryAdapter()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
GeometryAdapter<MeshConfig>::~GeometryAdapter()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const common::ArrayView<const math::DenseDMat<Real>, _1D, Uint> GeometryAdapter<
    MeshConfig>::compute_child_coords(PointSetTag const parent_elem_type,
                                      CellAdaptOpTag const adapt_op_tag,
                                      DofCoordinates<MeshConfig::GDIM> const &cell_coords)
{
  interpolator_key search_key(parent_elem_type, adapt_op_tag);

  map_iterator interp_mat_it = m_interp_matrices.find(search_key);

  if (interp_mat_it == m_interp_matrices.end())
  {
    interp_mat_it = build_interpolation_work_matrices(parent_elem_type, search_key, adapt_op_tag,
                                                      cell_coords.size());
  } // If interpolation matrices not found and had to be created

  std::vector<math::DenseDMat<Real>> &interp_matrices = interp_mat_it->second;
  map_iterator adapted_coord_it                       = m_adapted_coord_buffer.find(search_key);
  std::vector<math::DenseDMat<Real>> &adapted_coord_matrices = adapted_coord_it->second;

  math::DenseDMat<Real> &parent_coords = adapted_coord_matrices.back();

  for (Uint n = 0; n < cell_coords.size(); ++n)
  {
    parent_coords.insert_row(n, cell_coords.c(n));
  }

  for (Uint c = 0; c < interp_matrices.size(); ++c)
  {
    adapted_coord_matrices[c] = interp_matrices[c] * parent_coords;
  }

  // Return a proxy containing the computed child coordinates. Do not
  // include the last matrix, because that matrix contains the coordinates
  // of the parent cell
  const common::ArrayView<const math::DenseDMat<Real>, _1D, Uint> child_coords(
      adapted_coord_matrices.data(), adapted_coord_matrices.size() - 1);
  return child_coords;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const common::ArrayView<const math::DenseDMat<Real>, _1D, Uint> GeometryAdapter<
    MeshConfig>::compute_child_coords(PointSetTag const parent_elem_type,
                                      CellAdaptOpTag const adapt_op_tag,
                                      mesh::CellGeometry<MeshConfig::GDIM> const &cell_coords)
{
  interpolator_key search_key(parent_elem_type, adapt_op_tag);

  map_iterator interp_mat_it = m_interp_matrices.find(search_key);

  const mesh::StdRegion parent_std_reg(parent_elem_type);
  const Uint nb_nodes = parent_std_reg.get().nb_nodes();

  if (interp_mat_it == m_interp_matrices.end())
  {
    interp_mat_it =
        build_interpolation_work_matrices(parent_elem_type, search_key, adapt_op_tag, nb_nodes);
  } // If interpolation matrices not found and had to be created

  std::vector<math::DenseDMat<Real>> &interp_matrices = interp_mat_it->second;
  map_iterator adapted_coord_it                       = m_adapted_coord_buffer.find(search_key);
  std::vector<math::DenseDMat<Real>> &adapted_coord_matrices = adapted_coord_it->second;

  math::DenseDMat<Real> &parent_coords = adapted_coord_matrices.back();

  for (Uint n = 0; n < nb_nodes; ++n)
  {
    typename mesh::CellGeometry<MeshConfig::GDIM>::node_view_t node_c =
        cell_coords.const_node_view(n);
    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      parent_coords(n, d) = node_c[d];
    }
  }

  for (Uint c = 0; c < interp_matrices.size(); ++c)
  {
    adapted_coord_matrices[c] = interp_matrices[c] * parent_coords;
  }

  // Return a proxy containing the computed child coordinates. Do not
  // include the last matrix, because that matrix contains the coordinates
  // of the parent cell
  const common::ArrayView<const math::DenseDMat<Real>, _1D, Uint> child_coords(
      adapted_coord_matrices.data(), adapted_coord_matrices.size() - 1);
  return child_coords;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename GeometryAdapter<MeshConfig>::map_iterator GeometryAdapter<
    MeshConfig>::build_interpolation_work_matrices(PointSetTag const parent_elem_type,
                                                   interpolator_key const search_key,
                                                   CellAdaptOpTag const adapt_op_tag,
                                                   const Uint nb_nodes)
{
  // ************************************************************************
  // 1) Create the interpolation matrices that will compute coordinates
  //    of child elements from the coordinates of parent element
  // ************************************************************************

  std::pair<map_iterator, bool> out = m_interp_matrices.insert(
      std::pair<interpolator_key, std::vector<math::DenseDMat<Real>>>(search_key, {}));

  adapt::CellAdaptOp split_strategy(adapt_op_tag);

  // Here the child types will be changed to have the same polynomial
  // order and reference topology type as the parent:
  const Uint nb_children = split_strategy.get().nb_child_elems();

  std::vector<StdRegion> real_child_types(nb_children);
  for (Uint c = 0; c < nb_children; ++c)
  {
    const PointSetTag child_tag(split_strategy.get().child_elem_shape(c),
                                parent_elem_type.poly_order(), parent_elem_type.ref_topology());

    real_child_types[c].change_type(child_tag);
  }

  // Ask the local adaptation strategy to compute child coordinates
  std::vector<math::DenseDMat<Real>> child_coords;
  split_strategy.get().compute_child_coords(parent_elem_type, child_coords);

  const sf::SFTag parent_sf_tag(parent_elem_type.elem_shape(), SFunc::Lagrange,
                                parent_elem_type.poly_order(), ModalBasis::Modal);
  sf::ShapeFunction shape_func;
  shape_func.change_type(parent_elem_type, parent_sf_tag);

  // Update the iterator: now it does not point to the 'end' of the container
  map_iterator interp_mat_it                          = out.first;
  std::vector<math::DenseDMat<Real>> &interp_matrices = interp_mat_it->second;
  interp_matrices.resize(split_strategy.get().nb_child_elems());

  for (Uint c = 0; c < interp_matrices.size(); ++c)
  {
    // child_coord[c].rows() is the number of nodes in child element
    // Resize is done automatically when the shape function computes ref.
    // values interp_matrices[c].resize(child_coords[c].rows(),
    // child_coords[c].rows());
    shape_func.get().compute_ref_values(child_coords[c], interp_matrices[c]);
  }

  // ************************************************************************
  // 2) Create the matrices that will temporarily hold the coordinates
  //    of parent element and the computed coordinates of its children
  //    'Temporarily' means that every time this method is called, the
  //    matrices corresponding to given element type and adaptation
  //    operation in m_adapted_coords will be overwritten
  // ************************************************************************

  std::pair<map_iterator, bool> out2 = m_adapted_coord_buffer.insert(
      std::pair<interpolator_key, std::vector<math::DenseDMat<Real>>>(search_key, {}));

  map_iterator adapted_coord_it                              = out2.first;
  std::vector<math::DenseDMat<Real>> &adapted_coord_matrices = adapted_coord_it->second;
  // Convention: the first n matrices hold the coordinates of the children,
  // the last matrix has the coordinates of the original (parent) element
  adapted_coord_matrices.resize(split_strategy.get().nb_child_elems() + 1);

  for (Uint c = 0; c < real_child_types.size(); ++c)
  {
    adapted_coord_matrices[c].resize(real_child_types[c].get().nb_nodes(), MeshConfig::GDIM);
  }

  adapted_coord_matrices.back().resize(nb_nodes, MeshConfig::GDIM);

  return interp_mat_it;
}

// ----------------------------------------------------------------------------

} // namespace adapt

} // namespace mesh

} // namespace pdekit

#endif
