#ifndef PDEKIT_Interpolation_ElementAdaptInterpolator_hpp
#define PDEKIT_Interpolation_ElementAdaptInterpolator_hpp

#include "unordered_map"
#include "vector"

#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "mesh/DofCoordinates.hpp"
#include "mesh/adaptation/CellAdaptOp.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

class ElementAdaptInterpolator
{
  public:
  /// Default constructor
  ElementAdaptInterpolator();

  /// Deleted copy constructor
  ElementAdaptInterpolator(const ElementAdaptInterpolator &other) = delete;

  /// Default desctructor
  ~ElementAdaptInterpolator();

  /// Deleted assignment operator
  ElementAdaptInterpolator &operator=(const ElementAdaptInterpolator &rhs) = delete;

  /// Compute the coordinates of child elements after
  /// splitting given parent element
  /// @param parent_elem_type ... element type that is being split. This is to
  /// detect polynomial order
  ///                             and interpolation point distribution
  /// @param adapt_op_tag ... determines how is the element split (type of
  /// adaptation)
  /// @param cell_coordinates ... coordinates of parent element (before
  /// adaptation)
  void compute_child_values(mesh::PointSetTag const parent_elem_type,
                            mesh::adapt::CellAdaptOpTag const adapt_op_tag,
                            math::DenseDMat<Real> const &parent_values);

  /// Number of child blocks
  Uint nb_child_value_blocks() const;

  /// Get interpolated values of one child
  const math::DenseConstMatView<Real> interpolated_child_values(const Uint c) const;

  private:
  /// TYPES

  typedef std::pair<mesh::PointSetTag, mesh::adapt::CellAdaptOpTag> interpolator_key;

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
  map_iterator build_interpolation_work_matrices(mesh::PointSetTag const parent_elem_type,
                                                 interpolator_key const search_key,
                                                 mesh::adapt::CellAdaptOpTag const adapt_op_tag);

  /// Number of fields to interpolated
  Uint m_nb_fields;

  /// Number of children
  Uint m_nb_children;

  /// DATA
  /// Interpolation matrices are Vandermonde matrices of Lagrange shape
  /// functions of the parent element (including the correct polynomial order)
  /// evalated at the coordinates of each child element in the reference space
  interpolator_map m_interp_matrices;

  /// Buffer holds contatenated data values of all matrices with interpolated
  /// values
  std::vector<Real> m_adapted_coord_buffer;
  /// Buffer offsets delimit values belonging to each resulting matrix with
  /// interpolated values
  std::vector<Uint> m_buffer_offsets;
};

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
