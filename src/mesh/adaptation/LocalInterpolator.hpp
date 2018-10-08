#ifndef PDEKIT_Mesh_Adaptation_Local_Interpolator_hpp
#define PDEKIT_Mesh_Adaptation_Local_Interpolator_hpp

#include "unordered_map"
#include "vector"

#include "math/DenseConstMatView.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseMatView.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

template <Uint GeoDim>
class CellGeometry;

namespace adapt
{

// ----------------------------------------------------------------------------

class LocalInterpolator
{
  public:
  /// Default constructor
  LocalInterpolator() = default;

  /// Deleted copy constructor
  LocalInterpolator(const LocalInterpolator &other) = delete;

  /// Default desctructor
  ~LocalInterpolator() = default;

  /// Deleted assignment operator
  LocalInterpolator &operator=(const LocalInterpolator &other) = delete;

  /// Compute the values on element after interpolation
  /// @param std_reg_tag_in  ... element type that on input. We expect the
  ///                            data to be defined on this element type
  /// @param std_reg_tag_out ... element type on output. The data will be
  ///                            interpolated ONTO this element
  /// @param data_in         ... values to be interpolated. Each column represents
  ///                            a set of (scalar) values defined on the input
  ///                            element.
  /// @return                ... interpolated values
  const math::DenseConstMatView<Real> transfer_data(PointSetTag const std_reg_tag_in,
                                                    PointSetTag const std_reg_tag_out,
                                                    const math::DenseDMat<Real> &data_in);

  /// Compute the values on element after interpolation
  /// @param std_reg_tag_in  ... element type that on input. We expect the data
  ///                            to be defined on this element type
  /// @param std_reg_tag_out ... element type on output. The data will be
  ///                            interpolated ONTO this element
  /// @param data_in         ... values to be interpolated. Each column represents
  ///                            a set of (scalar) values defined on the input
  ///                            element.
  /// @return                ... interpolated values
  const math::DenseConstMatView<Real> transfer_data(PointSetTag const std_reg_tag_in,
                                                    PointSetTag const std_reg_tag_out,
                                                    const math::DenseConstMatView<Real> &data_in);

  /// Compute the values on element after interpolation
  /// @param std_reg_tag_in  ... element type that on input. We expect the
  ///                            data to be defined on this element type
  /// @param std_reg_tag_out ... element type on output. The data will be
  ///                            interpolated ONTO this element
  /// @param coords_in       ... values to be interpolated. Each column represents
  ///                            a set of (scalar) values defined on the input element.
  /// @return                ... interpolated values
  template <Uint GeoDim>
  const math::DenseConstMatView<Real> transfer_coords(PointSetTag const std_reg_tag_in,
                                                      PointSetTag const std_reg_tag_out,
                                                      const mesh::CellGeometry<GeoDim> &coords_in);

  private:
  /// TYPES
  using interpolator_key = std::pair<PointSetTag, PointSetTag>;

  struct InterpolatorKeyHasher
  {
    inline std::size_t operator()(const interpolator_key &key) const
    {
      return 1000 * key.first.store_value() ^ key.second.store_value();
    }
  };

  using interpolator_map =
      std::unordered_map<interpolator_key, math::DenseDMat<Real>, InterpolatorKeyHasher>;

  using map_iterator       = typename interpolator_map::iterator;
  using const_map_iterator = typename interpolator_map::const_iterator;

  /// METHODS
  map_iterator build_interpolation_work_matrix(interpolator_key const &search_key);

  /// DATA
  /// Interpolation matrices are Vandermonde matrices of Lagrange shape
  /// functions that transfer values from one point set to another
  interpolator_map m_interp_matrices;

  /// Work data to temporarily store the interpolated values
  std::vector<Real> m_work_data_out;
};

// ----------------------------------------------------------------------------

template <Uint GeoDim>
const math::DenseConstMatView<Real> LocalInterpolator::transfer_coords(
    PointSetTag const std_reg_tag_in, PointSetTag const std_reg_tag_out,
    const mesh::CellGeometry<GeoDim> &coords_in)
{
  interpolator_key search_key(std_reg_tag_in, std_reg_tag_out);
  map_iterator interp_mat_it = m_interp_matrices.find(search_key);

  if (interp_mat_it == m_interp_matrices.end())
  {
    interp_mat_it = build_interpolation_work_matrix(search_key);
  } // If interpolation matrices not found and had to be created

  const math::DenseDMat<Real> &I = interp_mat_it->second;
  const Uint nb_nodes_in         = coords_in.size();
  const Uint nb_nodes_out        = I.rows();

  m_work_data_out.resize(nb_nodes_out * GeoDim);
  m_work_data_out.assign(m_work_data_out.size(), 0.0);
  math::DenseMatView<Real> data_out(m_work_data_out.data(), GeoDim, nb_nodes_out, GeoDim);

  if (std_reg_tag_in == std_reg_tag_out)
  {
    for (Uint r = 0; r < nb_nodes_out; ++r)
    {
      const auto one_node_coord = coords_in.const_node_view(r);

      for (Uint dim = 0; dim < GeoDim; ++dim)
      {
        data_out(r, dim) = one_node_coord[dim];
      }
    }
  }
  else
  {
    for (Uint r = 0; r < nb_nodes_out; ++r)
    {
      for (Uint c = 0; c < nb_nodes_in; ++c)
      {
        const auto one_node_coord = coords_in.const_node_view(c);

        for (Uint dim = 0; dim < GeoDim; ++dim)
        {
          data_out(r, dim) += I(r, c) * one_node_coord[dim];
        }
      }
    }
  }

  math::DenseConstMatView<Real> data_out_const(m_work_data_out.data(), GeoDim, I.rows(), GeoDim);
  return data_out_const;
}

// ----------------------------------------------------------------------------

} // namespace adapt

} // namespace mesh

} // namespace pdekit

#endif
