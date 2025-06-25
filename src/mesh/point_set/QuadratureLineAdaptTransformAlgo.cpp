#include "mesh/point_set/QuadratureLineAdaptTransformAlgo.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

QuadLineToTwoSegmentsTransformAlgo::QuadLineToTwoSegmentsTransformAlgo()
{
}

// ----------------------------------------------------------------------------

QuadLineToTwoSegmentsTransformAlgo::~QuadLineToTwoSegmentsTransformAlgo()
{
}

// ----------------------------------------------------------------------------

void QuadLineToTwoSegmentsTransformAlgo::compute_transformed_coords(
    const math::DenseDMat<Real> &input_coords, const Uint local_id,
    math::DenseDMat<Real> &output_coords) const
{
  output_coords.resize(input_coords.rows(), input_coords.cols());

  const Real offset = (local_id == 0) ? -0.5 : 0.5;

  if ((local_id == 0) || (local_id == 1))
  {
    for (Uint r = 0; r < input_coords.rows(); ++r)
    {
      output_coords(r, X0) = 0.5 * input_coords(r, X0) + offset;
    }
  }
}

// ----------------------------------------------------------------------------

void QuadLineToTwoSegmentsTransformAlgo::compute_transformed_weights(
    const math::DenseDVec<Real> &input_weights, const Uint local_id,
    math::DenseDVec<Real> &output_weights) const
{
  output_weights.resize(input_weights.size());

  for (Uint i = 0; i < input_weights.size(); ++i)
  {
    output_weights[i] = 0.5 * input_weights[i];
  }
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
