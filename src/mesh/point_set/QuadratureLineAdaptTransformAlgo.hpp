#ifndef PDEKIT_Mesh_Line_Adapt_Transform_Algorithm_hpp
#define PDEKIT_Mesh_Line_Adapt_Transform_Algorithm_hpp

#include "mesh/point_set/QuadratureAdaptTransformAlgoBase.hpp"

namespace pdekit
{

namespace mesh
{

class QuadLineToTwoSegmentsTransformAlgo : public QuadratureAdaptTransformAlgoBase
{
  public:
  /// Default constructor
  QuadLineToTwoSegmentsTransformAlgo();

  /// Default destructor
  ~QuadLineToTwoSegmentsTransformAlgo() override;

  /// Compute the new (transformed coordinates)
  /// @param input_coords  ... original quadrature coordinates
  /// @param local_id      ... id that determines which sub-entity after
  ///                          partitioning do we consider
  /// @param output_coords ... final quadrature coordinates
  void compute_transformed_coords(const math::DenseDMat<Real> &input_coords, const Uint local_id,
                                  math::DenseDMat<Real> &output_coords) const override;

  /// Compute weights corresponding to the transformed coordinates
  /// @param input_weights  ... weights of the original quadrature prior to
  /// transformation
  /// @param local_id       ... id that determines which sub-entity after
  /// @param output_weights ... final quadrature weights

  void compute_transformed_weights(const math::DenseDVec<Real> &input_weights, const Uint local_id,
                                   math::DenseDVec<Real> &output_weights) const override;

  private:
};

} // namespace mesh

} // namespace pdekit

#endif
