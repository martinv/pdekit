#ifndef PDEKIT_Mesh_Transform_Algorithm_hpp
#define PDEKIT_Mesh_Transform_Algorithm_hpp

#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"

namespace pdekit
{

namespace mesh
{

class QuadratureAdaptTransformAlgoBase
{
  public:
  /// Default constructor
  QuadratureAdaptTransformAlgoBase();

  /// Default destructor
  virtual ~QuadratureAdaptTransformAlgoBase();

  static std::string type_name();

  /// Compute the new (transformed coordinates)
  /// @param input_coords  ... original quadrature coordinates
  /// @param local_id      ... id that determines which sub-entity after
  ///                          partitioning do we consider
  /// @param output_coords ... final quadrature coordinates
  virtual void compute_transformed_coords(const math::DenseDMat<Real> &input_coords,
                                          const Uint local_id,
                                          math::DenseDMat<Real> &output_coords) const = 0;

  /// Compute weights corresponding to the transformed coordinates
  /// @param input_weights  ... weights of the original quadrature prior to
  /// transformation
  /// @param local_id       ... id that determines which sub-entity after
  /// @param output_weights ... final quadrature weights

  virtual void compute_transformed_weights(const math::DenseDVec<Real> &input_weights,
                                           const Uint local_id,
                                           math::DenseDVec<Real> &output_weights) const = 0;

  private:
};

} // namespace mesh

} // namespace pdekit

#endif
