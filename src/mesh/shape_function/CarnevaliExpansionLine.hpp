#ifndef PDEKIT_Mesh_Shape_Function_Carnevali_Expansion_Line_hpp
#define PDEKIT_Mesh_Shape_Function_Carnevali_Expansion_Line_hpp

#include <cmath>

#include "common/Constants.hpp"
#include "mesh/shape_function/ModalExpansion.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

class CarnevaliExpansionLine : public ModalExpansion
{
  public:
  /// Default constructor
  CarnevaliExpansionLine();

  /// Constructor which sets the polynomial order
  CarnevaliExpansionLine(const Uint poly_order);

  /// Destructor
  ~CarnevaliExpansionLine() override;

  /// Return the number of modes
  Uint nb_modes() const override;

  /// Return the topological dimension
  Uint topo_dim() const override;

  /// Compute the values of prime basis modes at one point
  void evaluate_in_one_point(const math::DenseDVec<Real> &point,
                             math::DenseDVec<Real> &values) override;

  /// Evaluate the prime basis
  void Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                          math::DenseDMat<Real> &values) override;

  /// Compute the derivatives of the prime basis
  void evaluate_derivatives_in_one_point(const math::DenseDVec<Real> &point,
                                         const math::DenseDVec<Real> &values,
                                         math::DenseDMat<Real> &derivatives) override;

  /// Compute the derivatives of all modes in prime basis at a set of points
  void Vandermonde_matrix_derivatives(const math::DenseDMat<Real> &coordinates,
                                      std::vector<math::DenseDMat<Real>> &values) override;

  /// Determine for each mode, whether it is a leading mode, i.e. a mode that
  /// is present in expansion of order P, but NOT present in expansion of
  /// order (P-1) This makes sense for hierarchical modal bases. Leading terms
  /// cannot be defined for Lagrange bases. is_leading_term[i] is true for
  /// leading modes, false for the remaining modes
  void is_leading_expansion_term(math::DenseDVec<bool> &is_leading_term) override;

  /// For each expansion mode in the prime basis, list what is its polynomial
  /// degree
  void mode_poly_deg(math::DenseDVec<Uint> &poly_deg) override;

  private:
};

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit

#endif
