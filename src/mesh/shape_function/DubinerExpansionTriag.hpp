#ifndef PDEKIT_Mesh_Shape_Function_Orthonormal_Prime_Basis_Triag_hpp
#define PDEKIT_Mesh_Shape_Function_Orthonormal_Prime_Basis_Triag_hpp

#include <cmath>

#include "common/Constants.hpp"
#include "mesh/shape_function/ModalExpansion.hpp"

#define USE_SINGULARITY_FREE_EVALUATION_ON_TRIAG 1

namespace pdekit
{

namespace mesh
{

namespace sf
{

class DubinerExpansionTriag : public ModalExpansion
{
  public:
  /// Default constructor
  DubinerExpansionTriag();

  /// Constructor which sets the polynomial order
  DubinerExpansionTriag(const Uint poly_order);

  /// Destructor
  ~DubinerExpansionTriag() override;

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
  /// The modes indexed by (p,q) are put in an array - we need to assign
  /// consequent indexes to each pair (p,q). Note that the modes
  /// are such that p+q <= P, i.e. if we placed them matrix, only one
  /// triangle (lower or upper) + the diagonal would be used
  inline static Uint mode_id(const Uint p, const Uint q)
  {
    // return p*(P+1)+q - p*(p-1)/2;
    return (p + q) * (p + q + 1) / 2 + q;
  }

  /// Coefficient in the recursive formula for evaluation of Dubiner modes
  /// based on the paper of R. C. Kirby 'Singularity-free evaluation of
  /// collapsed-coordinate orthogonal polynomials'
  inline static Real a(const Uint n, const Real alpha, const Real beta)
  {
    return 0.5 * (2. * n + 1 + alpha + beta) * (2. * n + 2. + alpha + beta) /
           ((n + 1) * (n + 1. + alpha + beta));
  }

  /// Coefficient in the recursive formula for evaluation of Dubiner modes
  inline static Real b(const Uint n, const Real alpha, const Real beta)
  {
    return 0.5 * ((alpha * alpha - beta * beta) * (2. * n + 1 + alpha + beta)) /
           ((n + 1) * (2. * n + alpha + beta) * (n + 1. + alpha + beta));
  }

  /// Coefficient in the recursive formula for evaluation of Dubiner modes
  inline static Real c(const Uint n, const Real alpha, const Real beta)
  {
    return ((n + alpha) * (n + beta) * (2 * n + 2. + alpha + beta)) /
           ((n + 1) * (2. * n + alpha + beta) * (n + 1. + alpha + beta));
  }
};

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit

#endif
