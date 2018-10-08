#ifndef PDEKIT_Mesh_Shape_Function_Modal_Expansion_hpp
#define PDEKIT_Mesh_Shape_Function_Modal_Expansion_hpp

#include <memory>

#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

class ModalExpansion
{
  public:
  /// TYPEDEFS:
  using ptr       = std::shared_ptr<ModalExpansion>;
  using const_ptr = std::shared_ptr<ModalExpansion const>;

  /// Default constructor
  ModalExpansion();

  /// Constructor which sets the polynomial order
  ModalExpansion(const Uint poly_order);

  /// Destructor
  virtual ~ModalExpansion() = 0;

  static std::string type_name()
  {
    return "PrimeBasis";
  }

  /// Return the polynomial order
  inline Uint polynomial_order() const
  {
    return P;
  }

  /// Set the order of the expansion
  void set_polynomial_order(const Uint poly_order);

  /// Return the number of modes
  virtual Uint nb_modes() const = 0;

  /// Return the topological dimension
  virtual Uint topo_dim() const = 0;

  /// Compute the values of prime basis modes at one point
  virtual void evaluate_in_one_point(const math::DenseDVec<Real> &point,
                                     math::DenseDVec<Real> &values) = 0;

  /// Evaluate the prime basis
  virtual void Vandermonde_matrix(const math::DenseDMat<Real> &coordinates,
                                  math::DenseDMat<Real> &values) = 0;

  /// Compute the derivatives of the prime basis
  virtual void evaluate_derivatives_in_one_point(const math::DenseDVec<Real> &point,
                                                 const math::DenseDVec<Real> &values,
                                                 math::DenseDMat<Real> &derivatives) = 0;

  /// Compute the derivatives of all modes in prime basis at a set of points
  virtual void Vandermonde_matrix_derivatives(const math::DenseDMat<Real> &coordinates,
                                              std::vector<math::DenseDMat<Real>> &values) = 0;

  /// Determine for each mode, whether it is a leading mode, i.e. a mode that
  /// is present in expansion of order P, but NOT present in expansion of
  /// order (P-1) This makes sense for hierarchical modal bases. Leading terms
  /// cannot be defined for Lagrange bases. is_leading_term[i] is true for
  /// leading modes, false for the remaining modes
  virtual void is_leading_expansion_term(math::DenseDVec<bool> &is_leading_term) = 0;

  /// For each expansion mode in the prime basis, list what is its polynomial
  /// degree
  virtual void mode_poly_deg(math::DenseDVec<Uint> &poly_deg) = 0;

  protected:
  /// Polynomial order
  Uint P;
};

} // namespace sf

} // namespace mesh

} // namespace pdekit

#endif
