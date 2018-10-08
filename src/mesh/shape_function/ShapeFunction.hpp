#ifndef PDEKIT_Interpolation_Shape_Function_Instance_hpp
#define PDEKIT_Interpolation_Shape_Function_Instance_hpp

#include <memory>
#include <string>
#include <vector>

#include "common/Constants.hpp"
#include "common/Flyweight.hpp"

#include "math/DenseConstMatView.hpp"
#include "math/DenseConstVecView.hpp"

#include "mesh/shape_function/ModalExpansion.hpp"
#include "mesh/shape_function/SFTag.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

namespace detail
{

class ShapeFunctionInstance
{
  public:
  /// TYPEDEFS:
  using ptr       = std::shared_ptr<ShapeFunctionInstance>;
  using const_ptr = std::shared_ptr<ShapeFunctionInstance const>;
  using coord_t   = math::DenseDMat<Real>;

  /// Tuple to identify the shape function uniquely
  /// First entry:  element shape (e.g Line)
  /// Second entry: interpolant type (e.g Lagrange)
  /// Third entry:  interpolation point set (e.g Equidistant)
  /// Fourth entry: polynomial order
  /// Fifth entry:  prime basis type (e.g Modal)
  ///
  static const std::tuple<PointSetTag, SFTag> undefined;

  // typedef common::Tag5<ElemShape,SFunc,RefTopology,PolyOrder,PrimeBasis>
  // tag_type;

  /// Static method to fill a reference element
  static void construct(const std::tuple<PointSetTag, SFTag> sf_key,
                        ShapeFunctionInstance &sf_instance);

  /// Default constructor
  ShapeFunctionInstance();

  /// Constructor from shape function tag
  ShapeFunctionInstance(const std::tuple<PointSetTag, SFTag> sf_key, const SFTag sf_tag);

  /// Destructor
  ~ShapeFunctionInstance();

  /// Get the tag of this sf
  const SFTag &tag() const;

  static std::string type_name()
  {
    return "ShapeFunctionInstance";
  }

  /// Return the number of nodes (degrees of freedom) for this sf
  Uint nb_dof() const;

  /// Get the topological dimension of this shape function
  Uint topo_dim() const;

  /// Fill a matrix containing the reference coordinates
  const coord_t &ref_coords() const;

  /// Fill a matrix containing the reference coordinates
  void ref_coords(coord_t &coord) const;

  /// Compute the values of all shape functions at one point:
  void compute_ref_values(const math::DenseDMat<Real> &ref_coords,
                          math::DenseDMat<Real> &values) const;

  /// Compute the derivatives of all shape functions at one point:
  void compute_ref_derivatives(const math::DenseDMat<Real> &ref_coords,
                               std::vector<math::DenseDMat<Real>> &values) const;

  /// Compute the values of shape functions at a series of coordinates
  /// Store only values of shape functions associated with certain local
  /// entity of given standard region
  /// @param ref_coords      ... reference coordinates
  /// @param dim             ... dimension of sub-entity
  /// @param local_entity_id ... index of sub-entity
  /// @param values          ... matrix of computed values (output)
  void compute_ref_values(const math::DenseDMat<Real> &ref_coords, const Uint dim,
                          const Uint local_entity_id, math::DenseDMat<Real> &values) const;

  /// Compute the derivatives of shape functions at a series of coordinates
  /// Store only values of shape function derivatives associated with certain
  /// local entity of given standard region
  /// @param ref_coords      ... reference coordinates
  /// @param dim             ... dimension of sub-entity
  /// @param local_entity_id ... index of sub-entity
  /// @param values          ... vector of matrices of computed derivatives
  /// (output)
  ///                            each matrix holds derivatives with respect to
  ///                            one coordinate dimension
  void compute_ref_derivatives(const math::DenseDMat<Real> &ref_coords, const Uint dim,
                               const Uint local_entity_id,
                               std::vector<math::DenseDMat<Real>> &values) const;

  /// Return the information about which modes have the highest order in the
  /// expansion set
  /// @note The data is return in ConstVectorBlock to avoid returning
  ///       reference or pointer to member data
  const math::DenseConstVecView<bool> is_leading_expansion_term() const;

  /// Return the information about polynomial degree of each mode
  const math::DenseConstVecView<Uint> mode_poly_deg() const;

  private:
  /// Tuple to identify the shape function uniquely
  SFTag m_tag;

  /// Vandermonde matrix: prime basis evaluated at interpolation point set
  math::DenseDMat<Real> m_V;

  /// Inverse of Vandermonde matrix
  math::DenseDMat<Real> m_invV;

  /// Vector of boolean values indicating for each mode in the expansion
  /// whether the mode has the highest polynomial order in this shape
  /// function set (m_is_leading_expansion_term[i] == true)
  /// or not (m_is_leading_expantion_term[i] == false)
  math::DenseDVec<bool> m_is_leading_expansion_term;

  /// Polynomial degree of each mode in expansion (if this is modal
  /// expansion). If this is Lagrange set, then m_mode_poly_deg = P, where P
  /// is the polynomial degree of this shape function set (i.e. max P)
  math::DenseDVec<Uint> m_mode_poly_deg;

  /// Point coordinates in reference space
  coord_t m_prime_coords;

  /// Pointer to the prime basis
  ModalExpansion::ptr m_prime_basis;
};

struct ShapeFunctionFlyweightPolicy
{
  using key_type = std::tuple<PointSetTag, SFTag>;
};

} // namespace detail

using ShapeFunction =
    common::Flyweight<detail::ShapeFunctionInstance, detail::ShapeFunctionFlyweightPolicy>;

} // namespace sf

} // namespace mesh

} // namespace pdekit

#endif
