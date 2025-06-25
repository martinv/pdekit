#ifndef PDEKIT_Interpolation_FE_Values_hpp
#define PDEKIT_Interpolation_FE_Values_hpp

#include <memory>

#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "mesh/shape_function/SFTag.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace interpolation
{

class FEValues
{
  public:
  /// TYPEDEFS
  typedef std::shared_ptr<FEValues> ptr;
  typedef std::shared_ptr<FEValues const> const_ptr;

  typedef math::DenseDMat<Real> coord_t;

  /// METHODS
  static void copy(FEValues const &fe_source, FEValues &fe_target);

  /// Default constructor
  FEValues();

  /// Disabled copy constructor
  FEValues(const FEValues &other) = delete;

  /// Constructor
  FEValues(const mesh::PointSetTag std_region_tag, const mesh::sf::SFTag &sf_type);

  /// Destructor
  ~FEValues();

  /// Disabled assignment operator
  FEValues &operator=(const FEValues &rhs) = delete;

  /// Configure the FEValues
  void configure(const mesh::PointSetTag std_region_tag, const mesh::sf::SFTag &sf_type);

  /// Return the element type
  mesh::PointSetTag std_region_id() const;

  /// Return the shape function type
  mesh::sf::SFTag sf_type();

  /// Return the element type, const version
  mesh::sf::SFTag sf_type() const;

  /// Return the topological dimension
  Uint topo_dim() const;

  /// Return the number of nodes
  Uint nb_nodes() const;

  /// Return the number of quadrature points
  Uint nb_qd_pts() const;

  /// Return reference to the matrix of reference coordinates,
  /// of the element, const version. These are coords of the element dofs
  const coord_t &ref_coord() const;

  /// Fill the Vandermonde matrix of shape function values and shape function
  /// derivatives using quadrature of given order and type
  void fill_Vandermonde(const math::DenseDMat<Real> &point_set,
                        const math::DenseDVec<Real> &weights, const bool apply_filter = false);

  /// Return the Vandermonde matrix, ONLY const version
  math::DenseDMat<Real> const &Vandermonde() const;

  /// Return the Vandermonde matrix of partial derivatives
  /// ONLY const version
  math::DenseDMat<Real> const &deriv_Vandermonde(const Uint variable) const;

  /// Return the vector of quadrature weights, ONLY const version
  const math::DenseDVec<Real> &qw() const;

  /// Return the matrix of quadrature coordinates, ONLY const version
  const math::DenseDMat<Real> &qp() const;

  /// Print some information for debugging
  void print() const;

  private:
  /// METHODS
  void setup_type();

  /// DATA

  /// Element type
  mesh::PointSetTag m_std_region_tag;

  /// Shape function type
  mesh::sf::SFTag m_sf_type;

  /// Topological dimension (in reference space)
  Uint m_topo_dim;

  /// Coordinates of element in reference space
  coord_t m_ref_coord;

  /// Vandermonde Matrix
  math::DenseDMat<Real> m_V;

  /// Vandermonde matrices of derivatives
  std::vector<math::DenseDMat<Real>> m_dV;

  /// Coordinates of points in which the Vandermonde matrix
  /// was computed
  math::DenseDMat<Real> m_point_set;

  /// Weights of points in which the Vandermonde matrix and
  /// derivatives were evaluated
  math::DenseDVec<Real> m_pt_weights;
};

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
