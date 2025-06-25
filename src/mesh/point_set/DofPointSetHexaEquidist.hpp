#ifndef PDEKIT_Mesh_Dof_Point_Set_Hexa_Equidist_hpp
#define PDEKIT_Mesh_Dof_Point_Set_Hexa_Equidist_hpp

#include "mesh/point_set/StdPointSetBase.hpp"
#include "mesh/std_region/EquidistStdRegionHexa.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Equidistant point set for hexahedra
// ----------------------------------------------------------------------------

class DofPointSetHexaEquidist : public StdPointSetBase
{
  public:
  /// Default constructor
  DofPointSetHexaEquidist(const std::tuple<Uint> &p_order);

  /// Destructor
  ~DofPointSetHexaEquidist() override = default;

  static std::string type_name()
  {
    return "DofPointSetHexaEquidist";
  }

  std::string name() const override
  {
    return "DofPointSetHexaEquidist";
  }

  /// Order of polynomial which this quadrature integrates exactly
  Uint order() const override;

  /// Topological dimension of element for which this quadrature can be
  /// applied
  Uint dim() const override;

  /// Topological codimension of quadrature.
  Uint codim() const override;

  /// Return the number of local entities on which this quadrature
  /// has points defined
  Uint nb_local_entities() const override;

  /// Return the number of quadrature points
  Uint size(const Uint local_idx = 0) const override;

  /// Return a matrix containing the reference coordinates
  void reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0) const override;

  /// Compute weights
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
  using std_reg_t = mesh::EquidistStdRegionHexa;

  const Uint m_poly_order;
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
