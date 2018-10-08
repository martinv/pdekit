#ifndef PDEKIT_Mesh_Qd_Point_Set_Tetra_Gauss_hpp
#define PDEKIT_Mesh_Qd_Point_Set_Tetra_Gauss_hpp

#include "mesh/point_set/StdPointSetBase.hpp"

namespace pdekit
{

namespace mesh
{

// ============================================================================
// Quadrature for P1 polynomials on reference tetrahedron
// ============================================================================

class QdPointSetP1TetraGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP1TetraGauss();

  /// Destructor
  ~QdPointSetP1TetraGauss() override;

  static std::string type_name()
  {
    return "Tetra-P1-Gauss";
  }

  std::string name() const override
  {
    return "Tetra-P1-Gauss";
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

  /// Return a vector containing the reference coordinates
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of quadrature points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
};

// ============================================================================
// Quadrature for P2 polynomials on reference tetrahedron
// ============================================================================

class QdPointSetP2TetraGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP2TetraGauss();

  /// Destructor
  ~QdPointSetP2TetraGauss() override;

  static std::string type_name()
  {
    return "Tetra-P2-Gauss";
  }

  std::string name() const override
  {
    return "Tetra-P2-Gauss";
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

  /// Return a vector containing the reference coordinates
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of quadrature points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
};

// ============================================================================
// Quadrature for P3 polynomials on reference tetrahedron
// ============================================================================

class QdPointSetP3TetraGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP3TetraGauss();

  /// Destructor
  ~QdPointSetP3TetraGauss() override;

  static std::string type_name()
  {
    return "Tetra-P3-Gauss";
  }

  std::string name() const override
  {
    return "Tetra-P3-Gauss";
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

  /// Return a vector containing the reference coordinates
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of quadrature points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
};

// ============================================================================
// Quadrature for P4 polynomials on reference tetrahedron
// ============================================================================

class QdPointSetP4TetraGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP4TetraGauss();

  /// Destructor
  ~QdPointSetP4TetraGauss() override;

  static std::string type_name()
  {
    return "Tetra-P4-Gauss";
  }

  std::string name() const override
  {
    return "Tetra-P4-Gauss";
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

  /// Return a vector containing the reference coordinates
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of quadrature points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
};

// ============================================================================
// Quadrature for P5 polynomials on reference tetrahedron
// ============================================================================

class QdPointSetP5TetraGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP5TetraGauss();

  /// Destructor
  ~QdPointSetP5TetraGauss() override;

  static std::string type_name()
  {
    return "Tetra-P5-Gauss";
  }

  std::string name() const override
  {
    return "Tetra-P5-Gauss";
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

  /// Return a vector containing the reference coordinates
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of quadrature points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
};

// ============================================================================
// Quadrature for P6 polynomials on reference tetrahedron
// ============================================================================

class QdPointSetP6TetraGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP6TetraGauss();

  /// Destructor
  ~QdPointSetP6TetraGauss() override;

  static std::string type_name()
  {
    return "Tetra-P6-Gauss";
  }

  std::string name() const override
  {
    return "Tetra-P6-Gauss";
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

  /// Return a vector containing the reference coordinates
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of quadrature points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
};

// ============================================================================
// Quadrature for P7 polynomials on reference tetrahedron
// ============================================================================

class QdPointSetP7TetraGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP7TetraGauss();

  /// Destructor
  ~QdPointSetP7TetraGauss() override;

  static std::string type_name()
  {
    return "Tetra-P7-Gauss";
  }

  std::string name() const override
  {
    return "Tetra-P7-Gauss";
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

  /// Return a vector containing the reference coordinates
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of quadrature points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
};

// ============================================================================
// Quadrature for P8 polynomials on reference tetrahedron
// ============================================================================

class QdPointSetP8TetraGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP8TetraGauss();

  /// Destructor
  ~QdPointSetP8TetraGauss() override;

  static std::string type_name()
  {
    return "Tetra-P8-Gauss";
  }

  std::string name() const override
  {
    return "Tetra-P8-Gauss";
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

  /// Return a vector containing the reference coordinates
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of quadrature points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
};

// ============================================================================
// Quadrature for P9 polynomials on reference tetrahedron
// ============================================================================

// 5 negative weights, 12 points outside of the tetrahedron

class QdPointSetP9TetraGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP9TetraGauss();

  /// Destructor
  ~QdPointSetP9TetraGauss() override;

  static std::string type_name()
  {
    return "Tetra-P9-Gauss";
  }

  std::string name() const override
  {
    return "Tetra-P9-Gauss";
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

  /// Return a vector containing the reference coordinates
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of quadrature points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
};

// ============================================================================
// Quadrature for P10 polynomials on reference tetrahedron
// ============================================================================

// This is actually P11 tet. quadrature from Solin et. al used for P10
// polynomials

class QdPointSetP10TetraGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP10TetraGauss();

  /// Destructor
  ~QdPointSetP10TetraGauss() override;

  static std::string type_name()
  {
    return "Tetra-P10-Gauss";
  }

  std::string name() const override
  {
    return "Tetra-P10-Gauss";
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

  /// Return a vector containing the reference coordinates
  void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const override;

  /// Fill a vector which represents a permutation of quadrature points
  void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                   std::vector<Uint> &permutation_vec) override;

  private:
};

// ============================================================================

} // namespace mesh

} // namespace pdekit

#endif
