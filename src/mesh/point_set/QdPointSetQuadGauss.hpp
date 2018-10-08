#ifndef PDEKIT_Mesh_Qd_Point_Set_Quad_Gauss_hpp
#define PDEKIT_Mesh_Qd_Point_Set_Quad_Gauss_hpp

#include "mesh/point_set/StdPointSetBase.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Quad quadrature for P1 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP1QuadGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP1QuadGauss();

  /// Destructor
  ~QdPointSetP1QuadGauss() override;

  static std::string type_name()
  {
    return "Quad-P1-Gauss";
  }

  std::string name() const override
  {
    return "Quad-P1-Gauss";
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

// ----------------------------------------------------------------------------
// Quad quadrature for P2 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP2QuadGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP2QuadGauss();

  /// Destructor
  ~QdPointSetP2QuadGauss() override;

  static std::string type_name()
  {
    return "Quad-P2-Gauss";
  }

  std::string name() const override
  {
    return "Quad-P2-Gauss";
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

// ----------------------------------------------------------------------------
// Quad quadrature for P3 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP3QuadGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP3QuadGauss();

  /// Destructor
  ~QdPointSetP3QuadGauss() override;

  static std::string type_name()
  {
    return "Quad-P3-Gauss";
  }

  std::string name() const override
  {
    return "Quad-P3-Gauss";
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

// ----------------------------------------------------------------------------
// Quad quadrature for P4 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP4QuadGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP4QuadGauss();

  /// Destructor
  ~QdPointSetP4QuadGauss() override;

  static std::string type_name()
  {
    return "Quad-P4-Gauss";
  }

  std::string name() const override
  {
    return "Quad-P4-Gauss";
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

// ----------------------------------------------------------------------------
// Quad quadrature for P5 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP5QuadGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP5QuadGauss();

  /// Destructor
  ~QdPointSetP5QuadGauss() override;

  static std::string type_name()
  {
    return "Quad-P5-Gauss";
  }

  std::string name() const override
  {
    return "Quad-P5-Gauss";
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

// ----------------------------------------------------------------------------
// Quad quadrature for P6 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP6QuadGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP6QuadGauss();

  /// Destructor
  ~QdPointSetP6QuadGauss() override;

  static std::string type_name()
  {
    return "Quad-P6-Gauss";
  }

  std::string name() const override
  {
    return "Quad-P6-Gauss";
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

// ----------------------------------------------------------------------------
// Quad quadrature for P7 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP7QuadGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP7QuadGauss();

  /// Destructor
  ~QdPointSetP7QuadGauss() override;

  static std::string type_name()
  {
    return "Quad-P7-Gauss";
  }

  std::string name() const override
  {
    return "Quad-P7-Gauss";
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

// ----------------------------------------------------------------------------
// Quad quadrature for P8 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP8QuadGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP8QuadGauss();

  /// Destructor
  ~QdPointSetP8QuadGauss() override;

  static std::string type_name()
  {
    return "Quad-P8-Gauss";
  }

  std::string name() const override
  {
    return "Quad-P8-Gauss";
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

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
