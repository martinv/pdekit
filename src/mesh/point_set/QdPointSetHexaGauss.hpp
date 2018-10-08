#ifndef PDEKIT_Mesh_Qd_Point_Set_Hexa_Gauss_hpp
#define PDEKIT_Mesh_Qd_Point_Set_Hexa_Gauss_hpp

#include "mesh/point_set/StdPointSetBase.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Hexahedron - quadrature for P1 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP1HexaGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP1HexaGauss();

  /// Destructor
  ~QdPointSetP1HexaGauss() override;

  static std::string type_name()
  {
    return "Hexa-P1-Gauss";
  }

  std::string name() const override
  {
    return "Hexa-P1-Gauss";
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
// Hexahedron - quadrature for P2 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP2HexaGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP2HexaGauss();

  /// Destructor
  ~QdPointSetP2HexaGauss() override;

  static std::string type_name()
  {
    return "Hexa-P2-Gauss";
  }

  std::string name() const override
  {
    return "Hexa-P2-Gauss";
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
// Hexahedron - quadrature for P3 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP3HexaGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP3HexaGauss();

  /// Destructor
  ~QdPointSetP3HexaGauss() override;

  static std::string type_name()
  {
    return "Hexa-P3-Gauss";
  }

  std::string name() const override
  {
    return "Hexa-P3-Gauss";
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
// Hexahedron - quadrature for P4 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP4HexaGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP4HexaGauss();

  /// Destructor
  ~QdPointSetP4HexaGauss() override;

  static std::string type_name()
  {
    return "Hexa-P4-Gauss";
  }

  std::string name() const override
  {
    return "Hexa-P4-Gauss";
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
// Hexahedron - quadrature for P5 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP5HexaGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP5HexaGauss();

  /// Destructor
  ~QdPointSetP5HexaGauss() override;

  static std::string type_name()
  {
    return "Hexa-P5-Gauss";
  }

  std::string name() const override
  {
    return "Hexa-P5-Gauss";
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
// Hexahedron - quadrature for P6 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP6HexaGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP6HexaGauss();

  /// Destructor
  ~QdPointSetP6HexaGauss() override;

  static std::string type_name()
  {
    return "Hexa-P6-Gauss";
  }

  std::string name() const override
  {
    return "Hexa-P6-Gauss";
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
// Hexahedron - quadrature for P7 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP7HexaGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP7HexaGauss();

  /// Destructor
  ~QdPointSetP7HexaGauss() override;

  static std::string type_name()
  {
    return "Hexa-P7-Gauss";
  }

  std::string name() const override
  {
    return "Hexa-P7-Gauss";
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
