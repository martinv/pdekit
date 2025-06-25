#ifndef PDEKIT_Mesh_Qd_Point_Set_Line_Gauss_hpp
#define PDEKIT_Mesh_Qd_Point_Set_Line_Gauss_hpp

#include "mesh/point_set/StdPointSetBase.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Line quadrature for P1 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP1LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP1LineGauss();

  /// Destructor
  ~QdPointSetP1LineGauss() override;

  static std::string type_name()
  {
    return "Line-P1-Gauss";
  }

  std::string name() const override
  {
    return "Line-P1-Gauss";
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
  enum
  {
    N_QD_PTS = 2
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P2 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP2LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP2LineGauss();

  /// Destructor
  ~QdPointSetP2LineGauss() override;

  static std::string type_name()
  {
    return "Line-P2-Gauss";
  }

  std::string name() const override
  {
    return "Line-P2-Gauss";
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
  enum
  {
    N_QD_PTS = 2
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P3 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP3LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP3LineGauss();

  /// Destructor
  ~QdPointSetP3LineGauss() override;

  static std::string type_name()
  {
    return "Line-P3-Gauss";
  }

  std::string name() const override
  {
    return "Line-P3-Gauss";
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
  enum
  {
    N_QD_PTS = 2
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P4 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP4LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP4LineGauss();

  /// Destructor
  ~QdPointSetP4LineGauss() override;

  static std::string type_name()
  {
    return "Line-P4-Gauss";
  }

  std::string name() const override
  {
    return "Line-P4-Gauss";
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
  enum
  {
    N_QD_PTS = 3
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P5 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP5LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP5LineGauss();

  /// Destructor
  ~QdPointSetP5LineGauss() override;

  static std::string type_name()
  {
    return "Line-P5-Gauss";
  }

  std::string name() const override
  {
    return "Line-P5-Gauss";
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
  enum
  {
    N_QD_PTS = 3
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P6 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP6LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP6LineGauss();

  /// Destructor
  ~QdPointSetP6LineGauss() override;

  static std::string type_name()
  {
    return "Line-P6-Gauss";
  }

  std::string name() const override
  {
    return "Line-P6-Gauss";
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
  enum
  {
    N_QD_PTS = 4
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P7 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP7LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP7LineGauss();

  /// Destructor
  ~QdPointSetP7LineGauss() override;

  static std::string type_name()
  {
    return "Line-P7-Gauss";
  }

  std::string name() const override
  {
    return "Line-P7-Gauss";
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
  enum
  {
    N_QD_PTS = 4
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P8 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP8LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP8LineGauss();

  /// Destructor
  ~QdPointSetP8LineGauss() override;

  static std::string type_name()
  {
    return "Line-P8-Gauss";
  }

  std::string name() const override
  {
    return "Line-P8-Gauss";
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
  enum
  {
    N_QD_PTS = 5
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P9 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP9LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP9LineGauss();

  /// Destructor
  ~QdPointSetP9LineGauss() override;

  static std::string type_name()
  {
    return "Line-P9-Gauss";
  }

  std::string name() const override
  {
    return "Line-P9-Gauss";
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
  enum
  {
    N_QD_PTS = 5
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P10 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP10LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP10LineGauss();

  /// Destructor
  ~QdPointSetP10LineGauss() override;

  static std::string type_name()
  {
    return "Line-P10-Gauss";
  }

  std::string name() const override
  {
    return "Line-P10-Gauss";
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
  enum
  {
    N_QD_PTS = 6
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P11 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP11LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP11LineGauss();

  /// Destructor
  ~QdPointSetP11LineGauss() override;

  static std::string type_name()
  {
    return "Line-P11-Gauss";
  }

  std::string name() const override
  {
    return "Line-P11-Gauss";
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
  enum
  {
    N_QD_PTS = 6
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P12 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP12LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP12LineGauss();

  /// Destructor
  ~QdPointSetP12LineGauss() override;

  static std::string type_name()
  {
    return "Line-P12-Gauss";
  }

  std::string name() const override
  {
    return "Line-P12-Gauss";
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
  enum
  {
    N_QD_PTS = 7
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P13 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP13LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP13LineGauss();

  /// Destructor
  ~QdPointSetP13LineGauss() override;

  static std::string type_name()
  {
    return "Line-P13-Gauss";
  }

  std::string name() const override
  {
    return "Line-P13-Gauss";
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
  enum
  {
    N_QD_PTS = 7
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P14 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP14LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP14LineGauss();

  /// Destructor
  ~QdPointSetP14LineGauss() override;

  static std::string type_name()
  {
    return "Line-P14-Gauss";
  }

  std::string name() const override
  {
    return "Line-P14-Gauss";
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
  enum
  {
    N_QD_PTS = 8
  };
};

// ----------------------------------------------------------------------------
// Line quadrature for P15 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP15LineGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP15LineGauss();

  /// Destructor
  ~QdPointSetP15LineGauss() override;

  static std::string type_name()
  {
    return "Line-P15-Gauss";
  }

  std::string name() const override
  {
    return "Line-P15-Gauss";
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
  enum
  {
    N_QD_PTS = 8
  };
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
