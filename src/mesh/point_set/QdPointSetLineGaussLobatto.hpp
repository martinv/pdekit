#ifndef PDEKIT_Mesh_Qd_Point_Set_Line_Gauss_Lobatto_hpp
#define PDEKIT_Mesh_Qd_Point_Set_Line_Gauss_Lobatto_hpp

#include "mesh/point_set/StdPointSetBase.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Line quadrature for P1 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP1LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 2
  };

  /// Default constructor
  QdPointSetP1LineGaussLobatto();

  /// Destructor
  ~QdPointSetP1LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P1-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P1-Gauss-Lobatto";
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
// Line quadrature for P2 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP2LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 3
  };

  /// Default constructor
  QdPointSetP2LineGaussLobatto();

  /// Destructor
  ~QdPointSetP2LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P2-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P2-Gauss-Lobatto";
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
// Line quadrature for P3 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP3LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 3
  };

  /// Default constructor
  QdPointSetP3LineGaussLobatto();

  /// Destructor
  ~QdPointSetP3LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P3-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P3-Gauss-Lobatto";
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
// Line quadrature for P4 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP4LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 4
  };

  /// Default constructor
  QdPointSetP4LineGaussLobatto();

  /// Destructor
  ~QdPointSetP4LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P4-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P4-Gauss-Lobatto";
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
// Line quadrature for P5 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP5LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 4
  };

  /// Default constructor
  QdPointSetP5LineGaussLobatto();

  /// Destructor
  ~QdPointSetP5LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P5-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P5-Gauss-Lobatto";
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
// Line quadrature for P6 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP6LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 5
  };

  /// Default constructor
  QdPointSetP6LineGaussLobatto();

  /// Destructor
  ~QdPointSetP6LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P6-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P6-Gauss-Lobatto";
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
// Line quadrature for P7 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP7LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 5
  };

  /// Default constructor
  QdPointSetP7LineGaussLobatto();

  /// Destructor
  ~QdPointSetP7LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P7-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P7-Gauss-Lobatto";
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
// Line quadrature for P8 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP8LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 6
  };

  /// Default constructor
  QdPointSetP8LineGaussLobatto();

  /// Destructor
  ~QdPointSetP8LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P8-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P8-Gauss-Lobatto";
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
// Line quadrature for P9 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP9LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 6
  };

  /// Default constructor
  QdPointSetP9LineGaussLobatto();

  /// Destructor
  ~QdPointSetP9LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P9-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P9-Gauss-Lobatto";
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
// Line quadrature for P10 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP10LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 7
  };

  /// Default constructor
  QdPointSetP10LineGaussLobatto();

  /// Destructor
  ~QdPointSetP10LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P10-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P10-Gauss-Lobatto";
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
// Line quadrature for P11 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP11LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 7
  };

  /// Default constructor
  QdPointSetP11LineGaussLobatto();

  /// Destructor
  ~QdPointSetP11LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P11-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P11-Gauss-Lobatto";
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
// Line quadrature for P12 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP12LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 8
  };

  /// Default constructor
  QdPointSetP12LineGaussLobatto();

  /// Destructor
  ~QdPointSetP12LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P12-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P12-Gauss-Lobatto";
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
// Line quadrature for P13 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP13LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 8
  };

  /// Default constructor
  QdPointSetP13LineGaussLobatto();

  /// Destructor
  ~QdPointSetP13LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P13-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P13-Gauss-Lobatto";
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
// Line quadrature for P14 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP14LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 9
  };

  /// Default constructor
  QdPointSetP14LineGaussLobatto();

  /// Destructor
  ~QdPointSetP14LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P14-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P14-Gauss-Lobatto";
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
// Line quadrature for P15 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP15LineGaussLobatto : public StdPointSetBase
{
  public:
  enum
  {
    N_QD_PTS = 9
  };

  /// Default constructor
  QdPointSetP15LineGaussLobatto();

  /// Destructor
  ~QdPointSetP15LineGaussLobatto() override;

  static std::string type_name()
  {
    return "Line-P15-Gauss-Lobatto";
  }

  std::string name() const override
  {
    return "Line-P15-Gauss-Lobatto";
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
