#ifndef PDEKIT_Mesh_Qd_Point_Set_Triag_Gauss_hpp
#define PDEKIT_Mesh_Qd_Point_Set_Triag_Gauss_hpp

#include "mesh/point_set/StdPointSetBase.hpp"

namespace pdekit
{

namespace mesh
{

namespace detail
{
void fill_triag_quadrature_permutation(const Uint nb_qd_pts,
                                       const mesh::EntityRealignCode &permutation_code,
                                       const std::vector<Uint> &canonical_rot_permutation,
                                       const std::vector<Uint> &canonical_flip_permutation,
                                       std::vector<Uint> &permutation_vec);
}

// ----------------------------------------------------------------------------
// Triangle quadrature for P1 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP1TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP1TriagGauss();

  /// Destructor
  ~QdPointSetP1TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P1-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P1-Gauss";
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
    N_QD_PTS = 1
  };
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P2 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP2TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP2TriagGauss();

  /// Destructor
  ~QdPointSetP2TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P2-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P2-Gauss";
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

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P3 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP3TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP3TriagGauss();

  /// Destructor
  ~QdPointSetP3TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P3-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P3-Gauss";
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

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P4 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP4TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP4TriagGauss();

  /// Destructor
  ~QdPointSetP4TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P4-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P4-Gauss";
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

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P5 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP5TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP5TriagGauss();

  /// Destructor
  ~QdPointSetP5TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P5-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P5-Gauss";
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

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P6 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP6TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP6TriagGauss();

  /// Destructor
  ~QdPointSetP6TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P6-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P6-Gauss";
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
    N_QD_PTS = 12
  };

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P7 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP7TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP7TriagGauss();

  /// Destructor
  ~QdPointSetP7TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P7-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P7-Gauss";
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
    N_QD_PTS = 13
  };

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P8 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP8TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP8TriagGauss();

  /// Destructor
  ~QdPointSetP8TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P8-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P8-Gauss";
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
    N_QD_PTS = 16
  };

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P9 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP9TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP9TriagGauss();

  /// Destructor
  ~QdPointSetP9TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P9-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P9-Gauss";
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
    N_QD_PTS = 19
  };

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P10 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP10TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP10TriagGauss();

  /// Destructor
  ~QdPointSetP10TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P10-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P10-Gauss";
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
    N_QD_PTS = 25
  };

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P11 polynomials
// 0 negative weights, 3 points outside of the triangle
// ----------------------------------------------------------------------------

class QdPointSetP11TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP11TriagGauss();

  /// Destructor
  ~QdPointSetP11TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P11-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P11-Gauss";
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
    N_QD_PTS = 27
  };

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P12 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP12TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP12TriagGauss();

  /// Destructor
  ~QdPointSetP12TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P12-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P12-Gauss";
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
    N_QD_PTS = 33
  };

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P13 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP13TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP13TriagGauss();

  /// Destructor
  ~QdPointSetP13TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P13-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P13-Gauss";
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
    N_QD_PTS = 37
  };

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P14 polynomials
// ----------------------------------------------------------------------------

class QdPointSetP14TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP14TriagGauss();

  /// Destructor
  ~QdPointSetP14TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P14-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P14-Gauss";
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
    N_QD_PTS = 42
  };

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------
// Triangle quadrature for P15 polynomials
// 0 negative weights, 9 points outside of the triangle
// ----------------------------------------------------------------------------

class QdPointSetP15TriagGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetP15TriagGauss();

  /// Destructor
  ~QdPointSetP15TriagGauss() override;

  static std::string type_name()
  {
    return "Triag-P15-Gauss";
  }

  std::string name() const override
  {
    return "Triag-P15-Gauss";
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
    N_QD_PTS = 48
  };

  /// Compute quadrature point coordinates in reference space
  static void compute_reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0);

  /// Fill the vectors which determine how to reorder quadrature
  /// points to have one rotation/flip permutation
  static void initialize_canonical_permutations();

  /// Flag marking whether canonical permutations (i.e. vectors determining
  /// reordering of quadrature points) for rotation and flip have already been
  /// filled
  static bool canonical_permutations_initialized;

  /// Vector determining reordering of points for one rotation
  static std::vector<Uint> canonical_rot_permutation;

  /// Vector determining reordering of points for one flip
  static std::vector<Uint> canonical_flip_permutation;
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
