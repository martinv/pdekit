#ifndef PDEKIT_Mesh_Qd_Point_Set_Tetra_Face_Gauss_hpp
#define PDEKIT_Mesh_Qd_Point_Set_Tetra_Face_Gauss_hpp

#include "common/Constants.hpp"
#include "common/StringUtils.hpp"
#include "math/DenseSVec.hpp"
#include "mesh/point_set/QdPointSetTriagGauss.hpp"
#include "mesh/point_set/QuadratureTransformUtils.hpp"
#include "mesh/point_set/StdPointSetBase.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// HELPER FUNCTIONS
// ----------------------------------------------------------------------------

namespace detail
{

template <typename TriagQuadratureType>
void fill_tetra_face_quadrature_coords(math::DenseDMat<Real> &quadrature_coords,
                                       const Uint facet_id)
{
  TriagQuadratureType quad_triag_face;
  math::DenseDMat<Real> triag_face_quadrature_coords;
  quad_triag_face.reference_coords(triag_face_quadrature_coords);

  const math::DenseSVec<Real, 2> V0_tri = math::values_list(-1.)(-1.);
  const math::DenseSVec<Real, 2> V1_tri = math::values_list(1.)(-1.);
  const math::DenseSVec<Real, 2> V2_tri = math::values_list(-1.)(1.);

  const math::DenseSVec<Real, 3> V0_tet = math::values_list(-1.)(-1.)(-1.);
  const math::DenseSVec<Real, 3> V1_tet = math::values_list(1.)(-1.)(-1.);
  const math::DenseSVec<Real, 3> V2_tet = math::values_list(-1.)(1.)(-1.);
  const math::DenseSVec<Real, 3> V3_tet = math::values_list(-1.)(-1.)(1.);

  math::DenseDMat<Real> triag_face_bar_coords;
  detail::compute_triag_barycentric_coords(V0_tri, V1_tri, V2_tri, triag_face_quadrature_coords,
                                           triag_face_bar_coords);

  const Uint nb_qd_pts = triag_face_quadrature_coords.rows();
  quadrature_coords.resize(nb_qd_pts, _3D);

  switch (facet_id)
  {
    case (0):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        quadrature_coords(q, XI0) = V0_tet[XI0] * triag_face_bar_coords(q, 0) +
                                    V2_tet[XI0] * triag_face_bar_coords(q, 1) +
                                    V1_tet[XI0] * triag_face_bar_coords(q, 2);
        quadrature_coords(q, XI1) = V0_tet[XI1] * triag_face_bar_coords(q, 0) +
                                    V2_tet[XI1] * triag_face_bar_coords(q, 1) +
                                    V1_tet[XI1] * triag_face_bar_coords(q, 2);
        quadrature_coords(q, XI2) = V0_tet[XI2] * triag_face_bar_coords(q, 0) +
                                    V2_tet[XI2] * triag_face_bar_coords(q, 1) +
                                    V1_tet[XI2] * triag_face_bar_coords(q, 2);
      }
      break;
    }
    case (1):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        quadrature_coords(q, XI0) = V0_tet[XI0] * triag_face_bar_coords(q, 0) +
                                    V1_tet[XI0] * triag_face_bar_coords(q, 1) +
                                    V3_tet[XI0] * triag_face_bar_coords(q, 2);
        quadrature_coords(q, XI1) = V0_tet[XI1] * triag_face_bar_coords(q, 0) +
                                    V1_tet[XI1] * triag_face_bar_coords(q, 1) +
                                    V3_tet[XI1] * triag_face_bar_coords(q, 2);
        quadrature_coords(q, XI2) = V0_tet[XI2] * triag_face_bar_coords(q, 0) +
                                    V1_tet[XI2] * triag_face_bar_coords(q, 1) +
                                    V3_tet[XI2] * triag_face_bar_coords(q, 2);
      }
      break;
    }
    case (2):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        quadrature_coords(q, XI0) = V0_tet[XI0] * triag_face_bar_coords(q, 0) +
                                    V3_tet[XI0] * triag_face_bar_coords(q, 1) +
                                    V2_tet[XI0] * triag_face_bar_coords(q, 2);
        quadrature_coords(q, XI1) = V0_tet[XI1] * triag_face_bar_coords(q, 0) +
                                    V3_tet[XI1] * triag_face_bar_coords(q, 1) +
                                    V2_tet[XI1] * triag_face_bar_coords(q, 2);
        quadrature_coords(q, XI2) = V0_tet[XI2] * triag_face_bar_coords(q, 0) +
                                    V3_tet[XI2] * triag_face_bar_coords(q, 1) +
                                    V2_tet[XI2] * triag_face_bar_coords(q, 2);
      }
      break;
    }
    case (3):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        quadrature_coords(q, XI0) = V3_tet[XI0] * triag_face_bar_coords(q, 0) +
                                    V1_tet[XI0] * triag_face_bar_coords(q, 1) +
                                    V2_tet[XI0] * triag_face_bar_coords(q, 2);
        quadrature_coords(q, XI1) = V3_tet[XI1] * triag_face_bar_coords(q, 0) +
                                    V1_tet[XI1] * triag_face_bar_coords(q, 1) +
                                    V2_tet[XI1] * triag_face_bar_coords(q, 2);
        quadrature_coords(q, XI2) = V3_tet[XI2] * triag_face_bar_coords(q, 0) +
                                    V1_tet[XI2] * triag_face_bar_coords(q, 1) +
                                    V2_tet[XI2] * triag_face_bar_coords(q, 2);
      }
      break;
    }

  }; // switch
}

template <Uint QuadOrder>
struct SelectTetraFaceGaussQuadrature;

template <>
struct SelectTetraFaceGaussQuadrature<P1>
{
  typedef QdPointSetP1TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P2>
{
  typedef QdPointSetP2TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P3>
{
  typedef QdPointSetP3TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P4>
{
  typedef QdPointSetP4TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P5>
{
  typedef QdPointSetP5TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P6>
{
  typedef QdPointSetP6TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P7>
{
  typedef QdPointSetP7TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P8>
{
  typedef QdPointSetP8TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P9>
{
  typedef QdPointSetP9TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P10>
{
  typedef QdPointSetP10TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P11>
{
  typedef QdPointSetP11TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P12>
{
  typedef QdPointSetP12TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P13>
{
  typedef QdPointSetP13TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P14>
{
  typedef QdPointSetP14TriagGauss face_quad_type;
};

template <>
struct SelectTetraFaceGaussQuadrature<P15>
{
  typedef QdPointSetP15TriagGauss face_quad_type;
};

} // namespace detail

// ----------------------------------------------------------------------------
// Tetrahedra facet quadrature for P1 polynomials
// ----------------------------------------------------------------------------

template <Uint QuadOrder>
class QdPointSetTetraFaceGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetTetraFaceGauss();

  /// Destructor
  ~QdPointSetTetraFaceGauss() override;

  static std::string type_name()
  {
    return "Tetra-P" + common::StringUtils::to_string(QuadOrder) + "-FaceGauss";
  }

  std::string name() const override
  {
    return "Tetra-P" + common::StringUtils::to_string(QuadOrder) + "-FaceGauss";
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
  typedef typename detail::SelectTetraFaceGaussQuadrature<QuadOrder>::face_quad_type face_quad_type;
};

template <Uint QuadOrder>
QdPointSetTetraFaceGauss<QuadOrder>::QdPointSetTetraFaceGauss() : StdPointSetBase()
{
}

template <Uint QuadOrder>
QdPointSetTetraFaceGauss<QuadOrder>::~QdPointSetTetraFaceGauss()
{
}

template <Uint QuadOrder>
Uint QdPointSetTetraFaceGauss<QuadOrder>::order() const
{
  return QuadOrder;
}

template <Uint QuadOrder>
Uint QdPointSetTetraFaceGauss<QuadOrder>::dim() const
{
  return _3D;
}

template <Uint QuadOrder>
Uint QdPointSetTetraFaceGauss<QuadOrder>::codim() const
{
  return _1D;
}

template <Uint QuadOrder>
Uint QdPointSetTetraFaceGauss<QuadOrder>::nb_local_entities() const
{
  return 4u;
}

template <Uint QuadOrder>
Uint QdPointSetTetraFaceGauss<QuadOrder>::size(const Uint local_idx) const
{
  const face_quad_type face_quad;
  return face_quad.size();
}

template <Uint QuadOrder>
void QdPointSetTetraFaceGauss<QuadOrder>::reference_coords(math::DenseDMat<Real> &coords,
                                                           const Uint local_idx) const
{
  detail::fill_tetra_face_quadrature_coords<face_quad_type>(coords, local_idx);
}

template <Uint QuadOrder>
void QdPointSetTetraFaceGauss<QuadOrder>::weights(math::DenseDVec<Real> &wgt,
                                                  const Uint local_idx) const
{
  face_quad_type tri_face_quad;
  tri_face_quad.weights(wgt, local_idx);

  if (local_idx == 3)
  {
    for (Uint q = 0; q < wgt.size(); ++q)
    {
      wgt[q] *= std::sqrt(3.);
    }
  }
}

template <Uint QuadOrder>
void QdPointSetTetraFaceGauss<QuadOrder>::permutation(
    const Uint local_id, const mesh::EntityRealignCode &permutation_code,
    std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID
  face_quad_type tri_face_quad;
  tri_face_quad.permutation(local_id, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
