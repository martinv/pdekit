#ifndef PDEKIT_Mesh_Qd_Point_Set_Quad_Face_Gauss_hpp
#define PDEKIT_Mesh_Qd_Point_Set_Quad_Face_Gauss_hpp

#include "common/Constants.hpp"
#include "common/StringUtils.hpp"
#include "mesh/point_set/QdPointSetLineGauss.hpp"
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

template <typename LineQuadratureType>
void fill_quad_face_quadrature_coords(math::DenseDMat<Real> &quadrature_coords, const Uint facet_id)
{
  LineQuadratureType quad1d;
  math::DenseDMat<Real> line_quadrature_coords;
  quad1d.reference_coords(line_quadrature_coords);

  const Uint nb_qd_pts = line_quadrature_coords.rows();
  quadrature_coords.resize(nb_qd_pts, _2D);

  switch (facet_id)
  {
    case (0):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        quadrature_coords(q, XI0) = line_quadrature_coords(q, XI0);
        quadrature_coords(q, XI1) = -1.;
      }
      break;
    }
    case (1):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        quadrature_coords(q, XI0) = 1.;
        quadrature_coords(q, XI1) = line_quadrature_coords(q, XI0);
      }
      break;
    }
    case (2):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        quadrature_coords(q, XI0) = line_quadrature_coords(nb_qd_pts - q - 1, XI0);
        quadrature_coords(q, XI1) = 1.;
      }
      break;
    }
    case (3):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        quadrature_coords(q, XI0) = -1.;
        quadrature_coords(q, XI1) = line_quadrature_coords(nb_qd_pts - q - 1, XI0);
      }
      break;
    }
  }; // switch
}

template <typename LineQuadratureType>
void fill_quad_face_quadrature_weights(math::DenseDVec<Real> &wgt, const Uint facet_id)
{
  LineQuadratureType quad1d;
  math::DenseDVec<Real> line_quadrature_weights;
  quad1d.weights(line_quadrature_weights);

  const Uint nb_qd_pts = line_quadrature_weights.size();
  wgt.resize(nb_qd_pts);

  switch (facet_id)
  {
    case (0):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        wgt[q] = line_quadrature_weights[q];
      }
      break;
    }
    case (1):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        wgt[q] = line_quadrature_weights[q];
      }
      break;
    }
    case (2):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        wgt[q] = line_quadrature_weights[nb_qd_pts - q - 1];
      }
      break;
    }
    case (3):
    {
      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        wgt[q] = line_quadrature_weights[nb_qd_pts - q - 1];
      }
      break;
    }
  }; // switch
}

void fill_quad_face_quadrature_permutation(const Uint local_id, const Uint nb_qd_pts,
                                           const mesh::EntityRealignCode &permutation_code,
                                           std::vector<Uint> &permutation_vec);

template <Uint QuadOrder>
struct SelectQuadFaceGaussQuadrature;

template <>
struct SelectQuadFaceGaussQuadrature<P1>
{
  using face_quad_type = QdPointSetP1LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P2>
{
  using face_quad_type = QdPointSetP2LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P3>
{
  using face_quad_type = QdPointSetP3LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P4>
{
  using face_quad_type = QdPointSetP4LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P5>
{
  using face_quad_type = QdPointSetP5LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P6>
{
  using face_quad_type = QdPointSetP6LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P7>
{
  using face_quad_type = QdPointSetP7LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P8>
{
  using face_quad_type = QdPointSetP8LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P9>
{
  using face_quad_type = QdPointSetP9LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P10>
{
  using face_quad_type = QdPointSetP10LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P11>
{
  using face_quad_type = QdPointSetP11LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P12>
{
  using face_quad_type = QdPointSetP12LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P13>
{
  using face_quad_type = QdPointSetP13LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P14>
{
  using face_quad_type = QdPointSetP14LineGauss;
};

template <>
struct SelectQuadFaceGaussQuadrature<P15>
{
  using face_quad_type = QdPointSetP15LineGauss;
};

} // namespace detail

// ----------------------------------------------------------------------------
// Quadrilateral face quadrature for polynomials of order 'QuadOrder'
// ----------------------------------------------------------------------------

template <Uint QuadOrder>
class QdPointSetQuadFaceGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetQuadFaceGauss();

  /// Destructor
  ~QdPointSetQuadFaceGauss() override;

  static std::string type_name()
  {
    return "Quad-P" + common::StringUtils::to_string(QuadOrder) + "-FaceGauss";
  }

  std::string name() const override
  {
    return "Quad-P" + common::StringUtils::to_string(QuadOrder) + "-FaceGauss";
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
  typedef typename detail::SelectQuadFaceGaussQuadrature<QuadOrder>::face_quad_type face_quad_type;
};

template <Uint QuadOrder>
QdPointSetQuadFaceGauss<QuadOrder>::QdPointSetQuadFaceGauss() : StdPointSetBase()
{
}

template <Uint QuadOrder>
QdPointSetQuadFaceGauss<QuadOrder>::~QdPointSetQuadFaceGauss()
{
}

template <Uint QuadOrder>
Uint QdPointSetQuadFaceGauss<QuadOrder>::order() const
{
  return QuadOrder;
}

template <Uint QuadOrder>
Uint QdPointSetQuadFaceGauss<QuadOrder>::dim() const
{
  return _2D;
}

template <Uint QuadOrder>
Uint QdPointSetQuadFaceGauss<QuadOrder>::codim() const
{
  return _1D;
}

template <Uint QuadOrder>
Uint QdPointSetQuadFaceGauss<QuadOrder>::nb_local_entities() const
{
  return 4u;
}

template <Uint QuadOrder>
Uint QdPointSetQuadFaceGauss<QuadOrder>::size(const Uint local_idx) const
{
  const face_quad_type face_quad;
  return face_quad.size();
}

template <Uint QuadOrder>
void QdPointSetQuadFaceGauss<QuadOrder>::reference_coords(math::DenseDMat<Real> &coords,
                                                          const Uint local_idx) const
{
  detail::fill_quad_face_quadrature_coords<face_quad_type>(coords, local_idx);
}

template <Uint QuadOrder>
void QdPointSetQuadFaceGauss<QuadOrder>::weights(math::DenseDVec<Real> &wgt,
                                                 const Uint local_idx) const
{
  detail::fill_quad_face_quadrature_weights<face_quad_type>(wgt, local_idx);
}

template <Uint QuadOrder>
void QdPointSetQuadFaceGauss<QuadOrder>::permutation(
    const Uint local_id, const mesh::EntityRealignCode &permutation_code,
    std::vector<Uint> &permutation_vec)
{
  const face_quad_type face_quad;
  detail::fill_quad_face_quadrature_permutation(local_id, face_quad.size(), permutation_code,
                                                permutation_vec);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
