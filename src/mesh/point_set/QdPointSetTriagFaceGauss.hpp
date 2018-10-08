#ifndef PDEKIT_Mesh_Qd_Point_Set_Triag_Face_Gauss_hpp
#define PDEKIT_Mesh_Qd_Point_Set_Triag_Face_Gauss_hpp

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
void fill_triag_face_quadrature_coords(math::DenseDMat<Real> &quadrature_coords,
                                       const Uint facet_id)
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
        quadrature_coords(q, XI0) = line_quadrature_coords(nb_qd_pts - q - 1, XI0);
        quadrature_coords(q, XI1) = -quadrature_coords(q, XI0);
      }
      break;
    }
    case (2):
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
void fill_triag_face_quadrature_weights(math::DenseDVec<Real> &wgt, const Uint facet_id)
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
      const Real sqrt2 = std::sqrt(2.);

      for (Uint q = 0; q < nb_qd_pts; ++q)
      {
        wgt[q] = sqrt2 * line_quadrature_weights[nb_qd_pts - q - 1];
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
  }; // switch
}

void fill_triag_face_quadrature_permutation(const Uint local_id, const Uint nb_qd_pts,
                                            const mesh::EntityRealignCode &permutation_code,
                                            std::vector<Uint> &permutation_vec);

template <Uint QuadOrder>
struct SelectTriagFaceGaussQuadrature;

template <>
struct SelectTriagFaceGaussQuadrature<P1>
{
  typedef QdPointSetP1LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P2>
{
  typedef QdPointSetP2LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P3>
{
  typedef QdPointSetP3LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P4>
{
  typedef QdPointSetP4LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P5>
{
  typedef QdPointSetP5LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P6>
{
  typedef QdPointSetP6LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P7>
{
  typedef QdPointSetP7LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P8>
{
  typedef QdPointSetP8LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P9>
{
  typedef QdPointSetP9LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P10>
{
  typedef QdPointSetP10LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P11>
{
  typedef QdPointSetP11LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P12>
{
  typedef QdPointSetP12LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P13>
{
  typedef QdPointSetP13LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P14>
{
  typedef QdPointSetP14LineGauss face_quad_type;
};

template <>
struct SelectTriagFaceGaussQuadrature<P15>
{
  typedef QdPointSetP15LineGauss face_quad_type;
};

} // namespace detail

// ----------------------------------------------------------------------------
// Triangle face quadrature for polynomials of order 'QuadOrder'
// ----------------------------------------------------------------------------

template <Uint QuadOrder>
class QdPointSetTriagFaceGauss : public StdPointSetBase
{
  public:
  /// Default constructor
  QdPointSetTriagFaceGauss();

  /// Destructor
  ~QdPointSetTriagFaceGauss() override;

  static std::string type_name()
  {
    return "Triag-P" + common::StringUtils::to_string(QuadOrder) + "-FaceGauss";
  }

  std::string name() const override
  {
    return "Triag-P" + common::StringUtils::to_string(QuadOrder) + "-FaceGauss";
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
  typedef typename detail::SelectTriagFaceGaussQuadrature<QuadOrder>::face_quad_type face_quad_type;
};

template <Uint QuadOrder>
QdPointSetTriagFaceGauss<QuadOrder>::QdPointSetTriagFaceGauss() : StdPointSetBase()
{
}

template <Uint QuadOrder>
QdPointSetTriagFaceGauss<QuadOrder>::~QdPointSetTriagFaceGauss()
{
}

template <Uint QuadOrder>
Uint QdPointSetTriagFaceGauss<QuadOrder>::order() const
{
  return QuadOrder;
}

template <Uint QuadOrder>
Uint QdPointSetTriagFaceGauss<QuadOrder>::dim() const
{
  return _2D;
}

template <Uint QuadOrder>
Uint QdPointSetTriagFaceGauss<QuadOrder>::codim() const
{
  return _1D;
}

template <Uint QuadOrder>
Uint QdPointSetTriagFaceGauss<QuadOrder>::nb_local_entities() const
{
  return 3u;
}

template <Uint QuadOrder>
Uint QdPointSetTriagFaceGauss<QuadOrder>::size(const Uint local_idx) const
{
  const face_quad_type face_quad;
  return face_quad.size();
}

template <Uint QuadOrder>
void QdPointSetTriagFaceGauss<QuadOrder>::reference_coords(math::DenseDMat<Real> &coords,
                                                           const Uint local_idx) const
{
  detail::fill_triag_face_quadrature_coords<face_quad_type>(coords, local_idx);
}

template <Uint QuadOrder>
void QdPointSetTriagFaceGauss<QuadOrder>::weights(math::DenseDVec<Real> &wgt,
                                                  const Uint local_idx) const
{
  detail::fill_triag_face_quadrature_weights<face_quad_type>(wgt, local_idx);
}

template <Uint QuadOrder>
void QdPointSetTriagFaceGauss<QuadOrder>::permutation(
    const Uint local_id, const mesh::EntityRealignCode &permutation_code,
    std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID
  const face_quad_type face_quad;
  detail::fill_triag_face_quadrature_permutation(local_id, face_quad.size(), permutation_code,
                                                 permutation_vec);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
