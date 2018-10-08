#include "mesh/point_set/QdPointSetLineGaussLobatto.hpp"
#include "common/Constants.hpp"
#include "mesh/point_set/QuadratureTransformUtils.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Line quadrature for P1 polynomials
// ----------------------------------------------------------------------------

QdPointSetP1LineGaussLobatto::QdPointSetP1LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP1LineGaussLobatto::~QdPointSetP1LineGaussLobatto()
{
}

Uint QdPointSetP1LineGaussLobatto::order() const
{
  return P1;
}

Uint QdPointSetP1LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP1LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP1LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP1LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP1LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                    const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -1.0;
  coords(1, XI0) = 1.0;
}

void QdPointSetP1LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 1.0;
  wgt[1] = 1.0;
}

void QdPointSetP1LineGaussLobatto::permutation(const Uint local_id,
                                               const mesh::EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P2 polynomials
// ----------------------------------------------------------------------------

QdPointSetP2LineGaussLobatto::QdPointSetP2LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP2LineGaussLobatto::~QdPointSetP2LineGaussLobatto()
{
}

Uint QdPointSetP2LineGaussLobatto::order() const
{
  return P2;
}

Uint QdPointSetP2LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP2LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP2LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP2LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP2LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                    const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -1.0;
  coords(1, XI0) = 0.0;
  coords(2, XI0) = 1.0;
}

void QdPointSetP2LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 1. / 3.;
  wgt[1] = 4. / 3.;
  wgt[2] = 1. / 3.;
}

void QdPointSetP2LineGaussLobatto::permutation(const Uint local_id,
                                               const mesh::EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P3 polynomials
// ----------------------------------------------------------------------------

QdPointSetP3LineGaussLobatto::QdPointSetP3LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP3LineGaussLobatto::~QdPointSetP3LineGaussLobatto()
{
}

Uint QdPointSetP3LineGaussLobatto::order() const
{
  return P3;
}

Uint QdPointSetP3LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP3LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP3LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP3LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP3LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                    const Uint local_idx) const
{
  QdPointSetP2LineGaussLobatto quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP3LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP2LineGaussLobatto quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP3LineGaussLobatto::permutation(const Uint local_id,
                                               const mesh::EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P4 polynomials
// ----------------------------------------------------------------------------

QdPointSetP4LineGaussLobatto::QdPointSetP4LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP4LineGaussLobatto::~QdPointSetP4LineGaussLobatto()
{
}

Uint QdPointSetP4LineGaussLobatto::order() const
{
  return P4;
}

Uint QdPointSetP4LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP4LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP4LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP4LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP4LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                    const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -1.0;
  coords(1, XI0) = -0.4472135954999579392818;
  coords(2, XI0) = 0.4472135954999579392818;
  coords(3, XI0) = 1.0;
}

void QdPointSetP4LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.1666666666666666666667;
  wgt[1] = 0.833333333333333333333;
  wgt[2] = 0.833333333333333333333;
  wgt[3] = 0.1666666666666666666667;
}

void QdPointSetP4LineGaussLobatto::permutation(const Uint local_id,
                                               const mesh::EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P5 polynomials
// ----------------------------------------------------------------------------

QdPointSetP5LineGaussLobatto::QdPointSetP5LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP5LineGaussLobatto::~QdPointSetP5LineGaussLobatto()
{
}

Uint QdPointSetP5LineGaussLobatto::order() const
{
  return P5;
}

Uint QdPointSetP5LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP5LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP5LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP5LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP5LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                    const Uint local_idx) const
{
  QdPointSetP4LineGaussLobatto quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP5LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP4LineGaussLobatto quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP5LineGaussLobatto::permutation(const Uint local_id,
                                               const mesh::EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P6 polynomials
// ----------------------------------------------------------------------------

QdPointSetP6LineGaussLobatto::QdPointSetP6LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP6LineGaussLobatto::~QdPointSetP6LineGaussLobatto()
{
}

Uint QdPointSetP6LineGaussLobatto::order() const
{
  return P6;
}

Uint QdPointSetP6LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP6LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP6LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP6LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP6LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                    const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = 1.0;
  coords(1, XI0) = -0.6546536707079771437983;
  coords(2, XI0) = 0.0;
  coords(3, XI0) = 0.6546536707079771437983;
  coords(4, XI0) = 1.0;
}

void QdPointSetP6LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.1;
  wgt[1] = 0.5444444444444444444444;
  wgt[2] = 0.7111111111111111111111;
  wgt[3] = 0.5444444444444444444444;
  wgt[4] = 0.1;
}

void QdPointSetP6LineGaussLobatto::permutation(const Uint local_id,
                                               const mesh::EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P7 polynomials
// ----------------------------------------------------------------------------

QdPointSetP7LineGaussLobatto::QdPointSetP7LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP7LineGaussLobatto::~QdPointSetP7LineGaussLobatto()
{
}

Uint QdPointSetP7LineGaussLobatto::order() const
{
  return P7;
}

Uint QdPointSetP7LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP7LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP7LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP7LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP7LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                    const Uint local_idx) const
{
  QdPointSetP6LineGaussLobatto quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP7LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP6LineGaussLobatto quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP7LineGaussLobatto::permutation(const Uint local_id,
                                               const mesh::EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P8 polynomials
// ----------------------------------------------------------------------------

QdPointSetP8LineGaussLobatto::QdPointSetP8LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP8LineGaussLobatto::~QdPointSetP8LineGaussLobatto()
{
}

Uint QdPointSetP8LineGaussLobatto::order() const
{
  return P8;
}

Uint QdPointSetP8LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP8LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP8LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP8LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP8LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                    const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -1.0;
  coords(1, XI0) = -0.765055323929464692851;
  coords(2, XI0) = -0.2852315164806450963142;
  coords(3, XI0) = 0.2852315164806450963142;
  coords(4, XI0) = 0.765055323929464692851;
  coords(5, XI0) = 1.0;
}

void QdPointSetP8LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.06666666666666666666667;
  wgt[1] = 0.378474956297846980317;
  wgt[2] = 0.55485837703548635302;
  wgt[3] = 0.55485837703548635302;
  wgt[4] = 0.378474956297846980317;
  wgt[5] = 0.066666666666666666666;
}

void QdPointSetP8LineGaussLobatto::permutation(const Uint local_id,
                                               const mesh::EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P9 polynomials
// ----------------------------------------------------------------------------

QdPointSetP9LineGaussLobatto::QdPointSetP9LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP9LineGaussLobatto::~QdPointSetP9LineGaussLobatto()
{
}

Uint QdPointSetP9LineGaussLobatto::order() const
{
  return P9;
}

Uint QdPointSetP9LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP9LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP9LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP9LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP9LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                    const Uint local_idx) const
{
  QdPointSetP8LineGaussLobatto quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP9LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP8LineGaussLobatto quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP9LineGaussLobatto::permutation(const Uint local_id,
                                               const mesh::EntityRealignCode &permutation_code,
                                               std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P10 polynomials
// ----------------------------------------------------------------------------

QdPointSetP10LineGaussLobatto::QdPointSetP10LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP10LineGaussLobatto::~QdPointSetP10LineGaussLobatto()
{
}

Uint QdPointSetP10LineGaussLobatto::order() const
{
  return P10;
}

Uint QdPointSetP10LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP10LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP10LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP10LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP10LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                     const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -1.0;
  coords(1, XI0) = -0.830223896278566929872;
  coords(2, XI0) = -0.4688487934707142138038;
  coords(3, XI0) = 0.0;
  coords(4, XI0) = 0.4688487934707142138038;
  coords(5, XI0) = 0.830223896278566929872;
  coords(6, XI0) = 1.0;
}

void QdPointSetP10LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);

  wgt[0] = 0.04761904761904761904762;
  wgt[1] = 0.276826047361565948011;
  wgt[2] = 0.4317453812098626234179;
  wgt[3] = 0.4876190476190476190476;
  wgt[4] = 0.4317453812098626234179;
  wgt[5] = 0.276826047361565948011;
  wgt[6] = 0.04761904761904761904762;
}

void QdPointSetP10LineGaussLobatto::permutation(const Uint local_id,
                                                const mesh::EntityRealignCode &permutation_code,
                                                std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P11 polynomials
// ----------------------------------------------------------------------------

QdPointSetP11LineGaussLobatto::QdPointSetP11LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP11LineGaussLobatto::~QdPointSetP11LineGaussLobatto()
{
}

Uint QdPointSetP11LineGaussLobatto::order() const
{
  return P11;
}

Uint QdPointSetP11LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP11LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP11LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP11LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP11LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                     const Uint local_idx) const
{
  QdPointSetP10LineGaussLobatto quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP11LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP10LineGaussLobatto quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP11LineGaussLobatto::permutation(const Uint local_id,
                                                const mesh::EntityRealignCode &permutation_code,
                                                std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P12 polynomials
// ----------------------------------------------------------------------------

QdPointSetP12LineGaussLobatto::QdPointSetP12LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP12LineGaussLobatto::~QdPointSetP12LineGaussLobatto()
{
}

Uint QdPointSetP12LineGaussLobatto::order() const
{
  return P12;
}

Uint QdPointSetP12LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP12LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP12LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP12LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP12LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                     const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -1.0;
  coords(1, XI0) = -0.8717401485096066153374;
  coords(2, XI0) = -0.5917001814331423021445;
  coords(3, XI0) = -0.2092992179024788687687;
  coords(4, XI0) = 0.2092992179024788687687;
  coords(5, XI0) = 0.5917001814331423021445;
  coords(6, XI0) = 0.8717401485096066153374;
  coords(7, XI0) = 1.0;
}

void QdPointSetP12LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.03571428571428571428571;
  wgt[1] = 0.210704227143506039383;
  wgt[2] = 0.3411226924835043647642;
  wgt[3] = 0.412458794658703881567;
  wgt[4] = 0.412458794658703881567;
  wgt[5] = 0.3411226924835043647642;
  wgt[6] = 0.210704227143506039383;
  wgt[7] = 0.03571428571428571428571;
}

void QdPointSetP12LineGaussLobatto::permutation(const Uint local_id,
                                                const mesh::EntityRealignCode &permutation_code,
                                                std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P13 polynomials
// ----------------------------------------------------------------------------

QdPointSetP13LineGaussLobatto::QdPointSetP13LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP13LineGaussLobatto::~QdPointSetP13LineGaussLobatto()
{
}

Uint QdPointSetP13LineGaussLobatto::order() const
{
  return P13;
}

Uint QdPointSetP13LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP13LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP13LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP13LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP13LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                     const Uint local_idx) const
{
  QdPointSetP12LineGaussLobatto quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP13LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP12LineGaussLobatto quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP13LineGaussLobatto::permutation(const Uint local_id,
                                                const mesh::EntityRealignCode &permutation_code,
                                                std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P14 polynomials
// ----------------------------------------------------------------------------

QdPointSetP14LineGaussLobatto::QdPointSetP14LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP14LineGaussLobatto::~QdPointSetP14LineGaussLobatto()
{
}

Uint QdPointSetP14LineGaussLobatto::order() const
{
  return P14;
}

Uint QdPointSetP14LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP14LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP14LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP14LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP14LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                     const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -1.0;
  coords(1, XI0) = -0.8997579954114601573123;
  coords(2, XI0) = -0.6771862795107377534459;
  coords(3, XI0) = -0.3631174638261781587108;
  coords(4, XI0) = 0.0;
  coords(5, XI0) = 0.3631174638261781587108;
  coords(6, XI0) = 0.6771862795107377534459;
  coords(7, XI0) = 0.8997579954114601573123;
  coords(8, XI0) = 1.0;
}

void QdPointSetP14LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.02777777777777777777778;
  wgt[1] = 0.165495361560805525046;
  wgt[2] = 0.2745387125001617352807;
  wgt[3] = 0.346428510973046345115;
  wgt[4] = 0.3715192743764172335601;
  wgt[5] = 0.346428510973046345115;
  wgt[6] = 0.2745387125001617352807;
  wgt[7] = 0.165495361560805525046;
  wgt[8] = 0.02777777777777777777778;
}

void QdPointSetP14LineGaussLobatto::permutation(const Uint local_id,
                                                const mesh::EntityRealignCode &permutation_code,
                                                std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P15 polynomials
// ----------------------------------------------------------------------------

QdPointSetP15LineGaussLobatto::QdPointSetP15LineGaussLobatto() : StdPointSetBase()
{
}

QdPointSetP15LineGaussLobatto::~QdPointSetP15LineGaussLobatto()
{
}

Uint QdPointSetP15LineGaussLobatto::order() const
{
  return P15;
}

Uint QdPointSetP15LineGaussLobatto::dim() const
{
  return _1D;
}

Uint QdPointSetP15LineGaussLobatto::codim() const
{
  return _0D;
}

Uint QdPointSetP15LineGaussLobatto::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP15LineGaussLobatto::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP15LineGaussLobatto::reference_coords(math::DenseDMat<Real> &coords,
                                                     const Uint local_idx) const
{
  QdPointSetP14LineGaussLobatto quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP15LineGaussLobatto::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP14LineGaussLobatto quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP15LineGaussLobatto::permutation(const Uint local_id,
                                                const mesh::EntityRealignCode &permutation_code,
                                                std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
