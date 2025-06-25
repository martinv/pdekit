#include "mesh/point_set/QdPointSetLineGauss.hpp"
#include "common/Constants.hpp"
#include "mesh/point_set/QuadratureTransformUtils.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Line quadrature for P1 polynomials
// ----------------------------------------------------------------------------

QdPointSetP1LineGauss::QdPointSetP1LineGauss() : StdPointSetBase()
{
}

QdPointSetP1LineGauss::~QdPointSetP1LineGauss()
{
}

Uint QdPointSetP1LineGauss::order() const
{
  return P1;
}

Uint QdPointSetP1LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP1LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP1LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP1LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP1LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -0.5773502691896257645091488;
  coords(1, XI0) = 0.5773502691896257645091488;
}

void QdPointSetP1LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 1.0;
  wgt[1] = 1.0;
}

void QdPointSetP1LineGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P2 polynomials
// ----------------------------------------------------------------------------

QdPointSetP2LineGauss::QdPointSetP2LineGauss() : StdPointSetBase()
{
}

QdPointSetP2LineGauss::~QdPointSetP2LineGauss()
{
}

Uint QdPointSetP2LineGauss::order() const
{
  return P2;
}

Uint QdPointSetP2LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP2LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP2LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP2LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP2LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP1LineGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP2LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP1LineGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP2LineGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P3 polynomials
// ----------------------------------------------------------------------------

QdPointSetP3LineGauss::QdPointSetP3LineGauss() : StdPointSetBase()
{
}

QdPointSetP3LineGauss::~QdPointSetP3LineGauss()
{
}

Uint QdPointSetP3LineGauss::order() const
{
  return P3;
}

Uint QdPointSetP3LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP3LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP3LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP3LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP3LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP1LineGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP3LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP1LineGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP3LineGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P4 polynomials
// ----------------------------------------------------------------------------

QdPointSetP4LineGauss::QdPointSetP4LineGauss() : StdPointSetBase()
{
}

QdPointSetP4LineGauss::~QdPointSetP4LineGauss()
{
}

Uint QdPointSetP4LineGauss::order() const
{
  return P4;
}

Uint QdPointSetP4LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP4LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP4LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP4LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP4LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -0.7745966692414833770358531;
  coords(1, XI0) = 0.0;
  coords(2, XI0) = 0.7745966692414833770358531;
}

void QdPointSetP4LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 5. / 9.; // 0.555555555555555555
  wgt[1] = 8. / 9.; // 0.888888888888888888
  wgt[2] = 5. / 9.; // 0.555555555555555555
}

void QdPointSetP4LineGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P5 polynomials
// ----------------------------------------------------------------------------

QdPointSetP5LineGauss::QdPointSetP5LineGauss() : StdPointSetBase()
{
}

QdPointSetP5LineGauss::~QdPointSetP5LineGauss()
{
}

Uint QdPointSetP5LineGauss::order() const
{
  return P5;
}

Uint QdPointSetP5LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP5LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP5LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP5LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP5LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP4LineGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP5LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP4LineGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP5LineGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P6 polynomials
// ----------------------------------------------------------------------------

QdPointSetP6LineGauss::QdPointSetP6LineGauss() : StdPointSetBase()
{
}

QdPointSetP6LineGauss::~QdPointSetP6LineGauss()
{
}

Uint QdPointSetP6LineGauss::order() const
{
  return P6;
}

Uint QdPointSetP6LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP6LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP6LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP6LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP6LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -0.8611363115940525752239465;
  coords(1, XI0) = -0.3399810435848562648026658;
  coords(2, XI0) = 0.3399810435848562648026658;
  coords(3, XI0) = 0.8611363115940525752239465;
}

void QdPointSetP6LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.3478548451374538573730639;
  wgt[1] = 0.6521451548625461426269361;
  wgt[2] = 0.6521451548625461426269361;
  wgt[3] = 0.3478548451374538573730639;
}

void QdPointSetP6LineGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P7 polynomials
// ----------------------------------------------------------------------------

QdPointSetP7LineGauss::QdPointSetP7LineGauss() : StdPointSetBase()
{
}

QdPointSetP7LineGauss::~QdPointSetP7LineGauss()
{
}

Uint QdPointSetP7LineGauss::order() const
{
  return P7;
}

Uint QdPointSetP7LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP7LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP7LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP7LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP7LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP6LineGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP7LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP6LineGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP7LineGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P8 polynomials
// ----------------------------------------------------------------------------

QdPointSetP8LineGauss::QdPointSetP8LineGauss() : StdPointSetBase()
{
}

QdPointSetP8LineGauss::~QdPointSetP8LineGauss()
{
}

Uint QdPointSetP8LineGauss::order() const
{
  return P8;
}

Uint QdPointSetP8LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP8LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP8LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP8LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP8LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -0.9061798459386639927976269;
  coords(1, XI0) = -0.5384693101056830910363144;
  coords(2, XI0) = 0.0;
  coords(3, XI0) = 0.5384693101056830910363144;
  coords(4, XI0) = 0.9061798459386639927976269;
}

void QdPointSetP8LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.2369268850561890875142640;
  wgt[1] = 0.4786286704993664680412915;
  wgt[2] = 0.5688888888888888888888888;
  wgt[3] = 0.4786286704993664680412915;
  wgt[4] = 0.2369268850561890875142640;
}

void QdPointSetP8LineGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P9 polynomials
// ----------------------------------------------------------------------------

QdPointSetP9LineGauss::QdPointSetP9LineGauss() : StdPointSetBase()
{
}

QdPointSetP9LineGauss::~QdPointSetP9LineGauss()
{
}

Uint QdPointSetP9LineGauss::order() const
{
  return P9;
}

Uint QdPointSetP9LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP9LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP9LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP9LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP9LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP8LineGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP9LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP8LineGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP9LineGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P10 polynomials
// ----------------------------------------------------------------------------

QdPointSetP10LineGauss::QdPointSetP10LineGauss() : StdPointSetBase()
{
}

QdPointSetP10LineGauss::~QdPointSetP10LineGauss()
{
}

Uint QdPointSetP10LineGauss::order() const
{
  return P10;
}

Uint QdPointSetP10LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP10LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP10LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP10LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP10LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -0.9324695142031520278123016;
  coords(1, XI0) = -0.6612093864662645136613996;
  coords(2, XI0) = -0.2386191860831969086305017;
  coords(3, XI0) = 0.2386191860831969086305017;
  coords(4, XI0) = 0.6612093864662645136613996;
  coords(5, XI0) = 0.9324695142031520278123016;
}

void QdPointSetP10LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.1713244923791703450402961;
  wgt[1] = 0.3607615730481386075698335;
  wgt[2] = 0.4679139345726910473898703;
  wgt[3] = 0.4679139345726910473898703;
  wgt[4] = 0.3607615730481386075698335;
  wgt[5] = 0.1713244923791703450402961;
}

void QdPointSetP10LineGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P11 polynomials
// ----------------------------------------------------------------------------

QdPointSetP11LineGauss::QdPointSetP11LineGauss() : StdPointSetBase()
{
}

QdPointSetP11LineGauss::~QdPointSetP11LineGauss()
{
}

Uint QdPointSetP11LineGauss::order() const
{
  return P11;
}

Uint QdPointSetP11LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP11LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP11LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP11LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP11LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  QdPointSetP10LineGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP11LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP10LineGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP11LineGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P12 polynomials
// ----------------------------------------------------------------------------

QdPointSetP12LineGauss::QdPointSetP12LineGauss() : StdPointSetBase()
{
}

QdPointSetP12LineGauss::~QdPointSetP12LineGauss()
{
}

Uint QdPointSetP12LineGauss::order() const
{
  return P12;
}

Uint QdPointSetP12LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP12LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP12LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP12LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP12LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -0.9491079123427585245261897;
  coords(1, XI0) = -0.7415311855993944398638648;
  coords(2, XI0) = -0.4058451513773971669066064;
  coords(3, XI0) = 0.0;
  coords(4, XI0) = 0.4058451513773971669066064;
  coords(5, XI0) = 0.7415311855993944398638648;
  coords(6, XI0) = 0.9491079123427585245261897;
}

void QdPointSetP12LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.1294849661688696932706114;
  wgt[1] = 0.2797053914892766679014678;
  wgt[2] = 0.3818300505051189449503698;
  wgt[3] = 0.4179591836734693877551020;
  wgt[4] = 0.3818300505051189449503698;
  wgt[5] = 0.2797053914892766679014678;
  wgt[6] = 0.1294849661688696932706114;
}

void QdPointSetP12LineGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P13 polynomials
// ----------------------------------------------------------------------------

QdPointSetP13LineGauss::QdPointSetP13LineGauss() : StdPointSetBase()
{
}

QdPointSetP13LineGauss::~QdPointSetP13LineGauss()
{
}

Uint QdPointSetP13LineGauss::order() const
{
  return P13;
}

Uint QdPointSetP13LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP13LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP13LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP13LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP13LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  QdPointSetP12LineGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP13LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP12LineGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP13LineGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P14 polynomials
// ----------------------------------------------------------------------------

QdPointSetP14LineGauss::QdPointSetP14LineGauss() : StdPointSetBase()
{
}

QdPointSetP14LineGauss::~QdPointSetP14LineGauss()
{
}

Uint QdPointSetP14LineGauss::order() const
{
  return P14;
}

Uint QdPointSetP14LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP14LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP14LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP14LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP14LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(N_QD_PTS, 1);
  coords(0, XI0) = -0.9602898564975362316835609;
  coords(1, XI0) = -0.7966664774136267395915539;
  coords(2, XI0) = -0.5255324099163289858177390;
  coords(3, XI0) = -0.1834346424956498049394761;
  coords(4, XI0) = 0.1834346424956498049394761;
  coords(5, XI0) = 0.5255324099163289858177390;
  coords(6, XI0) = 0.7966664774136267395915539;
  coords(7, XI0) = 0.9602898564975362316835609;
}

void QdPointSetP14LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.1012285362903762591525314;
  wgt[1] = 0.2223810344533744705443560;
  wgt[2] = 0.3137066458778872873379622;
  wgt[3] = 0.3626837833783619829651504;
  wgt[4] = 0.3626837833783619829651504;
  wgt[5] = 0.3137066458778872873379622;
  wgt[6] = 0.2223810344533744705443560;
  wgt[7] = 0.1012285362903762591525314;
}

void QdPointSetP14LineGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------
// Line quadrature for P15 polynomials
// ----------------------------------------------------------------------------

QdPointSetP15LineGauss::QdPointSetP15LineGauss() : StdPointSetBase()
{
}

QdPointSetP15LineGauss::~QdPointSetP15LineGauss()
{
}

Uint QdPointSetP15LineGauss::order() const
{
  return P15;
}

Uint QdPointSetP15LineGauss::dim() const
{
  return _1D;
}

Uint QdPointSetP15LineGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP15LineGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP15LineGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP15LineGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  QdPointSetP14LineGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP15LineGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP14LineGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP15LineGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  detail::fill_line_quadrature_permutation(N_QD_PTS, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
