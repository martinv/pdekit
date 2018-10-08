#include "mesh/point_set/QdPointSetQuadGauss.hpp"
#include "common/Constants.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Quad quadrature for P1 polynomials
// Obtained from 1D Gauss quadrature (see section 4.1.4 in the book of Solin)
// ----------------------------------------------------------------------------

QdPointSetP1QuadGauss::QdPointSetP1QuadGauss() : StdPointSetBase()
{
}

QdPointSetP1QuadGauss::~QdPointSetP1QuadGauss()
{
}

Uint QdPointSetP1QuadGauss::order() const
{
  return P1;
}

Uint QdPointSetP1QuadGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP1QuadGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP1QuadGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP1QuadGauss::size(const Uint local_idx) const
{
  return 1u;
}

void QdPointSetP1QuadGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  /*
  coords.resize(4, 2);
  coords(0, XI0) = -0.57735026918962576451;
  coords(0, XI1) = -0.57735026918962576451;
  coords(1, XI0) = 0.57735026918962576451;
  coords(1, XI1) = -0.57735026918962576451;
  coords(2, XI0) = 0.57735026918962576451;
  coords(2, XI1) = 0.57735026918962576451;
  coords(3, XI0) = -0.57735026918962576451;
  coords(3, XI1) = 0.57735026918962576451;
  */

  coords.resize(1, 2);
  coords(0, XI0) = 0.0;
  coords(0, XI1) = 0.0;
}

void QdPointSetP1QuadGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  /*
  wgt.resize(4);
  wgt[0] = 1.0;
  wgt[1] = 1.0;
  wgt[2] = 1.0;
  wgt[3] = 1.0;
  */

  wgt.resize(1);
  wgt[0] = 4.0;
}

void QdPointSetP1QuadGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Quad quadrature for P2 polynomials
// Obtained from 1D Gauss quadrature (see section 4.1.4 in the book of Solin)
// ----------------------------------------------------------------------------

QdPointSetP2QuadGauss::QdPointSetP2QuadGauss() : StdPointSetBase()
{
}

QdPointSetP2QuadGauss::~QdPointSetP2QuadGauss()
{
}

Uint QdPointSetP2QuadGauss::order() const
{
  return P2;
}

Uint QdPointSetP2QuadGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP2QuadGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP2QuadGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP2QuadGauss::size(const Uint local_idx) const
{
  return 4u;
}

void QdPointSetP2QuadGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  /*
  coords.resize(9, 2);
  coords(0, XI0) = -0.77459666924148337704;
  coords(0, XI1) = -0.77459666924148337704;
  coords(1, XI0) = 0.0;
  coords(1, XI1) = -0.77459666924148337704;
  coords(2, XI0) = 0.77459666924148337704;
  coords(2, XI1) = -0.77459666924148337704;
  coords(3, XI0) = -0.77459666924148337704;
  coords(3, XI1) = 0.0;
  coords(4, XI0) = 0.0;
  coords(4, XI1) = 0.0;
  coords(5, XI0) = 0.77459666924148337704;
  coords(5, XI1) = 0.0;
  coords(6, XI0) = -0.77459666924148337704;
  coords(6, XI1) = 0.77459666924148337704;
  coords(7, XI0) = 0.0;
  coords(7, XI1) = 0.77459666924148337704;
  coords(8, XI0) = 0.77459666924148337704;
  coords(8, XI1) = 0.77459666924148337704;
  */

  coords.resize(4, 2);
  coords(0, XI0) = 0.577350269189626;
  coords(0, XI1) = 0.577350269189626;
  coords(1, XI0) = 0.577350269189626;
  coords(1, XI1) = -0.577350269189626;
  coords(2, XI0) = -0.577350269189626;
  coords(2, XI1) = 0.577350269189626;
  coords(3, XI0) = -0.577350269189626;
  coords(3, XI1) = -0.577350269189626;
}

void QdPointSetP2QuadGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  /*
  wgt.resize(9);
  wgt[0] = 10. / 27.;
  wgt[1] = 16. / 27.;
  wgt[2] = 10. / 27.;
  wgt[3] = 10. / 27.;
  wgt[4] = 16. / 27.;
  wgt[5] = 10. / 27.;
  wgt[6] = 10. / 27.;
  wgt[7] = 16. / 27.;
  wgt[8] = 10. / 27.;
  */

  wgt.resize(4);
  wgt[0] = 1.0;
  wgt[1] = 1.0;
  wgt[2] = 1.0;
  wgt[3] = 1.0;
}

void QdPointSetP2QuadGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Quad quadrature for P3 polynomials
// Obtained from 1D Gauss quadrature (see section 4.1.4 in the book of Solin)
// ----------------------------------------------------------------------------

QdPointSetP3QuadGauss::QdPointSetP3QuadGauss() : StdPointSetBase()
{
}

QdPointSetP3QuadGauss::~QdPointSetP3QuadGauss()
{
}

Uint QdPointSetP3QuadGauss::order() const
{
  return P3;
}

Uint QdPointSetP3QuadGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP3QuadGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP3QuadGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP3QuadGauss::size(const Uint local_idx) const
{
  return 4u;
}

void QdPointSetP3QuadGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP2QuadGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP3QuadGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP2QuadGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP3QuadGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Quad quadrature for P4 polynomials
// Obtained from 1D Gauss quadrature (see section 4.1.4 in the book of Solin)
// ----------------------------------------------------------------------------

QdPointSetP4QuadGauss::QdPointSetP4QuadGauss() : StdPointSetBase()
{
}

QdPointSetP4QuadGauss::~QdPointSetP4QuadGauss()
{
}

Uint QdPointSetP4QuadGauss::order() const
{
  return P4;
}

Uint QdPointSetP4QuadGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP4QuadGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP4QuadGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP4QuadGauss::size(const Uint local_idx) const
{
  return 8u; /*return 16u;*/
}

void QdPointSetP4QuadGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  /*
  coords.resize(16, 2);
  coords(0, XI0) = -0.861136311594953;
  coords(0, XI1) = -0.861136311594953;
  coords(1, XI0) = -0.339981043584856;
  coords(1, XI1) = -0.861136311594953;
  coords(2, XI0) = 0.339981043584856;
  coords(2, XI1) = -0.861136311594953;
  coords(3, XI0) = 0.861136311594953;
  coords(3, XI1) = -0.861136311594953;
  coords(4, XI0) = -0.861136311594953;
  coords(4, XI1) = -0.339981043584856;
  coords(5, XI0) = -0.339981043584856;
  coords(5, XI1) = -0.339981043584856;
  coords(6, XI0) = 0.339981043584856;
  coords(6, XI1) = -0.339981043584856;
  coords(7, XI0) = 0.861136311594953;
  coords(7, XI1) = -0.339981043584856;
  coords(8, XI0) = -0.861136311594953;
  coords(8, XI1) = 0.339981043584856;
  coords(9, XI0) = -0.339981043584856;
  coords(9, XI1) = 0.339981043584856;
  coords(10, XI0) = 0.339981043584856;
  coords(10, XI1) = 0.3399810435848563;
  coords(11, XI0) = 0.861136311594953;
  coords(11, XI1) = 0.339981043584856;
  coords(12, XI0) = -0.861136311594953;
  coords(12, XI1) = 0.861136311594953;
  coords(13, XI0) = -0.339981043584856;
  coords(13, XI1) = 0.861136311594953;
  coords(14, XI0) = 0.339981043584856;
  coords(14, XI1) = 0.861136311594953;
  coords(15, XI0) = 0.861136311594953;
  coords(15, XI1) = 0.861136311594953;
  */

  coords.resize(8, 2);
  coords(0, XI0) = 0.683130051063973;
  coords(0, XI1) = 0.0;

  coords(1, XI0) = -0.683130051063973;
  coords(1, XI1) = 0.0;

  coords(2, XI0) = 0.0;
  coords(2, XI1) = 0.683130051063973;

  coords(3, XI0) = 0.0;
  coords(3, XI1) = -0.683130051063973;

  coords(4, XI0) = 0.881917103688197;
  coords(4, XI1) = 0.881917103688197;

  coords(5, XI0) = 0.881917103688197;
  coords(5, XI1) = -0.881917103688197;

  coords(6, XI0) = -0.881917103688197;
  coords(6, XI1) = 0.881917103688197;

  coords(7, XI0) = -0.881917103688197;
  coords(7, XI1) = -0.881917103688197;
}

void QdPointSetP4QuadGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  /*
  wgt.resize(16);
  wgt[0] = 0.347854845137454 * 0.347854845137454;
  wgt[1] = 0.652145154862546 * 0.347854845137454;
  wgt[2] = 0.652145154862546 * 0.347854845137454;
  wgt[3] = 0.347854845137454 * 0.347854845137454;
  wgt[4] = 0.347854845137454 * 0.652145154862546;
  wgt[5] = 0.652145154862546 * 0.652145154862546;
  wgt[6] = 0.652145154862546 * 0.652145154862546;
  wgt[7] = 0.347854845137454 * 0.652145154862546;
  wgt[8] = 0.347854845137454 * 0.652145154862546;
  wgt[9] = 0.652145154862546 * 0.652145154862546;
  wgt[10] = 0.652145154862546 * 0.652145154862546;
  wgt[11] = 0.347854845137454 * 0.652145154862546;
  wgt[12] = 0.347854845137454 * 0.347854845137454;
  wgt[13] = 0.652145154862546 * 0.347854845137454;
  wgt[14] = 0.652145154862546 * 0.347854845137454;
  wgt[15] = 0.347854845137454 * 0.347854845137454;
  */

  wgt.resize(8);
  wgt[0] = 0.816326530612245;
  wgt[1] = 0.816326530612245;
  wgt[2] = 0.816326530612245;
  wgt[3] = 0.816326530612245;
  wgt[4] = 0.183673469387755;
  wgt[5] = 0.183673469387755;
  wgt[6] = 0.183673469387755;
  wgt[7] = 0.183673469387755;
}

void QdPointSetP4QuadGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Quad quadrature for P5 polynomials
// See section 4.2.3 in the book of Solin
// ----------------------------------------------------------------------------

QdPointSetP5QuadGauss::QdPointSetP5QuadGauss() : StdPointSetBase()
{
}

QdPointSetP5QuadGauss::~QdPointSetP5QuadGauss()
{
}

Uint QdPointSetP5QuadGauss::order() const
{
  return P5;
}

Uint QdPointSetP5QuadGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP5QuadGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP5QuadGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP5QuadGauss::size(const Uint local_idx) const
{
  return 8u;
}

void QdPointSetP5QuadGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP4QuadGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP5QuadGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP4QuadGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP5QuadGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Quad quadrature for P6 polynomials
// See section 4.1.4 in the book of Solin
// ----------------------------------------------------------------------------

QdPointSetP6QuadGauss::QdPointSetP6QuadGauss() : StdPointSetBase()
{
}

QdPointSetP6QuadGauss::~QdPointSetP6QuadGauss()
{
}

Uint QdPointSetP6QuadGauss::order() const
{
  return P6;
}

Uint QdPointSetP6QuadGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP6QuadGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP6QuadGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP6QuadGauss::size(const Uint local_idx) const
{
  return 12u;
}

void QdPointSetP6QuadGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  coords.resize(12, 2);
  coords(0, XI0) = 0.925820099772551;
  coords(0, XI1) = 0.0;

  coords(1, XI0) = -0.925820099772551;
  coords(1, XI1) = 0.0;

  coords(2, XI0) = 0.0;
  coords(2, XI1) = 0.925820099772551;

  coords(3, XI0) = 0.0;
  coords(3, XI1) = -0.925820099772551;

  coords(4, XI0) = 0.805979782918599;
  coords(4, XI1) = 0.805979782918599;

  coords(5, XI0) = 0.805979782918599;
  coords(5, XI1) = -0.805979782918599;

  coords(6, XI0) = -0.805979782918599;
  coords(6, XI1) = 0.805979782918599;

  coords(7, XI0) = -0.805979782918599;
  coords(7, XI1) = -0.805979782918599;

  coords(8, XI0) = 0.380554433208316;
  coords(8, XI1) = 0.380554433208316;

  coords(9, XI0) = 0.380554433208316;
  coords(9, XI1) = -0.380554433208316;

  coords(10, XI0) = -0.380554433208316;
  coords(10, XI1) = 0.380554433208316;

  coords(11, XI0) = -0.380554433208316;
  coords(11, XI1) = -0.380554433208316;
}

void QdPointSetP6QuadGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(12);
  wgt[0]  = 0.241975308641975;
  wgt[1]  = 0.241975308641975;
  wgt[2]  = 0.241975308641975;
  wgt[3]  = 0.241975308641975;
  wgt[4]  = 0.237431774690630;
  wgt[5]  = 0.237431774690630;
  wgt[6]  = 0.237431774690630;
  wgt[7]  = 0.237431774690630;
  wgt[8]  = 0.520592916667394;
  wgt[9]  = 0.520592916667394;
  wgt[10] = 0.520592916667394;
  wgt[11] = 0.520592916667394;
}

void QdPointSetP6QuadGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Quad quadrature for P7 polynomials
// See section 4.1.4 in the book of Solin
// ----------------------------------------------------------------------------

QdPointSetP7QuadGauss::QdPointSetP7QuadGauss() : StdPointSetBase()
{
}

QdPointSetP7QuadGauss::~QdPointSetP7QuadGauss()
{
}

Uint QdPointSetP7QuadGauss::order() const
{
  return P7;
}

Uint QdPointSetP7QuadGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP7QuadGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP7QuadGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP7QuadGauss::size(const Uint local_idx) const
{
  return 12u;
}

void QdPointSetP7QuadGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP6QuadGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP7QuadGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP6QuadGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP7QuadGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Quad quadrature for P7 polynomials
// Copied from quadrature tables supplied with the book of Solin
// ----------------------------------------------------------------------------

QdPointSetP8QuadGauss::QdPointSetP8QuadGauss() : StdPointSetBase()
{
}

QdPointSetP8QuadGauss::~QdPointSetP8QuadGauss()
{
}

Uint QdPointSetP8QuadGauss::order() const
{
  return P8;
}

Uint QdPointSetP8QuadGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP8QuadGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP8QuadGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP8QuadGauss::size(const Uint local_idx) const
{
  return 20u;
}

void QdPointSetP8QuadGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  coords.resize(20, 2);
  coords(0, XI0) = 1.121225763866564;
  coords(0, XI1) = 0.000000000000000;

  coords(1, XI0) = -1.121225763866564;
  coords(1, XI1) = 0.000000000000000;

  coords(2, XI0) = 0.000000000000000;
  coords(2, XI1) = 1.121225763866564;

  coords(3, XI0) = 0.000000000000000;
  coords(3, XI1) = -1.121225763866564;

  coords(4, XI0) = 0.451773049920657;
  coords(4, XI1) = 0.000000000000000;

  coords(5, XI0) = -0.451773049920657;
  coords(5, XI1) = 0.000000000000000;

  coords(6, XI0) = 0.000000000000000;
  coords(6, XI1) = 0.451773049920657;

  coords(7, XI0) = 0.000000000000000;
  coords(7, XI1) = -0.451773049920657;

  coords(8, XI0) = 0.891849420851512;
  coords(8, XI1) = 0.891849420851512;

  coords(9, XI0) = 0.891849420851512;
  coords(9, XI1) = -0.891849420851512;

  coords(10, XI0) = -0.891849420851512;
  coords(10, XI1) = 0.891849420851512;

  coords(11, XI0) = -0.891849420851512;
  coords(11, XI1) = -0.891849420851512;

  coords(12, XI0) = 0.824396370749276;
  coords(12, XI1) = 0.411623426336542;

  coords(13, XI0) = 0.824396370749276;
  coords(13, XI1) = -0.411623426336542;

  coords(14, XI0) = -0.824396370749276;
  coords(14, XI1) = 0.411623426336542;

  coords(15, XI0) = -0.824396370749276;
  coords(15, XI1) = -0.411623426336542;

  coords(16, XI0) = 0.411623426336542;
  coords(16, XI1) = 0.824396370749276;

  coords(17, XI0) = 0.411623426336542;
  coords(17, XI1) = -0.824396370749276;

  coords(18, XI0) = -0.411623426336542;
  coords(18, XI1) = 0.824396370749276;

  coords(19, XI0) = -0.411623426336542;
  coords(19, XI1) = -0.824396370749276;
}

void QdPointSetP8QuadGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(20);
  wgt[0]  = 0.018475842507491;
  wgt[1]  = 0.018475842507491;
  wgt[2]  = 0.018475842507491;
  wgt[3]  = 0.018475842507491;
  wgt[4]  = 0.390052939160735;
  wgt[5]  = 0.390052939160735;
  wgt[6]  = 0.390052939160735;
  wgt[7]  = 0.390052939160735;
  wgt[8]  = 0.083095178026482;
  wgt[9]  = 0.083095178026482;
  wgt[10] = 0.083095178026482;
  wgt[11] = 0.083095178026482;
  wgt[12] = 0.254188020152646;
  wgt[13] = 0.254188020152646;
  wgt[14] = 0.254188020152646;
  wgt[15] = 0.254188020152646;
  wgt[16] = 0.254188020152646;
  wgt[17] = 0.254188020152646;
  wgt[18] = 0.254188020152646;
  wgt[19] = 0.254188020152646;
}

void QdPointSetP8QuadGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
