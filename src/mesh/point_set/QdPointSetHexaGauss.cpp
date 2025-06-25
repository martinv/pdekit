#include "mesh/point_set/QdPointSetHexaGauss.hpp"
#include "common/Constants.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Hexa quadrature for P1 polynomials
// ----------------------------------------------------------------------------

QdPointSetP1HexaGauss::QdPointSetP1HexaGauss() : StdPointSetBase()
{
}

QdPointSetP1HexaGauss::~QdPointSetP1HexaGauss()
{
}

Uint QdPointSetP1HexaGauss::order() const
{
  return P1;
}

Uint QdPointSetP1HexaGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP1HexaGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP1HexaGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP1HexaGauss::size(const Uint local_idx) const
{
  return 1u;
}

void QdPointSetP1HexaGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  coords.resize(1, 3);
  coords(0, XI0) = 0.0;
  coords(0, XI1) = 0.0;
  coords(0, XI2) = 0.0;
}

void QdPointSetP1HexaGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(1);
  wgt[0] = 8.0;
}

void QdPointSetP1HexaGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Hexa quadrature for P2 polynomials
// ----------------------------------------------------------------------------

QdPointSetP2HexaGauss::QdPointSetP2HexaGauss() : StdPointSetBase()
{
}

QdPointSetP2HexaGauss::~QdPointSetP2HexaGauss()
{
}

Uint QdPointSetP2HexaGauss::order() const
{
  return P2;
}

Uint QdPointSetP2HexaGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP2HexaGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP2HexaGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP2HexaGauss::size(const Uint local_idx) const
{
  return 6u;
}

void QdPointSetP2HexaGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  coords.resize(6, 3);
  coords(0, XI0) = 1.0;
  coords(0, XI1) = 0.0;
  coords(0, XI2) = 0.0;

  coords(1, XI0) = -1.0;
  coords(1, XI1) = 0.0;
  coords(1, XI2) = 0.0;

  coords(2, XI0) = 0.0;
  coords(2, XI1) = 1.0;
  coords(2, XI2) = 0.0;

  coords(3, XI0) = 0.0;
  coords(3, XI1) = -1.0;
  coords(3, XI2) = 0.0;

  coords(4, XI0) = 0.0;
  coords(4, XI1) = 0.0;
  coords(4, XI2) = 1.0;

  coords(5, XI0) = 0.0;
  coords(5, XI1) = 0.0;
  coords(5, XI2) = -1.0;
}

void QdPointSetP2HexaGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(6);
  wgt[0] = 4. / 3.;
  wgt[1] = 4. / 3.;
  wgt[2] = 4. / 3.;
  wgt[3] = 4. / 3.;
  wgt[4] = 4. / 3.;
  wgt[5] = 4. / 3.;
}

void QdPointSetP2HexaGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Hexa quadrature for P3 polynomials
// ----------------------------------------------------------------------------

QdPointSetP3HexaGauss::QdPointSetP3HexaGauss() : StdPointSetBase()
{
}

QdPointSetP3HexaGauss::~QdPointSetP3HexaGauss()
{
}

Uint QdPointSetP3HexaGauss::order() const
{
  return P3;
}

Uint QdPointSetP3HexaGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP3HexaGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP3HexaGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP3HexaGauss::size(const Uint local_idx) const
{
  return 6u;
}

void QdPointSetP3HexaGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP2HexaGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP3HexaGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP2HexaGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP3HexaGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Hexa quadrature for P4 polynomials
// ----------------------------------------------------------------------------

QdPointSetP4HexaGauss::QdPointSetP4HexaGauss() : StdPointSetBase()
{
}

QdPointSetP4HexaGauss::~QdPointSetP4HexaGauss()
{
}

Uint QdPointSetP4HexaGauss::order() const
{
  return P4;
}

Uint QdPointSetP4HexaGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP4HexaGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP4HexaGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP4HexaGauss::size(const Uint local_idx) const
{
  return 14u;
}

void QdPointSetP4HexaGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  coords.resize(14, 3);
  coords(0, XI0) = 0.7958224257;
  coords(0, XI1) = 0.0;
  coords(0, XI2) = 0.0;

  coords(1, XI0) = -0.7958224257;
  coords(1, XI1) = 0.0;
  coords(1, XI2) = 0.0;

  coords(2, XI0) = 0.0;
  coords(2, XI1) = 0.7958224257;
  coords(2, XI2) = 0.0;

  coords(3, XI0) = 0.0;
  coords(3, XI1) = -0.7958224257;
  coords(3, XI2) = 0.0;

  coords(4, XI0) = 0.0;
  coords(4, XI1) = 0.0;
  coords(4, XI2) = 0.7958224257;

  coords(5, XI0) = 0.0;
  coords(5, XI1) = 0.0;
  coords(5, XI2) = -0.7958224257;

  coords(6, XI0) = 0.7587869106;
  coords(6, XI1) = 0.7587869106;
  coords(6, XI2) = 0.7587869106;

  coords(7, XI0) = 0.7587869106;
  coords(7, XI1) = -0.7587869106;
  coords(7, XI2) = 0.7587869106;

  coords(8, XI0) = 0.7587869106;
  coords(8, XI1) = 0.7587869106;
  coords(8, XI2) = -0.7587869106;

  coords(9, XI0) = 0.7587869106;
  coords(9, XI1) = -0.7587869106;
  coords(9, XI2) = -0.7587869106;

  coords(10, XI0) = -0.7587869106;
  coords(10, XI1) = 0.7587869106;
  coords(10, XI2) = 0.7587869106;

  coords(11, XI0) = -0.7587869106;
  coords(11, XI1) = -0.7587869106;
  coords(11, XI2) = 0.7587869106;

  coords(12, XI0) = -0.7587869106;
  coords(12, XI1) = 0.7587869106;
  coords(12, XI2) = -0.7587869106;

  coords(13, XI0) = -0.7587869106;
  coords(13, XI1) = -0.7587869106;
  coords(13, XI2) = -0.7587869106;
}

void QdPointSetP4HexaGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(14);
  wgt[0]  = 0.8864265927;
  wgt[1]  = 0.8864265927;
  wgt[2]  = 0.8864265927;
  wgt[3]  = 0.8864265927;
  wgt[4]  = 0.8864265927;
  wgt[5]  = 0.8864265927;
  wgt[6]  = 0.3351800554;
  wgt[7]  = 0.3351800554;
  wgt[8]  = 0.3351800554;
  wgt[9]  = 0.3351800554;
  wgt[10] = 0.3351800554;
  wgt[11] = 0.3351800554;
  wgt[12] = 0.3351800554;
  wgt[13] = 0.3351800554;
}

void QdPointSetP4HexaGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Hexa quadrature for P5 polynomials
// ----------------------------------------------------------------------------

QdPointSetP5HexaGauss::QdPointSetP5HexaGauss() : StdPointSetBase()
{
}

QdPointSetP5HexaGauss::~QdPointSetP5HexaGauss()
{
}

Uint QdPointSetP5HexaGauss::order() const
{
  return P5;
}

Uint QdPointSetP5HexaGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP5HexaGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP5HexaGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP5HexaGauss::size(const Uint local_idx) const
{
  return 14u;
}

void QdPointSetP5HexaGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP4HexaGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP5HexaGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP4HexaGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP5HexaGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Hexa quadrature for P6 polynomials
// ----------------------------------------------------------------------------

QdPointSetP6HexaGauss::QdPointSetP6HexaGauss() : StdPointSetBase()
{
}

QdPointSetP6HexaGauss::~QdPointSetP6HexaGauss()
{
}

Uint QdPointSetP6HexaGauss::order() const
{
  return P6;
}

Uint QdPointSetP6HexaGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP6HexaGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP6HexaGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP6HexaGauss::size(const Uint local_idx) const
{
  return 27u;
}

void QdPointSetP6HexaGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  coords.resize(27, 3);
  coords(0, XI0) = 0.0;
  coords(0, XI1) = 0.0;
  coords(0, XI2) = 0.0;

  coords(1, XI0) = 0.8484180014;
  coords(1, XI1) = 0.0;
  coords(1, XI2) = 0.0;

  coords(2, XI0) = -0.8484180014;
  coords(2, XI1) = 0.0;
  coords(2, XI2) = 0.0;

  coords(3, XI0) = 0.0;
  coords(3, XI1) = 0.8484180014;
  coords(3, XI2) = 0.0;

  coords(4, XI0) = 0.0;
  coords(4, XI1) = -0.8484180014;
  coords(4, XI2) = 0.0;

  coords(5, XI0) = 0.0;
  coords(5, XI1) = 0.0;
  coords(5, XI2) = 0.8484180014;

  coords(6, XI0) = 0.0;
  coords(6, XI1) = 0.0;
  coords(6, XI2) = -0.8484180014;

  coords(7, XI0) = 0.6528164721;
  coords(7, XI1) = 0.6528164721;
  coords(7, XI2) = 0.6528164721;

  coords(8, XI0) = 0.6528164721;
  coords(8, XI1) = -0.6528164721;
  coords(8, XI2) = 0.6528164721;

  coords(9, XI0) = 0.6528164721;
  coords(9, XI1) = 0.6528164721;
  coords(9, XI2) = -0.6528164721;

  coords(10, XI0) = 0.6528164721;
  coords(10, XI1) = -0.6528164721;
  coords(10, XI2) = -0.6528164721;

  coords(11, XI0) = -0.6528164721;
  coords(11, XI1) = 0.6528164721;
  coords(11, XI2) = 0.6528164721;

  coords(12, XI0) = -0.6528164721;
  coords(12, XI1) = -0.6528164721;
  coords(12, XI2) = 0.6528164721;

  coords(13, XI0) = -0.6528164721;
  coords(13, XI1) = 0.6528164721;
  coords(13, XI2) = -0.6528164721;

  coords(14, XI0) = -0.6528164721;
  coords(14, XI1) = -0.6528164721;
  coords(14, XI2) = -0.6528164721;

  coords(15, XI0) = 0.0;
  coords(15, XI1) = 1.1064128986;
  coords(15, XI2) = 1.1064128986;

  coords(16, XI0) = 0.0;
  coords(16, XI1) = -1.1064128986;
  coords(16, XI2) = 1.1064128986;

  coords(17, XI0) = 0.0;
  coords(17, XI1) = 1.1064128986;
  coords(17, XI2) = -1.1064128986;

  coords(18, XI0) = 0.0;
  coords(18, XI1) = -1.1064128986;
  coords(18, XI2) = -1.1064128986;

  coords(19, XI0) = 1.1064128986;
  coords(19, XI1) = 0.0;
  coords(19, XI2) = 1.1064128986;

  coords(20, XI0) = -1.1064128986;
  coords(20, XI1) = 0.0;
  coords(20, XI2) = 1.1064128986;

  coords(21, XI0) = 1.1064128986;
  coords(21, XI1) = 0.0;
  coords(21, XI2) = -1.1064128986;

  coords(22, XI0) = -1.1064128986;
  coords(22, XI1) = 0.0;
  coords(22, XI2) = -1.1064128986;

  coords(23, XI0) = 1.1064128986;
  coords(23, XI1) = 1.1064128986;
  coords(23, XI2) = 0.0;

  coords(24, XI0) = -1.1064128986;
  coords(24, XI1) = 1.1064128986;
  coords(24, XI2) = 0.0;

  coords(25, XI0) = 1.1064128986;
  coords(25, XI1) = -1.1064128986;
  coords(25, XI2) = 0.0;

  coords(26, XI0) = -1.1064128986;
  coords(26, XI1) = -1.1064128986;
  coords(26, XI2) = 0.0;
}

void QdPointSetP6HexaGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(27);
  wgt[0]  = 0.7880734827;
  wgt[1]  = 0.4993690023;
  wgt[2]  = 0.4993690023;
  wgt[3]  = 0.4993690023;
  wgt[4]  = 0.4993690023;
  wgt[5]  = 0.4993690023;
  wgt[6]  = 0.4993690023;
  wgt[7]  = 0.4785084494;
  wgt[8]  = 0.4785084494;
  wgt[9]  = 0.4785084494;
  wgt[10] = 0.4785084494;
  wgt[11] = 0.4785084494;
  wgt[12] = 0.4785084494;
  wgt[13] = 0.4785084494;
  wgt[14] = 0.4785084494;
  wgt[15] = 0.0323037423;
  wgt[16] = 0.0323037423;
  wgt[17] = 0.0323037423;
  wgt[18] = 0.0323037423;
  wgt[19] = 0.0323037423;
  wgt[20] = 0.0323037423;
  wgt[21] = 0.0323037423;
  wgt[22] = 0.0323037423;
  wgt[23] = 0.0323037423;
  wgt[24] = 0.0323037423;
  wgt[25] = 0.0323037423;
  wgt[26] = 0.0323037423;
}

void QdPointSetP6HexaGauss::permutation(const Uint local_id,
                                        const mesh::EntityRealignCode &permutation_code,
                                        std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ----------------------------------------------------------------------------
// Hexa quadrature for P7 polynomials
// ----------------------------------------------------------------------------

QdPointSetP7HexaGauss::QdPointSetP7HexaGauss() : StdPointSetBase()
{
}

QdPointSetP7HexaGauss::~QdPointSetP7HexaGauss()
{
}

Uint QdPointSetP7HexaGauss::order() const
{
  return P7;
}

Uint QdPointSetP7HexaGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP7HexaGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP7HexaGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP7HexaGauss::size(const Uint local_idx) const
{
  return 27u;
}

void QdPointSetP7HexaGauss::reference_coords(math::DenseDMat<Real> &coords,
                                             const Uint local_idx) const
{
  QdPointSetP6HexaGauss quad;
  quad.reference_coords(coords, local_idx);
}

void QdPointSetP7HexaGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  QdPointSetP6HexaGauss quad;
  quad.weights(wgt, local_idx);
}

void QdPointSetP7HexaGauss::permutation(const Uint local_id,
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
