#include "mesh/point_set/QdPointSetTetraGauss.hpp"
#include "common/Constants.hpp"

namespace pdekit
{

namespace mesh
{

// ============================================================================
// Tetrahedral quadrature for P1 polynomials
// ============================================================================

QdPointSetP1TetraGauss::QdPointSetP1TetraGauss() : StdPointSetBase()
{
}

QdPointSetP1TetraGauss::~QdPointSetP1TetraGauss()
{
}

Uint QdPointSetP1TetraGauss::order() const
{
  return P1;
}

Uint QdPointSetP1TetraGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP1TetraGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP1TetraGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP1TetraGauss::size(const Uint local_idx) const
{
  return 1u;
}

void QdPointSetP1TetraGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(1, 3);
  coords(0, XI0) = -0.5;
  coords(0, XI1) = -0.5;
  coords(0, XI2) = -0.5;
}

void QdPointSetP1TetraGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(1);
  wgt[0] = 4. / 3.;
}

void QdPointSetP1TetraGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ============================================================================
// Tetrahedral quadrature for P2 polynomials
// ============================================================================

QdPointSetP2TetraGauss::QdPointSetP2TetraGauss() : StdPointSetBase()
{
}

QdPointSetP2TetraGauss::~QdPointSetP2TetraGauss()
{
}

Uint QdPointSetP2TetraGauss::order() const
{
  return P2;
}

Uint QdPointSetP2TetraGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP2TetraGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP2TetraGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP2TetraGauss::size(const Uint local_idx) const
{
  return 4u;
}

void QdPointSetP2TetraGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(4, 3);
  coords(0, XI0) = -0.7236067977;
  coords(0, XI1) = -0.7236067977;
  coords(0, XI2) = -0.7236067977;

  coords(1, XI0) = 0.1708203932;
  coords(1, XI1) = -0.7236067977;
  coords(1, XI2) = -0.7236067977;

  coords(2, XI0) = -0.7236067977;
  coords(2, XI1) = 0.1708203932;
  coords(2, XI2) = -0.7236067977;

  coords(3, XI0) = -0.7236067977;
  coords(3, XI1) = -0.7236067977;
  coords(3, XI2) = 0.1708203932;
}

void QdPointSetP2TetraGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(4);
  wgt[0] = 1. / 3.;
  wgt[1] = 1. / 3.;
  wgt[2] = 1. / 3.;
  wgt[3] = 1. / 3.;
}

void QdPointSetP2TetraGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ============================================================================
// Tetrahedral quadrature for P3 polynomials
// ============================================================================

QdPointSetP3TetraGauss::QdPointSetP3TetraGauss() : StdPointSetBase()
{
}

QdPointSetP3TetraGauss::~QdPointSetP3TetraGauss()
{
}

Uint QdPointSetP3TetraGauss::order() const
{
  return P3;
}

Uint QdPointSetP3TetraGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP3TetraGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP3TetraGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP3TetraGauss::size(const Uint local_idx) const
{
  return 5u;
}

void QdPointSetP3TetraGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(5, 3);
  coords(0, XI0) = -0.5;
  coords(0, XI1) = -0.5;
  coords(0, XI2) = -0.5;

  coords(1, XI0) = -2. / 3.;
  coords(1, XI1) = -2. / 3.;
  coords(1, XI2) = -2. / 3.;

  coords(2, XI0) = -2. / 3.;
  coords(2, XI1) = -2. / 3.;
  coords(2, XI2) = 0.0;

  coords(3, XI0) = -2. / 3.;
  coords(3, XI1) = 0.0;
  coords(3, XI2) = -2. / 3.;

  coords(4, XI0) = 0.0;
  coords(4, XI1) = -2. / 3.;
  coords(4, XI2) = -2. / 3.;
}

void QdPointSetP3TetraGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(5);
  wgt[0] = -16. / 15.; // -1.0666666666666666
  wgt[1] = 0.6;
  wgt[2] = 0.6;
  wgt[3] = 0.6;
  wgt[4] = 0.6;
}

void QdPointSetP3TetraGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ============================================================================
// Tetrahedral quadrature for P4 polynomials
// ============================================================================

QdPointSetP4TetraGauss::QdPointSetP4TetraGauss() : StdPointSetBase()
{
}

QdPointSetP4TetraGauss::~QdPointSetP4TetraGauss()
{
}

Uint QdPointSetP4TetraGauss::order() const
{
  return P4;
}

Uint QdPointSetP4TetraGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP4TetraGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP4TetraGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP4TetraGauss::size(const Uint local_idx) const
{
  return 11u;
}

void QdPointSetP4TetraGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(11, 3);
  coords(0, XI0) = -0.5;
  coords(0, XI1) = -0.5;
  coords(0, XI2) = -0.5;

  coords(1, XI0) = -0.8571428571;
  coords(1, XI1) = -0.8571428571;
  coords(1, XI2) = -0.8571428571;

  coords(2, XI0) = -0.8571428571;
  coords(2, XI1) = -0.8571428571;
  coords(2, XI2) = 0.5714285714;

  coords(3, XI0) = -0.8571428571;
  coords(3, XI1) = 0.5714285714;
  coords(3, XI2) = -0.8571428571;

  coords(4, XI0) = 0.5714285714;
  coords(4, XI1) = -0.8571428571;
  coords(4, XI2) = -0.8571428571;

  coords(5, XI0) = -0.2011928476;
  coords(5, XI1) = -0.2011928476;
  coords(5, XI2) = -0.7988071523;

  coords(6, XI0) = -0.2011928476;
  coords(6, XI1) = -0.7988071523;
  coords(6, XI2) = -0.2011928476;

  coords(7, XI0) = -0.7988071523;
  coords(7, XI1) = -0.2011928476;
  coords(7, XI2) = -0.2011928476;

  coords(8, XI0) = -0.2011928476;
  coords(8, XI1) = -0.7988071523;
  coords(8, XI2) = -0.7988071523;

  coords(9, XI0) = -0.7988071523;
  coords(9, XI1) = -0.2011928476;
  coords(9, XI2) = -0.7988071523;

  coords(10, XI0) = -0.7988071523;
  coords(10, XI1) = -0.7988071523;
  coords(10, XI2) = -0.2011928476;
}

void QdPointSetP4TetraGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(11);
  wgt[0]  = -0.1052444444;
  wgt[1]  = 0.0609777777;
  wgt[2]  = 0.0609777777;
  wgt[3]  = 0.0609777777;
  wgt[4]  = 0.0609777777;
  wgt[5]  = 0.1991111111;
  wgt[6]  = 0.1991111111;
  wgt[7]  = 0.1991111111;
  wgt[8]  = 0.1991111111;
  wgt[9]  = 0.1991111111;
  wgt[10] = 0.1991111111;
}

void QdPointSetP4TetraGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ============================================================================
// Tetrahedral quadrature for P5 polynomials
// ============================================================================

QdPointSetP5TetraGauss::QdPointSetP5TetraGauss() : StdPointSetBase()
{
}

QdPointSetP5TetraGauss::~QdPointSetP5TetraGauss()
{
}

Uint QdPointSetP5TetraGauss::order() const
{
  return P5;
}

Uint QdPointSetP5TetraGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP5TetraGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP5TetraGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP5TetraGauss::size(const Uint local_idx) const
{
  return 14u;
}

void QdPointSetP5TetraGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(14, 3);
  coords(0, XI0) = -0.8145294993;
  coords(0, XI1) = -0.8145294993;
  coords(0, XI2) = -0.8145294993;

  coords(1, XI0) = 0.4435884981;
  coords(1, XI1) = -0.8145294993;
  coords(1, XI2) = -0.8145294993;

  coords(2, XI0) = -0.8145294993;
  coords(2, XI1) = 0.4435884981;
  coords(2, XI2) = -0.8145294993;

  coords(3, XI0) = -0.8145294993;
  coords(3, XI1) = -0.8145294993;
  coords(3, XI2) = 0.4435884981;

  coords(4, XI0) = -0.3782281614;
  coords(4, XI1) = -0.3782281614;
  coords(4, XI2) = -0.3782281614;

  coords(5, XI0) = -0.8653155155;
  coords(5, XI1) = -0.3782281614;
  coords(5, XI2) = -0.3782281614;

  coords(6, XI0) = -0.3782281614;
  coords(6, XI1) = -0.8653155155;
  coords(6, XI2) = -0.3782281614;

  coords(7, XI0) = -0.3782281614;
  coords(7, XI1) = -0.3782281614;
  coords(7, XI2) = -0.8653155155;

  coords(8, XI0) = -0.0910074082;
  coords(8, XI1) = -0.0910074082;
  coords(8, XI2) = -0.9089925917;

  coords(9, XI0) = -0.0910074082;
  coords(9, XI1) = -0.9089925917;
  coords(9, XI2) = -0.0910074082;

  coords(10, XI0) = -0.9089925917;
  coords(10, XI1) = -0.0910074082;
  coords(10, XI2) = -0.0910074082;

  coords(11, XI0) = -0.0910074082;
  coords(11, XI1) = -0.9089925917;
  coords(11, XI2) = -0.9089925917;

  coords(12, XI0) = -0.9089925917;
  coords(12, XI1) = -0.0910074082;
  coords(12, XI2) = -0.9089925917;

  coords(13, XI0) = -0.9089925917;
  coords(13, XI1) = -0.9089925917;
  coords(13, XI2) = -0.0910074082;
}

void QdPointSetP5TetraGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(14);
  wgt[0]  = 0.0979907241;
  wgt[1]  = 0.0979907241;
  wgt[2]  = 0.0979907241;
  wgt[3]  = 0.0979907241;
  wgt[4]  = 0.1502505676;
  wgt[5]  = 0.1502505676;
  wgt[6]  = 0.1502505676;
  wgt[7]  = 0.1502505676;
  wgt[8]  = 0.0567280277;
  wgt[9]  = 0.0567280277;
  wgt[10] = 0.0567280277;
  wgt[11] = 0.0567280277;
  wgt[12] = 0.0567280277;
  wgt[13] = 0.0567280277;
}

void QdPointSetP5TetraGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ============================================================================
// Tetrahedral quadrature for P6 polynomials
// ============================================================================

QdPointSetP6TetraGauss::QdPointSetP6TetraGauss() : StdPointSetBase()
{
}

QdPointSetP6TetraGauss::~QdPointSetP6TetraGauss()
{
}

Uint QdPointSetP6TetraGauss::order() const
{
  return P6;
}

Uint QdPointSetP6TetraGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP6TetraGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP6TetraGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP6TetraGauss::size(const Uint local_idx) const
{
  return 24u;
}

void QdPointSetP6TetraGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(24, 3);
  coords(0, XI0) = -0.5707942574;
  coords(0, XI1) = -0.5707942574;
  coords(0, XI2) = -0.5707942574;

  coords(1, XI0) = -0.2876172275;
  coords(1, XI1) = -0.5707942574;
  coords(1, XI2) = -0.5707942574;

  coords(2, XI0) = -0.5707942574;
  coords(2, XI1) = -0.2876172275;
  coords(2, XI2) = -0.5707942574;

  coords(3, XI0) = -0.5707942574;
  coords(3, XI1) = -0.5707942574;
  coords(3, XI2) = -0.2876172275;

  coords(4, XI0) = -0.9186520829;
  coords(4, XI1) = -0.9186520829;
  coords(4, XI2) = -0.9186520829;

  coords(5, XI0) = 0.7559562487;
  coords(5, XI1) = -0.9186520829;
  coords(5, XI2) = -0.9186520829;

  coords(6, XI0) = -0.9186520829;
  coords(6, XI1) = 0.7559562487;
  coords(6, XI2) = -0.9186520829;

  coords(7, XI0) = -0.9186520829;
  coords(7, XI1) = -0.9186520829;
  coords(7, XI2) = 0.7559562487;

  coords(8, XI0) = -0.3553242197;
  coords(8, XI1) = -0.3553242197;
  coords(8, XI2) = -0.3553242197;

  coords(9, XI0) = -0.9340273408;
  coords(9, XI1) = -0.3553242197;
  coords(9, XI2) = -0.3553242197;

  coords(10, XI0) = -0.3553242197;
  coords(10, XI1) = -0.9340273408;
  coords(10, XI2) = -0.3553242197;

  coords(11, XI0) = -0.3553242197;
  coords(11, XI1) = -0.3553242197;
  coords(11, XI2) = -0.9340273408;

  coords(12, XI0) = -0.8726779962;
  coords(12, XI1) = -0.8726779962;
  coords(12, XI2) = -0.4606553370;

  coords(13, XI0) = -0.8726779962;
  coords(13, XI1) = -0.4606553370;
  coords(13, XI2) = -0.8726779962;

  coords(14, XI0) = -0.8726779962;
  coords(14, XI1) = -0.8726779962;
  coords(14, XI2) = 0.2060113295;

  coords(15, XI0) = -0.8726779962;
  coords(15, XI1) = 0.2060113295;
  coords(15, XI2) = -0.8726779962;

  coords(16, XI0) = -0.8726779962;
  coords(16, XI1) = -0.4606553370;
  coords(16, XI2) = 0.2060113295;

  coords(17, XI0) = -0.8726779962;
  coords(17, XI1) = 0.2060113295;
  coords(17, XI2) = -0.4606553370;

  coords(18, XI0) = -0.4606553370;
  coords(18, XI1) = -0.8726779962;
  coords(18, XI2) = -0.8726779962;

  coords(19, XI0) = -0.4606553370;
  coords(19, XI1) = -0.8726779962;
  coords(19, XI2) = 0.2060113295;

  coords(20, XI0) = -0.4606553370;
  coords(20, XI1) = 0.2060113295;
  coords(20, XI2) = -0.8726779962;

  coords(21, XI0) = 0.2060113295;
  coords(21, XI1) = -0.8726779962;
  coords(21, XI2) = -0.4606553370;

  coords(22, XI0) = 0.2060113295;
  coords(22, XI1) = -0.8726779962;
  coords(22, XI2) = -0.8726779962;

  coords(23, XI0) = 0.2060113295;
  coords(23, XI1) = -0.4606553370;
  coords(23, XI2) = -0.8726779962;
}

void QdPointSetP6TetraGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(24);
  wgt[0]  = 0.0532303336;
  wgt[1]  = 0.0532303336;
  wgt[2]  = 0.0532303336;
  wgt[3]  = 0.0532303336;
  wgt[4]  = 0.0134362814;
  wgt[5]  = 0.0134362814;
  wgt[6]  = 0.0134362814;
  wgt[7]  = 0.0134362814;
  wgt[8]  = 0.0738095753;
  wgt[9]  = 0.0738095753;
  wgt[10] = 0.0738095753;
  wgt[11] = 0.0738095753;
  wgt[12] = 0.0642857142;
  wgt[13] = 0.0642857142;
  wgt[14] = 0.0642857142;
  wgt[15] = 0.0642857142;
  wgt[16] = 0.0642857142;
  wgt[17] = 0.0642857142;
  wgt[18] = 0.0642857142;
  wgt[19] = 0.0642857142;
  wgt[20] = 0.0642857142;
  wgt[21] = 0.0642857142;
  wgt[22] = 0.0642857142;
  wgt[23] = 0.0642857142;
}

void QdPointSetP6TetraGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ============================================================================
// Tetrahedral quadrature for P7 polynomials
// ============================================================================

QdPointSetP7TetraGauss::QdPointSetP7TetraGauss() : StdPointSetBase()
{
}

QdPointSetP7TetraGauss::~QdPointSetP7TetraGauss()
{
}

Uint QdPointSetP7TetraGauss::order() const
{
  return P7;
}

Uint QdPointSetP7TetraGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP7TetraGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP7TetraGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP7TetraGauss::size(const Uint local_idx) const
{
  return 31u;
}

void QdPointSetP7TetraGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(31, 3);
  coords(0, XI0) = 0.000000000000000;
  coords(0, XI1) = 0.000000000000000;
  coords(0, XI2) = -1.000000000000000;

  coords(1, XI0) = 0.000000000000000;
  coords(1, XI1) = -1.000000000000000;
  coords(1, XI2) = 0.000000000000000;

  coords(2, XI0) = -1.000000000000000;
  coords(2, XI1) = 0.000000000000000;
  coords(2, XI2) = 0.000000000000000;

  coords(3, XI0) = -1.000000000000000;
  coords(3, XI1) = -1.000000000000000;
  coords(3, XI2) = 0.000000000000000;

  coords(4, XI0) = -1.000000000000000;
  coords(4, XI1) = 0.000000000000000;
  coords(4, XI2) = -1.000000000000000;

  coords(5, XI0) = 0.000000000000000;
  coords(5, XI1) = -1.000000000000000;
  coords(5, XI2) = -1.000000000000000;

  coords(6, XI0) = -0.500000000000000;
  coords(6, XI1) = -0.500000000000000;
  coords(6, XI2) = -0.500000000000000;

  coords(7, XI0) = -0.843573615339400;
  coords(7, XI1) = -0.843573615339400;
  coords(7, XI2) = -0.843573615339400;

  coords(8, XI0) = -0.843573615339400;
  coords(8, XI1) = -0.843573615339400;
  coords(8, XI2) = 0.530720846018000;

  coords(9, XI0) = -0.843573615339400;
  coords(9, XI1) = 0.530720846018000;
  coords(9, XI2) = -0.843573615339400;

  coords(10, XI0) = 0.530720846018000;
  coords(10, XI1) = -0.843573615339400;
  coords(10, XI2) = -0.843573615339400;

  coords(11, XI0) = -0.756313566672000;
  coords(11, XI1) = -0.756313566672000;
  coords(11, XI2) = -0.756313566672000;

  coords(12, XI0) = -0.756313566672000;
  coords(12, XI1) = -0.756313566672000;
  coords(12, XI2) = 0.268940700016000;

  coords(13, XI0) = -0.756313566672000;
  coords(13, XI1) = 0.268940700016000;
  coords(13, XI2) = -0.756313566672000;

  coords(14, XI0) = 0.268940700016000;
  coords(14, XI1) = -0.756313566672000;
  coords(14, XI2) = -0.756313566672000;

  coords(15, XI0) = -0.334921671108000;
  coords(15, XI1) = -0.334921671108000;
  coords(15, XI2) = -0.334921671108000;

  coords(16, XI0) = -0.334921671108000;
  coords(16, XI1) = -0.334921671108000;
  coords(16, XI2) = -0.995234986678520;

  coords(17, XI0) = -0.334921671108000;
  coords(17, XI1) = -0.995234986678520;
  coords(17, XI2) = -0.334921671108000;

  coords(18, XI0) = -0.995234986678520;
  coords(18, XI1) = -0.334921671108000;
  coords(18, XI2) = -0.334921671108000;

  coords(19, XI0) = -0.800000000000000;
  coords(19, XI1) = -0.800000000000000;
  coords(19, XI2) = -0.600000000000000;

  coords(20, XI0) = -0.800000000000000;
  coords(20, XI1) = -0.600000000000000;
  coords(20, XI2) = -0.800000000000000;

  coords(21, XI0) = -0.800000000000000;
  coords(21, XI1) = -0.800000000000000;
  coords(21, XI2) = 0.200000000000000;

  coords(22, XI0) = -0.800000000000000;
  coords(22, XI1) = 0.200000000000000;
  coords(22, XI2) = -0.800000000000000;

  coords(23, XI0) = -0.800000000000000;
  coords(23, XI1) = -0.600000000000000;
  coords(23, XI2) = 0.200000000000000;

  coords(24, XI0) = -0.800000000000000;
  coords(24, XI1) = 0.200000000000000;
  coords(24, XI2) = -0.600000000000000;

  coords(25, XI0) = -0.600000000000000;
  coords(25, XI1) = -0.800000000000000;
  coords(25, XI2) = -0.800000000000000;

  coords(26, XI0) = -0.600000000000000;
  coords(26, XI1) = -0.800000000000000;
  coords(26, XI2) = 0.200000000000000;

  coords(27, XI0) = -0.600000000000000;
  coords(27, XI1) = 0.200000000000000;
  coords(27, XI2) = -0.800000000000000;

  coords(28, XI0) = 0.200000000000000;
  coords(28, XI1) = -0.800000000000000;
  coords(28, XI2) = -0.600000000000000;

  coords(29, XI0) = 0.200000000000000;
  coords(29, XI1) = -0.800000000000000;
  coords(29, XI2) = -0.800000000000000;

  coords(30, XI0) = 0.200000000000000;
  coords(30, XI1) = -0.600000000000000;
  coords(30, XI2) = -0.800000000000000;
}

void QdPointSetP7TetraGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(31);
  wgt[0]  = 0.007760141093480;
  wgt[1]  = 0.007760141093480;
  wgt[2]  = 0.007760141093480;
  wgt[3]  = 0.007760141093480;
  wgt[4]  = 0.007760141093480;
  wgt[5]  = 0.007760141093480;
  wgt[6]  = 0.146113787729336;
  wgt[7]  = 0.084799532195336;
  wgt[8]  = 0.084799532195336;
  wgt[9]  = 0.084799532195336;
  wgt[10] = 0.084799532195336;
  wgt[11] = -0.500141920914664;
  wgt[12] = -0.500141920914664;
  wgt[13] = -0.500141920914664;
  wgt[14] = -0.500141920914664;
  wgt[15] = 0.039131402104536;
  wgt[16] = 0.039131402104536;
  wgt[17] = 0.039131402104536;
  wgt[18] = 0.039131402104536;
  wgt[19] = 0.220458553792000;
  wgt[20] = 0.220458553792000;
  wgt[21] = 0.220458553792000;
  wgt[22] = 0.220458553792000;
  wgt[23] = 0.220458553792000;
  wgt[24] = 0.220458553792000;
  wgt[25] = 0.220458553792000;
  wgt[26] = 0.220458553792000;
  wgt[27] = 0.220458553792000;
  wgt[28] = 0.220458553792000;
  wgt[29] = 0.220458553792000;
  wgt[30] = 0.220458553792000;
}

void QdPointSetP7TetraGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ============================================================================
// Tetrahedral quadrature for P8 polynomials
// ============================================================================

QdPointSetP8TetraGauss::QdPointSetP8TetraGauss() : StdPointSetBase()
{
}

QdPointSetP8TetraGauss::~QdPointSetP8TetraGauss()
{
}

Uint QdPointSetP8TetraGauss::order() const
{
  return P8;
}

Uint QdPointSetP8TetraGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP8TetraGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP8TetraGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP8TetraGauss::size(const Uint local_idx) const
{
  return 43u;
}

void QdPointSetP8TetraGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(43, 3);
  coords(0, XI0) = -0.500000000000000;
  coords(0, XI1) = -0.500000000000000;
  coords(0, XI2) = -0.500000000000000;

  coords(1, XI0) = -0.586340136778000;
  coords(1, XI1) = -0.586340136778000;
  coords(1, XI2) = -0.586340136778000;

  coords(2, XI0) = -0.586340136778000;
  coords(2, XI1) = -0.586340136778000;
  coords(2, XI2) = -0.240979589664000;

  coords(3, XI0) = -0.586340136778000;
  coords(3, XI1) = -0.240979589664000;
  coords(3, XI2) = -0.586340136778000;

  coords(4, XI0) = -0.240979589664000;
  coords(4, XI1) = -0.586340136778000;
  coords(4, XI2) = -0.586340136778000;

  coords(5, XI0) = -0.835792823379000;
  coords(5, XI1) = -0.835792823379000;
  coords(5, XI2) = -0.835792823379000;

  coords(6, XI0) = -0.835792823379000;
  coords(6, XI1) = -0.835792823379000;
  coords(6, XI2) = 0.507378470136000;

  coords(7, XI0) = -0.835792823379000;
  coords(7, XI1) = 0.507378470136000;
  coords(7, XI2) = -0.835792823379000;

  coords(8, XI0) = 0.507378470136000;
  coords(8, XI1) = -0.835792823379000;
  coords(8, XI2) = -0.835792823379000;

  coords(9, XI0) = -0.988436098989600;
  coords(9, XI1) = -0.988436098989600;
  coords(9, XI2) = -0.988436098989600;

  coords(10, XI0) = -0.988436098989600;
  coords(10, XI1) = -0.988436098989600;
  coords(10, XI2) = 0.965308296968000;

  coords(11, XI0) = -0.988436098989600;
  coords(11, XI1) = 0.965308296968000;
  coords(11, XI2) = -0.988436098989600;

  coords(12, XI0) = 0.965308296968000;
  coords(12, XI1) = -0.988436098989600;
  coords(12, XI2) = -0.988436098989600;

  coords(13, XI0) = -0.898934519962200;
  coords(13, XI1) = -0.898934519962200;
  coords(13, XI2) = -0.101065480038000;

  coords(14, XI0) = -0.898934519962200;
  coords(14, XI1) = -0.101065480038000;
  coords(14, XI2) = -0.898934519962200;

  coords(15, XI0) = -0.101065480038000;
  coords(15, XI1) = -0.898934519962200;
  coords(15, XI2) = -0.898934519962200;

  coords(16, XI0) = -0.898934519962200;
  coords(16, XI1) = -0.101065480038000;
  coords(16, XI2) = -0.101065480038000;

  coords(17, XI0) = -0.101065480038000;
  coords(17, XI1) = -0.898934519962200;
  coords(17, XI2) = -0.101065480038000;

  coords(18, XI0) = -0.101065480038000;
  coords(18, XI1) = -0.101065480038000;
  coords(18, XI2) = -0.898934519962200;

  coords(19, XI0) = -0.541866927766000;
  coords(19, XI1) = -0.541866927766000;
  coords(19, XI2) = -0.928720834423000;

  coords(20, XI0) = -0.541866927766000;
  coords(20, XI1) = -0.928720834423000;
  coords(20, XI2) = -0.541866927766000;

  coords(21, XI0) = -0.541866927766000;
  coords(21, XI1) = -0.541866927766000;
  coords(21, XI2) = 0.012454689956000;

  coords(22, XI0) = -0.541866927766000;
  coords(22, XI1) = 0.012454689956000;
  coords(22, XI2) = -0.541866927766000;

  coords(23, XI0) = -0.541866927766000;
  coords(23, XI1) = -0.928720834423000;
  coords(23, XI2) = 0.012454689956000;

  coords(24, XI0) = -0.541866927766000;
  coords(24, XI1) = 0.012454689956000;
  coords(24, XI2) = -0.928720834423000;

  coords(25, XI0) = -0.928720834423000;
  coords(25, XI1) = -0.541866927766000;
  coords(25, XI2) = -0.541866927766000;

  coords(26, XI0) = -0.928720834423000;
  coords(26, XI1) = -0.541866927766000;
  coords(26, XI2) = 0.012454689956000;

  coords(27, XI0) = -0.928720834423000;
  coords(27, XI1) = 0.012454689956000;
  coords(27, XI2) = -0.541866927766000;

  coords(28, XI0) = 0.012454689956000;
  coords(28, XI1) = -0.541866927766000;
  coords(28, XI2) = -0.928720834423000;

  coords(29, XI0) = 0.012454689956000;
  coords(29, XI1) = -0.541866927766000;
  coords(29, XI2) = -0.541866927766000;

  coords(30, XI0) = 0.012454689956000;
  coords(30, XI1) = -0.928720834423000;
  coords(30, XI2) = -0.541866927766000;

  coords(31, XI0) = -0.926784500893600;
  coords(31, XI1) = -0.926784500893600;
  coords(31, XI2) = -0.619027916130000;

  coords(32, XI0) = -0.926784500893600;
  coords(32, XI1) = -0.619027916130000;
  coords(32, XI2) = -0.926784500893600;

  coords(33, XI0) = -0.926784500893600;
  coords(33, XI1) = -0.926784500893600;
  coords(33, XI2) = 0.472596917918000;

  coords(34, XI0) = -0.926784500893600;
  coords(34, XI1) = 0.472596917918000;
  coords(34, XI2) = -0.926784500893600;

  coords(35, XI0) = -0.926784500893600;
  coords(35, XI1) = -0.619027916130000;
  coords(35, XI2) = 0.472596917918000;

  coords(36, XI0) = -0.926784500893600;
  coords(36, XI1) = 0.472596917918000;
  coords(36, XI2) = -0.619027916130000;

  coords(37, XI0) = -0.619027916130000;
  coords(37, XI1) = -0.926784500893600;
  coords(37, XI2) = -0.926784500893600;

  coords(38, XI0) = -0.619027916130000;
  coords(38, XI1) = -0.926784500893600;
  coords(38, XI2) = 0.472596917918000;

  coords(39, XI0) = -0.619027916130000;
  coords(39, XI1) = 0.472596917918000;
  coords(39, XI2) = -0.926784500893600;

  coords(40, XI0) = 0.472596917918000;
  coords(40, XI1) = -0.926784500893600;
  coords(40, XI2) = -0.619027916130000;

  coords(41, XI0) = 0.472596917918000;
  coords(41, XI1) = -0.926784500893600;
  coords(41, XI2) = -0.926784500893600;

  coords(42, XI0) = 0.472596917918000;
  coords(42, XI1) = -0.619027916130000;
  coords(42, XI2) = -0.926784500893600;
}

void QdPointSetP8TetraGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(43);
  wgt[0]  = -0.164001509269336;
  wgt[1]  = 0.114002446582936;
  wgt[2]  = 0.114002446582936;
  wgt[3]  = 0.114002446582936;
  wgt[4]  = 0.114002446582936;
  wgt[5]  = 0.015736266505064;
  wgt[6]  = 0.015736266505064;
  wgt[7]  = 0.015736266505064;
  wgt[8]  = 0.015736266505064;
  wgt[9]  = 0.001358672872744;
  wgt[10] = 0.001358672872744;
  wgt[11] = 0.001358672872744;
  wgt[12] = 0.001358672872744;
  wgt[13] = 0.036637470595736;
  wgt[14] = 0.036637470595736;
  wgt[15] = 0.036637470595736;
  wgt[16] = 0.036637470595736;
  wgt[17] = 0.036637470595736;
  wgt[18] = 0.036637470595736;
  wgt[19] = 0.045635886469464;
  wgt[20] = 0.045635886469464;
  wgt[21] = 0.045635886469464;
  wgt[22] = 0.045635886469464;
  wgt[23] = 0.045635886469464;
  wgt[24] = 0.045635886469464;
  wgt[25] = 0.045635886469464;
  wgt[26] = 0.045635886469464;
  wgt[27] = 0.045635886469464;
  wgt[28] = 0.045635886469464;
  wgt[29] = 0.045635886469464;
  wgt[30] = 0.045635886469464;
  wgt[31] = 0.017124153129336;
  wgt[32] = 0.017124153129336;
  wgt[33] = 0.017124153129336;
  wgt[34] = 0.017124153129336;
  wgt[35] = 0.017124153129336;
  wgt[36] = 0.017124153129336;
  wgt[37] = 0.017124153129336;
  wgt[38] = 0.017124153129336;
  wgt[39] = 0.017124153129336;
  wgt[40] = 0.017124153129336;
  wgt[41] = 0.017124153129336;
  wgt[42] = 0.017124153129336;
}

void QdPointSetP8TetraGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ============================================================================
// Tetrahedral quadrature for P9 polynomials
// ============================================================================

QdPointSetP9TetraGauss::QdPointSetP9TetraGauss() : StdPointSetBase()
{
}

QdPointSetP9TetraGauss::~QdPointSetP9TetraGauss()
{
}

Uint QdPointSetP9TetraGauss::order() const
{
  return P9;
}

Uint QdPointSetP9TetraGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP9TetraGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP9TetraGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP9TetraGauss::size(const Uint local_idx) const
{
  return 53u;
}

void QdPointSetP9TetraGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  // 5 negative weights, 12 points outside of the tetrahedron
  coords.resize(53, 3);
  coords(0, XI0) = -0.500000000000000;
  coords(0, XI1) = -0.500000000000000;
  coords(0, XI2) = -0.500000000000000;

  coords(1, XI0) = -0.903297922900600;
  coords(1, XI1) = -0.903297922900600;
  coords(1, XI2) = -0.903297922900600;

  coords(2, XI0) = -0.903297922900600;
  coords(2, XI1) = -0.903297922900600;
  coords(2, XI2) = 0.709893768702000;

  coords(3, XI0) = -0.903297922900600;
  coords(3, XI1) = 0.709893768702000;
  coords(3, XI2) = -0.903297922900600;

  coords(4, XI0) = 0.709893768702000;
  coords(4, XI1) = -0.903297922900600;
  coords(4, XI2) = -0.903297922900600;

  coords(5, XI0) = -0.350841439764000;
  coords(5, XI1) = -0.350841439764000;
  coords(5, XI2) = -0.350841439764000;

  coords(6, XI0) = -0.350841439764000;
  coords(6, XI1) = -0.350841439764000;
  coords(6, XI2) = -0.947475680707200;

  coords(7, XI0) = -0.350841439764000;
  coords(7, XI1) = -0.947475680707200;
  coords(7, XI2) = -0.350841439764000;

  coords(8, XI0) = -0.947475680707200;
  coords(8, XI1) = -0.350841439764000;
  coords(8, XI2) = -0.350841439764000;

  coords(9, XI0) = -0.770766919552000;
  coords(9, XI1) = -0.770766919552000;
  coords(9, XI2) = -0.770766919552000;

  coords(10, XI0) = -0.770766919552000;
  coords(10, XI1) = -0.770766919552000;
  coords(10, XI2) = 0.312300758656000;

  coords(11, XI0) = -0.770766919552000;
  coords(11, XI1) = 0.312300758656000;
  coords(11, XI2) = -0.770766919552000;

  coords(12, XI0) = 0.312300758656000;
  coords(12, XI1) = -0.770766919552000;
  coords(12, XI2) = -0.770766919552000;

  coords(13, XI0) = -0.549020096176000;
  coords(13, XI1) = -0.549020096176000;
  coords(13, XI2) = -0.549020096176000;

  coords(14, XI0) = -0.549020096176000;
  coords(14, XI1) = -0.549020096176000;
  coords(14, XI2) = -0.352939711470000;

  coords(15, XI0) = -0.549020096176000;
  coords(15, XI1) = -0.352939711470000;
  coords(15, XI2) = -0.549020096176000;

  coords(16, XI0) = -0.352939711470000;
  coords(16, XI1) = -0.549020096176000;
  coords(16, XI2) = -0.549020096176000;

  coords(17, XI0) = -0.736744381506000;
  coords(17, XI1) = -0.736744381506000;
  coords(17, XI2) = -0.832670596765600;

  coords(18, XI0) = -0.736744381506000;
  coords(18, XI1) = -0.832670596765600;
  coords(18, XI2) = -0.736744381506000;

  coords(19, XI0) = -0.736744381506000;
  coords(19, XI1) = -0.736744381506000;
  coords(19, XI2) = 0.306159359778000;

  coords(20, XI0) = -0.736744381506000;
  coords(20, XI1) = 0.306159359778000;
  coords(20, XI2) = -0.736744381506000;

  coords(21, XI0) = -0.736744381506000;
  coords(21, XI1) = -0.832670596765600;
  coords(21, XI2) = 0.306159359778000;

  coords(22, XI0) = -0.736744381506000;
  coords(22, XI1) = 0.306159359778000;
  coords(22, XI2) = -0.832670596765600;

  coords(23, XI0) = -0.832670596765600;
  coords(23, XI1) = -0.736744381506000;
  coords(23, XI2) = -0.736744381506000;

  coords(24, XI0) = -0.832670596765600;
  coords(24, XI1) = -0.736744381506000;
  coords(24, XI2) = 0.306159359778000;

  coords(25, XI0) = -0.832670596765600;
  coords(25, XI1) = 0.306159359778000;
  coords(25, XI2) = -0.736744381506000;

  coords(26, XI0) = 0.306159359778000;
  coords(26, XI1) = -0.736744381506000;
  coords(26, XI2) = -0.832670596765600;

  coords(27, XI0) = 0.306159359778000;
  coords(27, XI1) = -0.736744381506000;
  coords(27, XI2) = -0.736744381506000;

  coords(28, XI0) = 0.306159359778000;
  coords(28, XI1) = -0.832670596765600;
  coords(28, XI2) = -0.736744381506000;

  coords(29, XI0) = -0.132097077178000;
  coords(29, XI1) = -0.132097077178000;
  coords(29, XI2) = -0.784460280902000;

  coords(30, XI0) = -0.132097077178000;
  coords(30, XI1) = -0.784460280902000;
  coords(30, XI2) = -0.132097077178000;

  coords(31, XI0) = -0.132097077178000;
  coords(31, XI1) = -0.132097077178000;
  coords(31, XI2) = -0.951345564744400;

  coords(32, XI0) = -0.132097077178000;
  coords(32, XI1) = -0.951345564744400;
  coords(32, XI2) = -0.132097077178000;

  coords(33, XI0) = -0.132097077178000;
  coords(33, XI1) = -0.784460280902000;
  coords(33, XI2) = -0.951345564744400;

  coords(34, XI0) = -0.132097077178000;
  coords(34, XI1) = -0.951345564744400;
  coords(34, XI2) = -0.784460280902000;

  coords(35, XI0) = -0.784460280902000;
  coords(35, XI1) = -0.132097077178000;
  coords(35, XI2) = -0.132097077178000;

  coords(36, XI0) = -0.784460280902000;
  coords(36, XI1) = -0.132097077178000;
  coords(36, XI2) = -0.951345564744400;

  coords(37, XI0) = -0.784460280902000;
  coords(37, XI1) = -0.951345564744400;
  coords(37, XI2) = -0.132097077178000;

  coords(38, XI0) = -0.951345564744400;
  coords(38, XI1) = -0.132097077178000;
  coords(38, XI2) = -0.784460280902000;

  coords(39, XI0) = -0.951345564744400;
  coords(39, XI1) = -0.132097077178000;
  coords(39, XI2) = -0.132097077178000;

  coords(40, XI0) = -0.951345564744400;
  coords(40, XI1) = -0.784460280902000;
  coords(40, XI2) = -0.132097077178000;

  coords(41, XI0) = -1.002752554636280;
  coords(41, XI1) = -1.002752554636280;
  coords(41, XI2) = -0.446893054726000;

  coords(42, XI0) = -1.002752554636280;
  coords(42, XI1) = -0.446893054726000;
  coords(42, XI2) = -1.002752554636280;

  coords(43, XI0) = -1.002752554636280;
  coords(43, XI1) = -1.002752554636280;
  coords(43, XI2) = 0.452398163998000;

  coords(44, XI0) = -1.002752554636280;
  coords(44, XI1) = 0.452398163998000;
  coords(44, XI2) = -1.002752554636280;

  coords(45, XI0) = -1.002752554636280;
  coords(45, XI1) = -0.446893054726000;
  coords(45, XI2) = 0.452398163998000;

  coords(46, XI0) = -1.002752554636280;
  coords(46, XI1) = 0.452398163998000;
  coords(46, XI2) = -0.446893054726000;

  coords(47, XI0) = -0.446893054726000;
  coords(47, XI1) = -1.002752554636280;
  coords(47, XI2) = -1.002752554636280;

  coords(48, XI0) = -0.446893054726000;
  coords(48, XI1) = -1.002752554636280;
  coords(48, XI2) = 0.452398163998000;

  coords(49, XI0) = -0.446893054726000;
  coords(49, XI1) = 0.452398163998000;
  coords(49, XI2) = -1.002752554636280;

  coords(50, XI0) = 0.452398163998000;
  coords(50, XI1) = -1.002752554636280;
  coords(50, XI2) = -0.446893054726000;

  coords(51, XI0) = 0.452398163998000;
  coords(51, XI1) = -1.002752554636280;
  coords(51, XI2) = -1.002752554636280;

  coords(52, XI0) = 0.452398163998000;
  coords(52, XI1) = -0.446893054726000;
  coords(52, XI2) = -1.002752554636280;
}

void QdPointSetP9TetraGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(53);
  wgt[0]  = -1.102392306609336;
  wgt[1]  = 0.014922692552664;
  wgt[2]  = 0.014922692552664;
  wgt[3]  = 0.014922692552664;
  wgt[4]  = 0.014922692552664;
  wgt[5]  = 0.034475391756000;
  wgt[6]  = 0.034475391756000;
  wgt[7]  = 0.034475391756000;
  wgt[8]  = 0.034475391756000;
  wgt[9]  = -0.721478131849336;
  wgt[10] = -0.721478131849336;
  wgt[11] = -0.721478131849336;
  wgt[12] = -0.721478131849336;
  wgt[13] = 0.357380609620000;
  wgt[14] = 0.357380609620000;
  wgt[15] = 0.357380609620000;
  wgt[16] = 0.357380609620000;
  wgt[17] = 0.277603247076000;
  wgt[18] = 0.277603247076000;
  wgt[19] = 0.277603247076000;
  wgt[20] = 0.277603247076000;
  wgt[21] = 0.277603247076000;
  wgt[22] = 0.277603247076000;
  wgt[23] = 0.277603247076000;
  wgt[24] = 0.277603247076000;
  wgt[25] = 0.277603247076000;
  wgt[26] = 0.277603247076000;
  wgt[27] = 0.277603247076000;
  wgt[28] = 0.277603247076000;
  wgt[29] = 0.026820671221336;
  wgt[30] = 0.026820671221336;
  wgt[31] = 0.026820671221336;
  wgt[32] = 0.026820671221336;
  wgt[33] = 0.026820671221336;
  wgt[34] = 0.026820671221336;
  wgt[35] = 0.026820671221336;
  wgt[36] = 0.026820671221336;
  wgt[37] = 0.026820671221336;
  wgt[38] = 0.026820671221336;
  wgt[39] = 0.026820671221336;
  wgt[40] = 0.026820671221336;
  wgt[41] = 0.003453031004456;
  wgt[42] = 0.003453031004456;
  wgt[43] = 0.003453031004456;
  wgt[44] = 0.003453031004456;
  wgt[45] = 0.003453031004456;
  wgt[46] = 0.003453031004456;
  wgt[47] = 0.003453031004456;
  wgt[48] = 0.003453031004456;
  wgt[49] = 0.003453031004456;
  wgt[50] = 0.003453031004456;
  wgt[51] = 0.003453031004456;
  wgt[52] = 0.003453031004456;
}

void QdPointSetP9TetraGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ============================================================================
// Tetrahedral quadrature for P10 polynomials
// ============================================================================

QdPointSetP10TetraGauss::QdPointSetP10TetraGauss() : StdPointSetBase()
{
}

QdPointSetP10TetraGauss::~QdPointSetP10TetraGauss()
{
}

Uint QdPointSetP10TetraGauss::order() const
{
  return P10;
}

Uint QdPointSetP10TetraGauss::dim() const
{
  return _3D;
}

Uint QdPointSetP10TetraGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP10TetraGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP10TetraGauss::size(const Uint local_idx) const
{
  return 126u;
}

void QdPointSetP10TetraGauss::reference_coords(math::DenseDMat<Real> &coords,
                                               const Uint local_idx) const
{
  // This is actually P11 tet. quadrature from Solin et. al used for P10
  // polynomials
  coords.resize(126, 3);
  coords(0, XI0) = -0.857142857142800;
  coords(0, XI1) = -0.857142857142800;
  coords(0, XI2) = 0.571428571428000;

  coords(1, XI0) = -0.857142857142800;
  coords(1, XI1) = -0.571428571428000;
  coords(1, XI2) = 0.285714285714000;

  coords(2, XI0) = -0.857142857142800;
  coords(2, XI1) = -0.285714285714000;
  coords(2, XI2) = 0.000000000000000;

  coords(3, XI0) = -0.857142857142800;
  coords(3, XI1) = 0.000000000000000;
  coords(3, XI2) = -0.285714285714000;

  coords(4, XI0) = -0.857142857142800;
  coords(4, XI1) = 0.285714285714000;
  coords(4, XI2) = -0.571428571428000;

  coords(5, XI0) = -0.857142857142800;
  coords(5, XI1) = 0.571428571428000;
  coords(5, XI2) = -0.857142857142800;

  coords(6, XI0) = -0.571428571428000;
  coords(6, XI1) = -0.857142857142800;
  coords(6, XI2) = 0.285714285714000;

  coords(7, XI0) = -0.571428571428000;
  coords(7, XI1) = -0.571428571428000;
  coords(7, XI2) = 0.000000000000000;

  coords(8, XI0) = -0.571428571428000;
  coords(8, XI1) = -0.285714285714000;
  coords(8, XI2) = -0.285714285714000;

  coords(9, XI0) = -0.571428571428000;
  coords(9, XI1) = 0.000000000000000;
  coords(9, XI2) = -0.571428571428000;

  coords(10, XI0) = -0.571428571428000;
  coords(10, XI1) = 0.285714285714000;
  coords(10, XI2) = -0.857142857142800;

  coords(11, XI0) = -0.285714285714000;
  coords(11, XI1) = -0.857142857142800;
  coords(11, XI2) = 0.000000000000000;

  coords(12, XI0) = -0.285714285714000;
  coords(12, XI1) = -0.571428571428000;
  coords(12, XI2) = -0.285714285714000;

  coords(13, XI0) = -0.285714285714000;
  coords(13, XI1) = -0.285714285714000;
  coords(13, XI2) = -0.571428571428000;

  coords(14, XI0) = -0.285714285714000;
  coords(14, XI1) = 0.000000000000000;
  coords(14, XI2) = -0.857142857142800;

  coords(15, XI0) = 0.000000000000000;
  coords(15, XI1) = -0.857142857142800;
  coords(15, XI2) = -0.285714285714000;

  coords(16, XI0) = 0.000000000000000;
  coords(16, XI1) = -0.571428571428000;
  coords(16, XI2) = -0.571428571428000;

  coords(17, XI0) = 0.000000000000000;
  coords(17, XI1) = -0.285714285714000;
  coords(17, XI2) = -0.857142857142800;

  coords(18, XI0) = 0.285714285714000;
  coords(18, XI1) = -0.857142857142800;
  coords(18, XI2) = -0.571428571428000;

  coords(19, XI0) = 0.285714285714000;
  coords(19, XI1) = -0.571428571428000;
  coords(19, XI2) = -0.857142857142800;

  coords(20, XI0) = 0.571428571428000;
  coords(20, XI1) = -0.857142857142800;
  coords(20, XI2) = -0.857142857142800;

  coords(21, XI0) = -0.857142857142800;
  coords(21, XI1) = -0.857142857142800;
  coords(21, XI2) = 0.285714285714000;

  coords(22, XI0) = -0.857142857142800;
  coords(22, XI1) = -0.571428571428000;
  coords(22, XI2) = 0.000000000000000;

  coords(23, XI0) = -0.857142857142800;
  coords(23, XI1) = -0.285714285714000;
  coords(23, XI2) = -0.285714285714000;

  coords(24, XI0) = -0.857142857142800;
  coords(24, XI1) = 0.000000000000000;
  coords(24, XI2) = -0.571428571428000;

  coords(25, XI0) = -0.857142857142800;
  coords(25, XI1) = 0.285714285714000;
  coords(25, XI2) = -0.857142857142800;

  coords(26, XI0) = -0.571428571428000;
  coords(26, XI1) = -0.857142857142800;
  coords(26, XI2) = 0.000000000000000;

  coords(27, XI0) = -0.571428571428000;
  coords(27, XI1) = -0.571428571428000;
  coords(27, XI2) = -0.285714285714000;

  coords(28, XI0) = -0.571428571428000;
  coords(28, XI1) = -0.285714285714000;
  coords(28, XI2) = -0.571428571428000;

  coords(29, XI0) = -0.571428571428000;
  coords(29, XI1) = 0.000000000000000;
  coords(29, XI2) = -0.857142857142800;

  coords(30, XI0) = -0.285714285714000;
  coords(30, XI1) = -0.857142857142800;
  coords(30, XI2) = -0.285714285714000;

  coords(31, XI0) = -0.285714285714000;
  coords(31, XI1) = -0.571428571428000;
  coords(31, XI2) = -0.571428571428000;

  coords(32, XI0) = -0.285714285714000;
  coords(32, XI1) = -0.285714285714000;
  coords(32, XI2) = -0.857142857142800;

  coords(33, XI0) = 0.000000000000000;
  coords(33, XI1) = -0.857142857142800;
  coords(33, XI2) = -0.571428571428000;

  coords(34, XI0) = 0.000000000000000;
  coords(34, XI1) = -0.571428571428000;
  coords(34, XI2) = -0.857142857142800;

  coords(35, XI0) = 0.285714285714000;
  coords(35, XI1) = -0.857142857142800;
  coords(35, XI2) = -0.857142857142800;

  coords(36, XI0) = -0.857142857142800;
  coords(36, XI1) = -0.857142857142800;
  coords(36, XI2) = 0.000000000000000;

  coords(37, XI0) = -0.857142857142800;
  coords(37, XI1) = -0.571428571428000;
  coords(37, XI2) = -0.285714285714000;

  coords(38, XI0) = -0.857142857142800;
  coords(38, XI1) = -0.285714285714000;
  coords(38, XI2) = -0.571428571428000;

  coords(39, XI0) = -0.857142857142800;
  coords(39, XI1) = 0.000000000000000;
  coords(39, XI2) = -0.857142857142800;

  coords(40, XI0) = -0.571428571428000;
  coords(40, XI1) = -0.857142857142800;
  coords(40, XI2) = -0.285714285714000;

  coords(41, XI0) = -0.571428571428000;
  coords(41, XI1) = -0.571428571428000;
  coords(41, XI2) = -0.571428571428000;

  coords(42, XI0) = -0.571428571428000;
  coords(42, XI1) = -0.285714285714000;
  coords(42, XI2) = -0.857142857142800;

  coords(43, XI0) = -0.285714285714000;
  coords(43, XI1) = -0.857142857142800;
  coords(43, XI2) = -0.571428571428000;

  coords(44, XI0) = -0.285714285714000;
  coords(44, XI1) = -0.571428571428000;
  coords(44, XI2) = -0.857142857142800;

  coords(45, XI0) = 0.000000000000000;
  coords(45, XI1) = -0.857142857142800;
  coords(45, XI2) = -0.857142857142800;

  coords(46, XI0) = -0.857142857142800;
  coords(46, XI1) = -0.857142857142800;
  coords(46, XI2) = -0.285714285714000;

  coords(47, XI0) = -0.857142857142800;
  coords(47, XI1) = -0.571428571428000;
  coords(47, XI2) = -0.571428571428000;

  coords(48, XI0) = -0.857142857142800;
  coords(48, XI1) = -0.285714285714000;
  coords(48, XI2) = -0.857142857142800;

  coords(49, XI0) = -0.571428571428000;
  coords(49, XI1) = -0.857142857142800;
  coords(49, XI2) = -0.571428571428000;

  coords(50, XI0) = -0.571428571428000;
  coords(50, XI1) = -0.571428571428000;
  coords(50, XI2) = -0.857142857142800;

  coords(51, XI0) = -0.285714285714000;
  coords(51, XI1) = -0.857142857142800;
  coords(51, XI2) = -0.857142857142800;

  coords(52, XI0) = -0.857142857142800;
  coords(52, XI1) = -0.857142857142800;
  coords(52, XI2) = -0.571428571428000;

  coords(53, XI0) = -0.857142857142800;
  coords(53, XI1) = -0.571428571428000;
  coords(53, XI2) = -0.857142857142800;

  coords(54, XI0) = -0.571428571428000;
  coords(54, XI1) = -0.857142857142800;
  coords(54, XI2) = -0.857142857142800;

  coords(55, XI0) = -0.857142857142800;
  coords(55, XI1) = -0.857142857142800;
  coords(55, XI2) = -0.857142857142800;

  coords(56, XI0) = -0.833333333333400;
  coords(56, XI1) = -0.833333333333400;
  coords(56, XI2) = 0.500000000000000;

  coords(57, XI0) = -0.833333333333400;
  coords(57, XI1) = -0.500000000000000;
  coords(57, XI2) = 0.166666666666000;

  coords(58, XI0) = -0.833333333333400;
  coords(58, XI1) = -0.166666666666000;
  coords(58, XI2) = -0.166666666666000;

  coords(59, XI0) = -0.833333333333400;
  coords(59, XI1) = 0.166666666666000;
  coords(59, XI2) = -0.500000000000000;

  coords(60, XI0) = -0.833333333333400;
  coords(60, XI1) = 0.500000000000000;
  coords(60, XI2) = -0.833333333333400;

  coords(61, XI0) = -0.500000000000000;
  coords(61, XI1) = -0.833333333333400;
  coords(61, XI2) = 0.166666666666000;

  coords(62, XI0) = -0.500000000000000;
  coords(62, XI1) = -0.500000000000000;
  coords(62, XI2) = -0.166666666666000;

  coords(63, XI0) = -0.500000000000000;
  coords(63, XI1) = -0.166666666666000;
  coords(63, XI2) = -0.500000000000000;

  coords(64, XI0) = -0.500000000000000;
  coords(64, XI1) = 0.166666666666000;
  coords(64, XI2) = -0.833333333333400;

  coords(65, XI0) = -0.166666666666000;
  coords(65, XI1) = -0.833333333333400;
  coords(65, XI2) = -0.166666666666000;

  coords(66, XI0) = -0.166666666666000;
  coords(66, XI1) = -0.500000000000000;
  coords(66, XI2) = -0.500000000000000;

  coords(67, XI0) = -0.166666666666000;
  coords(67, XI1) = -0.166666666666000;
  coords(67, XI2) = -0.833333333333400;

  coords(68, XI0) = 0.166666666666000;
  coords(68, XI1) = -0.833333333333400;
  coords(68, XI2) = -0.500000000000000;

  coords(69, XI0) = 0.166666666666000;
  coords(69, XI1) = -0.500000000000000;
  coords(69, XI2) = -0.833333333333400;

  coords(70, XI0) = 0.500000000000000;
  coords(70, XI1) = -0.833333333333400;
  coords(70, XI2) = -0.833333333333400;

  coords(71, XI0) = -0.833333333333400;
  coords(71, XI1) = -0.833333333333400;
  coords(71, XI2) = 0.166666666666000;

  coords(72, XI0) = -0.833333333333400;
  coords(72, XI1) = -0.500000000000000;
  coords(72, XI2) = -0.166666666666000;

  coords(73, XI0) = -0.833333333333400;
  coords(73, XI1) = -0.166666666666000;
  coords(73, XI2) = -0.500000000000000;

  coords(74, XI0) = -0.833333333333400;
  coords(74, XI1) = 0.166666666666000;
  coords(74, XI2) = -0.833333333333400;

  coords(75, XI0) = -0.500000000000000;
  coords(75, XI1) = -0.833333333333400;
  coords(75, XI2) = -0.166666666666000;

  coords(76, XI0) = -0.500000000000000;
  coords(76, XI1) = -0.500000000000000;
  coords(76, XI2) = -0.500000000000000;

  coords(77, XI0) = -0.500000000000000;
  coords(77, XI1) = -0.166666666666000;
  coords(77, XI2) = -0.833333333333400;

  coords(78, XI0) = -0.166666666666000;
  coords(78, XI1) = -0.833333333333400;
  coords(78, XI2) = -0.500000000000000;

  coords(79, XI0) = -0.166666666666000;
  coords(79, XI1) = -0.500000000000000;
  coords(79, XI2) = -0.833333333333400;

  coords(80, XI0) = 0.166666666666000;
  coords(80, XI1) = -0.833333333333400;
  coords(80, XI2) = -0.833333333333400;

  coords(81, XI0) = -0.833333333333400;
  coords(81, XI1) = -0.833333333333400;
  coords(81, XI2) = -0.166666666666000;

  coords(82, XI0) = -0.833333333333400;
  coords(82, XI1) = -0.500000000000000;
  coords(82, XI2) = -0.500000000000000;

  coords(83, XI0) = -0.833333333333400;
  coords(83, XI1) = -0.166666666666000;
  coords(83, XI2) = -0.833333333333400;

  coords(84, XI0) = -0.500000000000000;
  coords(84, XI1) = -0.833333333333400;
  coords(84, XI2) = -0.500000000000000;

  coords(85, XI0) = -0.500000000000000;
  coords(85, XI1) = -0.500000000000000;
  coords(85, XI2) = -0.833333333333400;

  coords(86, XI0) = -0.166666666666000;
  coords(86, XI1) = -0.833333333333400;
  coords(86, XI2) = -0.833333333333400;

  coords(87, XI0) = -0.833333333333400;
  coords(87, XI1) = -0.833333333333400;
  coords(87, XI2) = -0.500000000000000;

  coords(88, XI0) = -0.833333333333400;
  coords(88, XI1) = -0.500000000000000;
  coords(88, XI2) = -0.833333333333400;

  coords(89, XI0) = -0.500000000000000;
  coords(89, XI1) = -0.833333333333400;
  coords(89, XI2) = -0.833333333333400;

  coords(90, XI0) = -0.833333333333400;
  coords(90, XI1) = -0.833333333333400;
  coords(90, XI2) = -0.833333333333400;

  coords(91, XI0) = -0.800000000000000;
  coords(91, XI1) = -0.800000000000000;
  coords(91, XI2) = 0.400000000000000;

  coords(92, XI0) = -0.800000000000000;
  coords(92, XI1) = -0.400000000000000;
  coords(92, XI2) = 0.000000000000000;

  coords(93, XI0) = -0.800000000000000;
  coords(93, XI1) = 0.000000000000000;
  coords(93, XI2) = -0.400000000000000;

  coords(94, XI0) = -0.800000000000000;
  coords(94, XI1) = 0.400000000000000;
  coords(94, XI2) = -0.800000000000000;

  coords(95, XI0) = -0.400000000000000;
  coords(95, XI1) = -0.800000000000000;
  coords(95, XI2) = 0.000000000000000;

  coords(96, XI0) = -0.400000000000000;
  coords(96, XI1) = -0.400000000000000;
  coords(96, XI2) = -0.400000000000000;

  coords(97, XI0) = -0.400000000000000;
  coords(97, XI1) = 0.000000000000000;
  coords(97, XI2) = -0.800000000000000;

  coords(98, XI0) = 0.000000000000000;
  coords(98, XI1) = -0.800000000000000;
  coords(98, XI2) = -0.400000000000000;

  coords(99, XI0) = 0.000000000000000;
  coords(99, XI1) = -0.400000000000000;
  coords(99, XI2) = -0.800000000000000;

  coords(100, XI0) = 0.400000000000000;
  coords(100, XI1) = -0.800000000000000;
  coords(100, XI2) = -0.800000000000000;

  coords(101, XI0) = -0.800000000000000;
  coords(101, XI1) = -0.800000000000000;
  coords(101, XI2) = 0.000000000000000;

  coords(102, XI0) = -0.800000000000000;
  coords(102, XI1) = -0.400000000000000;
  coords(102, XI2) = -0.400000000000000;

  coords(103, XI0) = -0.800000000000000;
  coords(103, XI1) = 0.000000000000000;
  coords(103, XI2) = -0.800000000000000;

  coords(104, XI0) = -0.400000000000000;
  coords(104, XI1) = -0.800000000000000;
  coords(104, XI2) = -0.400000000000000;

  coords(105, XI0) = -0.400000000000000;
  coords(105, XI1) = -0.400000000000000;
  coords(105, XI2) = -0.800000000000000;

  coords(106, XI0) = 0.000000000000000;
  coords(106, XI1) = -0.800000000000000;
  coords(106, XI2) = -0.800000000000000;

  coords(107, XI0) = -0.800000000000000;
  coords(107, XI1) = -0.800000000000000;
  coords(107, XI2) = -0.400000000000000;

  coords(108, XI0) = -0.800000000000000;
  coords(108, XI1) = -0.400000000000000;
  coords(108, XI2) = -0.800000000000000;

  coords(109, XI0) = -0.400000000000000;
  coords(109, XI1) = -0.800000000000000;
  coords(109, XI2) = -0.800000000000000;

  coords(110, XI0) = -0.800000000000000;
  coords(110, XI1) = -0.800000000000000;
  coords(110, XI2) = -0.800000000000000;

  coords(111, XI0) = -0.750000000000000;
  coords(111, XI1) = -0.750000000000000;
  coords(111, XI2) = 0.250000000000000;

  coords(112, XI0) = -0.750000000000000;
  coords(112, XI1) = -0.250000000000000;
  coords(112, XI2) = -0.250000000000000;

  coords(113, XI0) = -0.750000000000000;
  coords(113, XI1) = 0.250000000000000;
  coords(113, XI2) = -0.750000000000000;

  coords(114, XI0) = -0.250000000000000;
  coords(114, XI1) = -0.750000000000000;
  coords(114, XI2) = -0.250000000000000;

  coords(115, XI0) = -0.250000000000000;
  coords(115, XI1) = -0.250000000000000;
  coords(115, XI2) = -0.750000000000000;

  coords(116, XI0) = 0.250000000000000;
  coords(116, XI1) = -0.750000000000000;
  coords(116, XI2) = -0.750000000000000;

  coords(117, XI0) = -0.750000000000000;
  coords(117, XI1) = -0.750000000000000;
  coords(117, XI2) = -0.250000000000000;

  coords(118, XI0) = -0.750000000000000;
  coords(118, XI1) = -0.250000000000000;
  coords(118, XI2) = -0.750000000000000;

  coords(119, XI0) = -0.250000000000000;
  coords(119, XI1) = -0.750000000000000;
  coords(119, XI2) = -0.750000000000000;

  coords(120, XI0) = -0.750000000000000;
  coords(120, XI1) = -0.750000000000000;
  coords(120, XI2) = -0.750000000000000;

  coords(121, XI0) = -0.666666666666000;
  coords(121, XI1) = -0.666666666666000;
  coords(121, XI2) = 0.000000000000000;

  coords(122, XI0) = -0.666666666666000;
  coords(122, XI1) = 0.000000000000000;
  coords(122, XI2) = -0.666666666666000;

  coords(123, XI0) = 0.000000000000000;
  coords(123, XI1) = -0.666666666666000;
  coords(123, XI2) = -0.666666666666000;

  coords(124, XI0) = -0.666666666666000;
  coords(124, XI1) = -0.666666666666000;
  coords(124, XI2) = -0.666666666666000;

  coords(125, XI0) = -0.500000000000000;
  coords(125, XI1) = -0.500000000000000;
  coords(125, XI2) = -0.500000000000000;
}

void QdPointSetP10TetraGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(126);
  wgt[0]   = 0.362902592520000;
  wgt[1]   = 0.362902592520000;
  wgt[2]   = 0.362902592520000;
  wgt[3]   = 0.362902592520000;
  wgt[4]   = 0.362902592520000;
  wgt[5]   = 0.362902592520000;
  wgt[6]   = 0.362902592520000;
  wgt[7]   = 0.362902592520000;
  wgt[8]   = 0.362902592520000;
  wgt[9]   = 0.362902592520000;
  wgt[10]  = 0.362902592520000;
  wgt[11]  = 0.362902592520000;
  wgt[12]  = 0.362902592520000;
  wgt[13]  = 0.362902592520000;
  wgt[14]  = 0.362902592520000;
  wgt[15]  = 0.362902592520000;
  wgt[16]  = 0.362902592520000;
  wgt[17]  = 0.362902592520000;
  wgt[18]  = 0.362902592520000;
  wgt[19]  = 0.362902592520000;
  wgt[20]  = 0.362902592520000;
  wgt[21]  = 0.362902592520000;
  wgt[22]  = 0.362902592520000;
  wgt[23]  = 0.362902592520000;
  wgt[24]  = 0.362902592520000;
  wgt[25]  = 0.362902592520000;
  wgt[26]  = 0.362902592520000;
  wgt[27]  = 0.362902592520000;
  wgt[28]  = 0.362902592520000;
  wgt[29]  = 0.362902592520000;
  wgt[30]  = 0.362902592520000;
  wgt[31]  = 0.362902592520000;
  wgt[32]  = 0.362902592520000;
  wgt[33]  = 0.362902592520000;
  wgt[34]  = 0.362902592520000;
  wgt[35]  = 0.362902592520000;
  wgt[36]  = 0.362902592520000;
  wgt[37]  = 0.362902592520000;
  wgt[38]  = 0.362902592520000;
  wgt[39]  = 0.362902592520000;
  wgt[40]  = 0.362902592520000;
  wgt[41]  = 0.362902592520000;
  wgt[42]  = 0.362902592520000;
  wgt[43]  = 0.362902592520000;
  wgt[44]  = 0.362902592520000;
  wgt[45]  = 0.362902592520000;
  wgt[46]  = 0.362902592520000;
  wgt[47]  = 0.362902592520000;
  wgt[48]  = 0.362902592520000;
  wgt[49]  = 0.362902592520000;
  wgt[50]  = 0.362902592520000;
  wgt[51]  = 0.362902592520000;
  wgt[52]  = 0.362902592520000;
  wgt[53]  = 0.362902592520000;
  wgt[54]  = 0.362902592520000;
  wgt[55]  = 0.362902592520000;
  wgt[56]  = -0.932187812188000;
  wgt[57]  = -0.932187812188000;
  wgt[58]  = -0.932187812188000;
  wgt[59]  = -0.932187812188000;
  wgt[60]  = -0.932187812188000;
  wgt[61]  = -0.932187812188000;
  wgt[62]  = -0.932187812188000;
  wgt[63]  = -0.932187812188000;
  wgt[64]  = -0.932187812188000;
  wgt[65]  = -0.932187812188000;
  wgt[66]  = -0.932187812188000;
  wgt[67]  = -0.932187812188000;
  wgt[68]  = -0.932187812188000;
  wgt[69]  = -0.932187812188000;
  wgt[70]  = -0.932187812188000;
  wgt[71]  = -0.932187812188000;
  wgt[72]  = -0.932187812188000;
  wgt[73]  = -0.932187812188000;
  wgt[74]  = -0.932187812188000;
  wgt[75]  = -0.932187812188000;
  wgt[76]  = -0.932187812188000;
  wgt[77]  = -0.932187812188000;
  wgt[78]  = -0.932187812188000;
  wgt[79]  = -0.932187812188000;
  wgt[80]  = -0.932187812188000;
  wgt[81]  = -0.932187812188000;
  wgt[82]  = -0.932187812188000;
  wgt[83]  = -0.932187812188000;
  wgt[84]  = -0.932187812188000;
  wgt[85]  = -0.932187812188000;
  wgt[86]  = -0.932187812188000;
  wgt[87]  = -0.932187812188000;
  wgt[88]  = -0.932187812188000;
  wgt[89]  = -0.932187812188000;
  wgt[90]  = -0.932187812188000;
  wgt[91]  = 0.815498319838664;
  wgt[92]  = 0.815498319838664;
  wgt[93]  = 0.815498319838664;
  wgt[94]  = 0.815498319838664;
  wgt[95]  = 0.815498319838664;
  wgt[96]  = 0.815498319838664;
  wgt[97]  = 0.815498319838664;
  wgt[98]  = 0.815498319838664;
  wgt[99]  = 0.815498319838664;
  wgt[100] = 0.815498319838664;
  wgt[101] = 0.815498319838664;
  wgt[102] = 0.815498319838664;
  wgt[103] = 0.815498319838664;
  wgt[104] = 0.815498319838664;
  wgt[105] = 0.815498319838664;
  wgt[106] = 0.815498319838664;
  wgt[107] = 0.815498319838664;
  wgt[108] = 0.815498319838664;
  wgt[109] = 0.815498319838664;
  wgt[110] = 0.815498319838664;
  wgt[111] = -0.280203089092000;
  wgt[112] = -0.280203089092000;
  wgt[113] = -0.280203089092000;
  wgt[114] = -0.280203089092000;
  wgt[115] = -0.280203089092000;
  wgt[116] = -0.280203089092000;
  wgt[117] = -0.280203089092000;
  wgt[118] = -0.280203089092000;
  wgt[119] = -0.280203089092000;
  wgt[120] = -0.280203089092000;
  wgt[121] = 0.032544642857200;
  wgt[122] = 0.032544642857200;
  wgt[123] = 0.032544642857200;
  wgt[124] = 0.032544642857200;
  wgt[125] = -0.000752498530272;
}

void QdPointSetP10TetraGauss::permutation(const Uint local_id,
                                          const mesh::EntityRealignCode &permutation_code,
                                          std::vector<Uint> &permutation_vec)
{
  // Do not do anything, since this is a 'volume' quadrature
  // Quadrature permutations should only be allowed for quadrature rules
  // that place nodes on element faces
  permutation_vec.resize(0);
}

// ============================================================================

} // namespace mesh

} // namespace pdekit
