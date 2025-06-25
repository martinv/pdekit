#include "mesh/point_set/QdPointSetTriagGauss.hpp"
#include "common/Constants.hpp"
#include "mesh/point_set/QuadratureTransformUtils.hpp"

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
                                       std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  // First set permutation_vec as identity permutation: p[i] = i;

  permutation_vec.resize(nb_qd_pts);
  for (Uint i = 0; i < nb_qd_pts; ++i)
  {
    permutation_vec[i] = i;
  }

  // If this permutation is identity, we are done
  if (permutation_code.is_identity(ElemShape::Triag))
  {
    return;
  }

  // Temporary vector
  std::vector<Uint> old_permutation(nb_qd_pts);

  // Apply all flips first
  for (Uint i = 0; i < permutation_code.nb_flips(); ++i)
  {
    old_permutation.swap(permutation_vec);

    // Compute the flip permutation
    for (Uint j = 0; j < nb_qd_pts; ++j)
    {
      permutation_vec[j] = old_permutation[canonical_flip_permutation[j]];
    }
  }

  // Then apply all rotations
  for (Uint i = 0; i < permutation_code.nb_rotations(); ++i)
  {
    old_permutation.swap(permutation_vec);

    // Compute rotation permutation:
    for (Uint j = 0; j < nb_qd_pts; ++j)
    {
      permutation_vec[j] = old_permutation[canonical_rot_permutation[j]];
    }
  }

} // function 'fill_triag_quadrature_permutation'

} // namespace detail

// ----------------------------------------------------------------------------
// Triangle quadrature for P1 polynomials
// ----------------------------------------------------------------------------

QdPointSetP1TriagGauss::QdPointSetP1TriagGauss() : StdPointSetBase()
{
}

QdPointSetP1TriagGauss::~QdPointSetP1TriagGauss()
{
}

Uint QdPointSetP1TriagGauss::order() const
{
  return P1;
}

Uint QdPointSetP1TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP1TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP1TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP1TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP1TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  coords.resize(1, 2);
  coords(0, XI0) = -1. / 3.;
  coords(0, XI1) = -1. / 3.;
}

void QdPointSetP1TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(1);
  wgt[0] = 2.0;
}

void QdPointSetP1TriagGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // Whatever is the rotation/flip, there is only one quadrature point -
  // nothing to permute
  permutation_vec.resize(1);
  permutation_vec[0] = 0;
}

// ----------------------------------------------------------------------------
// Triangle quadrature for P2 polynomials
// ----------------------------------------------------------------------------

QdPointSetP2TriagGauss::QdPointSetP2TriagGauss() : StdPointSetBase()
{
}

QdPointSetP2TriagGauss::~QdPointSetP2TriagGauss()
{
}

Uint QdPointSetP2TriagGauss::order() const
{
  return P2;
}

Uint QdPointSetP2TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP2TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP2TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP2TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP2TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP2TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 2. / 3.;
  wgt[1] = 2. / 3.;
  wgt[2] = 2. / 3.;
}

void QdPointSetP2TriagGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP2TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                      const Uint local_idx)
{
  coords.resize(3, 2);
  coords(0, XI0) = -2. / 3.;
  coords(0, XI1) = -2. / 3.;

  coords(1, XI0) = -2. / 3.;
  coords(1, XI1) = 1. / 3.;

  coords(2, XI0) = 1. / 3.;
  coords(2, XI1) = -2. / 3.;
}

void QdPointSetP2TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP2TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP2TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP2TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P3 polynomials
// ----------------------------------------------------------------------------

QdPointSetP3TriagGauss::QdPointSetP3TriagGauss() : StdPointSetBase()
{
}

QdPointSetP3TriagGauss::~QdPointSetP3TriagGauss()
{
}

Uint QdPointSetP3TriagGauss::order() const
{
  return P3;
}

Uint QdPointSetP3TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP3TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP3TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP3TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP3TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP3TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = -27. / 24.;
  wgt[1] = 25. / 24.;
  wgt[2] = 25. / 24.;
  wgt[3] = 25. / 24.;
}

void QdPointSetP3TriagGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP3TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                      const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -1. / 3.;
  coords(0, XI1) = -1. / 3.;

  coords(1, XI0) = -0.6;
  coords(1, XI1) = -0.6;

  coords(2, XI0) = 0.2;
  coords(2, XI1) = -0.6;

  coords(3, XI0) = -0.6;
  coords(3, XI1) = 0.2;
}

void QdPointSetP3TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP3TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP3TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP3TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P4 polynomials
// ----------------------------------------------------------------------------

QdPointSetP4TriagGauss::QdPointSetP4TriagGauss() : StdPointSetBase()
{
}

QdPointSetP4TriagGauss::~QdPointSetP4TriagGauss()
{
}

Uint QdPointSetP4TriagGauss::order() const
{
  return P4;
}

Uint QdPointSetP4TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP4TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP4TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP4TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP4TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP4TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.446763179356022;
  wgt[1] = 0.446763179356022;
  wgt[2] = 0.446763179356022;
  wgt[3] = 0.219903487310644;
  wgt[4] = 0.219903487310644;
  wgt[5] = 0.219903487310644;
}

void QdPointSetP4TriagGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP4TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                      const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -0.108103018168070;
  coords(0, XI1) = -0.108103018168070;

  coords(1, XI0) = -0.108103018168070;
  coords(1, XI1) = -0.783793963663860;

  coords(2, XI0) = -0.783793963663860;
  coords(2, XI1) = -0.108103018168070;

  coords(3, XI0) = -0.816847572980458;
  coords(3, XI1) = -0.816847572980458;

  coords(4, XI0) = -0.816847572980458;
  coords(4, XI1) = 0.633695145960918;

  coords(5, XI0) = 0.633695145960918;
  coords(5, XI1) = -0.816847572980458;
}

void QdPointSetP4TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP4TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP4TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP4TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P5 polynomials
// ----------------------------------------------------------------------------

QdPointSetP5TriagGauss::QdPointSetP5TriagGauss() : StdPointSetBase()
{
}

QdPointSetP5TriagGauss::~QdPointSetP5TriagGauss()
{
}

Uint QdPointSetP5TriagGauss::order() const
{
  return P5;
}

Uint QdPointSetP5TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP5TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP5TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP5TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP5TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP5TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0] = 0.450000000000000;
  wgt[1] = 0.264788305577012;
  wgt[2] = 0.264788305577012;
  wgt[3] = 0.264788305577012;
  wgt[4] = 0.251878361089654;
  wgt[5] = 0.251878361089654;
  wgt[6] = 0.251878361089654;
}

void QdPointSetP5TriagGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP5TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                      const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -0.333333333333333;
  coords(0, XI1) = -0.333333333333333;

  coords(1, XI0) = -0.059715871789770;
  coords(1, XI1) = -0.059715871789770;

  coords(2, XI0) = -0.059715871789770;
  coords(2, XI1) = -0.880568256420460;

  coords(3, XI0) = -0.880568256420460;
  coords(3, XI1) = -0.059715871789770;

  coords(4, XI0) = -0.797426985353088;
  coords(4, XI1) = -0.797426985353088;

  coords(5, XI0) = -0.797426985353088;
  coords(5, XI1) = 0.594853970706174;

  coords(6, XI0) = 0.594853970706174;
  coords(6, XI1) = -0.797426985353088;
}

void QdPointSetP5TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP5TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP5TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP5TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P6 polynomials
// ----------------------------------------------------------------------------

QdPointSetP6TriagGauss::QdPointSetP6TriagGauss() : StdPointSetBase()
{
}

QdPointSetP6TriagGauss::~QdPointSetP6TriagGauss()
{
}

Uint QdPointSetP6TriagGauss::order() const
{
  return P6;
}

Uint QdPointSetP6TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP6TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP6TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP6TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP6TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP6TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0]  = 0.233572551452758;
  wgt[1]  = 0.233572551452758;
  wgt[2]  = 0.233572551452758;
  wgt[3]  = 0.101689812740414;
  wgt[4]  = 0.101689812740414;
  wgt[5]  = 0.101689812740414;
  wgt[6]  = 0.165702151236748;
  wgt[7]  = 0.165702151236748;
  wgt[8]  = 0.165702151236748;
  wgt[9]  = 0.165702151236748;
  wgt[10] = 0.165702151236748;
  wgt[11] = 0.165702151236748;
}

void QdPointSetP6TriagGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP6TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                      const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -0.501426509658180;
  coords(0, XI1) = -0.501426509658180;

  coords(1, XI0) = -0.501426509658180;
  coords(1, XI1) = 0.002853019316358;

  coords(2, XI0) = 0.002853019316358;
  coords(2, XI1) = -0.501426509658180;

  coords(3, XI0) = -0.873821971016996;
  coords(3, XI1) = -0.873821971016996;

  coords(4, XI0) = -0.873821971016996;
  coords(4, XI1) = 0.747643942033992;

  coords(5, XI0) = 0.747643942033992;
  coords(5, XI1) = -0.873821971016996;

  coords(6, XI0) = -0.379295097932432;
  coords(6, XI1) = 0.273004998242798;

  coords(7, XI0) = 0.273004998242798;
  coords(7, XI1) = -0.893709900310366;

  coords(8, XI0) = -0.893709900310366;
  coords(8, XI1) = -0.379295097932432;

  coords(9, XI0) = -0.379295097932432;
  coords(9, XI1) = -0.893709900310366;

  coords(10, XI0) = 0.273004998242798;
  coords(10, XI1) = -0.379295097932432;

  coords(11, XI0) = -0.893709900310366;
  coords(11, XI1) = 0.273004998242798;
}

void QdPointSetP6TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP6TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP6TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP6TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P7 polynomials
// ----------------------------------------------------------------------------

QdPointSetP7TriagGauss::QdPointSetP7TriagGauss() : StdPointSetBase()
{
}

QdPointSetP7TriagGauss::~QdPointSetP7TriagGauss()
{
}

Uint QdPointSetP7TriagGauss::order() const
{
  return P7;
}

Uint QdPointSetP7TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP7TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP7TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP7TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP7TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP7TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0]  = -0.299140088935364;
  wgt[1]  = 0.351230514866416;
  wgt[2]  = 0.351230514866416;
  wgt[3]  = 0.351230514866416;
  wgt[4]  = 0.106694471217676;
  wgt[5]  = 0.106694471217676;
  wgt[6]  = 0.106694471217676;
  wgt[7]  = 0.154227521780514;
  wgt[8]  = 0.154227521780514;
  wgt[9]  = 0.154227521780514;
  wgt[10] = 0.154227521780514;
  wgt[11] = 0.154227521780514;
  wgt[12] = 0.154227521780514;
}

void QdPointSetP7TriagGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP7TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                      const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -0.333333333333333;
  coords(0, XI1) = -0.333333333333333;

  coords(1, XI0) = -0.479308067841920;
  coords(1, XI1) = -0.479308067841920;

  coords(2, XI0) = -0.479308067841920;
  coords(2, XI1) = -0.041383864316160;

  coords(3, XI0) = -0.041383864316160;
  coords(3, XI1) = -0.479308067841920;

  coords(4, XI0) = -0.869739794195568;
  coords(4, XI1) = -0.869739794195568;

  coords(5, XI0) = -0.869739794195568;
  coords(5, XI1) = 0.739479588391136;

  coords(6, XI0) = 0.739479588391136;
  coords(6, XI1) = -0.869739794195568;

  coords(7, XI0) = -0.374269007990252;
  coords(7, XI1) = 0.276888377139620;

  coords(8, XI0) = 0.276888377139620;
  coords(8, XI1) = -0.902619369149368;

  coords(9, XI0) = -0.902619369149368;
  coords(9, XI1) = -0.374269007990252;

  coords(10, XI0) = -0.374269007990252;
  coords(10, XI1) = -0.902619369149368;

  coords(11, XI0) = 0.276888377139620;
  coords(11, XI1) = -0.374269007990252;

  coords(12, XI0) = -0.902619369149368;
  coords(12, XI1) = 0.276888377139620;
}

void QdPointSetP7TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP7TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP7TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP7TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P8 polynomials
// ----------------------------------------------------------------------------

QdPointSetP8TriagGauss::QdPointSetP8TriagGauss() : StdPointSetBase()
{
}

QdPointSetP8TriagGauss::~QdPointSetP8TriagGauss()
{
}

Uint QdPointSetP8TriagGauss::order() const
{
  return P8;
}

Uint QdPointSetP8TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP8TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP8TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP8TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP8TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP8TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0]  = 0.288631215355576;
  wgt[1]  = 0.190183268534572;
  wgt[2]  = 0.190183268534572;
  wgt[3]  = 0.190183268534572;
  wgt[4]  = 0.206434741069436;
  wgt[5]  = 0.206434741069436;
  wgt[6]  = 0.206434741069436;
  wgt[7]  = 0.064916995246396;
  wgt[8]  = 0.064916995246396;
  wgt[9]  = 0.064916995246396;
  wgt[10] = 0.054460628348868;
  wgt[11] = 0.054460628348868;
  wgt[12] = 0.054460628348868;
  wgt[13] = 0.054460628348868;
  wgt[14] = 0.054460628348868;
  wgt[15] = 0.054460628348868;
}

void QdPointSetP8TriagGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP8TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                      const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -0.333333333333334;
  coords(0, XI1) = -0.333333333333334;

  coords(1, XI0) = -0.081414823414554;
  coords(1, XI1) = -0.081414823414554;

  coords(2, XI0) = -0.081414823414554;
  coords(2, XI1) = -0.837170353170892;

  coords(3, XI0) = -0.837170353170892;
  coords(3, XI1) = -0.081414823414554;

  coords(4, XI0) = -0.658861384496480;
  coords(4, XI1) = -0.658861384496480;

  coords(5, XI0) = -0.658861384496480;
  coords(5, XI1) = 0.317722768992960;

  coords(6, XI0) = 0.317722768992960;
  coords(6, XI1) = -0.658861384496480;

  coords(7, XI0) = -0.898905543365938;
  coords(7, XI1) = -0.898905543365938;

  coords(8, XI0) = -0.898905543365938;
  coords(8, XI1) = 0.797811086731876;

  coords(9, XI0) = 0.797811086731876;
  coords(9, XI1) = -0.898905543365938;

  coords(10, XI0) = -0.473774340730724;
  coords(10, XI1) = 0.456984785910808;

  coords(11, XI0) = 0.456984785910808;
  coords(11, XI1) = -0.983210445180084;

  coords(12, XI0) = -0.983210445180084;
  coords(12, XI1) = -0.473774340730724;

  coords(13, XI0) = -0.473774340730724;
  coords(13, XI1) = -0.983210445180084;

  coords(14, XI0) = 0.456984785910808;
  coords(14, XI1) = -0.473774340730724;

  coords(15, XI0) = -0.983210445180084;
  coords(15, XI1) = 0.456984785910808;
}

void QdPointSetP8TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP8TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP8TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP8TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P9 polynomials
// ----------------------------------------------------------------------------

QdPointSetP9TriagGauss::QdPointSetP9TriagGauss() : StdPointSetBase()
{
}

QdPointSetP9TriagGauss::~QdPointSetP9TriagGauss()
{
}

Uint QdPointSetP9TriagGauss::order() const
{
  return P9;
}

Uint QdPointSetP9TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP9TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP9TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP9TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP9TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                              const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP9TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0]  = 0.194271592565600;
  wgt[1]  = 0.062669400454200;
  wgt[2]  = 0.062669400454200;
  wgt[3]  = 0.062669400454200;
  wgt[4]  = 0.155655082009600;
  wgt[5]  = 0.155655082009600;
  wgt[6]  = 0.155655082009600;
  wgt[7]  = 0.159295477854400;
  wgt[8]  = 0.159295477854400;
  wgt[9]  = 0.159295477854400;
  wgt[10] = 0.051155351317400;
  wgt[11] = 0.051155351317400;
  wgt[12] = 0.051155351317400;
  wgt[13] = 0.086567078754600;
  wgt[14] = 0.086567078754600;
  wgt[15] = 0.086567078754600;
  wgt[16] = 0.086567078754600;
  wgt[17] = 0.086567078754600;
  wgt[18] = 0.086567078754600;
}

void QdPointSetP9TriagGauss::permutation(const Uint local_id,
                                         const mesh::EntityRealignCode &permutation_code,
                                         std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP9TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                      const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -0.333333333334000;
  coords(0, XI1) = -0.333333333333400;

  coords(1, XI0) = -0.020634961602000;
  coords(1, XI1) = -0.020634961602000;

  coords(2, XI0) = -0.020634961602000;
  coords(2, XI1) = -0.958730076795000;

  coords(3, XI0) = -0.958730076795000;
  coords(3, XI1) = -0.020634961602000;

  coords(4, XI0) = -0.125820817014000;
  coords(4, XI1) = -0.125820817014000;

  coords(5, XI0) = -0.125820817014000;
  coords(5, XI1) = -0.748358365972000;

  coords(6, XI0) = -0.748358365972000;
  coords(6, XI1) = -0.125820817014000;

  coords(7, XI0) = -0.623592928762000;
  coords(7, XI1) = -0.623592928762000;

  coords(8, XI0) = -0.623592928762000;
  coords(8, XI1) = 0.247185857524000;

  coords(9, XI0) = 0.247185857524000;
  coords(9, XI1) = -0.623592928762000;

  coords(10, XI0) = -0.910540973211000;
  coords(10, XI1) = -0.910540973211000;

  coords(11, XI0) = -0.910540973211000;
  coords(11, XI1) = 0.821081946422000;

  coords(12, XI0) = 0.821081946422000;
  coords(12, XI1) = -0.910540973211000;

  coords(13, XI0) = -0.556074021678000;
  coords(13, XI1) = 0.482397197568000;

  coords(14, XI0) = 0.482397197568000;
  coords(14, XI1) = -0.926323175890600;

  coords(15, XI0) = -0.926323175890600;
  coords(15, XI1) = -0.556074021678000;

  coords(16, XI0) = -0.556074021678000;
  coords(16, XI1) = -0.926323175890600;

  coords(17, XI0) = 0.482397197568000;
  coords(17, XI1) = -0.556074021678000;

  coords(18, XI0) = -0.926323175890600;
  coords(18, XI1) = 0.482397197568000;
}

void QdPointSetP9TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP9TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP9TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP9TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P10 polynomials
// ----------------------------------------------------------------------------

QdPointSetP10TriagGauss::QdPointSetP10TriagGauss() : StdPointSetBase()
{
}

QdPointSetP10TriagGauss::~QdPointSetP10TriagGauss()
{
}

Uint QdPointSetP10TriagGauss::order() const
{
  return P10;
}

Uint QdPointSetP10TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP10TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP10TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP10TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP10TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                               const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP10TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0]  = 0.181635980765600;
  wgt[1]  = 0.073451915513000;
  wgt[2]  = 0.073451915513000;
  wgt[3]  = 0.073451915513000;
  wgt[4]  = 0.090642118871000;
  wgt[5]  = 0.090642118871000;
  wgt[6]  = 0.090642118871000;
  wgt[7]  = 0.145515833690800;
  wgt[8]  = 0.145515833690800;
  wgt[9]  = 0.145515833690800;
  wgt[10] = 0.145515833690800;
  wgt[11] = 0.145515833690800;
  wgt[12] = 0.145515833690800;
  wgt[13] = 0.056654485062200;
  wgt[14] = 0.056654485062200;
  wgt[15] = 0.056654485062200;
  wgt[16] = 0.056654485062200;
  wgt[17] = 0.056654485062200;
  wgt[18] = 0.056654485062200;
  wgt[19] = 0.018843333927400;
  wgt[20] = 0.018843333927400;
  wgt[21] = 0.018843333927400;
  wgt[22] = 0.018843333927400;
  wgt[23] = 0.018843333927400;
  wgt[24] = 0.018843333927400;
}

void QdPointSetP10TriagGauss::permutation(const Uint local_id,
                                          const mesh::EntityRealignCode &permutation_code,
                                          std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP10TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                       const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -0.333333333334000;
  coords(0, XI1) = -0.333333333334000;

  coords(1, XI0) = -0.028844733232000;
  coords(1, XI1) = -0.028844733232000;

  coords(2, XI0) = -0.028844733232000;
  coords(2, XI1) = -0.942310533534600;

  coords(3, XI0) = -0.942310533534600;
  coords(3, XI1) = -0.028844733232000;

  coords(4, XI0) = -0.781036849030000;
  coords(4, XI1) = -0.781036849030000;

  coords(5, XI0) = -0.781036849030000;
  coords(5, XI1) = 0.562073698060000;

  coords(6, XI0) = 0.562073698060000;
  coords(6, XI1) = -0.781036849030000;

  coords(7, XI0) = -0.384120322472000;
  coords(7, XI1) = 0.100705883642000;

  coords(8, XI0) = 0.100705883642000;
  coords(8, XI1) = -0.716585561170000;

  coords(9, XI0) = -0.716585561170000;
  coords(9, XI1) = -0.384120322472000;

  coords(10, XI0) = -0.384120322472000;
  coords(10, XI1) = -0.716585561170000;

  coords(11, XI0) = 0.100705883642000;
  coords(11, XI1) = -0.384120322472000;

  coords(12, XI0) = -0.716585561170000;
  coords(12, XI1) = 0.100705883642000;

  coords(13, XI0) = -0.506654878720000;
  coords(13, XI1) = 0.456647809194000;

  coords(14, XI0) = 0.456647809194000;
  coords(14, XI1) = -0.949992930474600;

  coords(15, XI0) = -0.949992930474600;
  coords(15, XI1) = -0.506654878720000;

  coords(16, XI0) = -0.506654878720000;
  coords(16, XI1) = -0.949992930474600;

  coords(17, XI0) = 0.456647809194000;
  coords(17, XI1) = -0.506654878720000;

  coords(18, XI0) = -0.949992930474600;
  coords(18, XI1) = 0.456647809194000;

  coords(19, XI0) = -0.866393497975600;
  coords(19, XI1) = 0.847311867174000;

  coords(20, XI0) = 0.847311867174000;
  coords(20, XI1) = -0.980918369199400;

  coords(21, XI0) = -0.980918369199400;
  coords(21, XI1) = -0.866393497975600;

  coords(22, XI0) = -0.866393497975600;
  coords(22, XI1) = -0.980918369199400;

  coords(23, XI0) = 0.847311867174000;
  coords(23, XI1) = -0.866393497975600;

  coords(24, XI0) = -0.980918369199400;
  coords(24, XI1) = 0.847311867174000;
}

void QdPointSetP10TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP10TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP10TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP10TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P11 polynomials
// ----------------------------------------------------------------------------

QdPointSetP11TriagGauss::QdPointSetP11TriagGauss() : StdPointSetBase()
{
}

QdPointSetP11TriagGauss::~QdPointSetP11TriagGauss()
{
}

Uint QdPointSetP11TriagGauss::order() const
{
  return P11;
}

Uint QdPointSetP11TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP11TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP11TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP11TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP11TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                               const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP11TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0]  = 0.001854012657920;
  wgt[1]  = 0.001854012657920;
  wgt[2]  = 0.001854012657920;
  wgt[3]  = 0.154299069829600;
  wgt[4]  = 0.154299069829600;
  wgt[5]  = 0.154299069829600;
  wgt[6]  = 0.118645954761600;
  wgt[7]  = 0.118645954761600;
  wgt[8]  = 0.118645954761600;
  wgt[9]  = 0.072369081006800;
  wgt[10] = 0.072369081006800;
  wgt[11] = 0.072369081006800;
  wgt[12] = 0.027319462005400;
  wgt[13] = 0.027319462005400;
  wgt[14] = 0.027319462005400;
  wgt[15] = 0.104674223924400;
  wgt[16] = 0.104674223924400;
  wgt[17] = 0.104674223924400;
  wgt[18] = 0.104674223924400;
  wgt[19] = 0.104674223924400;
  wgt[20] = 0.104674223924400;
  wgt[21] = 0.041415319278200;
  wgt[22] = 0.041415319278200;
  wgt[23] = 0.041415319278200;
  wgt[24] = 0.041415319278200;
  wgt[25] = 0.041415319278200;
  wgt[26] = 0.041415319278200;
}

void QdPointSetP11TriagGauss::permutation(const Uint local_id,
                                          const mesh::EntityRealignCode &permutation_code,
                                          std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP11TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                       const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = 0.069222096542000;
  coords(0, XI1) = 0.069222096542000;

  coords(1, XI0) = -1.138444193083000;
  coords(1, XI1) = 0.069222096542000;

  coords(2, XI0) = 0.069222096542000;
  coords(2, XI1) = -1.138444193083000;

  coords(3, XI0) = -0.202061394068000;
  coords(3, XI1) = -0.202061394068000;

  coords(4, XI0) = -0.595877211864000;
  coords(4, XI1) = -0.202061394068000;

  coords(5, XI0) = -0.202061394068000;
  coords(5, XI1) = -0.595877211864000;

  coords(6, XI0) = -0.593380199138000;
  coords(6, XI1) = -0.593380199138000;

  coords(7, XI0) = 0.186760398274000;
  coords(7, XI1) = -0.593380199138000;

  coords(8, XI0) = -0.593380199138000;
  coords(8, XI1) = 0.186760398274000;

  coords(9, XI0) = -0.761298175434000;
  coords(9, XI1) = -0.761298175434000;

  coords(10, XI0) = 0.522596350870000;
  coords(10, XI1) = -0.761298175434000;

  coords(11, XI0) = -0.761298175434000;
  coords(11, XI1) = 0.522596350870000;

  coords(12, XI0) = -0.935270103777400;
  coords(12, XI1) = -0.935270103777400;

  coords(13, XI0) = 0.870540207554000;
  coords(13, XI1) = -0.935270103777400;

  coords(14, XI0) = -0.935270103777400;
  coords(14, XI1) = 0.870540207554000;

  coords(15, XI0) = 0.186402426856000;
  coords(15, XI1) = -0.286758703478000;

  coords(16, XI0) = -0.899643723379000;
  coords(16, XI1) = 0.186402426856000;

  coords(17, XI0) = -0.286758703478000;
  coords(17, XI1) = -0.899643723379000;

  coords(18, XI0) = -0.899643723379000;
  coords(18, XI1) = -0.286758703478000;

  coords(19, XI0) = -0.286758703478000;
  coords(19, XI1) = 0.186402426856000;

  coords(20, XI0) = 0.186402426856000;
  coords(20, XI1) = -0.899643723379000;

  coords(21, XI0) = 0.614978006320000;
  coords(21, XI1) = -0.657022039392000;

  coords(22, XI0) = -0.957955966927600;
  coords(22, XI1) = 0.614978006320000;

  coords(23, XI0) = -0.657022039392000;
  coords(23, XI1) = -0.957955966927600;

  coords(24, XI0) = -0.957955966927600;
  coords(24, XI1) = -0.657022039392000;

  coords(25, XI0) = -0.657022039392000;
  coords(25, XI1) = 0.614978006320000;

  coords(26, XI0) = 0.614978006320000;
  coords(26, XI1) = -0.957955966927600;
}

void QdPointSetP11TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP11TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP11TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP11TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P12 polynomials
// ----------------------------------------------------------------------------

QdPointSetP12TriagGauss::QdPointSetP12TriagGauss() : StdPointSetBase()
{
}

QdPointSetP12TriagGauss::~QdPointSetP12TriagGauss()
{
}

Uint QdPointSetP12TriagGauss::order() const
{
  return P12;
}

Uint QdPointSetP12TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP12TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP12TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP12TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP12TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                               const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP12TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0]  = 0.051462132881000;
  wgt[1]  = 0.051462132881000;
  wgt[2]  = 0.051462132881000;
  wgt[3]  = 0.087385089076000;
  wgt[4]  = 0.087385089076000;
  wgt[5]  = 0.087385089076000;
  wgt[6]  = 0.125716448435800;
  wgt[7]  = 0.125716448435800;
  wgt[8]  = 0.125716448435800;
  wgt[9]  = 0.069592225861400;
  wgt[10] = 0.069592225861400;
  wgt[11] = 0.069592225861400;
  wgt[12] = 0.012332522103120;
  wgt[13] = 0.012332522103120;
  wgt[14] = 0.012332522103120;
  wgt[15] = 0.080743115532800;
  wgt[16] = 0.080743115532800;
  wgt[17] = 0.080743115532800;
  wgt[18] = 0.080743115532800;
  wgt[19] = 0.080743115532800;
  wgt[20] = 0.080743115532800;
  wgt[21] = 0.044713546404600;
  wgt[22] = 0.044713546404600;
  wgt[23] = 0.044713546404600;
  wgt[24] = 0.044713546404600;
  wgt[25] = 0.044713546404600;
  wgt[26] = 0.044713546404600;
  wgt[27] = 0.034632462217400;
  wgt[28] = 0.034632462217400;
  wgt[29] = 0.034632462217400;
  wgt[30] = 0.034632462217400;
  wgt[31] = 0.034632462217400;
  wgt[32] = 0.034632462217400;
}

void QdPointSetP12TriagGauss::permutation(const Uint local_id,
                                          const mesh::EntityRealignCode &permutation_code,
                                          std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP12TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                       const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -0.023565220452000;
  coords(0, XI1) = -0.023565220452000;

  coords(1, XI0) = -0.023565220452000;
  coords(1, XI1) = -0.952869559095200;

  coords(2, XI0) = -0.952869559095200;
  coords(2, XI1) = -0.023565220452000;

  coords(3, XI0) = -0.120551215412000;
  coords(3, XI1) = -0.120551215412000;

  coords(4, XI0) = -0.120551215412000;
  coords(4, XI1) = -0.758897569178000;

  coords(5, XI0) = -0.758897569178000;
  coords(5, XI1) = -0.120551215412000;

  coords(6, XI0) = -0.457579229976000;
  coords(6, XI1) = -0.457579229976000;

  coords(7, XI0) = -0.457579229976000;
  coords(7, XI1) = -0.084841540048000;

  coords(8, XI0) = -0.084841540048000;
  coords(8, XI1) = -0.457579229976000;

  coords(9, XI0) = -0.744847708916000;
  coords(9, XI1) = -0.744847708916000;

  coords(10, XI0) = -0.744847708916000;
  coords(10, XI1) = 0.489695417834000;

  coords(11, XI0) = 0.489695417834000;
  coords(11, XI1) = -0.744847708916000;

  coords(12, XI0) = -0.957365299093600;
  coords(12, XI1) = -0.957365299093600;

  coords(13, XI0) = -0.957365299093600;
  coords(13, XI1) = 0.914730598188000;

  coords(14, XI0) = 0.914730598188000;
  coords(14, XI1) = -0.957365299093600;

  coords(15, XI0) = -0.448573460628000;
  coords(15, XI1) = 0.217886471560000;

  coords(16, XI0) = 0.217886471560000;
  coords(16, XI1) = -0.769313010930000;

  coords(17, XI0) = -0.769313010930000;
  coords(17, XI1) = -0.448573460628000;

  coords(18, XI0) = -0.448573460628000;
  coords(18, XI1) = -0.769313010930000;

  coords(19, XI0) = 0.217886471560000;
  coords(19, XI1) = -0.448573460628000;

  coords(20, XI0) = -0.769313010930000;
  coords(20, XI1) = 0.217886471560000;

  coords(21, XI0) = -0.437348838020000;
  coords(21, XI1) = 0.391672173576000;

  coords(22, XI0) = 0.391672173576000;
  coords(22, XI1) = -0.954323335555400;

  coords(23, XI0) = -0.954323335555400;
  coords(23, XI1) = -0.437348838020000;

  coords(24, XI0) = -0.437348838020000;
  coords(24, XI1) = -0.954323335555400;

  coords(25, XI0) = 0.391672173576000;
  coords(25, XI1) = -0.437348838020000;

  coords(26, XI0) = -0.954323335555400;
  coords(26, XI1) = 0.391672173576000;

  coords(27, XI0) = -0.767496168184000;
  coords(27, XI1) = 0.716028067088000;

  coords(28, XI0) = 0.716028067088000;
  coords(28, XI1) = -0.948531898903400;

  coords(29, XI0) = -0.948531898903400;
  coords(29, XI1) = -0.767496168184000;

  coords(30, XI0) = -0.767496168184000;
  coords(30, XI1) = -0.948531898903400;

  coords(31, XI0) = 0.716028067088000;
  coords(31, XI1) = -0.767496168184000;

  coords(32, XI0) = -0.948531898903400;
  coords(32, XI1) = 0.716028067088000;
}

void QdPointSetP12TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP12TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP12TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP12TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P13 polynomials
// ----------------------------------------------------------------------------

QdPointSetP13TriagGauss::QdPointSetP13TriagGauss() : StdPointSetBase()
{
}

QdPointSetP13TriagGauss::~QdPointSetP13TriagGauss()
{
}

Uint QdPointSetP13TriagGauss::order() const
{
  return P13;
}

Uint QdPointSetP13TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP13TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP13TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP13TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP13TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                               const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP13TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0]  = 0.105041846801600;
  wgt[1]  = 0.022560290418600;
  wgt[2]  = 0.022560290418600;
  wgt[3]  = 0.022560290418600;
  wgt[4]  = 0.062847036725000;
  wgt[5]  = 0.062847036725000;
  wgt[6]  = 0.062847036725000;
  wgt[7]  = 0.094145005008400;
  wgt[8]  = 0.094145005008400;
  wgt[9]  = 0.094145005008400;
  wgt[10] = 0.094727173072800;
  wgt[11] = 0.094727173072800;
  wgt[12] = 0.094727173072800;
  wgt[13] = 0.062335058091600;
  wgt[14] = 0.062335058091600;
  wgt[15] = 0.062335058091600;
  wgt[16] = 0.015951542930140;
  wgt[17] = 0.015951542930140;
  wgt[18] = 0.015951542930140;
  wgt[19] = 0.073696805457400;
  wgt[20] = 0.073696805457400;
  wgt[21] = 0.073696805457400;
  wgt[22] = 0.073696805457400;
  wgt[23] = 0.073696805457400;
  wgt[24] = 0.073696805457400;
  wgt[25] = 0.034802926607600;
  wgt[26] = 0.034802926607600;
  wgt[27] = 0.034802926607600;
  wgt[28] = 0.034802926607600;
  wgt[29] = 0.034802926607600;
  wgt[30] = 0.034802926607600;
  wgt[31] = 0.031043573678000;
  wgt[32] = 0.031043573678000;
  wgt[33] = 0.031043573678000;
  wgt[34] = 0.031043573678000;
  wgt[35] = 0.031043573678000;
  wgt[36] = 0.031043573678000;
}

void QdPointSetP13TriagGauss::permutation(const Uint local_id,
                                          const mesh::EntityRealignCode &permutation_code,
                                          std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP13TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                       const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -0.333333333333400;
  coords(0, XI1) = -0.333333333333400;

  coords(1, XI0) = -0.009903630120000;
  coords(1, XI1) = -0.009903630120000;

  coords(2, XI0) = -0.980192739758820;
  coords(2, XI1) = -0.009903630120000;

  coords(3, XI0) = -0.009903630120000;
  coords(3, XI1) = -0.980192739758820;

  coords(4, XI0) = -0.062566729780000;
  coords(4, XI1) = -0.062566729780000;

  coords(5, XI0) = -0.874866540438200;
  coords(5, XI1) = -0.062566729780000;

  coords(6, XI0) = -0.062566729780000;
  coords(6, XI1) = -0.874866540438200;

  coords(7, XI0) = -0.170957326398000;
  coords(7, XI1) = -0.170957326398000;

  coords(8, XI0) = -0.658085347206000;
  coords(8, XI1) = -0.170957326398000;

  coords(9, XI0) = -0.170957326398000;
  coords(9, XI1) = -0.658085347206000;

  coords(10, XI0) = -0.541200855914000;
  coords(10, XI1) = -0.541200855914000;

  coords(11, XI0) = 0.082401711828000;
  coords(11, XI1) = -0.541200855914000;

  coords(12, XI0) = -0.541200855914000;
  coords(12, XI1) = 0.082401711828000;

  coords(13, XI0) = -0.771151009608000;
  coords(13, XI1) = -0.771151009608000;

  coords(14, XI0) = 0.542302019214000;
  coords(14, XI1) = -0.771151009608000;

  coords(15, XI0) = -0.771151009608000;
  coords(15, XI1) = 0.542302019214000;

  coords(16, XI0) = -0.950377217273000;
  coords(16, XI1) = -0.950377217273000;

  coords(17, XI0) = 0.900754434546000;
  coords(17, XI1) = -0.950377217273000;

  coords(18, XI0) = -0.950377217273000;
  coords(18, XI1) = 0.900754434546000;

  coords(19, XI0) = 0.272702349124000;
  coords(19, XI1) = -0.462410005882000;

  coords(20, XI0) = -0.810292343240800;
  coords(20, XI1) = 0.272702349124000;

  coords(21, XI0) = -0.462410005882000;
  coords(21, XI1) = -0.810292343240800;

  coords(22, XI0) = -0.810292343240800;
  coords(22, XI1) = -0.462410005882000;

  coords(23, XI0) = -0.462410005882000;
  coords(23, XI1) = 0.272702349124000;

  coords(24, XI0) = 0.272702349124000;
  coords(24, XI1) = -0.810292343240800;

  coords(25, XI0) = 0.380338319974000;
  coords(25, XI1) = -0.416539866532000;

  coords(26, XI0) = -0.963798453442400;
  coords(26, XI1) = 0.380338319974000;

  coords(27, XI0) = -0.416539866532000;
  coords(27, XI1) = -0.963798453442400;

  coords(28, XI0) = -0.963798453442400;
  coords(28, XI1) = -0.416539866532000;

  coords(29, XI0) = -0.416539866532000;
  coords(29, XI1) = 0.380338319974000;

  coords(30, XI0) = 0.380338319974000;
  coords(30, XI1) = -0.963798453442400;

  coords(31, XI0) = 0.702819075668000;
  coords(31, XI1) = -0.747285229016000;

  coords(32, XI0) = -0.955533846651800;
  coords(32, XI1) = 0.702819075668000;

  coords(33, XI0) = -0.747285229016000;
  coords(33, XI1) = -0.955533846651800;

  coords(34, XI0) = -0.955533846651800;
  coords(34, XI1) = -0.747285229016000;

  coords(35, XI0) = -0.747285229016000;
  coords(35, XI1) = 0.702819075668000;

  coords(36, XI0) = 0.702819075668000;
  coords(36, XI1) = -0.955533846651800;
}

void QdPointSetP13TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP13TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP13TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP13TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P14 polynomials
// ----------------------------------------------------------------------------

QdPointSetP14TriagGauss::QdPointSetP14TriagGauss() : StdPointSetBase()
{
}

QdPointSetP14TriagGauss::~QdPointSetP14TriagGauss()
{
}

Uint QdPointSetP14TriagGauss::order() const
{
  return P14;
}

Uint QdPointSetP14TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP14TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP14TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP14TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP14TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                               const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP14TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0]  = 0.043767162738800;
  wgt[1]  = 0.043767162738800;
  wgt[2]  = 0.043767162738800;
  wgt[3]  = 0.065576707088200;
  wgt[4]  = 0.065576707088200;
  wgt[5]  = 0.065576707088200;
  wgt[6]  = 0.103548209014600;
  wgt[7]  = 0.103548209014600;
  wgt[8]  = 0.103548209014600;
  wgt[9]  = 0.084325177474000;
  wgt[10] = 0.084325177474000;
  wgt[11] = 0.084325177474000;
  wgt[12] = 0.028867399339600;
  wgt[13] = 0.028867399339600;
  wgt[14] = 0.028867399339600;
  wgt[15] = 0.009846807204800;
  wgt[16] = 0.009846807204800;
  wgt[17] = 0.009846807204800;
  wgt[18] = 0.049331506425200;
  wgt[19] = 0.049331506425200;
  wgt[20] = 0.049331506425200;
  wgt[21] = 0.049331506425200;
  wgt[22] = 0.049331506425200;
  wgt[23] = 0.049331506425200;
  wgt[24] = 0.077143021574200;
  wgt[25] = 0.077143021574200;
  wgt[26] = 0.077143021574200;
  wgt[27] = 0.077143021574200;
  wgt[28] = 0.077143021574200;
  wgt[29] = 0.077143021574200;
  wgt[30] = 0.028872616227000;
  wgt[31] = 0.028872616227000;
  wgt[32] = 0.028872616227000;
  wgt[33] = 0.028872616227000;
  wgt[34] = 0.028872616227000;
  wgt[35] = 0.028872616227000;
  wgt[36] = 0.010020457677000;
  wgt[37] = 0.010020457677000;
  wgt[38] = 0.010020457677000;
  wgt[39] = 0.010020457677000;
  wgt[40] = 0.010020457677000;
  wgt[41] = 0.010020457677000;
}

void QdPointSetP14TriagGauss::permutation(const Uint local_id,
                                          const mesh::EntityRealignCode &permutation_code,
                                          std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP14TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                       const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = -0.022072179276000;
  coords(0, XI1) = -0.022072179276000;

  coords(1, XI0) = -0.955855641448800;
  coords(1, XI1) = -0.022072179276000;

  coords(2, XI0) = -0.022072179276000;
  coords(2, XI1) = -0.955855641448800;

  coords(3, XI0) = -0.164710561320000;
  coords(3, XI1) = -0.164710561320000;

  coords(4, XI0) = -0.670578877362000;
  coords(4, XI1) = -0.164710561320000;

  coords(5, XI0) = -0.164710561320000;
  coords(5, XI1) = -0.670578877362000;

  coords(6, XI0) = -0.453044943382000;
  coords(6, XI1) = -0.453044943382000;

  coords(7, XI0) = -0.093910113236000;
  coords(7, XI1) = -0.453044943382000;

  coords(8, XI0) = -0.453044943382000;
  coords(8, XI1) = -0.093910113236000;

  coords(9, XI0) = -0.645588935174000;
  coords(9, XI1) = -0.645588935174000;

  coords(10, XI0) = 0.291177870350000;
  coords(10, XI1) = -0.645588935174000;

  coords(11, XI0) = -0.645588935174000;
  coords(11, XI1) = 0.291177870350000;

  coords(12, XI0) = -0.876400233818200;
  coords(12, XI1) = -0.876400233818200;

  coords(13, XI0) = 0.752800467636000;
  coords(13, XI1) = -0.876400233818200;

  coords(14, XI0) = -0.876400233818200;
  coords(14, XI1) = 0.752800467636000;

  coords(15, XI0) = -0.961218077502600;
  coords(15, XI1) = -0.961218077502600;

  coords(16, XI0) = 0.922436155006000;
  coords(16, XI1) = -0.961218077502600;

  coords(17, XI0) = -0.961218077502600;
  coords(17, XI1) = 0.922436155006000;

  coords(18, XI0) = 0.541217109550000;
  coords(18, XI1) = -0.655466624358000;

  coords(19, XI0) = -0.885750485192800;
  coords(19, XI1) = 0.541217109550000;

  coords(20, XI0) = -0.655466624358000;
  coords(20, XI1) = -0.885750485192800;

  coords(21, XI0) = -0.885750485192800;
  coords(21, XI1) = -0.655466624358000;

  coords(22, XI0) = -0.655466624358000;
  coords(22, XI1) = 0.541217109550000;

  coords(23, XI0) = 0.541217109550000;
  coords(23, XI1) = -0.885750485192800;

  coords(24, XI0) = 0.140444581694000;
  coords(24, XI1) = -0.326277080408000;

  coords(25, XI0) = -0.814167501286000;
  coords(25, XI1) = 0.140444581694000;

  coords(26, XI0) = -0.326277080408000;
  coords(26, XI1) = -0.814167501286000;

  coords(27, XI0) = -0.814167501286000;
  coords(27, XI1) = -0.326277080408000;

  coords(28, XI0) = -0.326277080408000;
  coords(28, XI1) = 0.140444581694000;

  coords(29, XI0) = 0.140444581694000;
  coords(29, XI1) = -0.814167501286000;

  coords(30, XI0) = 0.373960335616000;
  coords(30, XI1) = -0.403254235728000;

  coords(31, XI0) = -0.970706099888600;
  coords(31, XI1) = 0.373960335616000;

  coords(32, XI0) = -0.403254235728000;
  coords(32, XI1) = -0.970706099888600;

  coords(33, XI0) = -0.970706099888600;
  coords(33, XI1) = -0.403254235728000;

  coords(34, XI0) = -0.403254235728000;
  coords(34, XI1) = 0.373960335616000;

  coords(35, XI0) = 0.373960335616000;
  coords(35, XI1) = -0.970706099888600;

  coords(36, XI0) = 0.759514342740000;
  coords(36, XI1) = -0.762051004606000;

  coords(37, XI0) = -0.997463338134260;
  coords(37, XI1) = 0.759514342740000;

  coords(38, XI0) = -0.762051004606000;
  coords(38, XI1) = -0.997463338134260;

  coords(39, XI0) = -0.997463338134260;
  coords(39, XI1) = -0.762051004606000;

  coords(40, XI0) = -0.762051004606000;
  coords(40, XI1) = 0.759514342740000;

  coords(41, XI0) = 0.759514342740000;
  coords(41, XI1) = -0.997463338134260;
}

void QdPointSetP14TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP14TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP14TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP14TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------
// Triangle quadrature for P15 polynomials
// ----------------------------------------------------------------------------

QdPointSetP15TriagGauss::QdPointSetP15TriagGauss() : StdPointSetBase()
{
}

QdPointSetP15TriagGauss::~QdPointSetP15TriagGauss()
{
}

Uint QdPointSetP15TriagGauss::order() const
{
  return P15;
}

Uint QdPointSetP15TriagGauss::dim() const
{
  return _2D;
}

Uint QdPointSetP15TriagGauss::codim() const
{
  return _0D;
}

Uint QdPointSetP15TriagGauss::nb_local_entities() const
{
  return 1u;
}

Uint QdPointSetP15TriagGauss::size(const Uint local_idx) const
{
  return N_QD_PTS;
}

void QdPointSetP15TriagGauss::reference_coords(math::DenseDMat<Real> &coords,
                                               const Uint local_idx) const
{
  compute_reference_coords(coords, local_idx);
}

void QdPointSetP15TriagGauss::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(N_QD_PTS);
  wgt[0]  = 0.003833751285700;
  wgt[1]  = 0.003833751285700;
  wgt[2]  = 0.003833751285700;
  wgt[3]  = 0.088498054542200;
  wgt[4]  = 0.088498054542200;
  wgt[5]  = 0.088498054542200;
  wgt[6]  = 0.102373097437800;
  wgt[7]  = 0.102373097437800;
  wgt[8]  = 0.102373097437800;
  wgt[9]  = 0.047375471741400;
  wgt[10] = 0.047375471741400;
  wgt[11] = 0.047375471741400;
  wgt[12] = 0.026579551380000;
  wgt[13] = 0.026579551380000;
  wgt[14] = 0.026579551380000;
  wgt[15] = 0.009497833216380;
  wgt[16] = 0.009497833216380;
  wgt[17] = 0.009497833216380;
  wgt[18] = 0.077100145199200;
  wgt[19] = 0.077100145199200;
  wgt[20] = 0.077100145199200;
  wgt[21] = 0.077100145199200;
  wgt[22] = 0.077100145199200;
  wgt[23] = 0.077100145199200;
  wgt[24] = 0.054431628641200;
  wgt[25] = 0.054431628641200;
  wgt[26] = 0.054431628641200;
  wgt[27] = 0.054431628641200;
  wgt[28] = 0.054431628641200;
  wgt[29] = 0.054431628641200;
  wgt[30] = 0.004364154733600;
  wgt[31] = 0.004364154733600;
  wgt[32] = 0.004364154733600;
  wgt[33] = 0.004364154733600;
  wgt[34] = 0.004364154733600;
  wgt[35] = 0.004364154733600;
  wgt[36] = 0.043010639695400;
  wgt[37] = 0.043010639695400;
  wgt[38] = 0.043010639695400;
  wgt[39] = 0.043010639695400;
  wgt[40] = 0.043010639695400;
  wgt[41] = 0.043010639695400;
  wgt[42] = 0.015347885262100;
  wgt[43] = 0.015347885262100;
  wgt[44] = 0.015347885262100;
  wgt[45] = 0.015347885262100;
  wgt[46] = 0.015347885262100;
  wgt[47] = 0.015347885262100;
}

void QdPointSetP15TriagGauss::permutation(const Uint local_id,
                                          const mesh::EntityRealignCode &permutation_code,
                                          std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  if (!canonical_permutations_initialized)
  {
    initialize_canonical_permutations();
  }
  detail::fill_triag_quadrature_permutation(N_QD_PTS, permutation_code, canonical_rot_permutation,
                                            canonical_flip_permutation, permutation_vec);
}

void QdPointSetP15TriagGauss::compute_reference_coords(math::DenseDMat<Real> &coords,
                                                       const Uint local_idx)
{
  coords.resize(N_QD_PTS, 2);
  coords(0, XI0) = 0.013945833716000;
  coords(0, XI1) = 0.013945833716000;

  coords(1, XI0) = -1.027891667433000;
  coords(1, XI1) = 0.013945833716000;

  coords(2, XI0) = 0.013945833716000;
  coords(2, XI1) = -1.027891667433000;

  coords(3, XI0) = -0.137187291434000;
  coords(3, XI1) = -0.137187291434000;

  coords(4, XI0) = -0.725625417132000;
  coords(4, XI1) = -0.137187291434000;

  coords(5, XI0) = -0.137187291434000;
  coords(5, XI1) = -0.725625417132000;

  coords(6, XI0) = -0.444612710306000;
  coords(6, XI1) = -0.444612710306000;

  coords(7, XI0) = -0.110774579388000;
  coords(7, XI1) = -0.444612710306000;

  coords(8, XI0) = -0.444612710306000;
  coords(8, XI1) = -0.110774579388000;

  coords(9, XI0) = -0.747070217918000;
  coords(9, XI1) = -0.747070217918000;

  coords(10, XI0) = 0.494140435834000;
  coords(10, XI1) = -0.747070217918000;

  coords(11, XI0) = -0.747070217918000;
  coords(11, XI1) = 0.494140435834000;

  coords(12, XI0) = -0.858383228050600;
  coords(12, XI1) = -0.858383228050600;

  coords(13, XI0) = 0.716766456102000;
  coords(13, XI1) = -0.858383228050600;

  coords(14, XI0) = -0.858383228050600;
  coords(14, XI1) = 0.716766456102000;

  coords(15, XI0) = -0.962069659517800;
  coords(15, XI1) = -0.962069659517800;

  coords(16, XI0) = 0.924139319036000;
  coords(16, XI1) = -0.962069659517800;

  coords(17, XI0) = -0.962069659517800;
  coords(17, XI1) = 0.924139319036000;

  coords(18, XI0) = 0.209908933786000;
  coords(18, XI1) = -0.477377257720000;

  coords(19, XI0) = -0.732531676066000;
  coords(19, XI1) = 0.209908933786000;

  coords(20, XI0) = -0.477377257720000;
  coords(20, XI1) = -0.732531676066000;

  coords(21, XI0) = -0.732531676066000;
  coords(21, XI1) = -0.477377257720000;

  coords(22, XI0) = -0.477377257720000;
  coords(22, XI1) = 0.209908933786000;

  coords(23, XI0) = 0.209908933786000;
  coords(23, XI1) = -0.732531676066000;

  coords(24, XI0) = 0.151173111026000;
  coords(24, XI1) = -0.223906465820000;

  coords(25, XI0) = -0.927266645206200;
  coords(25, XI1) = 0.151173111026000;

  coords(26, XI0) = -0.223906465820000;
  coords(26, XI1) = -0.927266645206200;

  coords(27, XI0) = -0.927266645206200;
  coords(27, XI1) = -0.223906465820000;

  coords(28, XI0) = -0.223906465820000;
  coords(28, XI1) = 0.151173111026000;

  coords(29, XI0) = 0.151173111026000;
  coords(29, XI1) = -0.927266645206200;

  coords(30, XI0) = 0.448925326154000;
  coords(30, XI1) = -0.428575559900000;

  coords(31, XI0) = -1.020349766253200;
  coords(31, XI1) = 0.448925326154000;

  coords(32, XI0) = -0.428575559900000;
  coords(32, XI1) = -1.020349766253200;

  coords(33, XI0) = -1.020349766253200;
  coords(33, XI1) = -0.428575559900000;

  coords(34, XI0) = -0.428575559900000;
  coords(34, XI1) = 0.448925326154000;

  coords(35, XI0) = 0.448925326154000;
  coords(35, XI1) = -1.020349766253200;

  coords(36, XI0) = 0.495112932104000;
  coords(36, XI1) = -0.568800671856000;

  coords(37, XI0) = -0.926312260248200;
  coords(37, XI1) = 0.495112932104000;

  coords(38, XI0) = -0.568800671856000;
  coords(38, XI1) = -0.926312260248200;

  coords(39, XI0) = -0.926312260248200;
  coords(39, XI1) = -0.568800671856000;

  coords(40, XI0) = -0.568800671856000;
  coords(40, XI1) = 0.495112932104000;

  coords(41, XI0) = 0.495112932104000;
  coords(41, XI1) = -0.926312260248200;

  coords(42, XI0) = 0.767929148184000;
  coords(42, XI1) = -0.792848766848000;

  coords(43, XI0) = -0.975080381337600;
  coords(43, XI1) = 0.767929148184000;

  coords(44, XI0) = -0.792848766848000;
  coords(44, XI1) = -0.975080381337600;

  coords(45, XI0) = -0.975080381337600;
  coords(45, XI1) = -0.792848766848000;

  coords(46, XI0) = -0.792848766848000;
  coords(46, XI1) = 0.767929148184000;

  coords(47, XI0) = 0.767929148184000;
  coords(47, XI1) = -0.975080381337600;
}

void QdPointSetP15TriagGauss::initialize_canonical_permutations()
{
  math::DenseDMat<Real> coords;
  compute_reference_coords(coords, 0);

  detail::canonical_triag_quadrature_rotation_permutation(coords, canonical_rot_permutation);
  detail::canonical_triag_quadrature_flip_permutation(coords, canonical_flip_permutation);
  canonical_permutations_initialized = true;
}

bool QdPointSetP15TriagGauss::canonical_permutations_initialized = false;

std::vector<Uint> QdPointSetP15TriagGauss::canonical_rot_permutation = {};

std::vector<Uint> QdPointSetP15TriagGauss::canonical_flip_permutation = {};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
