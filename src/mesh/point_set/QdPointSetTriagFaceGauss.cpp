#include "mesh/point_set/QdPointSetTriagFaceGauss.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// HELPER FUNCTIONS
// ----------------------------------------------------------------------------

namespace detail
{

void fill_triag_face_quadrature_permutation(const Uint local_id, const Uint nb_qd_pts,
                                            const mesh::EntityRealignCode &permutation_code,
                                            std::vector<Uint> &permutation_vec)
{
  // NOTE THAT THE PERMUTATION WILL BE THE SAME FOR ALL FACES (LOCAL ENTITIES)
  // REGARDLESS OF THEIR LOCAL ID

  // First set permutation_vec as identity permutation: p[i] = i;

  permutation_vec.resize(nb_qd_pts);
  for (Uint i = 0; i < permutation_vec.size(); ++i)
  {
    permutation_vec[i] = i;
  }

  // If this permutation is identity, we are done
  if (permutation_code.is_identity(ElemShape::Line))
  {
    return;
  }

  // Temporary vector
  std::vector<Uint> old_permutation(nb_qd_pts);
  std::vector<Uint> new_permutation(nb_qd_pts);

  // Apply all flips first
  for (Uint i = 0; i < permutation_code.nb_flips(); ++i)
  {
    old_permutation.swap(permutation_vec);

    // Fill flip permutation:
    new_permutation.resize(nb_qd_pts);
    for (Uint q = 0; q < nb_qd_pts; ++q)
    {
      new_permutation[q] = nb_qd_pts - q - 1;
    }

    for (Uint j = 0; j < new_permutation.size(); ++j)
    {
      permutation_vec[j] = old_permutation[new_permutation[j]];
    }
  }

  // Then apply all rotations
  for (Uint i = 0; i < permutation_code.nb_rotations(); ++i)
  {
    old_permutation.swap(permutation_vec);

    // Fill rotation permutation:
    new_permutation.resize(nb_qd_pts);
    for (Uint q = 0; q < nb_qd_pts; ++q)
    {
      new_permutation[q] = nb_qd_pts - q - 1;
    }

    for (Uint j = 0; j < new_permutation.size(); ++j)
    {
      permutation_vec[j] = old_permutation[new_permutation[j]];
    }
  }
}

} // namespace detail

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
