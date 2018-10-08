#include "mesh/io/MeshCreator.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

MeshCreator::MeshCreator()
{
}

// ----------------------------------------------------------------------------

MeshCreator::~MeshCreator()
{
}

// ----------------------------------------------------------------------------

void MeshCreator::generate_unit_triangle_coords(std::vector<Real> &coords, const Uint N)
{
  // We suppose the triangle is defined as T = { (x,y) | -1< = x < = 1 && -1<
  // = y < = 1 && x+y < = 0 } In other words, the limits in x and y are in the
  // interval (-1,1), not (0,1)

  const Uint tot_nb_pts = N * (N + 1) / 2;
  const Real delta      = 2.0 / (N - 1);

  coords.resize(tot_nb_pts * _2D);

  Uint pt_idx = 0;

  for (Uint j = 0; j < N; ++j)
  {
    for (Uint i = 0; i < (N - j); ++i)
    {
      const Real xi0 = -1. + i * delta;
      const Real xi1 = -1. + j * delta;

      coords[pt_idx * _2D + XI0] = xi0;
      coords[pt_idx * _2D + XI1] = xi1;
      pt_idx++;
    }
  }
}

// ----------------------------------------------------------------------------

void MeshCreator::generate_unit_quad_coords(std::vector<Real> &coords, const Uint N)
{
  // We suppose the quad is defined as Q = { (x,y) | -1< = x < = 1 && -1< = y
  // < = 1 } In other words, the limits in x and y are in the domain (-1,1) x
  // (-1,1)

  const Uint tot_nb_pts = N * N;
  const Real delta      = 2.0 / (N - 1);

  coords.resize(tot_nb_pts * _2D);
  // std::vector<Real> one_node_coord(_2D);

  Uint pt_idx = 0;

  for (Uint j = 0; j < N; ++j)
  {
    for (Uint i = 0; i < N; ++i)
    {
      const Real xi0 = -1. + i * delta;
      const Real xi1 = -1. + j * delta;

      coords[pt_idx * _2D + XI0] = xi0;
      coords[pt_idx * _2D + XI1] = xi1;

      pt_idx++;
    }
  }
  // std::cout << "Mesh geometry for unit quad" << std::endl;
  // std::cout << coords << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
