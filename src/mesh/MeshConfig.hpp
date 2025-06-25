#ifndef PDEKIT_Mesh_Mesh_Config_hpp
#define PDEKIT_Mesh_Mesh_Config_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace mesh
{

template <typename MeshConfig>
class Tria;
template <typename MeshConfig>
class MeshBoundarySet;
template <typename MeshConfig>
class DofMap;

// ----------------------------------------------------------------------------

// Configuration for 2D Cartesian mesh
struct Cart2D
{
  // Geometrical dimension of this mesh (number of coordinates)
  enum
  {
    GDIM = DIM_2D
  };

  // Topological dimension of this mesh (surface mesh has TDIM = 2,
  //                                     volume mesh has TDIM = 3)
  enum
  {
    TDIM = DIM_2D
  };

  // Topological dimension of edges in Finite Element mesh
  enum
  {
    EDGE_DIM = DIM_1D
  };

  // Topological dimension of faces in Finite Element mesh
  enum
  {
    FACET_DIM = DIM_1D
  };
};

// ----------------------------------------------------------------------------

// Configuration for 3D Cartesian mesh
struct Cart3D
{
  // Geometrical dimension of this mesh (number of coordinates)
  enum
  {
    GDIM = DIM_3D
  };
  // Topological dimension of this mesh (surface mesh has TDIM = 2,
  //                                     volume mesh has TDIM = 3)
  enum
  {
    TDIM = DIM_3D
  };

  // Topological dimension of edges in Finite Element mesh
  enum
  {
    EDGE_DIM = DIM_1D
  };

  // Topological dimension of faces in Finite Element mesh
  enum
  {
    FACET_DIM = DIM_2D
  };
};

// ----------------------------------------------------------------------------

// Configuration for 2D surface mesh in 3D space
struct Surf2D3D
{
  // Geometrical dimension of this mesh (number of coordinates)
  enum
  {
    GDIM = DIM_3D
  };
  // Topological dimension of this mesh (surface mesh has TDIM = 2,
  //                                     volume mesh has TDIM = 3)
  enum
  {
    TDIM = DIM_2D
  };

  // Topological dimension of edges in Finite Element mesh
  enum
  {
    EDGE_DIM = DIM_1D
  };

  // Topological dimension of faces in Finite Element mesh
  enum
  {
    FACET_DIM = DIM_1D
  };
};

// ----------------------------------------------------------------------------

// Classes to label different types of entities in mesh

struct PointEntity
{
};

struct EdgeEntity
{
};

struct FacetEntity
{
};

struct VolumeEntity
{
};

// ----------------------------------------------------------------------------

} // namespace mesh

namespace result_of
{

// ----------------------
// Triangulation type
// ----------------------

template <typename MeshConfig>
using tria_t = mesh::Tria<MeshConfig>;

// ----------------------
// Mesh boundary set type
// ----------------------

template <typename MeshConfig>
using mesh_boundary_set_t = mesh::MeshBoundarySet<MeshConfig>;

// -----------------
// DOF storage type
// -----------------

template <typename MeshConfig>
using dof_map_t = mesh::DofMap<MeshConfig>;

} // Namespace result_of

// ----------------------------------------------------------------------------

} // namespace pdekit

#endif
