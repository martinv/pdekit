#ifndef PDEKIT_Mesh_Mesh_Creator_hpp
#define PDEKIT_Mesh_Mesh_Creator_hpp

#include <vector>

#include "mesh/CellBuffer.hpp"
#include "mesh/Tria.hpp"

namespace pdekit
{

namespace mesh
{

class MeshCreator
{

  public:
  /// Default constructor
  MeshCreator();

  /// Destructor
  ~MeshCreator();

  /// Generate a unit triangle
  template <typename MeshConfig>
  static void make_unit_triangle(Tria<MeshConfig> &mesh, const std::string &geo_dofs_name,
                                 const Uint N);

  /// Generate a unit quadrilateral
  template <typename MeshConfig>
  static void make_unit_quad(Tria<MeshConfig> &mesh, const std::string &geo_dofs_name, const Uint N,
                             const bool split_into_triags = false);

  private:
  /// Generate data needed for creation of a mesh on unit triangle

  /// Generate data needed for creation of a mesh on unit triangle
  static void generate_unit_triangle_coords(std::vector<Real> &coords, const Uint N);

  /// Generate connectivities for creation of a mesh on unit triangle
  template <typename MeshConfig>
  static void generate_unit_triangle_connectivity(
      Tria<MeshConfig> &interior_cells, typename result_of::dof_map_t<MeshConfig> &cell_dofs,
      MeshBoundarySet<MeshConfig> &boundaries, const Uint N);

  /// Generate data needed for creation of a mesh on unit quad

  static void generate_unit_quad_coords(std::vector<Real> &coords, const Uint N);

  /// Generate connectivities for creation of a mesh on unit quad. The 2D
  /// elements in
  /// connectivities are QUADS
  template <typename MeshConfig>
  static void generate_unit_quad_quad_connectivity(
      Tria<MeshConfig> &interior_cells, typename result_of::dof_map_t<MeshConfig> &cell_dofs,
      MeshBoundarySet<MeshConfig> &boundaries, const Uint N);

  /// Generate connectivities for creation of a mesh on unit quad. The 2D
  /// elements in
  /// connectivities are TRIANGLES
  template <typename MeshConfig>
  static void generate_unit_quad_triag_connectivity(
      Tria<MeshConfig> &interior_cells, typename result_of::dof_map_t<MeshConfig> &cell_dofs,
      MeshBoundarySet<MeshConfig> &boundaries, const Uint N);
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshCreator::make_unit_triangle(Tria<MeshConfig> &mesh, const std::string &geo_dofs_name,
                                     const Uint N)
{
  mesh.create_dof_storage(geo_dofs_name);

  generate_unit_triangle_connectivity(mesh, *mesh.dof_storage(geo_dofs_name), mesh.all_boundaries(),
                                      N);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshCreator::make_unit_quad(Tria<MeshConfig> &mesh, const std::string &geo_dofs_name,
                                 const Uint N, const bool split_into_triags)
{
  mesh.create_dof_storage(geo_dofs_name);

  if (!split_into_triags)
  {
    generate_unit_quad_quad_connectivity(mesh, *mesh.dof_storage(geo_dofs_name),
                                         mesh.all_boundaries(), N);
  }
  else
  {
    generate_unit_quad_triag_connectivity(mesh, *mesh.dof_storage(geo_dofs_name),
                                          mesh.all_boundaries(), N);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshCreator::generate_unit_triangle_connectivity(
    Tria<MeshConfig> &interior_cells, typename result_of::dof_map_t<MeshConfig> &cell_dofs,
    MeshBoundarySet<MeshConfig> &boundaries, const Uint N)
{
  // Create coordinates
  std::vector<Real> coords;
  generate_unit_triangle_coords(coords, N);

  // One triangle will have 3 * 2 = 6 coordinate values
  std::vector<Real> one_cell_coords(3 * _2D);

  // -----------------------------------------------------
  // N is the number of points on one edge of the triangle
  // -----------------------------------------------------
  CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> buffer;
  const PointSetTag elem_type = PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist);

  const Uint nb_triag = (N - 1) * (N - 1);

  // Connectivity of one triangle
  std::vector<Uint> connections(3);
  buffer.reserve(3 * nb_triag, nb_triag);
  Uint triag_idx = 0;

  /*
  //  ^    |  \
  //  |    |   \                  \
  //  |    |    \|    \            \
  //  |    N----N+1---N+2---     --2N-2
  //  |    |\    |\    |\            \
  //  | j  | \ b | \ d | \ f |        \
  //  |    |  \  |  \  |  \  |         \
  //  |    | a \ | c \ | e \ |          \
  //  |    |    \|    \|    \|           \
  //  |    0-----1-----2-----3     ---   N-1
  //  |
  //  |          i
  //  x----------------------->
  */

  // For each pair of triangles that form a rectangle together,
  // lower_left_triag_node is the bottom left index of this rectangle
  // Example: for the rectangle formed by triags c and d,
  // lower_left_triag_node = 1
  Uint lower_left_triag_node = 0;

  for (Uint j = 0; j < (N - 1); ++j)
  {
    for (Uint i = 0; i < (N - j - 1); ++i)
    {
      const Uint a = lower_left_triag_node;
      const Uint b = a + 1;
      const Uint c = lower_left_triag_node + (N - j);
      const Uint d = c + 1;

      // std::cout << "Lower left triag idx = " << lower_left_triag_node
      //          << " (i,j) = (" << i << "," << j << ")" << std::endl;
      // std::cout << "  [" << a << "," << b << "," << c << "]" <<
      // std::endl;
      connections[0] = a;
      connections[1] = b;
      connections[2] = c;
      // std::cout << "[" << connections[3*triag_idx] << "," <<
      // connections[3*triag_idx+1]
      //          << "," << connections[3*triag_idx+2] << "]" <<
      //          std::endl;

      // Fill the coordinates of one element
      for (Uint n = 0; n < connections.size(); ++n)
      {
        one_cell_coords[_2D * n]     = coords[_2D * connections[n]];
        one_cell_coords[_2D * n + 1] = coords[_2D * connections[n] + 1];
      }

      buffer.push_back_cell(triag_idx, elem_type, connections, 0u, one_cell_coords);
      triag_idx++;

      if (i < (N - j - 2))
      {
        // std::cout << "  [" << b << "," << d << "," << c << "]" <<
        // std::endl;
        connections[0] = b;
        connections[1] = d;
        connections[2] = c;
        // std::cout << "[" << connections[3*triag_idx] << "," <<
        // connections[3*triag_idx+1]
        //          << "," << connections[3*triag_idx+2] << "]" <<
        //          std::endl;

        // Fill the coordinates of one element
        for (Uint n = 0; n < connections.size(); ++n)
        {
          one_cell_coords[_2D * n]     = coords[_2D * connections[n]];
          one_cell_coords[_2D * n + 1] = coords[_2D * connections[n] + 1];
        }

        buffer.push_back_cell(triag_idx, elem_type, connections, 0u, one_cell_coords);
        triag_idx++;
      }
      else
      {
        lower_left_triag_node++;
      }
      lower_left_triag_node++;
    }
  }

  // Create the topology
  interior_cells.create_from_cells(buffer);

  // Create the dof storage
  cell_dofs.create_from_cells(buffer);

  // Tag all cells to create a domain
  std::vector<Uint> domain_tag(nb_triag, 1U);

  std::vector<std::pair<Uint, std::string>> cell_tag_names;
  cell_tag_names.push_back(std::pair<Uint, std::string>(1U, "Interior"));

  cell_dofs.tag_all_active_cells(cell_tag_names, domain_tag);

  // ******
  // Edge 0
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge0 =
      boundaries.create(_1D, "Edge0");

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data0(new std::vector<IncidenceEntry>());
  cell_id_data0->resize(N - 1);
  for (Uint j = 0; j < cell_id_data0->size(); ++j)
  {
    (*cell_id_data0)[j].cell_idx = 2 * j;
    (*cell_id_data0)[j].local_id = 0u;
  }

  (*edge0).emplace_bdry_cell_ids(std::move(cell_id_data0), 2U);

  // ******
  // Edge 1
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge1 =
      boundaries.create(_1D, "Edge1");

  triag_idx = 2 * (N - 2);

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data1(new std::vector<IncidenceEntry>());
  cell_id_data1->resize(N - 1);
  for (Uint j = 0; j < cell_id_data1->size(); ++j)
  {
    (*cell_id_data1)[j].cell_idx = triag_idx;
    (*cell_id_data1)[j].local_id = 1u;

    triag_idx += 2 * (N - (j + 1) - 2) + 1;
  }

  (*edge1).emplace_bdry_cell_ids(std::move(cell_id_data1), 3U);

  // ******
  // Edge 2
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge2 =
      boundaries.create(_1D, "Edge2");

  triag_idx = 0;

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data2(new std::vector<IncidenceEntry>());
  cell_id_data2->resize(N - 1);
  for (Uint j = 0; j < cell_id_data2->size(); ++j)
  {
    (*cell_id_data2)[j].cell_idx = triag_idx;
    (*cell_id_data2)[j].local_id = 2u;

    triag_idx += 2 * (N - j - 2) + 1;
  }

  (*edge2).emplace_bdry_cell_ids(std::move(cell_id_data2), 4U);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshCreator::generate_unit_quad_quad_connectivity(
    Tria<MeshConfig> &interior_cells, typename result_of::dof_map_t<MeshConfig> &cell_dofs,
    MeshBoundarySet<MeshConfig> &boundaries, const Uint N)
{
  // Create the coordinates
  std::vector<Real> coords;
  generate_unit_quad_coords(coords, N);

  // One quadrilateral will have 4 * 2 = 8 coordinate values
  std::vector<Real> one_cell_coords(4 * _2D);

  // ----------------------------------------------------------
  // N is the number of points on one edge of the quadrilateral
  // ----------------------------------------------------------
  CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> buffer;
  const PointSetTag elem_type = PointSetTag(ElemShape::Quad, P1, PointSetID::Equidist);

  const Uint nb_quad = (N - 1) * (N - 1);

  std::vector<Uint> connections(4);

  buffer.reserve(4 * nb_quad, nb_quad);

  Uint quad_idx = 0;

  /*
  //  ^    |     .     .
  //  |    |     .     .             |       |
  //  |    |     |     |             |       |
  //  |    N----N+1---N+2---   ... -2N-2----2N-1
  //  |    |     |     |             |       |
  //  | j  |     |     |     |       |       |
  //  |    |  a  |  b  |  c  |       |       |
  //  |    |     |     |     |       |       |
  //  |    |     |     |     |               |
  //  |    0-----1-----2-----3  ...  N-2----N-1
  //  |
  //  |          i
  //  x----------------------->
  */

  // lower_left_quad_node is the bottom left index of each rectangle
  // Example: for a rectangle 'c', lower_left_quad_node = 2
  Uint lower_left_quad_node = 0;

  for (Uint j = 0; j < (N - 1); ++j)
  {
    for (Uint i = 0; i < (N - 1); ++i)
    {
      connections[0] = lower_left_quad_node;
      connections[1] = lower_left_quad_node + 1;
      connections[2] = lower_left_quad_node + N + 1;
      connections[3] = lower_left_quad_node + N;

      // std::cout << "[" << connections[4*quad_idx]   << "," <<
      // connections[4*quad_idx+1]
      //          << "," << connections[4*quad_idx+2] << "," <<
      // connections[4*quad_idx+3] << "]" << std::endl;

      // Fill the coordinates of one element
      for (Uint n = 0; n < connections.size(); ++n)
      {
        one_cell_coords[_2D * n]     = coords[_2D * connections[n]];
        one_cell_coords[_2D * n + 1] = coords[_2D * connections[n] + 1];
      }

      buffer.push_back_cell(quad_idx, elem_type, connections, 0u, one_cell_coords);
      quad_idx++;
      lower_left_quad_node++;

      // If we are at the end of one row, we move one row up:
      if (i == (N - 2))
      {
        lower_left_quad_node++;
      }
    }
  }

  // Create the topology
  interior_cells.create_from_cells(buffer);

  // Create the dof storage
  cell_dofs.create_from_cells(buffer);

  // Tag all cells to create a domain
  std::vector<Uint> domain_tag(nb_quad, 1U);

  std::vector<std::pair<Uint, std::string>> cell_tag_names;
  cell_tag_names.push_back(std::pair<Uint, std::string>(1U, "Interior"));

  cell_dofs.tag_all_active_cells(cell_tag_names, domain_tag);

  // ******
  // Edge 0
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge0 =
      boundaries.create(_1D, "Edge0");

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data0(new std::vector<IncidenceEntry>());
  cell_id_data0->resize(N - 1);
  for (Uint i = 0; i < cell_id_data0->size(); ++i)
  {
    (*cell_id_data0)[i].cell_idx = i;
    (*cell_id_data0)[i].local_id = 0u;
  }

  (*edge0).emplace_bdry_cell_ids(std::move(cell_id_data0), 2U);

  // ******
  // Edge 1
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge1 =
      boundaries.create(_1D, "Edge1");

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data1(new std::vector<IncidenceEntry>());
  cell_id_data1->resize(N - 1);
  for (Uint i = 0; i < cell_id_data1->size(); ++i)
  {
    (*cell_id_data1)[i].cell_idx = (i + 1) * (N - 1) - 1;
    (*cell_id_data1)[i].local_id = 1u;
  }

  (*edge1).emplace_bdry_cell_ids(std::move(cell_id_data1), 3U);

  // ******
  // Edge 2
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge2 =
      boundaries.create(_1D, "Edge2");

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data2(new std::vector<IncidenceEntry>());
  cell_id_data2->resize(N - 1);
  for (Uint i = 0; i < cell_id_data2->size(); ++i)
  {
    (*cell_id_data2)[i].cell_idx = (N - 1) * (N - 1) - 1 - i;
    (*cell_id_data2)[i].local_id = 2u;
  }

  (*edge2).emplace_bdry_cell_ids(std::move(cell_id_data2), 4U);

  // ******
  // Edge 3
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge3 =
      boundaries.create(_1D, "Edge3");

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data3(new std::vector<IncidenceEntry>());
  cell_id_data3->resize(N - 1);
  for (Uint i = 0; i < cell_id_data3->size(); ++i)
  {
    (*cell_id_data3)[i].cell_idx = (N - 1) * (N - 2 - i);
    (*cell_id_data3)[i].local_id = 3u;
  }

  (*edge3).emplace_bdry_cell_ids(std::move(cell_id_data3), 5U);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshCreator::generate_unit_quad_triag_connectivity(
    Tria<MeshConfig> &interior_cells, typename result_of::dof_map_t<MeshConfig> &cell_dofs,
    MeshBoundarySet<MeshConfig> &boundaries, const Uint N)
{
  // Create coordinates
  std::vector<Real> coords;
  generate_unit_quad_coords(coords, N);

  // One triangle will have 3 * 2 = 6 coordinate values
  std::vector<Real> one_cell_coords(3 * _2D);

  // ----------------------------------------------------------
  // N is the number of points on one edge of the quadrilateral
  // ----------------------------------------------------------
  CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> buffer;
  const PointSetTag elem_type = PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist);

  const Uint nb_triag = 2 * (N - 1) * (N - 1);

  // Connectivity of one triangle
  std::vector<Uint> connections(3);
  buffer.reserve(3 * nb_triag, nb_triag);
  Uint triag_idx = 0;

  /*
  //  ^    |     .     .
  //  |    |     .     .             |     |
  //  |    |     |     |             |     |
  //  |    N----N+1---N+2---   ... -2N-2---2N-1
  //  |    |    /|    /|    /        |    /|
  //  | j  | b / | d / | f / |       |   / |
  //  |    |  /  |  /  |  /  |       |  /  |
  //  |    | / a | / c | / e | g     | /   |
  //  |    |/    |/    |/    |       |/    |
  //  |    0-----1-----2-----3  ...  N-2----N-1
  //  |
  //  |          i
  //  x----------------------->
  */

  // We have triangle pairs forming quad: (a,b), (c,d), (e,f) etc.
  // For each pair of triangles that form a rectangle together,
  // lower_left_triag_node is the bottom left index of each rectangle
  // Example: for the triangle 'c', lower_left_quad_node = 1
  Uint lower_left_triag_node = 0;

  for (Uint j = 0; j < (N - 1); ++j)
  {
    for (Uint i = 0; i < (N - 1); ++i)
    {
      // First triangle in each triangle pair
      connections[0] = lower_left_triag_node;
      connections[1] = lower_left_triag_node + 1;
      connections[2] = lower_left_triag_node + N + 1;

      // std::cout << "[" << connections[3*triag_idx]   << "," <<
      // connections[3*triag_idx+1]
      //          << "," << connections[3*triag_idx+2] << "]" <<
      //          std::endl;

      // Fill the coordinates of one element
      for (Uint n = 0; n < connections.size(); ++n)
      {
        one_cell_coords[_2D * n]     = coords[_2D * connections[n]];
        one_cell_coords[_2D * n + 1] = coords[_2D * connections[n] + 1];
      }

      buffer.push_back_cell(triag_idx, elem_type, connections, 0u, one_cell_coords);
      triag_idx++;

      // Second triangle in each triangle pair
      connections[0] = lower_left_triag_node + N + 1;
      connections[1] = lower_left_triag_node + N;
      connections[2] = lower_left_triag_node;

      // std::cout << "[" << connections[3*triag_idx]   << "," <<
      // connections[3*triag_idx+1]
      //          << "," << connections[3*triag_idx+2] << "]" <<
      //          std::endl;

      // Fill the coordinates of one element
      for (Uint n = 0; n < connections.size(); ++n)
      {
        one_cell_coords[_2D * n]     = coords[_2D * connections[n]];
        one_cell_coords[_2D * n + 1] = coords[_2D * connections[n] + 1];
      }

      buffer.push_back_cell(triag_idx, elem_type, connections, 0u, one_cell_coords);
      triag_idx++;

      lower_left_triag_node++;

      // If we are at the end of one row, we move one row up:
      if (i == (N - 2))
      {
        lower_left_triag_node++;
      }
    }
  }

  // Create the topology
  interior_cells.create_from_cells(buffer);

  // Create the dof storage
  cell_dofs.create_from_cells(buffer);

  // Tag all cells to create a domain
  std::vector<Uint> domain_tag(nb_triag, 1U);

  std::vector<std::pair<Uint, std::string>> cell_tag_names;
  cell_tag_names.push_back(std::pair<Uint, std::string>(1U, "Interior"));

  cell_dofs.tag_all_active_cells(cell_tag_names, domain_tag);

  // ******
  // Edge 0
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge0 =
      boundaries.create(_1D, "Edge0");

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data0(new std::vector<IncidenceEntry>());
  cell_id_data0->resize(N - 1);
  for (Uint i = 0; i < cell_id_data0->size(); ++i)
  {
    (*cell_id_data0)[i].cell_idx = 2 * i;
    (*cell_id_data0)[i].local_id = 0u;
  }

  (*edge0).emplace_bdry_cell_ids(std::move(cell_id_data0), 2U);

  // ******
  // Edge 1
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge1 =
      boundaries.create(_1D, "Edge1");

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data1(new std::vector<IncidenceEntry>());
  cell_id_data1->resize(N - 1);
  for (Uint i = 0; i < cell_id_data1->size(); ++i)
  {
    (*cell_id_data1)[i].cell_idx = 2 * (i + 1) * (N - 1) - 2;
    (*cell_id_data1)[i].local_id = 1u;
  }

  (*edge1).emplace_bdry_cell_ids(std::move(cell_id_data1), 3U);

  // ******
  // Edge 2
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge2 =
      boundaries.create(_1D, "Edge2");

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data2(new std::vector<IncidenceEntry>());
  cell_id_data2->resize(N - 1);
  for (Uint i = 0; i < cell_id_data2->size(); ++i)
  {
    (*cell_id_data2)[i].cell_idx = 2 * (N - 1) * (N - 1) - 1 - 2 * i;
    (*cell_id_data2)[i].local_id = 0u;
  }

  (*edge2).emplace_bdry_cell_ids(std::move(cell_id_data2), 4U);

  // ******
  // Edge 3
  typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr edge3 =
      boundaries.create(_1D, "Edge3");

  std::unique_ptr<std::vector<IncidenceEntry>> cell_id_data3(new std::vector<IncidenceEntry>());
  cell_id_data3->resize(N - 1);
  for (Uint i = 0; i < cell_id_data3->size(); ++i)
  {
    (*cell_id_data3)[i].cell_idx = 2 * (N - 1) * (N - 2 - i) + 1;
    (*cell_id_data3)[i].local_id = 1u;
  }

  (*edge3).emplace_bdry_cell_ids(std::move(cell_id_data3), 5U);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
