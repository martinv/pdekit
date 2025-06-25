#ifndef PDEKIT_Solver_FE_Helmholtz_Solver_CG_HDG_hpp
#define PDEKIT_Solver_FE_Helmholtz_Solver_CG_HDG_hpp

/// Standard template library headers
#include <array>
#include <cmath>
#include <ctime>
#include <forward_list>
#include <iostream>
#include <map>

/// PDEKIT headers
#include "common/MPI/MPIEnv.hpp"
#include "common/StringUtils.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "math/MathConstants.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "solver/fe/InteriorSolverCGHDG.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class HelmholtzSolverCGHDG
{
  public:
  /// TYPEDEFS
  using tria_t    = typename mesh::Tria<MeshConfig>;
  using dof_map_t = typename pdekit::result_of::dof_map_t<MeshConfig>;

  /// Default constructor
  HelmholtzSolverCGHDG();

  /// Destructor
  ~HelmholtzSolverCGHDG() = default;

  void setup(const tria_t &tria, const dof_map_t &dofs, const std::vector<Uint> &cell_to_part_ids);

  template <typename RHS>
  void assemble(const RHS &rhs);

  template <typename DirichletBC>
  void add_weak_bc(const DirichletBC &dirichlet_bc);

  void interior_solve();

  private:
  enum Dimensions : Uint
  {
    FACET = MeshConfig::TDIM - 1,
    CELL  = MeshConfig::TDIM
  };

  using interior_solver_t = InteriorSolverCGHDG<MeshConfig>;

  // Multi-edge key is a pair of partition indices [left, right]
  // All edges adjacent to cell such that the left incident cell has partition
  // index 'left' and the right adjacent cell has partition index 'right'
  // are stored in a multi-edge with the same key
  using multi_edge_key_t = std::tuple<Uint, Uint, Uint>;

  struct MultiEdgeKeyHash
  {
    std::size_t operator()(const multi_edge_key_t &tuple) const noexcept
    {
      return std::get<0>(tuple) && (std::get<1>(tuple) << 1) && (std::get<2>(tuple) << 2);
    }
  };

  using multi_edge_t = mesh::internal::TriaFacets<MeshConfig>;
  using multi_edge_map_t =
      std::unordered_map<multi_edge_key_t, std::shared_ptr<multi_edge_t>, MultiEdgeKeyHash>;

  using cell_buffer_t = mesh::CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM>;

  using cell_geo_metric_type =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             MeshConfig::GDIM>::cellwise_metric;
  using facet_geo_metric_type =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             MeshConfig::GDIM - 1>::cellwise_metric;

  void detect_internal_multiedges(const tria_t &tria_in, const std::vector<Uint> &cell_to_part_ids,
                                  multi_edge_map_t &multi_edges);

  void detect_boundary_multiedges(const tria_t &tria_in, const std::vector<Uint> &cell_to_part_ids,
                                  multi_edge_map_t &multi_edges);

  void write_multi_edge_to_gmsh(const tria_t &tria, const multi_edge_t &m_edge,
                                const std::string &filename, const std::string &gmsh_phys_id);

  multi_edge_map_t m_multi_edges;
  std::vector<std::vector<std::string>> m_interior_solver_boundaries;
  std::vector<std::unique_ptr<interior_solver_t>> m_interior_solvers;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
HelmholtzSolverCGHDG<MeshConfig>::HelmholtzSolverCGHDG()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void HelmholtzSolverCGHDG<MeshConfig>::setup(const tria_t &tria, const dof_map_t &dofs,
                                             const std::vector<Uint> &cell_to_part_ids)
{
  // -----------------------
  // Detect all multi-edges
  // -----------------------

  // Find multi-edges in mesh interior
  detect_internal_multiedges(tria, cell_to_part_ids, m_multi_edges);
  // Find multi-edges on mesh boundary
  detect_boundary_multiedges(tria, cell_to_part_ids, m_multi_edges);

  /*
  std::cout << "Total number of multi-edges: " << multi_edges.size() << std::endl;

  for (const auto &entry : multi_edges)
  {
    const auto key = std::get<0>(entry);
    std::cout << "Multi edge [" << std::get<0>(key) << "," << std::get<1>(key) << "("
              << std::get<2>(key) << ")]" << std::endl;
    const multi_edge_t &m_edge = *std::get<1>(entry);

    // std::cout << m_edge << std::endl;

    const std::string edge_uid = common::StringUtils::to_string(std::get<0>(key)) +
                                 common::StringUtils::to_string(std::get<1>(key)) +
                                 common::StringUtils::to_string(std::get<2>(key));

    const std::string edge_name = "MEdge_" + edge_uid + ".msh";
    write_multi_edge_to_gmsh(tria, m_edge, edge_name, edge_uid);
  }
  */

  // ------------------------------------------
  // Count the number of elements in each patch
  // ------------------------------------------
  std::unordered_map<Uint, Uint> patch_sizes;

  for (Uint ac = 0; ac < cell_to_part_ids.size(); ++ac)
  {
    patch_sizes[cell_to_part_ids[ac]]++;
  }

  const Uint nb_parts = patch_sizes.size();

  std::cout << "TOTAL " << nb_parts << " PATCHES" << std::endl;
  for (const auto &entry : patch_sizes)
  {
    std::cout << "Partition " << std::get<0>(entry) << " has " << std::get<1>(entry) << " elements"
              << std::endl;
  }

  std::vector<std::vector<mesh::ActiveIdx>> partition_to_cell_id(nb_parts);
  for (Uint p = 0; p < nb_parts; ++p)
  {
    std::vector<mesh::ActiveIdx> &cell_ids_in_part = partition_to_cell_id[p];
    const Uint nb_cells_in_part                    = patch_sizes[p];
    cell_ids_in_part.reserve(nb_cells_in_part);
  }

  for (Uint ac = 0; ac < cell_to_part_ids.size(); ++ac)
  {
    partition_to_cell_id[cell_to_part_ids[ac]].push_back(mesh::ActiveIdx(ac));
  }

  /*
  for (Uint p = 0; p < nb_parts; ++p)
  {
    std::cout << "PART " << p << std::endl;
    for (const auto val : partition_to_cell_id[p])
    {
      std::cout << val << " ";
    }
    std::cout << std::endl;
  }
  */

  m_interior_solvers.clear();
  m_interior_solvers.resize(0);
  m_interior_solver_boundaries.clear();
  m_interior_solver_boundaries.resize(nb_parts);

  for (Uint p = 0; p < nb_parts; ++p)
  {
    std::unique_ptr<interior_solver_t> int_solver_ptr(new interior_solver_t());
    int_solver_ptr->setup(tria, dofs, partition_to_cell_id[p]);

    for (const auto &entry : m_multi_edges)
    {
      const auto key             = std::get<0>(entry);
      const multi_edge_t &m_edge = *std::get<1>(entry);

      const Uint part_L = std::get<0>(key);
      const Uint part_R = std::get<1>(key);

#if 1
      // If this is a boundary multiedge or multiedge such that the partition
      // lies on its left side ...
      if (part_L == p)
      {
        const std::string edge_uid = common::StringUtils::to_string(part_L) +
                                     common::StringUtils::to_string(part_R) +
                                     common::StringUtils::to_string(std::get<2>(key));

        const std::string edge_name = "MEdge_" + edge_uid + ".msh";
        int_solver_ptr->add_boundary(tria, m_edge, LeftRightOrientation::LEFT, edge_name);
        m_interior_solver_boundaries[p].push_back(edge_name);
      }
      else if ((part_L != p) && (part_R == p))
      {
        const std::string edge_uid = common::StringUtils::to_string(part_L) +
                                     common::StringUtils::to_string(part_R) +
                                     common::StringUtils::to_string(std::get<2>(key));

        const std::string edge_name = "MEdge_" + edge_uid + ".msh";
        int_solver_ptr->add_boundary(tria, m_edge, LeftRightOrientation::RIGHT, edge_name);
        m_interior_solver_boundaries[p].push_back(edge_name);
      }
#endif
    } // Loop over all multi-edges

    m_interior_solvers.push_back(std::move(int_solver_ptr));
  }

  std::cout << "Boundary names for interior solvers:" << std::endl;
  for (Uint p = 0; p < nb_parts; ++p)
  {
    std::cout << "PATCH " << p << std::endl;
    for (const auto &name : m_interior_solver_boundaries[p])
    {
      std::cout << "    " << name << std::endl;
    }

    m_interior_solvers[p]->write_to_gmsh("interior_patch_" + common::StringUtils::to_string(p) +
                                         ".msh");
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename RHS>
void HelmholtzSolverCGHDG<MeshConfig>::assemble(const RHS &rhs)
{
  for (auto &int_solver : m_interior_solvers)
  {
    int_solver->assemble(rhs);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename DirichletBC>
void HelmholtzSolverCGHDG<MeshConfig>::add_weak_bc(const DirichletBC &dirichlet_bc)
{
  for (Uint is = 0; is < m_interior_solvers.size(); ++is)
  {
    const auto &bc_names = m_interior_solver_boundaries[is];
    for (const auto &name : bc_names)
    {
      m_interior_solvers[is]->add_weak_bc(name, dirichlet_bc);
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void HelmholtzSolverCGHDG<MeshConfig>::interior_solve()
{
  for (auto &is : m_interior_solvers)
  {
    is->solve();
  }

  for (Uint is = 0; is < m_interior_solvers.size(); ++is)
  {
    const std::string filename = "interior_solve_00" + common::StringUtils::to_string(is) + ".msh";
    m_interior_solvers[is]->write_to_gmsh(filename);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void HelmholtzSolverCGHDG<MeshConfig>::detect_internal_multiedges(
    const tria_t &tria, const std::vector<Uint> &cell_to_part_ids, multi_edge_map_t &multi_edges)
{
  std::array<mesh::IncidenceEntry, 2> inc_entry_buffer;
  std::array<mesh::EntityDofRealign, 2> perm_buffer;

  // Get all boundaries
  mesh::MeshBoundarySet<MeshConfig> const &boundaries = tria.all_boundaries();
  const Uint num_boundaries                           = boundaries.nb_domains();

  for (typename tria_t::const_skeleton_iterator facet_it = tria.cbegin_skeleton(Dimensions::FACET);
       facet_it != tria.cend_skeleton(Dimensions::FACET); ++facet_it)
  {
    const auto facet_incidences = facet_it->incidences();

    if (facet_incidences.size() == 2)
    {
      const mesh::CellTopologyView<MeshConfig> tcell_L =
          tria.cell(mesh::FlatIdx(facet_incidences.cell_id(0)));
      const mesh::CellTopologyView<MeshConfig> tcell_R =
          tria.cell(mesh::FlatIdx(facet_incidences.cell_id(1)));

      const Uint part_id_L = cell_to_part_ids[tcell_L.active_idx().id()];
      const Uint part_id_R = cell_to_part_ids[tcell_R.active_idx().id()];

      // If this is an edge between two different partitions ...
      if (part_id_L != part_id_R)
      {
        const bool flip_facet = part_id_L > part_id_R;
        const multi_edge_key_t key =
            flip_facet ? multi_edge_key_t(part_id_R, part_id_L, num_boundaries + 1)
                       : multi_edge_key_t(part_id_L, part_id_R, num_boundaries + 1);

        typename multi_edge_map_t::iterator it = multi_edges.find(key);

        if (it == multi_edges.end())
        {
          multi_edges[key] = std::make_shared<multi_edge_t>();
          multi_edges[key]->set_dim(Dimensions::FACET);
          it = multi_edges.find(key);
        }
        multi_edge_t &m_edge = *std::get<1>(*it);

        if (!flip_facet)
        {
          m_edge.add_facet(facet_incidences.incidences(), facet_incidences.permutations());
        }
        else
        {
          inc_entry_buffer[0] = facet_incidences.incidence_entry(1);
          inc_entry_buffer[1] = facet_incidences.incidence_entry(0);

          perm_buffer[0] = facet_incidences.permutation(1);
          perm_buffer[1] = facet_incidences.permutation(0);

          const common::ArrayView<const mesh::IncidenceEntry, _1D, Uint> incidences(
              inc_entry_buffer.data(), inc_entry_buffer.size());
          const common::ArrayView<const mesh::EntityDofRealign, _1D, Uint> facet_perms(
              perm_buffer.data(), perm_buffer.size());
          m_edge.add_facet(incidences, facet_perms);
        }
      } // if part_id_L != part_id_R

    } // if facet_incidences.size() == 2

  } // Loop over mesh skeleton

  for (const auto &entry : multi_edges)
  {
    multi_edge_t &m_edge = *std::get<1>(entry);
    m_edge.set_status(mesh::EntityStatus::Active);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void HelmholtzSolverCGHDG<MeshConfig>::detect_boundary_multiedges(
    const tria_t &tria, const std::vector<Uint> &cell_to_part_ids, multi_edge_map_t &multi_edges)
{
  // Get all boundaries
  mesh::MeshBoundarySet<MeshConfig> const &boundaries = tria.all_boundaries();
  const Uint num_boundaries                           = boundaries.nb_domains();

  const std::vector<typename mesh::MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr>
      &all_boundaries = boundaries.all_domains();

  std::vector<std::string> boundary_ids(num_boundaries);
  Uint bdry_idx = 0;

  for (const auto &bdry_ptr : all_boundaries)
  {
    boundary_ids[bdry_idx++] = bdry_ptr->name();
  }

  bdry_idx = 0;

  for (const auto &bdry_ptr : all_boundaries)
  {
    const mesh::BoundaryFacets<MeshConfig, Dimensions::FACET> &bdry = *bdry_ptr;

    for (auto bdry_facet_it = bdry.cbegin(); bdry_facet_it != bdry.cend(); ++bdry_facet_it)
    {
      const mesh::TraceIncidences bdry_incidences = bdry_facet_it->incidences();
      const mesh::FlatIdx adj_tcell_idx(bdry_incidences.cell_id(0));
      const mesh::CellTopologyView<MeshConfig> adjacent_tcell_view = tria.cell(adj_tcell_idx);

      const Uint part_id_L = cell_to_part_ids[adjacent_tcell_view.active_idx().id()];
      const multi_edge_key_t key(part_id_L, part_id_L, bdry_idx);

      typename multi_edge_map_t::iterator it = multi_edges.find(key);

      if (it == multi_edges.end())
      {
        multi_edges[key] = std::make_shared<multi_edge_t>();
        multi_edges[key]->set_dim(Dimensions::FACET);
        it = multi_edges.find(key);
      }
      multi_edge_t &m_edge = *std::get<1>(*it);
      m_edge.add_facet(bdry_incidences.incidences(), bdry_incidences.permutations());
    }
    bdry_idx++;
  }

  for (const auto &entry : multi_edges)
  {
    multi_edge_t &m_edge = *std::get<1>(entry);
    m_edge.set_status(mesh::EntityStatus::Active);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void HelmholtzSolverCGHDG<MeshConfig>::write_multi_edge_to_gmsh(const tria_t &tria,
                                                                const multi_edge_t &m_edge,
                                                                const std::string &filename,
                                                                const std::string &gmsh_phys_id)
{
  std::ofstream out_stream;
  out_stream.open(filename.c_str());

  const auto it_end = m_edge.cend();

  Uint nb_nodes  = 0;
  nb_nodes       = 0;
  Uint nb_facets = 0;

  for (auto it = m_edge.cbegin(); it != it_end; ++it)
  {
    const mesh::TraceIncidences trace_inc = it->incidences();
    const mesh::FlatIdx cell_idx(trace_inc.cell_id(0));
    const Uint local_id = trace_inc.local_id(0);

    const mesh::CellTopologyView<MeshConfig> tcell_view = tria.cell(cell_idx);

    const auto facet_coord = tcell_view.coordinates(Dimensions::FACET, local_id);
    nb_nodes += facet_coord.size();
    nb_facets += facet_coord.size() - 1;
  }

  out_stream << "$MeshFormat\n2.2 0 8\n$EndMeshFormat" << std::endl;
  out_stream << "$Nodes\n" << nb_nodes << std::endl;

  nb_nodes = 0;

  for (auto it = m_edge.cbegin(); it != it_end; ++it)
  {
    const mesh::TraceIncidences trace_inc = it->incidences();
    const mesh::FlatIdx cell_idx(trace_inc.cell_id(0));
    const Uint local_id = trace_inc.local_id(0);

    const mesh::CellTopologyView<MeshConfig> tcell_view = tria.cell(cell_idx);

    const auto facet_coord = tcell_view.coordinates(Dimensions::FACET, local_id);
    for (Uint n = 0; n < facet_coord.size(); ++n)
    {
      out_stream << nb_nodes + 1 << " " << facet_coord.const_node_view(n);

      if (facet_coord.dim() < _3D)
      {
        out_stream << " 0.0";
      }
      out_stream << std::endl;
      nb_nodes++;
    }
  }

  nb_nodes = 0;
  out_stream << "$EndNodes\n$Elements\n" << nb_facets << std::endl;
  nb_facets = 0;

  for (auto it = m_edge.cbegin(); it != it_end; ++it)
  {
    const mesh::TraceIncidences trace_inc = it->incidences();
    const mesh::FlatIdx cell_idx(trace_inc.cell_id(0));
    const Uint local_id = trace_inc.local_id(0);

    const mesh::CellTopologyView<MeshConfig> tcell_view = tria.cell(cell_idx);

    const auto facet_coord = tcell_view.coordinates(Dimensions::FACET, local_id);
    for (Uint n = 0; (n + 1) < facet_coord.size(); ++n)
    {
      out_stream << nb_facets + 1 << " 1 2 " << gmsh_phys_id << " " << gmsh_phys_id << " "
                 << nb_nodes + 1 << " " << nb_nodes + 2 << std::endl;
      nb_facets++;
      nb_nodes++;
    }
    nb_nodes++;
  }

  out_stream << "$EndElements" << std::endl;

  out_stream.close();
}
// ----------------------------------------------------------------------------

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
