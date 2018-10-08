#ifndef PDEKIT_Mesh_Mesh_Manipulator_hpp
#define PDEKIT_Mesh_Mesh_Manipulator_hpp

#include <vector>

#include "math/DenseSVec.hpp"
#include "mesh/Tria.hpp"

namespace pdekit
{

namespace mesh
{

class MeshManipulator
{

  public:
  /// Default constructor
  MeshManipulator();

  /// Destructor
  ~MeshManipulator();

  /// Scale a mesh by given factor
  template <typename MeshConfig>
  static void scale_mesh(Tria<MeshConfig> &mesh, const std::string &dof_handler_name,
                         const Real factor);

  /// Move a mesh by given factor
  /// @param mesh ... mesh that should be moved
  /// @param vec  ... vector defining the translation
  template <typename MeshConfig>
  static void translate_mesh(Tria<MeshConfig> &mesh, const std::string &dof_handler_name,
                             const math::DenseSVec<Real, MeshConfig::GDIM> &vec);

  private:
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshManipulator::scale_mesh(Tria<MeshConfig> &mesh, const std::string &dof_handler_name,
                                 const Real factor)
{
#if 0
  auto scale_transf = [factor](const math::DenseConstVecView<Real> &old_coord,
                               math::DenseVecView<Real> &new_coord) {
    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      new_coord[d] = factor * old_coord[d];
    }
  };

  mesh.geo_transform(scale_transf);

#else

  typename result_of::dof_map_t<MeshConfig> &cell_dofs = *(mesh.dof_storage(dof_handler_name));

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
    const MeshEntity cell                         = cell_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::CellGeometry<MeshConfig::GDIM> cell_geometry = tcell_view.coordinates();

    auto scale_transf = [factor, cell_geometry](const Uint n,
                                                const math::DenseConstVecView<Real> &old_coord,
                                                math::DenseVecView<Real> &new_coord) {
      const math::DenseVecView<const Real> old_geo_node = cell_geometry.const_node_view(n);
      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        new_coord[d] = factor * old_geo_node[d];
      }
    };

    mesh.geo_transform(mesh::ActiveIdx(c), scale_transf);
  }
#endif
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshManipulator::translate_mesh(Tria<MeshConfig> &mesh, const std::string &dof_handler_name,
                                     const math::DenseSVec<Real, MeshConfig::GDIM> &vec)
{
#if 0
  auto translate_transf = [vec](const math::DenseConstVecView<Real> &old_coord,
                                math::DenseVecView<Real> &new_coord) {
    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      new_coord[d] = old_coord[d] + vec[d];
    }
  };

  mesh.geo_transform(translate_transf);

#else

  typename result_of::dof_map_t<MeshConfig> &cell_dofs = *(mesh.dof_storage(dof_handler_name));

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
    const MeshEntity cell                         = cell_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::CellGeometry<MeshConfig::GDIM> cell_geometry = tcell_view.coordinates();

    auto translate_transf = [vec, cell_geometry](const Uint n,
                                                 const math::DenseConstVecView<Real> &old_coord,
                                                 math::DenseVecView<Real> &new_coord) {
      const math::DenseVecView<const Real> old_geo_node = cell_geometry.const_node_view(n);
      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        new_coord[d] = old_geo_node[d] + vec[d];
      }
    };

    mesh.geo_transform(mesh::ActiveIdx(c), translate_transf);
  }
#endif
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
