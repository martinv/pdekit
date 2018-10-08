/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE geometry_cache_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include "interpolation/GeometryCache.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"

// ----------------------------------------------------------------------------

using namespace pdekit;

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class DiscreteElemKeyGenerator
{
  public:
  DiscreteElemKeyGenerator(const mesh::DofMap<MeshConfig> &dofs,
                           const common::Range1D<Uint> &cell_range, SFunc sf_type,
                           PointSetID quad_type, Uint quad_order)
      : m_dofs(dofs), m_cell_range(cell_range), m_sf_type(sf_type), m_quad_type(quad_type),
        m_quad_order(quad_order), m_current(cell_range.lbound())
  {
  }

  DiscreteElemKeyGenerator(const DiscreteElemKeyGenerator &other) = delete;

  ~DiscreteElemKeyGenerator() = default;

  DiscreteElemKeyGenerator &operator=(const DiscreteElemKeyGenerator &other) = delete;

  void init()
  {
    m_current = m_cell_range.lbound();
  }

  mesh::DiscreteElemKey make_next()
  {
    const mesh::MeshEntity cell      = m_dofs.active_cell(mesh::ActiveIdx(m_current));
    const mesh::PointSetTag cell_tag = cell.pt_set_id();

    const ElemShape cell_shape = cell_tag.elem_shape();
    const Uint cell_order      = cell_tag.poly_order();

    const mesh::PointSetTagExt cell_tag_ext(cell_tag, P0, mesh::CellTransform::NO_TRANS, 0u);

    const mesh::sf::SFTag cell_basis(cell_shape, m_sf_type, cell_order, ModalBasis::Modal);

    const mesh::PointSetTag cell_quad(cell_shape, m_quad_order, m_quad_type);
    const mesh::PointSetTagExt cell_quad_ext(cell_quad, P0, mesh::CellTransform::NO_TRANS, 0u);
    m_current++;

    return mesh::DiscreteElemKey{cell_tag_ext, cell_basis, cell_quad_ext};
  }

  bool reached_end() const
  {
    return m_current == (m_cell_range.ubound() + 1);
  }

  private:
  const mesh::DofMap<MeshConfig> &m_dofs;
  const common::Range1D<Uint> m_cell_range;
  const SFunc m_sf_type;
  const PointSetID m_quad_type;
  const Uint m_quad_order;
  Uint m_current;
};

// ----------------------------------------------------------------------------

using MeshType2D = mesh::Tria<mesh::Cart2D>;
using MeshType3D = mesh::Tria<mesh::Cart3D>;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(geo_cache_2D_utest)
{
  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_mixed_p2.msh", *mesh2d, "geo_dofs");

  const auto geo_dofs = mesh2d->dof_storage("geo_dofs");
  const common::Range1D<Uint> cell_range(0, (*geo_dofs).nb_active_cells() - 1);

  DiscreteElemKeyGenerator<mesh::Cart2D> key_generator(*geo_dofs, cell_range, SFunc::Lagrange,
                                                       PointSetID::Gauss, P2);

  interpolation::GeometryCache<mesh::Cart2D::GDIM, interpolation::CacheInsertAutomatic> geo_cache;
  geo_cache.record_insertion_pattern(key_generator, cell_range.size());
  // geo_cache.print_contents();
  geo_cache.print_types();
}

// ----------------------------------------------------------------------------
