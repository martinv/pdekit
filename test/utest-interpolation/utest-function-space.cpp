/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE function_space_utest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>

#include "common/PDEKit.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

struct ComputeQuadOrder
{
  Uint operator()(const mesh::PointSetTag &geo_tag, const mesh::PointSetTag &sol_tag) const
  {
    return std::max(1u, std::max(2 * geo_tag.poly_order(), 2 * sol_tag.poly_order()));
  }
};

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(function_space_find_fe_pairs_test)
{
  typedef result_of::dof_map_t<mesh::Cart2D> dof_map_type;

  using map_type = common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>;

  mesh::Tria<mesh::Cart2D> mesh2d("mesh2d");
  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_mixed_p2.msh", mesh2d, "geo_dofs");

  common::PtrHandle<dof_map_type> geo_dofs = mesh2d.dof_storage("geo_dofs");
  common::PtrHandle<dof_map_type> sol_dofs = mesh2d.create_dof_storage("sol_dofs");

  dof_map_type::clone_continuous(mesh2d, *geo_dofs, *sol_dofs, P3, PointSetID::Warpblend);

  map_type geo_fe_values;
  map_type sol_fe_values;

  std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> tag_combinations;

  auto sf_generator = [=](const ElemShape shape, const Uint order,
                          const PointSetID rt_type) -> mesh::sf::SFTag {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto quad_generator = [=](const ElemShape shape, const Uint geo_elem_order,
                            const Uint sol_elem_order) -> mesh::PointSetTag {
    return mesh::PointSetTag(shape, std::max(1u, std::max(2 * geo_elem_order, 2 * sol_elem_order)),
                             PointSetID::Gauss);
  };

  interpolation::fill_fe_value_pairs<mesh::Cart2D>(*geo_dofs, *sol_dofs, sf_generator,
                                                   quad_generator, geo_fe_values, sol_fe_values,
                                                   tag_combinations);

  for (map_type::const_iterator it = geo_fe_values.cbegin(); it != geo_fe_values.cend(); ++it)
  {
    const common::PtrHandle<interpolation::FEValues const> geo_fe_val_ptr = it.data_ptr();

    const interpolation::FEValues &fe_geo = (*geo_fe_val_ptr);

    std::cout << "****************************************************" << std::endl;
    std::cout << "[ " << fe_geo.std_region_id().as_string() << " ]" << std::endl;

    std::cout << "GEOMETRY:" << std::endl;
    fe_geo.print();
  }

  for (map_type::const_iterator it = sol_fe_values.cbegin(); it != sol_fe_values.cend(); ++it)
  {
    const common::PtrHandle<interpolation::FEValues const> sol_fe_val_ptr =
        sol_fe_values.std_region_data(it.key_value());

    const interpolation::FEValues &fe_sol = (*sol_fe_val_ptr);

    std::cout << "****************************************************" << std::endl;
    std::cout << "[ " << fe_sol.std_region_id().as_string() << " ]" << std::endl;

    std::cout << "SOLUTION:" << std::endl;
    fe_sol.print();
  }

  std::cout << "Detected tag pairs:" << std::endl;
  for (auto tag_pair : tag_combinations)
  {
    std::cout << "{ " << tag_pair.first << " - " << tag_pair.second << " }" << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(function_space_find_facet_fe_pairs_test)
{
  typedef result_of::dof_map_t<mesh::Cart2D> dof_map_type;

  typedef common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> map_type;

  mesh::Tria<mesh::Cart2D> mesh2d("mesh2d");
  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_mixed_p2.msh", mesh2d, "geo_dofs");

  common::PtrHandle<dof_map_type> geo_dofs = mesh2d.dof_storage("geo_dofs");
  common::PtrHandle<dof_map_type> sol_dofs = mesh2d.create_dof_storage("sol_dofs");

  dof_map_type::clone_continuous(mesh2d, *geo_dofs, *sol_dofs, P3, PointSetID::Warpblend);

  map_type geo_fe_values;
  map_type sol_fe_values;

  std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> tag_combinations;

  ComputeQuadOrder quad_order_rule;

  interpolation::fill_fe_value_pairs<mesh::Cart2D, ComputeQuadOrder, _1D>(
      mesh2d, *geo_dofs, *sol_dofs, SFunc::Lagrange, PointSetID::Gauss, quad_order_rule,
      geo_fe_values, sol_fe_values, tag_combinations);

  for (map_type::const_iterator it = geo_fe_values.cbegin(); it != geo_fe_values.cend(); ++it)
  {
    const common::PtrHandle<interpolation::FEValues const> geo_fe_val_ptr = it.data_ptr();

    const interpolation::FEValues &fe_geo = (*geo_fe_val_ptr);

    std::cout << "****************************************************" << std::endl;
    std::cout << "[ " << fe_geo.std_region_id().as_string() << " ]" << std::endl;

    std::cout << "GEOMETRY:" << std::endl;
    fe_geo.print();
  }

  for (map_type::const_iterator it = sol_fe_values.cbegin(); it != sol_fe_values.cend(); ++it)
  {
    const common::PtrHandle<interpolation::FEValues const> sol_fe_val_ptr =
        sol_fe_values.std_region_data(it.key_value());

    const interpolation::FEValues &fe_sol = (*sol_fe_val_ptr);

    std::cout << "****************************************************" << std::endl;
    std::cout << "[ " << fe_sol.std_region_id().as_string() << " ]" << std::endl;

    std::cout << "SOLUTION:" << std::endl;
    fe_sol.print();
  }

  std::cout << "Detected tag pairs on skeleton:" << std::endl;
  for (auto tag_pair : tag_combinations)
  {
    std::cout << "{ " << tag_pair.first << " - " << tag_pair.second << " }" << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(function_space_find_facet_fe_pairs_on_boundary_test)
{
  typedef result_of::dof_map_t<mesh::Cart2D> dof_map_type;

  typedef common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> map_type;

  mesh::Tria<mesh::Cart2D> mesh2d("mesh2d");
  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_mixed_p2.msh", mesh2d, "geo_dofs");

  common::PtrHandle<dof_map_type> geo_dofs = mesh2d.dof_storage("geo_dofs");
  common::PtrHandle<dof_map_type> sol_dofs = mesh2d.create_dof_storage("sol_dofs");

  dof_map_type::clone_continuous(mesh2d, *geo_dofs, *sol_dofs, P3, PointSetID::Warpblend);

  map_type geo_fe_values;
  map_type sol_fe_values;

  std::vector<std::pair<mesh::PointSetTagExt, mesh::PointSetTagExt>> tag_combinations;

  ComputeQuadOrder quad_order_rule;

  const mesh::MeshBoundarySet<mesh::Cart2D>::bdry_facets_shared_ptr bottom =
      mesh2d.all_boundaries().domain("bottom");

  interpolation::fill_fe_value_pairs<mesh::Cart2D, ComputeQuadOrder>(
      *bottom, *geo_dofs, *sol_dofs, SFunc::Lagrange, PointSetID::Gauss, quad_order_rule,
      geo_fe_values, sol_fe_values, tag_combinations);

  for (map_type::const_iterator it = geo_fe_values.cbegin(); it != geo_fe_values.cend(); ++it)
  {
    const common::PtrHandle<interpolation::FEValues const> geo_fe_val_ptr = it.data_ptr();

    const interpolation::FEValues &fe_geo = (*geo_fe_val_ptr);

    std::cout << "****************************************************" << std::endl;
    std::cout << "[ " << fe_geo.std_region_id().as_string() << " ]" << std::endl;

    std::cout << "GEOMETRY:" << std::endl;
    fe_geo.print();
  }

  for (map_type::const_iterator it = sol_fe_values.cbegin(); it != sol_fe_values.cend(); ++it)
  {
    const common::PtrHandle<interpolation::FEValues const> sol_fe_val_ptr =
        sol_fe_values.std_region_data(it.key_value());

    const interpolation::FEValues &fe_sol = (*sol_fe_val_ptr);

    std::cout << "****************************************************" << std::endl;
    std::cout << "[ " << fe_sol.std_region_id().as_string() << " ]" << std::endl;

    std::cout << "SOLUTION:" << std::endl;
    fe_sol.print();
  }

  std::cout << "Detected tag pairs on boundary:" << std::endl;
  for (auto tag_pair : tag_combinations)
  {
    std::cout << "{ " << tag_pair.first << " - " << tag_pair.second << " }" << std::endl;
  }
}

// ----------------------------------------------------------------------------
