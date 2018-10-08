#ifndef PDEKIT_Solver_Setup_Algorithm_hpp
#define PDEKIT_Solver_Setup_Algorithm_hpp

#include <unordered_map>

#include "common/DataMap.hpp"
#include "common/IteratorRange.hpp"
#include "interpolation/FEValues.hpp"
#include "mesh/DiscreteElemKey.hpp"
#include "mesh/MeshEntity.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "mesh/std_region/PointSetTagExt.hpp"

namespace pdekit
{

namespace solver
{

// ----------------------------------------------------------------------------

class SolverSetupAlgorithm
{
  public:
  SolverSetupAlgorithm() = default;

  ~SolverSetupAlgorithm() = default;

  template <typename GeoDofIterType, typename SolDofIterType, typename PolySelectRule>
  static void generate_fe_pairs(
      const GeoDofIterType &geo_dofs_begin,
      const std::vector<common::IteratorRange<SolDofIterType>> &sol_dofs,
      const PolySelectRule &rule, std::vector<mesh::DiscreteElemKey> &elem_types_geo,
      common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> &elem_types_sol);

  template <typename PolySelectRule>
  static mesh::DiscreteElemKey geo_key(const mesh::PointSetTag geo_cell_tag,
                                       const mesh::PointSetTag sol_cell_tag,
                                       const PolySelectRule &rule);
};

// ----------------------------------------------------------------------------

template <typename GeoDofIterType, typename SolDofIterType, typename PolySelectRule>
void SolverSetupAlgorithm::generate_fe_pairs(
    const GeoDofIterType &geo_dofs_begin,
    const std::vector<common::IteratorRange<SolDofIterType>> &sol_dofs, const PolySelectRule &rule,
    std::vector<mesh::DiscreteElemKey> &elem_types_geo,
    common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> &elem_types_sol)
{
  mesh::StdPointSet quad;

  clock_t start, end;
  Real elapsed;

  start = clock();

  elem_types_geo.clear();

  // --------------------------------------------------------------------------
  // Loop over all cells in the solution mesh
  // Fill maps with FeValues in geometry and solution
  // --------------------------------------------------------------------------

  GeoDofIterType geo_iterator = geo_dofs_begin;

  for (const auto &range : sol_dofs)
  {
    for (auto sol_iterator = range.begin(); sol_iterator != range.end(); ++sol_iterator)
    {
      mesh::synchronize_dof_iterators(sol_iterator, geo_iterator);

      const mesh::MeshEntity geo_cell = geo_iterator->mesh_entity();
      const mesh::MeshEntity sol_cell = sol_iterator->mesh_entity();

      const mesh::PointSetTag geo_cell_tag = geo_cell.pt_set_id();
      const ElemShape eshape_geo           = geo_cell_tag.elem_shape();
      const Uint poly_order_geo            = geo_cell_tag.poly_order();

      const mesh::PointSetTag sol_cell_tag = sol_cell.pt_set_id();

      const mesh::sf::SFTag space_tag = rule.sf_tag(geo_cell_tag, sol_cell_tag);

      const Uint quad_order = rule.quadrature_order(geo_cell_tag, sol_cell_tag);

      const mesh::PointSetTagExt geo_cell_tag_ext(geo_cell_tag, P0, mesh::CellTransform::NO_TRANS,
                                                  0);
      const mesh::sf::SFTag geo_sf_tag(eshape_geo, SFunc::Lagrange, poly_order_geo,
                                       ModalBasis::Modal);
      const mesh::PointSetTag quad_tag_geo(eshape_geo, quad_order, PointSetID::Gauss);
      const mesh::PointSetTagExt quad_tag_geo_ext(quad_tag_geo, P0, mesh::CellTransform::NO_TRANS,
                                                  0);
      const mesh::DiscreteElemKey geo_key(geo_cell_tag_ext, geo_sf_tag, quad_tag_geo_ext);

      mesh::add_unique_discr_elem_key(elem_types_geo, geo_key);

      quad.change_type(quad_tag_geo);

      common::PtrHandle<interpolation::FEValues> fe_values_sol_ptr =
          elem_types_sol.std_region_data(mesh::PointSetTagExt(sol_cell_tag, quad_order));

      if (fe_values_sol_ptr.is_null())
      {
        fe_values_sol_ptr = elem_types_sol.create(mesh::PointSetTagExt(sol_cell_tag, quad_order));
        (*fe_values_sol_ptr).configure(sol_cell_tag, space_tag);
        (*fe_values_sol_ptr).fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
      }

    } // Loop over cells of one range
  }   // Loop over ranges

  // --------------

  end = clock();

  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  // std::cout << "Elemental matrix builder: mass matrices built in " <<
  // elapsed
  // << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename PolySelectRule>
mesh::DiscreteElemKey SolverSetupAlgorithm::geo_key(const mesh::PointSetTag geo_cell_tag,
                                                    const mesh::PointSetTag sol_cell_tag,
                                                    const PolySelectRule &rule)
{
  const ElemShape eshape_geo = geo_cell_tag.elem_shape();
  const Uint poly_order_geo  = geo_cell_tag.poly_order();

  const Uint quad_order = rule.quadrature_order(geo_cell_tag, sol_cell_tag);

  const mesh::PointSetTagExt geo_cell_tag_ext(geo_cell_tag, P0, mesh::CellTransform::NO_TRANS, 0);
  const mesh::sf::SFTag geo_sf_tag(eshape_geo, SFunc::Lagrange, poly_order_geo, ModalBasis::Modal);
  const mesh::PointSetTag quad_tag_geo(eshape_geo, quad_order, PointSetID::Gauss);
  const mesh::PointSetTagExt quad_tag_geo_ext(quad_tag_geo, P0, mesh::CellTransform::NO_TRANS, 0);
  const mesh::DiscreteElemKey geo_key(geo_cell_tag_ext, geo_sf_tag, quad_tag_geo_ext);

  return geo_key;
}

// ----------------------------------------------------------------------------

} // namespace solver

} // namespace pdekit

#endif
