#ifndef PDEKIT_Interpolation_Function_Space_Cells_hpp
#define PDEKIT_Interpolation_Function_Space_Cells_hpp

#include <functional>
#include <memory>

#include "common/DataMap.hpp"
#include "common/IteratorRange.hpp"
#include "interpolation/FEValues.hpp"
#include "mesh/DiscreteElemKey.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/MeshEntity.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "mesh/view/CellTopologyView.hpp"

namespace pdekit
{

namespace interpolation
{

template <typename MeshConfig>
class FunctionSpaceCells
{
  public:
  /// TYPEDEFS:
  using ref_elem_map = common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>;
  // One discrete element has
  // 1) A tag of point set (reference element) where it is defined
  // 2) A tag of polynomial (Finite Element basis) used in this element
  // 3) A tag of point set where the basis is evaluated (quadrature points)
  using discrete_elem_store = std::vector<mesh::DiscreteElemKey>;
  using sf_generator_fcn = std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)>;
  using quad_generator_fcn =
      std::function<mesh::PointSetTag(const ElemShape, const Uint elem_order)>;

  /// Constructor
  FunctionSpaceCells();

  /// Deleted copy constructor
  FunctionSpaceCells(const FunctionSpaceCells &other) = delete;

  /// Destructor
  ~FunctionSpaceCells();

  /// Deleted assignment operator
  FunctionSpaceCells &operator=(const FunctionSpaceCells &other) = delete;

  const sf_generator_fcn &sf_generator() const;

  const quad_generator_fcn &quad_generator() const;

  /// Configure the reference elements for all element types on this mesh
  template <typename DofIterator>
  void set_reference_fe_values(
      const common::IteratorRange<DofIterator> &dofs,
      const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
      const std::function<mesh::PointSetTag(const ElemShape shape, const Uint elem_order)>
          &quad_generator,
      const bool use_filter = false);

  /// Configure the reference elements for all element types on this mesh
  /*
  template <typename DofIterator, typename SFGenerator, typename
  QuadGenerator> void set_reference_fe_values(const
  common::IteratorRange<DofIterator> &dofs, const SFGenerator &sf_generator,
  const QuadGenerator &quad_generator, const bool use_filter = false);
  */

  private:
  /// METHODS
  inline static bool discrete_elem_found(const discrete_elem_store &elems,
                                         const mesh::DiscreteElemKey candidate)
  {
    for (const auto &elem : elems)
    {
      if (candidate == elem)
      {
        return true;
      }
    }
    return false;
  }

  protected:
  /// DATA
  /// Function that generates suitable shape functions
  sf_generator_fcn m_sf_generator;

  /// Function that generates appropriate quadratures
  quad_generator_fcn m_quad_generator;

  /// Shape functions for all reference elements in this function space
  ref_elem_map m_elements;

  discrete_elem_store m_discrete_elems;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
FunctionSpaceCells<MeshConfig>::FunctionSpaceCells()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
FunctionSpaceCells<MeshConfig>::~FunctionSpaceCells()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const typename FunctionSpaceCells<MeshConfig>::sf_generator_fcn &FunctionSpaceCells<
    MeshConfig>::sf_generator() const
{
  return m_sf_generator;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const typename FunctionSpaceCells<MeshConfig>::quad_generator_fcn &FunctionSpaceCells<
    MeshConfig>::quad_generator() const
{
  return m_quad_generator;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename DofIterator>
void FunctionSpaceCells<MeshConfig>::set_reference_fe_values(
    const common::IteratorRange<DofIterator> &dofs,
    const std::function<mesh::sf::SFTag(const ElemShape shape, const Uint order)> &sf_generator,
    const std::function<mesh::PointSetTag(const ElemShape shape, const Uint elem_order)>
        &quad_generator,
    const bool use_filter)
{
  m_elements.clear();
  m_discrete_elems.clear();

  m_sf_generator   = sf_generator;
  m_quad_generator = quad_generator;

  mesh::StdPointSet quadrature;

  for (DofIterator dof_iter = dofs.begin(); dof_iter != dofs.end(); ++dof_iter)
  {
    const mesh::PointSetTag dof_std_reg_tag = dof_iter->pt_set_id();
    const ElemShape shape                   = dof_std_reg_tag.elem_shape();
    const Uint generator_poly_deg           = dof_std_reg_tag.poly_order();

    mesh::sf::SFTag const sf_type_tag = sf_generator(shape, generator_poly_deg);

    const mesh::PointSetTag quad_tag = quad_generator(shape, generator_poly_deg);

    quadrature.change_type(quad_tag);

    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      const mesh::PointSetTagExt ext_tag(dof_std_reg_tag, P0, mesh::CellTransform::NO_TRANS, ie);
      common::PtrHandle<FEValues> elem_ptr = m_elements.std_region_data(ext_tag);

      if (elem_ptr.is_null())
      {
        elem_ptr = m_elements.create(ext_tag);
        (*elem_ptr).configure(dof_std_reg_tag, sf_type_tag);

        (*elem_ptr).fill_Vandermonde(quadrature.get().coordinates(ie), quadrature.get().weights(ie),
                                     use_filter);
      }
    }

    // Loop that generates new discrete elements
    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      const mesh::PointSetTagExt ext_std_reg_tag(dof_std_reg_tag, P0, mesh::CellTransform::NO_TRANS,
                                                 ie);
      const mesh::PointSetTagExt ext_quad_tag(quad_tag, P0, mesh::CellTransform::NO_TRANS, ie);
      const mesh::DiscreteElemKey candidate(ext_std_reg_tag, sf_type_tag, ext_quad_tag);
      if (!discrete_elem_found(m_discrete_elems, candidate))
      {
        m_discrete_elems.push_back(candidate);
      }
    }

  } // Iterator loop over dofs
}

// ----------------------------------------------------------------------------

#if 0
template <typename MeshConfig>
template <typename DofIterator, typename SFGenerator, typename QuadGenerator>
void FunctionSpaceCells<MeshConfig>::set_reference_fe_values(
  const common::IteratorRange<DofIterator> &dofs, const SFGenerator &sf_generator,
  const QuadGenerator &quad_generator, const bool use_filter)
{
  m_elements.clear();
  m_discrete_elems.clear();

  mesh::StdPointSet quadrature;

  for (DofIterator dof_iter = dofs.begin(); dof_iter != dofs.end(); ++dof_iter)
  {
    const mesh::PointSetTag dof_std_reg_tag = dof_iter->cell_type();
    const ElemShape shape = dof_std_reg_tag.elem_shape();
    const Uint generator_poly_deg = dof_std_reg_tag.poly_order();

    mesh::sf::SFTag const sf_type_tag = sf_generator(shape, generator_poly_deg);

    const mesh::PointSetTag quad_tag = quad_generator(shape, generator_poly_deg);

    quadrature.change_type(quad_tag);

    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      const mesh::PointSetTagExt ext_tag(dof_std_reg_tag, P0, mesh::CellTransform::NO_TRANS, ie);
      common::PtrHandle<FEValues> elem_ptr = m_elements.std_region_data(ext_tag);

      if (elem_ptr.is_null())
      {
        elem_ptr = m_elements.create(ext_tag);
        (*elem_ptr).configure(dof_std_reg_tag, sf_type_tag);

        (*elem_ptr).fill_Vandermonde(quadrature.get().coordinates(ie), quadrature.get().weights(ie),
                                     use_filter);
      }
    }

    // Loop that generates new discrete elements
    for (Uint ie = 0; ie < quadrature.get().nb_local_entities(); ++ie)
    {
      const mesh::PointSetTagExt ext_std_reg_tag(dof_std_reg_tag, P0, mesh::CellTransform::NO_TRANS,
                                                 ie);
      const mesh::PointSetTagExt ext_quad_tag(quad_tag, P0, mesh::CellTransform::NO_TRANS, ie);
      const mesh::DiscreteElemKey candidate(ext_std_reg_tag, sf_type_tag, ext_quad_tag);
      if (!discrete_elem_found(m_discrete_elems, candidate))
      {
        m_discrete_elems.push_back(candidate);
      }
    }

  } // Iterator loop over dofs
}
#endif

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // namespace pdekit

#endif
