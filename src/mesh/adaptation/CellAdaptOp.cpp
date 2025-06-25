#include <iostream>

#include "mesh/CellTransform.hpp"
#include "mesh/adaptation/CellAdaptOp.hpp"
#include "mesh/adaptation/CellAdaptOpFactory.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

// ----------------------------------------------------------------------------

const CellAdaptOpTag CellAdaptOpInstance::undefined =
    CellAdaptOpTag(ElemShape::Undefined, CellTransform::NO_TRANS);

// ----------------------------------------------------------------------------

void CellAdaptOpInstance::construct(const CellAdaptOpTag key, CellAdaptOpInstance &concrete_op)
{
  CellAdaptOpFactory::instance_type &adapt_op_factory = CellAdaptOpFactory::instance();

  const CellAdaptOpFactory::instance_type::const_product_base_ptr adapt_op_base =
      adapt_op_factory.create(key);

  adapt_op_base->set_facet_adapt_op_ids(concrete_op.m_cell_adapt_op_tag,
                                        concrete_op.m_facet_adapt_op_id);
  concrete_op.m_nb_parent_facets = adapt_op_base->nb_parent_facets();
  concrete_op.m_nb_child_elems   = adapt_op_base->nb_child_elems();

  adapt_op_base->fill_parent_ref_coords(concrete_op.m_parent_ref_coords);
  adapt_op_base->set_child_elem_shapes(concrete_op.m_child_elem_shapes);
  adapt_op_base->fill_local_child_incidences(concrete_op.m_local_internal_child_incidences,
                                             concrete_op.m_local_internal_child_permutations,
                                             concrete_op.m_local_facet_offsets);

  concrete_op.m_parent_child_incidences_on_facets.resize(adapt_op_base->nb_parent_facets());
  concrete_op.m_parent_child_permutations_on_facets.resize(adapt_op_base->nb_parent_facets());

  for (Uint f = 0; f < adapt_op_base->nb_parent_facets(); ++f)
  {
    adapt_op_base->fill_parent_child_incidences_on_facet(
        f, concrete_op.m_parent_child_incidences_on_facets[f],
        concrete_op.m_parent_child_permutations_on_facets[f]);
  }

  adapt_op_base->fill_child_coord_transformers(concrete_op.m_coord_transformers);
}

// ----------------------------------------------------------------------------

CellAdaptOpInstance::CellAdaptOpInstance()
{
  m_cell_adapt_op_tag = undefined;
  m_nb_parent_facets  = 0u;
  m_nb_child_elems    = 0u;
  m_parent_ref_coords.resize(0, 0);
  m_facet_adapt_op_id.resize(0);
  m_child_elem_shapes.resize(0);
  m_local_internal_child_incidences.resize(0);
  m_local_internal_child_permutations.resize(0);
  m_local_facet_offsets.resize(0);
  m_parent_child_incidences_on_facets.resize(0);
  m_parent_child_permutations_on_facets.resize(0);
  m_coord_transformers.resize(0);
}

// ----------------------------------------------------------------------------

CellAdaptOpInstance::CellAdaptOpInstance(CellAdaptOpTag &key)
{
  construct(key, *this);
}

// ----------------------------------------------------------------------------

CellAdaptOpInstance::CellAdaptOpInstance(const CellAdaptOpInstance &other_split_strategy)
{
  m_cell_adapt_op_tag = other_split_strategy.m_cell_adapt_op_tag;
  m_nb_parent_facets  = other_split_strategy.m_nb_parent_facets;
  m_nb_child_elems    = other_split_strategy.m_nb_child_elems;
  m_parent_ref_coords.resize(other_split_strategy.m_parent_ref_coords.rows(),
                             other_split_strategy.m_parent_ref_coords.cols());
  m_parent_ref_coords = other_split_strategy.m_parent_ref_coords;

  m_facet_adapt_op_id                 = other_split_strategy.m_facet_adapt_op_id;
  m_child_elem_shapes                 = other_split_strategy.m_child_elem_shapes;
  m_local_internal_child_incidences   = other_split_strategy.m_local_internal_child_incidences;
  m_local_internal_child_permutations = other_split_strategy.m_local_internal_child_permutations;
  m_local_facet_offsets               = other_split_strategy.m_local_facet_offsets;

  m_parent_child_incidences_on_facets.resize(
      other_split_strategy.m_parent_child_incidences_on_facets.size());
  m_parent_child_permutations_on_facets.resize(
      other_split_strategy.m_parent_child_permutations_on_facets.size());

  for (Uint f = 0; f < m_nb_parent_facets; ++f)
  {
    m_parent_child_incidences_on_facets[f] =
        other_split_strategy.m_parent_child_incidences_on_facets[f];
    m_parent_child_permutations_on_facets[f] =
        other_split_strategy.m_parent_child_permutations_on_facets[f];
  }

  m_coord_transformers = other_split_strategy.m_coord_transformers;
}

// ----------------------------------------------------------------------------

CellAdaptOpInstance &CellAdaptOpInstance::operator=(const CellAdaptOpInstance &other_split_strategy)
{
  m_cell_adapt_op_tag = other_split_strategy.m_cell_adapt_op_tag;
  m_nb_parent_facets  = other_split_strategy.m_nb_parent_facets;
  m_nb_child_elems    = other_split_strategy.m_nb_child_elems;

  m_parent_ref_coords.resize(other_split_strategy.m_parent_ref_coords.rows(),
                             other_split_strategy.m_parent_ref_coords.cols());
  m_parent_ref_coords = other_split_strategy.m_parent_ref_coords;

  m_facet_adapt_op_id                 = other_split_strategy.m_facet_adapt_op_id;
  m_child_elem_shapes                 = other_split_strategy.m_child_elem_shapes;
  m_local_internal_child_incidences   = other_split_strategy.m_local_internal_child_incidences;
  m_local_internal_child_permutations = other_split_strategy.m_local_internal_child_permutations;
  m_local_facet_offsets               = other_split_strategy.m_local_facet_offsets;

  m_parent_child_incidences_on_facets.resize(
      other_split_strategy.m_parent_child_incidences_on_facets.size());
  m_parent_child_permutations_on_facets.resize(
      other_split_strategy.m_parent_child_permutations_on_facets.size());

  for (Uint f = 0; f < m_nb_parent_facets; ++f)
  {
    m_parent_child_incidences_on_facets[f] =
        other_split_strategy.m_parent_child_incidences_on_facets[f];
    m_parent_child_permutations_on_facets[f] =
        other_split_strategy.m_parent_child_permutations_on_facets[f];
  }

  m_coord_transformers = other_split_strategy.m_coord_transformers;

  return *this;
}

// ----------------------------------------------------------------------------

CellAdaptOpInstance::~CellAdaptOpInstance()
{
}

// ----------------------------------------------------------------------------

const std::tuple<SUint, SUint> CellAdaptOpInstance::containing_parent_facet_id(
    const Uint child_id, const Uint child_facet_id) const
{
  for (Uint parent_f = 0; parent_f < m_parent_child_incidences_on_facets.size(); ++parent_f)
  {
    const std::vector<IncidenceEntry> &facet_incidences =
        m_parent_child_incidences_on_facets[parent_f];
    for (Uint i = 0; i < facet_incidences.size(); ++i)
    {
      if ((facet_incidences[i].cell_idx == child_id) &&
          (facet_incidences[i].local_id == child_facet_id))
      {
        return std::tuple<SUint, SUint>(parent_f, i);
      }
    }
  }
  return std::tuple<SUint, SUint>(INVALID_LOC_ENTITY_ID, INVALID_LOC_ENTITY_ID);
}

// ----------------------------------------------------------------------------

Uint CellAdaptOpInstance::containing_parent_interior_facet(const Uint child_id,
                                                           const Uint child_facet_id) const
{
  for (Uint f = 0; f < m_local_internal_child_incidences.size(); ++f)
  {
    if ((m_local_internal_child_incidences[f].cell_idx == child_id) &&
        (m_local_internal_child_incidences[f].local_id == child_facet_id))
    {
      return f;
    }
  }

  return INVALID_LOC_ENTITY_ID;
}

// ----------------------------------------------------------------------------

void CellAdaptOpInstance::compute_child_coords(
    const PointSetTag parent_type, std::vector<math::DenseDMat<Real>> &child_coords) const
{
  StdRegion parent_region_no_transform(parent_type);
  math::DenseDMat<Real> const &parent_elem_coords = parent_region_no_transform.get().coordinates();

  child_coords.resize(m_nb_child_elems);
  const Uint nb_nodes = parent_elem_coords.rows();
  const Uint dim      = parent_elem_coords.cols();

  for (Uint c = 0; c < m_nb_child_elems; ++c)
  {
    child_coords[c].resize(nb_nodes, dim);
  }

  for (Uint n = 0; n < nb_nodes; ++n)
  {
    const math::DenseConstVecView<Real> node_coord = parent_elem_coords.const_row_transp(n);
    for (Uint c = 0; c < m_nb_child_elems; ++c)
    {
      math::DenseVecView<Real> child_elem_coord = child_coords[c].row_transp(n);
      m_coord_transformers[c](m_parent_ref_coords, node_coord, child_elem_coord);
    }
  }
}

// ----------------------------------------------------------------------------

void CellAdaptOpInstance::transform_coords(const math::DenseDMat<Real> &parent_ref_coords,
                                           const math::DenseDMat<Real> &coords_in,
                                           const Uint child_id,
                                           math::DenseDMat<Real> &coords_out) const
{
  coords_out.resize(coords_in.rows(), coords_in.cols());

  for (Uint r = 0; r < coords_in.rows(); ++r)
  {
    const math::DenseConstVecView<Real> node_coord_in = coords_in.const_row_transp(r);
    math::DenseVecView<Real> node_coord_out           = coords_out.row_transp(r);
    m_coord_transformers[child_id](parent_ref_coords, node_coord_in, node_coord_out);
  }
}

// ----------------------------------------------------------------------------

void CellAdaptOpInstance::print() const
{
  std::cout << m_cell_adapt_op_tag.as_string() << std::endl;
  std::cout << "Number of parent facets: " << m_nb_parent_facets << std::endl;
  std::cout << "Number of child elements: " << m_nb_child_elems << std::endl;
  std::cout << "Parent reference coordinates:" << std::endl;
  std::cout << m_parent_ref_coords << std::endl;

  std::cout << "Children:" << std::endl;

  for (Uint s = 0; s < nb_child_elems(); ++s)
  {
    std::cout << ElemShapeInfo::name(m_child_elem_shapes[s]) << std::endl;
  }

  std::cout << "Parent-child incidences on facets of parent:" << std::endl;
  for (Uint f = 0; f < m_nb_parent_facets; ++f)
  {
    auto parent_child_inc = parent_child_incidences_on_facet(f);
    std::cout << " Parent facet " << f << std::endl;
    std::cout << " " << parent_child_inc << std::endl;
  }

  const Uint nb_internal_facets =
      m_local_facet_offsets.size() <= 1 ? 0 : m_local_facet_offsets.size() - 1;
  std::cout << "Number of local internal child facets: " << nb_internal_facets << std::endl;

  for (Uint f = 0; f < nb_internal_facets; ++f)
  {
    auto internal_child_inc = internal_child_facet(f);
    std::cout << " Child facet " << f << std::endl;
    std::cout << " " << internal_child_inc << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace adapt

} // namespace mesh

} // namespace pdekit
