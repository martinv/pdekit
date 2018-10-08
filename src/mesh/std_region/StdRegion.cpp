#include "mesh/std_region/StdRegion.hpp"
#include "mesh/std_region/StdRegionFactory.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Class 'StdRegionInstance'
// ----------------------------------------------------------------------------

const PointSetTag StdRegionInstance::undefined =
    PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined);

// ----------------------------------------------------------------------------

void StdRegionInstance::construct(const PointSetTag type_id, StdRegionInstance &std_reg_instance)
{
  // Get the interpolation point set from factory to fill the local
  // connectivity and coordinates
  StdRegionFactory::instance_type &factory = StdRegionFactory::instance();
  const StdRegionFactory::instance_type::const_product_base_ptr std_region_builder =
      factory.create(type_id);

  std_reg_instance.m_type_id     = type_id;
  std_reg_instance.m_topo_dim    = std_region_builder->topo_dim();
  std_reg_instance.m_nb_p1_nodes = std_region_builder->nb_p1_dof();

  std_region_builder->fill_reference_topology(std_reg_instance.m_incidences,
                                              std_reg_instance.m_entity_storage);
  std_region_builder->facet_normals(std_reg_instance.m_facet_normals);
}

// ----------------------------------------------------------------------------

StdRegionInstance::StdRegionInstance() : m_type_id(StdRegionInstance::undefined), m_topo_dim(_0D)
{
}

// ----------------------------------------------------------------------------

StdRegionInstance::StdRegionInstance(const PointSetTag type_id)
{
  /*
  // Get the interpolation point set from factory to fill the local
  connectivity and
  // coordinates
  ReferenceTopologyFactory::instance_type& factory =
  ReferenceTopologyFactory::instance();
  InterpolationPointSet::ptr ips = factory.create(type_id);

  m_type_id = type_id;
  m_topo_dim = ips->topo_dim();
  ips->fill_reference_topology(m_incidences,m_entity_storage);
  */
  construct(type_id, *this);
}

// ----------------------------------------------------------------------------

StdRegionInstance::StdRegionInstance(const StdRegionInstance &std_reg_instance)
{
  m_type_id        = std_reg_instance.m_type_id;
  m_topo_dim       = std_reg_instance.m_topo_dim;
  m_nb_p1_nodes    = std_reg_instance.m_nb_p1_nodes;
  m_incidences     = std_reg_instance.m_incidences;
  m_entity_storage = std_reg_instance.m_entity_storage;
}

// ----------------------------------------------------------------------------

StdRegionInstance &StdRegionInstance::operator=(const StdRegionInstance &std_reg_instance)
{
  m_type_id        = std_reg_instance.m_type_id;
  m_topo_dim       = std_reg_instance.m_topo_dim;
  m_nb_p1_nodes    = std_reg_instance.m_nb_p1_nodes;
  m_incidences     = std_reg_instance.m_incidences;
  m_entity_storage = std_reg_instance.m_entity_storage;
  return *this;
}

// ----------------------------------------------------------------------------

StdRegionInstance::~StdRegionInstance()
{
}

// ----------------------------------------------------------------------------

Uint StdRegionInstance::nb_entities(const Uint dim) const
{
  return m_entity_storage.const_block(dim).size();
}

// ----------------------------------------------------------------------------

StdRegionInstance::coordinates_type const &StdRegionInstance::coordinates() const
{
  return m_entity_storage.const_block(m_topo_dim)[0]->coordinates();
}

// ----------------------------------------------------------------------------

std::shared_ptr<StdRegionEntity const> StdRegionInstance::elem_entity(const Uint dim,
                                                                      const Uint id) const
{
  const common::ArrayView<const std::shared_ptr<StdRegionEntity>, _1D, Uint> entities_in_dim =
      m_entity_storage.const_block(dim);
  std::shared_ptr<StdRegionEntity const> entity_ptr(entities_in_dim[id]);
  return entity_ptr;
}

// ----------------------------------------------------------------------------

StdRegionInstance::incidence_list_type StdRegionInstance::incident_entities(
    const Uint my_dim, const Uint my_id, const Uint other_dim) const
{
  const incidence_table_type &itable = m_incidences[my_dim * (_3D + 1) + other_dim];
  return itable.const_block(my_id);
}

// ----------------------------------------------------------------------------

const math::DenseConstVecView<Real> StdRegionInstance::facet_normal(const Uint facet_id) const
{
  return m_facet_normals.const_row_transp(facet_id);
}

// ----------------------------------------------------------------------------

const std::string StdRegionInstance::type_name() const
{
  return m_type_id.as_string();
}

// ----------------------------------------------------------------------------

void StdRegionInstance::print_complete_topology() const
{
  std::cout << "Complete topology of reference element " << type_name() << ": " << std::endl;

  for (Uint d = 1; d <= m_topo_dim; ++d)
  {
    std::cout << " " << d << "D entities:" << std::endl;
    for (Uint e = 0; e < nb_entities(d); ++e)
    {
      std::cout << "   " << *(elem_entity(d, e)) << std::endl;
    }
    std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
