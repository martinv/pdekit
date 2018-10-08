#include <cassert>

#include "mesh/MeshConstants.hpp"
#include "mesh/std_region/StdRegionEntity.hpp"
#include "mesh/std_region/StdRegionFactory.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Methods of the class StdRegionEntity
// ----------------------------------------------------------------------------

StdRegionEntity::StdRegionEntity()
    : m_idx(INVALID_REF_ENTITY_ID), m_dim(0),
      m_std_region_id(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined))
{
  m_vert.resize(0);
  m_p1_vert_flag.resize(0);
}

// ----------------------------------------------------------------------------

// StdRegionEntity::StdRegionEntity( const Uint size, const Uint idx, const Uint
// topo_dim, const Uint type_id ) :
//  m_size(size),
//  m_idx(idx),
//  m_dim(topo_dim),
//  m_type(type_id)
//{
//  m_vert = new Uint[m_size];
//}

// ----------------------------------------------------------------------------

StdRegionEntity::StdRegionEntity(const StdRegionEntity &other)
    : m_idx(other.m_idx), m_dim(other.m_dim), m_std_region_id(other.m_std_region_id)
{
  m_coordinates.resize(other.m_coordinates.rows(), other.m_coordinates.cols());
  m_coordinates = other.m_coordinates;

  m_vert.resize(other.m_vert.size());
  std::copy(other.m_vert.begin(), other.m_vert.end(), m_vert.begin());

  m_p1_vert_flag.resize(other.m_p1_vert_flag.size());
  std::copy(other.m_p1_vert_flag.begin(), other.m_p1_vert_flag.end(), m_p1_vert_flag.begin());

  assert(m_p1_vert_flag.size() == m_vert.size());
}

// ----------------------------------------------------------------------------

StdRegionEntity &StdRegionEntity::operator=(const StdRegionEntity &other)
{
  m_vert.resize(other.m_vert.size());
  std::copy(other.m_vert.begin(), other.m_vert.end(), m_vert.begin());

  m_p1_vert_flag.resize(other.m_p1_vert_flag.size());
  std::copy(other.m_p1_vert_flag.begin(), other.m_p1_vert_flag.end(), m_p1_vert_flag.begin());

  assert(m_p1_vert_flag.size() == m_vert.size());

  m_idx = other.m_idx;
  m_dim = other.m_dim;

  if (m_std_region_id != other.m_std_region_id)
  {
    m_coordinates.resize(other.m_coordinates.rows(), other.m_coordinates.cols());
    m_coordinates = other.m_coordinates;
  }

  m_std_region_id = other.m_std_region_id;
  return *this;
}

// ----------------------------------------------------------------------------

StdRegionEntity::~StdRegionEntity()
{
  m_idx           = INVALID_REF_ENTITY_ID;
  m_dim           = 0;
  m_std_region_id = PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined);
  m_vert.clear();
  m_p1_vert_flag.clear();
}

// ----------------------------------------------------------------------------

void StdRegionEntity::construct(const PointSetTag pt_set_id,
                                const common::ArrayView<const Uint, _1D, Uint> &vertices,
                                const common::ArrayView<const Uint, _1D, Uint> &p1_flags,
                                const SUint entity_id)
{
  m_std_region_id = pt_set_id;

  StdRegionFactory::instance_type &std_region_factory = StdRegionFactory::instance();

  if (std_region_factory.key_is_registered(pt_set_id))
  {
    StdRegionBuilder::const_ptr std_reg_builder = std_region_factory.create(pt_set_id);
    std_reg_builder->coordinates(m_coordinates);

    const Uint nb_vert = std_reg_builder->nb_dof();

    if ((nb_vert != vertices.size()) || (nb_vert != p1_flags.size()))
    {
      std::cerr << "StdRegionEntity: Trying to construct " << pt_set_id.as_string()
                << " with input lengths " << vertices.size() << " and " << p1_flags.size()
                << std::endl;
      std::cerr << "StdRegionEntity: can't construct, wrong length of input arrays!" << std::endl;
      return;
    }

    m_vert.resize(nb_vert);
    // m_vert.assign(nb_vert, INVALID_NODE_ID);

    m_p1_vert_flag.resize(nb_vert);
    // m_p1_vert_flag.assign(nb_vert, false);

    m_dim = std_reg_builder->topo_dim();

    for (Uint i = 0; i < vertices.size(); ++i)
    {
      m_vert[i] = vertices[i];
    }
    for (Uint i = 0; i < p1_flags.size(); ++i)
    {
      if (p1_flags[i] != 0)
      {
        m_p1_vert_flag[i] = true;
      }
      else
      {
        m_p1_vert_flag[i] = false;
      }
    }

    m_idx = entity_id;
  }
  else
  {
    m_vert.resize(0u);
    m_p1_vert_flag.resize(0u);

    m_dim = 0;
  }
}

// ----------------------------------------------------------------------------

void StdRegionEntity::change_point_set_type(const PointSetID pt_set_id)
{
  const PointSetTag old_tag = m_std_region_id;
  m_std_region_id           = PointSetTag(old_tag.elem_shape(), old_tag.poly_order(), pt_set_id);

  StdRegionFactory::instance_type &std_region_factory = StdRegionFactory::instance();

  if (std_region_factory.key_is_registered(m_std_region_id))
  {
    StdRegionBuilder::const_ptr std_reg_builder = std_region_factory.create(m_std_region_id);
    std_reg_builder->coordinates(m_coordinates);
  }
}

// ----------------------------------------------------------------------------

// void StdRegionEntity::resize(const Uint size)
//{
//  m_vert.resize(size);
//  m_vert.assign(size,INVALID_NODE_ID);

//  m_p1_vert_flag.resize(size);
//  m_p1_vert_flag.assign(size,false);
//}

// ----------------------------------------------------------------------------

/*
void StdRegionEntity::set_vertex(const Uint vertex_pos, const Uint vertex_value)
{
  m_vert[vertex_pos] = vertex_value;
}

// ----------------------------------------------------------------------------

void StdRegionEntity::set_p1_vert_flag(const Uint p1_vertex_pos, const bool flag)
{
  m_p1_vert_flag[p1_vertex_pos] = flag;
}

// ----------------------------------------------------------------------------

void StdRegionEntity::set_id(const SUint my_id)
{
  m_idx = my_id;
}
*/

// ----------------------------------------------------------------------------

// void StdRegionEntity::set_dim(const SUint dim)
//{
//  m_dim = dim;
//}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const StdRegionEntity &entity)
{
  //   os.clear();
  if (entity.nb_vert() == 0)
    return os;
  for (Uint n = 0; n < (entity.nb_vert() - 1); ++n)
    os << (entity.vertex(n)) << " ";
  os << (entity.vertex(entity.nb_vert() - 1));

  return os;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
