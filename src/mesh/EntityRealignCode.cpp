#include "mesh/EntityRealignCode.hpp"
#include "common/StringUtils.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

const EntityRealignCode EntityRealignCode::identity(const ElemShape eshape)
{
  EntityRealignCode code(eshape, CellTransform::NO_TRANS, 0, eshape, 0, 0);
  return code;
}

// ----------------------------------------------------------------------------

const EntityRealignCode EntityRealignCode::single_flip(const ElemShape eshape)
{
  EntityRealignCode code(eshape, CellTransform::NO_TRANS, 0, eshape, 1, 0);
  return code;
}

// ----------------------------------------------------------------------------

const EntityRealignCode EntityRealignCode::single_rotation(const ElemShape eshape)
{
  EntityRealignCode code(eshape, CellTransform::NO_TRANS, 0, eshape, 0, 1);
  return code;
}

// ----------------------------------------------------------------------------

EntityRealignCode::EntityRealignCode()
{
}

// ----------------------------------------------------------------------------

EntityRealignCode::EntityRealignCode(pack_type::store_type store_value) : m_tag_data(store_value)
{
}

// ----------------------------------------------------------------------------

EntityRealignCode::EntityRealignCode(const ElemShape elem_shape, const CellTransform adapt_op_id,
                                     const Uint pos_in_parent, const ElemShape parent_shape,
                                     const Uint nb_flips, const Uint nb_rotations)
{
  m_tag_data.set_field<0>(static_cast<std::underlying_type<ElemShape>::type>(elem_shape));
  m_tag_data.set_field<1>(static_cast<std::underlying_type<CellTransform>::type>(adapt_op_id));
  m_tag_data.set_field<2>(pos_in_parent);
  m_tag_data.set_field<3>(static_cast<std::underlying_type<ElemShape>::type>(parent_shape));
  m_tag_data.set_field<4>(nb_flips);
  m_tag_data.set_field<5>(nb_rotations);
}

// ----------------------------------------------------------------------------

void EntityRealignCode::reset()
{
  m_tag_data = pack_type::store_type();
}

// ----------------------------------------------------------------------------

bool EntityRealignCode::is_identity(const ElemShape eshape) const
{
  const EntityRealignCode identity_code =
      EntityRealignCode(eshape, CellTransform::NO_TRANS, 0, eshape, 0, 0);
  return (m_tag_data.store_value() == identity_code.m_tag_data.store_value());
}

// ----------------------------------------------------------------------------

void EntityRealignCode::set_nb_flips(const Uint nb_flips)
{
  m_tag_data.set_field<4>(nb_flips);
}

// ----------------------------------------------------------------------------

void EntityRealignCode::set_nb_rotations(const Uint nb_rotations)
{
  m_tag_data.set_field<5>(nb_rotations);
}

// ----------------------------------------------------------------------------

void EntityRealignCode::add_flip()
{
  m_tag_data.set_field<4>(m_tag_data.field<4>() + 1);
}

// ----------------------------------------------------------------------------

void EntityRealignCode::add_rotation()
{
  m_tag_data.set_field<5>(m_tag_data.field<5>() + 1);
}

// ----------------------------------------------------------------------------

const std::string EntityRealignCode::as_string() const
{
  return ElemShapeInfo::name(elem_shape()) + "-" + CellTransformName::value[m_tag_data.field<1>()] +
         "-" + common::StringUtils::to_string(m_tag_data.field<2>()) + "-" +
         ElemShapeInfo::name(parent_shape()) + "_(parent)-" +
         common::StringUtils::to_string(m_tag_data.field<4>()) + "-" +
         common::StringUtils::to_string(m_tag_data.field<5>());
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, EntityRealignCode const &code)
{
  if (code.nb_rotations() == 0 && code.nb_flips() == 0)
  {
    os << "i";
    return os;
  }

  for (Uint f = 0; f < code.nb_flips(); ++f)
  {
    os << "f";
  }

  for (Uint r = 0; r < code.nb_rotations(); ++r)
  {
    os << "r";
  }

  return os;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
