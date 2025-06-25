#ifndef PDEKIT_Mesh_Mesh_Predicates_hpp
#define PDEKIT_Mesh_Mesh_Predicates_hpp

#include <vector>

#include "common/Meta.hpp"
#include "mesh/MeshEntity.hpp"

namespace pdekit
{

namespace mesh
{

// ============================================================================

class AllCells
{
  public:
  typedef common::NullType ConfigParamT1;

  AllCells();

  ~AllCells();

  protected:
  inline bool is_valid(const MeshEntity &entity)
  {
    return true;
  }
};

// ============================================================================

class AllEntities
{
  public:
  typedef common::NullType ConfigParamT1;

  AllEntities();

  ~AllEntities();

  protected:
  inline bool is_valid(const MeshEntity &entity)
  {
    return true;
  }
};

// ============================================================================

class CellGroup
{
  public:
  typedef PointSetTag ConfigParamT1;

  CellGroup();

  ~CellGroup();

  inline void configure(const ConfigParamT1 etype)
  {
    m_reference_etype = etype;
  }

  inline bool is_valid(const MeshEntity &entity)
  {
    return m_reference_etype == entity.pt_set_id();
  }

  void print() const
  {
    std::cout << "GroupFilter: the reference type id = " << m_reference_etype << std::endl;
  }

  protected:
  // Id of the group for which we're looking
  ConfigParamT1 m_reference_etype;

  ConfigParamT1 m_current_etype;
};

// ============================================================================

} // namespace mesh

} // namespace pdekit

#endif
