#ifndef PDEKIT_RDM_BC_Factory_hpp
#define PDEKIT_RDM_BC_Factory_hpp

#include <string>

#include "common/FactoryPool.hpp"
#include "common/FactoryT.hpp"
#include "common/Singleton.hpp"

#include "physics/PhysModelT.hpp"

#include "solver/rdm/bc/StrongDirichletBC.hpp"
#include "solver/rdm/bc/WeakDirichletBC.hpp"
#include "solver/rdm/bc/WeakFarfield.hpp"
#include "solver/rdm/bc/WeakSubInlet.hpp"
#include "solver/rdm/bc/WeakSubOutlet.hpp"
#include "solver/rdm/bc/WeakSuperInlet.hpp"
#include "solver/rdm/bc/WeakSymmetry.hpp"
#include "solver/rdm/bc/WeakWall.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

namespace detail
{

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM, typename PhysicsClass>
class RDMBcFactorySetup;

// ----------------------------------------------------------------------------
// Specialization of Boundary Condition Factory Setup for BC factory for
// scalar problems
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
class RDMBcFactorySetup<MeshConfig, Physics, BcDIM, physics::PhysicsClassScalar>
{
  public:
  static void setup_instance(
      common::FactoryT<RDMBCBase<MeshConfig, Physics, BcDIM>, std::string> &factory);
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBcFactorySetup<MeshConfig, Physics, BcDIM, physics::PhysicsClassScalar>::setup_instance(
    common::FactoryT<RDMBCBase<MeshConfig, Physics, BcDIM>, std::string> &factory)
{
  factory.template register_builder<StrongDirichletBC<MeshConfig, Physics, BcDIM>>(
      "StrongDirichlet");
  factory.template register_builder<WeakDirichletBC<MeshConfig, Physics, BcDIM>>("WeakDirichlet");
}

// ----------------------------------------------------------------------------
// Specialization of Boundary Condition Factory Setup for BC factory for
// Euler
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
class RDMBcFactorySetup<MeshConfig, Physics, BcDIM, physics::PhysicsClassEuler>
{
  public:
  static void setup_instance(
      common::FactoryT<RDMBCBase<MeshConfig, Physics, BcDIM>, std::string> &factory);
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBcFactorySetup<MeshConfig, Physics, BcDIM, physics::PhysicsClassEuler>::setup_instance(
    common::FactoryT<RDMBCBase<MeshConfig, Physics, BcDIM>, std::string> &factory)
{
  factory.template register_builder<StrongDirichletBC<MeshConfig, Physics, BcDIM>>(
      "StrongDirichlet");
  factory.template register_builder<WeakDirichletBC<MeshConfig, Physics, BcDIM>>("WeakDirichlet");
  factory.template register_builder<WeakWall<MeshConfig, Physics, BcDIM>>("WeakWall");
  factory.template register_builder<WeakFarfield<MeshConfig, Physics, BcDIM>>("WeakFarfield");
  factory.template register_builder<WeakSubInlet<MeshConfig, Physics, BcDIM>>("WeakSubInlet");
  factory.template register_builder<WeakSubOutlet<MeshConfig, Physics, BcDIM>>("WeakSubOutlet");
  factory.template register_builder<WeakSuperInlet<MeshConfig, Physics, BcDIM>>("WeakSuperInlet");
  factory.template register_builder<WeakSymmetry<MeshConfig, Physics, BcDIM>>("WeakSymmetry");
}

} // namespace detail

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM = Physics::DIM - 1>
class RDMBCManager
{
  public:
  typedef RDMBCBase<MeshConfig, Physics, BcDIM> bc_base;
  typedef typename std::map<std::string, std::shared_ptr<bc_base>>::const_iterator const_iterator;

  /// Default constructor, no parameters
  RDMBCManager();

  /// Default destructor
  ~RDMBCManager();

  /// Clear all boundary conditions
  void clear_all_bcs();

  /// Create a new boundary condition with given name
  std::shared_ptr<bc_base> create_bc(const std::string &bc_type_name, const std::string &name);

  /// Get a boundary condition with given name
  /// @param bc_name ... name of the boundary condition we want to retrieve
  std::shared_ptr<bc_base> get_bc(const std::string &bc_name) const;

  /// Remove a boundary condition with given name
  void remove_bc(const std::string &bc_name);

  /// Return const iterator to the beginning of the map container
  const_iterator cbegin();

  /// Return const iterator to the beginning of the map container, const
  /// version
  const_iterator cbegin() const;

  /// Return const iterator to the end of the map container
  const_iterator cend();

  /// Return const iterator to the end of the map container, const version
  const_iterator cend() const;

  /// Print all boundary conditions
  void print_all_bcs() const;

  private:
  /// TYPEDEFS
  typedef common::Singleton<
      common::FactoryT<bc_base, std::string>,
      detail::RDMBcFactorySetup<MeshConfig, Physics, BcDIM, typename Physics::physics_class>>
      bc_factory_type;

  /// Map of all boundary conditions
  std::map<std::string, std::shared_ptr<bc_base>> m_bc_map;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
RDMBCManager<MeshConfig, Physics, BcDIM>::RDMBCManager()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
RDMBCManager<MeshConfig, Physics, BcDIM>::~RDMBCManager()
{
  clear_all_bcs();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCManager<MeshConfig, Physics, BcDIM>::clear_all_bcs()
{
  m_bc_map.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
std::shared_ptr<typename RDMBCManager<MeshConfig, Physics, BcDIM>::bc_base> RDMBCManager<
    MeshConfig, Physics, BcDIM>::create_bc(const std::string &bc_type_name, const std::string &name)
{
  // First check if a boundary condition with the same name does not exist
  // already
  typename std::map<std::string, std::shared_ptr<bc_base>>::const_iterator it = m_bc_map.find(name);
  if (it != m_bc_map.end())
  {
    std::cerr << "RDMBCManager::create_bc::BC with name \"" << name << "\" already exists. ";
    std::cerr << "Not creating a new BC." << std::endl;
    return std::shared_ptr<bc_base>();
  }

  // Use factory to create the bc
  typename bc_factory_type::instance_type &bc_factory = bc_factory_type::instance();

  // Here we are converting a unique_ptr returned by the factor to a
  // shared_ptr needed by 'm_bc_map'
  std::shared_ptr<bc_base> boundary_cond = bc_factory.create(bc_type_name);
  boundary_cond->set_name(name);

  // Register it in local database - should RDMBCManager be a singleton ???
  m_bc_map.insert(std::pair<std::string, std::shared_ptr<bc_base>>(name, boundary_cond));

  // Return the new boundary condition
  return boundary_cond;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
std::shared_ptr<typename RDMBCManager<MeshConfig, Physics, BcDIM>::bc_base> RDMBCManager<
    MeshConfig, Physics, BcDIM>::get_bc(const std::string &bc_name) const
{
  typename std::map<std::string, std::shared_ptr<bc_base>>::const_iterator it =
      m_bc_map.find(bc_name);
  if (it != m_bc_map.end())
  {
    return it->second;
  }
  return std::shared_ptr<bc_base>();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCManager<MeshConfig, Physics, BcDIM>::remove_bc(const std::string &bc_name)
{
  m_bc_map.erase(bc_name);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
typename RDMBCManager<MeshConfig, Physics, BcDIM>::const_iterator RDMBCManager<MeshConfig, Physics,
                                                                               BcDIM>::cbegin()
{
  return m_bc_map.cbegin();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
typename RDMBCManager<MeshConfig, Physics, BcDIM>::const_iterator RDMBCManager<
    MeshConfig, Physics, BcDIM>::cbegin() const
{
  return m_bc_map.cbegin();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
typename RDMBCManager<MeshConfig, Physics, BcDIM>::const_iterator RDMBCManager<MeshConfig, Physics,
                                                                               BcDIM>::cend()
{
  return m_bc_map.cend();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
typename RDMBCManager<MeshConfig, Physics, BcDIM>::const_iterator RDMBCManager<MeshConfig, Physics,
                                                                               BcDIM>::cend() const
{
  return m_bc_map.cend();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDIM>
void RDMBCManager<MeshConfig, Physics, BcDIM>::print_all_bcs() const
{
  Uint counter = 1u;
  std::cout << "BCManager: total " << m_bc_map.size() << " boundary conditions" << std::endl;
  for (typename std::map<std::string, std::shared_ptr<bc_base>>::const_iterator it =
           m_bc_map.begin();
       it != m_bc_map.end(); ++it)
  {
    std::cout << counter << ") " << it->first << std::endl;
    counter++;
  }
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
