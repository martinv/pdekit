#ifndef PDEKIT_Common_Factory_Pool_hpp
#define PDEKIT_Common_Factory_Pool_hpp

#include <map>
#include <memory>

#include "common/AbstractFactory.hpp"
#include "common/Singleton.hpp"

namespace pdekit
{

namespace common
{

class FactoryPoolImpl // : public Singleton<FactoryPoolImpl>
{
  public:
  /// Storage type for all factories
  typedef std::map<std::string, std::shared_ptr<AbstractFactoryBase>> factory_storage_type;
  typedef factory_storage_type::iterator factory_iterator;
  typedef factory_storage_type::const_iterator const_factory_iterator;

  /// Register a new factory
  void register_factory(const std::string &name, std::shared_ptr<AbstractFactoryBase> f_ptr);

  /// List all factories registered
  void list_factories() const;

  static std::string type_name()
  {
    return "FactoryPool";
  }

  private:
  factory_storage_type m_factories;
};

// ============================================================================

typedef Singleton<FactoryPoolImpl> FactoryPool;

} // Namespace common

} // Namespace pdekit

#endif
