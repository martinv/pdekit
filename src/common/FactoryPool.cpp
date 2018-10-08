#include "common/FactoryPool.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

void FactoryPoolImpl::register_factory(const std::string &name,
                                       std::shared_ptr<AbstractFactoryBase> f_ptr)
{
  const_factory_iterator fact_iter;
  fact_iter = m_factories.find(name);

  if (fact_iter == m_factories.end())
  {
    m_factories[name] = f_ptr;
  }
}

// ----------------------------------------------------------------------------

void FactoryPoolImpl::list_factories() const
{
  std::cout << "Factory pool -> listing factories:" << std::endl;

  for (const_factory_iterator fact_iter = m_factories.begin(); fact_iter != m_factories.end();
       ++fact_iter)
  {
    std::cout << "  " << (*fact_iter).second->product_type_name() << std::endl;
  }
}

// ----------------------------------------------------------------------------

template class Singleton<FactoryPoolImpl>;

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit
