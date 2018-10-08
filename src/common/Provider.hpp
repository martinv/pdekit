#ifndef PDEKIT_Common_Provider_hpp
#define PDEKIT_Common_Provider_hpp

namespace pdekit
{

namespace common
{

template <typename FactoryT, typename ConcreteType>
class Provider
{
  public:
  typedef typename FactoryT::instance_type factory_instance_type;
  typedef typename factory_instance_type::key_type key_type;

  Provider(const key_type &key)
  {
    factory_instance_type &factory = FactoryT::instance();
    factory.template register_builder<ConcreteType>(key);
  }

  private:
};

} // Namespace common

} // Namespace pdekit

#endif
