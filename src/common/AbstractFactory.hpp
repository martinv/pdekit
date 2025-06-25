#ifndef PDEKIT_Common_Abstract_Factory_hpp
#define PDEKIT_Common_Abstract_Factory_hpp

#include <string>
//#include <boost/shared_ptr.hpp>

namespace pdekit
{

namespace common
{

// template<typename ProductBaseClass, typename ProductKey=std::string> class
// FactoryT;

class AbstractFactoryBase
{
  public:
  AbstractFactoryBase();

  virtual ~AbstractFactoryBase();

  /// Type name of this class
  static std::string type_name()
  {
    return "AbstractFactoryBase";
  }

  /// Get the name of the parent class name of generated objects
  virtual std::string product_type_name() const = 0;
};

} // Namespace common

} // Namespace pdekit

#endif
