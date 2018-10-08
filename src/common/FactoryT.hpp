#ifndef PDEKIT_Common_FactoryT_hpp
#define PDEKIT_Common_FactoryT_hpp

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>

#include "common/FactoryPool.hpp"
#include "common/IOStreamUtils.hpp"
#include "common/PtrHandle.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

namespace detail
{
template <typename ProductBaseClass, typename ProductKeyType>
class ProductBuilderBase
{

  public:
  ProductBuilderBase(const ProductKeyType &key) : m_key(key)
  {
  }

  virtual ~ProductBuilderBase()
  {
  }

  virtual typename std::unique_ptr<ProductBaseClass> create() const = 0;

  virtual ProductKeyType const &key() const
  {
    return m_key;
  }

  protected:
  ProductKeyType const m_key;
};

// ----------------------------------------------------------------------------
// CONCRETE PRODUCT BUILDER WHICH CALLS DEFAULT CONSTRUCTOR (NO ARGUMENTS)
// OF THE CONCRETE PRODUCT
// ----------------------------------------------------------------------------

template <typename ProductBaseClass, typename ConcreteProduct, typename ProductKeyType>
class ConcreteProductBuilder : public ProductBuilderBase<ProductBaseClass, ProductKeyType>
{
  private:
  using builder_base_type = ProductBuilderBase<ProductBaseClass, ProductKeyType>;

  public:
  /// Constructor
  explicit ConcreteProductBuilder(ProductKeyType const &key) : builder_base_type(key)
  {
  }

  /// Destructor
  ~ConcreteProductBuilder() override
  {
  }

  /// Create a concrete type
  typename std::unique_ptr<ProductBaseClass> create() const override
  {
    typename std::unique_ptr<ProductBaseClass> ptr(new ConcreteProduct());
    return ptr;
  }
};

// ----------------------------------------------------------------------------
// CONCRETE PRODUCT BUILDER WHICH CALLS CONSTRUCTOR
// OF THE CONCRETE PRODUCT WITH EXTRA ARGUMENTS
// ----------------------------------------------------------------------------

template <typename ProductBaseClass, typename ConcreteProduct, typename ProductKeyType,
          typename... Args>
class ConcreteProductBuilderWParams : public ProductBuilderBase<ProductBaseClass, ProductKeyType>
{
  private:
  using builder_base_type = ProductBuilderBase<ProductBaseClass, ProductKeyType>;

  public:
  /// Constructor
  explicit ConcreteProductBuilderWParams(ProductKeyType const &key, const Args &... args)
      : builder_base_type(key), m_construction_params(std::make_tuple(args...))
  {
  }

  /// Destructor
  ~ConcreteProductBuilderWParams() override
  {
  }

  /// Create a concrete type
  typename std::unique_ptr<ProductBaseClass> create() const override
  {
    typename std::unique_ptr<ProductBaseClass> ptr(new ConcreteProduct(m_construction_params));
    return ptr;
  }

  private:
  std::tuple<Args...> m_construction_params;
};

} // namespace detail

// ----------------------------------------------------------------------------

template <typename ProductBaseClass, typename ProductKeyType = std::string>
class FactoryT : public AbstractFactoryBase
{
  private:
  /// FORWARD DECLARATIONS

  public:
  /// TYPEDEFS

  using product_base_ptr       = std::unique_ptr<ProductBaseClass>;
  using const_product_base_ptr = std::unique_ptr<ProductBaseClass>;
  using key_type               = ProductKeyType;

  /// Type name of this class
  static std::string type_name()
  {
    return "FactoryT<" + ProductBaseClass::type_name() + ">";
  }

  /// Constructor
  FactoryT();

  /// Destructor
  ~FactoryT() override;

  /// Get the name of the parent class name of generated objects
  std::string product_type_name() const override
  {
    return ProductBaseClass::type_name();
  }

  /// Register a new concrete builder
  template <typename ConcreteType>
  void register_builder(const key_type &key);

  template <typename ConcreteType, typename... Args>
  void register_builder(const key_type &key, Args const &... args);

  /// Create a concrete object with given description
  product_base_ptr create(const key_type &key);

  /// Check if key is registered
  bool key_is_registered(const key_type &key) const;

  /// List all builders available
  void print_builders() const;

  private:
  /// TYPEDEFS
  using builder_base = typename detail::ProductBuilderBase<ProductBaseClass, ProductKeyType>;
  using builder_storage_type = std::map<key_type, builder_base *>;

  /// MEMBER VARIABLES

  /// A storage for all available builders in this factory
  builder_storage_type m_builders;
};

// ----------------------------------------------------------------------------
// METHODS OF FactoryT
// ----------------------------------------------------------------------------

template <typename ProductBaseClass, typename ProductKeyType>
FactoryT<ProductBaseClass, ProductKeyType>::FactoryT() : AbstractFactoryBase()
{
}

// ----------------------------------------------------------------------------

template <typename ProductBaseClass, typename ProductKeyType>
FactoryT<ProductBaseClass, ProductKeyType>::~FactoryT()
{
  typename builder_storage_type::iterator iter;
  for (iter = m_builders.begin(); iter != m_builders.end(); ++iter)
  {
    delete iter->second;
  }
}

// ----------------------------------------------------------------------------

template <typename ProductBaseClass, typename ProductKeyType>
template <typename ConcreteType>
void FactoryT<ProductBaseClass, ProductKeyType>::register_builder(const key_type &key)
{
  if (m_builders.size() > 0)
  {
    if (key_is_registered(key))
    {
      std::cerr << "Error, there's already a builder with the key " << key << " registered"
                << std::endl;
      return;
    }
  }

  using concrete_builder =
      typename detail::ConcreteProductBuilder<ProductBaseClass, ConcreteType, ProductKeyType>;
  builder_base *builder_ptr = new concrete_builder(key);
  m_builders.insert(std::pair<key_type, builder_base *>(key, builder_ptr));
}

// ----------------------------------------------------------------------------

template <typename ProductBaseClass, typename ProductKeyType>
template <typename ConcreteType, typename... Args>
void FactoryT<ProductBaseClass, ProductKeyType>::register_builder(const key_type &key,
                                                                  Args const &... args)
{
  if (m_builders.size() > 0)
  {
    if (key_is_registered(key))
    {
      std::cerr << "Error, there's already a builder with the key " << key << " registered"
                << std::endl;
      return;
    }
  }

  using concrete_builder =
      typename detail::ConcreteProductBuilderWParams<ProductBaseClass, ConcreteType, ProductKeyType,
                                                     Args...>;
  builder_base *builder_ptr = new concrete_builder(key, std::forward<Args>(args)...);
  // builder_base *builder_ptr = new concrete_builder(key, args...);
  m_builders.insert(std::pair<key_type, builder_base *>(key, builder_ptr));
}

// ----------------------------------------------------------------------------

template <typename ProductBaseClass, typename ProductKeyType>
typename FactoryT<ProductBaseClass, ProductKeyType>::product_base_ptr FactoryT<
    ProductBaseClass, ProductKeyType>::create(const ProductKeyType &key)
{
  typename builder_storage_type::const_iterator builder_iter;

  builder_iter = m_builders.find(key);
  if (builder_iter == m_builders.end())
  {
    std::cerr << type_name() << " : didn't find a builder for class with key " << key << std::endl;
    product_base_ptr null_ptr;
    return null_ptr;
  }

  return builder_iter->second->create();
}

// ----------------------------------------------------------------------------

template <typename ProductBaseClass, typename ProductKeyType>
bool FactoryT<ProductBaseClass, ProductKeyType>::key_is_registered(const key_type &key) const
{
  typename builder_storage_type::const_iterator builder_iter;
  builder_iter = m_builders.find(key);
  return (builder_iter != m_builders.end());
}

// ----------------------------------------------------------------------------

template <typename ProductBaseClass, typename ProductKeyType>
void FactoryT<ProductBaseClass, ProductKeyType>::print_builders() const
{
  typename builder_storage_type::const_iterator builder_iter;
  for (builder_iter = m_builders.begin(); builder_iter != m_builders.end(); ++builder_iter)
  {
    std::cout << (*builder_iter).second->key() << std::endl;
  }
  std::cout << "Total number of builders = " << m_builders.size() << std::endl;
}

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit

#endif
