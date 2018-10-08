#ifndef PDEKIT_Common_Component_hpp
#define PDEKIT_Common_Component_hpp

#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

#include "common/URI.hpp"

namespace pdekit
{

namespace common
{

class Component : public std::enable_shared_from_this<Component>
{

  /// type for storing the sub components
  using comp_storage_type = std::vector<std::shared_ptr<Component>>;

  public:
  /// Get the class name
  static std::string type_name()
  {
    return "Component";
  }

  /// Contructor
  /// @param name of the component
  Component(const std::string &name);

  /// Virtual destructor
  virtual ~Component();

  /// Access the name of the component
  const std::string &name() const;

  /// Rename the component
  // void rename(const std::string &name);

  /// Construct the full path
  const URI uri() const;

  /// Resolves relative elements within a path to complete it.
  /// The path may be relative to this component or absolute.
  /// @param path to a component
  /// @pre path must point to an existing component
  /// @post path statisfies URI::is_complete()
  /// @post path statisfies URI::is_absolute()
  // void complete_path(URI &path) const;

  /// Looks for a component via its path
  /// @param path to the component
  /// @return handle to component or null if it doesn't exist
  /// @warning the return type is non-const!!! ( same reasoning as for
  /// parent()
  /// )
  std::shared_ptr<Component> access_component(const URI &path);

  /// @returns handle to this component
  std::shared_ptr<Component> handle();

  /// @returns handle to this component, const version
  std::shared_ptr<Component const> handle() const;

  /// @returns handle to this component cast to type ComponentT
  template <typename ComponentT>
  std::shared_ptr<ComponentT> handle();

  /// @returns handle to this component cast to type ComponentT
  template <typename ComponentT>
  std::shared_ptr<ComponentT const> handle() const;

  /// @returns the handle to the parent component, which
  /// can be null if there is no parent
  std::shared_ptr<Component> parent();

  /// @returns the upper-most component in the tree,
  /// or self if there is no parent
  /// This is const version
  std::shared_ptr<Component const> root() const;

  /// @returns the upper-most component in the tree,
  /// or self if there is no parent
  std::shared_ptr<Component> root();

  /// Gets the named child component from the list of direct subcomponents,
  /// const version
  /// @return handle to the component. Empty shared_ptr if not found.
  std::shared_ptr<Component const> get_child(const std::string &name) const;

  /// Gets the named child component from the list of direct subcomponents.
  /// @return handle to the component. Empty shared_ptr if not found.
  std::shared_ptr<Component> get_child(const std::string &name);

  /// Create a (sub)component of this component automatically cast to the
  /// specified type
  template <typename T>
  std::shared_ptr<T> create_component(const std::string &name);

  /// Add the passed component as a subcomponent
  Component &add_component(const std::shared_ptr<Component> &subcomp);

  /// @returns a string representation of the tree below this component
  /// @param [in] depth       defines how many recursions should maximally be
  /// performed
  ///                         (default value depth=0 means full tree)
  /// @param [in] level       recursion parameter, should not be touched
  std::string tree(Uint depth = 0, Uint recursion_level = 0) const;

  /// @return Returns the number of children this component has.
  size_t count_children() const;

  /// @return Returns the type name of the subclass, according to
  /// @c cf3::common::TypeInfo
  virtual std::string derived_type_name() const = 0;

  private:
  /// METHODS

  template <typename TrueType, typename T>
  struct HandleConverter;

  template <typename T>
  struct HandleConverter<std::true_type, T>
  {
    static std::shared_ptr<T> apply(const std::shared_ptr<Component> &ptr)
    {
      // std::cout << "TRUE TYPE" << std::endl;
      // Don't use 'static_cast' here, but static_pointer_cast, because
      // we're casting pointers that are trying to manage the same
      // resource!
      return std::shared_ptr<T>(std::static_pointer_cast<T>(ptr));
    }
  };

  template <typename T>
  struct HandleConverter<std::false_type, T>
  {
    static std::shared_ptr<T> apply(const std::shared_ptr<Component> &ptr)
    {
      // std::cout << "FALSE TYPE" << std::endl;
      // Don't use 'dynamic_cast' here, but dynamic_pointer_cast, because
      // we're casting pointers that are trying to manage the same
      // resource!
      return std::shared_ptr<T>(std::dynamic_pointer_cast<T>(ptr));
    }
  };

  /// DATA
  /// Name of this component
  std::string m_name;

  /// Naked pointer to parent
  std::weak_ptr<Component> m_parent;

  /// Storage of children
  comp_storage_type m_child_components;
};

// ----------------------------------------------------------------------------
//             FREE FUNCTION TO CREATE COMPONENT
// ----------------------------------------------------------------------------

template <typename ComponentT>
std::shared_ptr<ComponentT> make_component(const std::string &name)
{
  std::shared_ptr<ComponentT> new_component = std::make_shared<ComponentT>(name);
  return new_component;
}

// ----------------------------------------------------------------------------

template <typename ComponentT>
std::shared_ptr<ComponentT> Component::handle()
{
  return HandleConverter<typename std::is_base_of<Component, ComponentT>::type, ComponentT>::apply(
      shared_from_this());
  // return
  // std::shared_ptr<ComponentT>(std::dynamic_pointer_cast<ComponentT>(shared_from_this()));
}

// ----------------------------------------------------------------------------

template <typename ComponentT>
std::shared_ptr<ComponentT const> Component::handle() const
{

  return HandleConverter<typename std::is_base_of<Component, ComponentT>::type,
                         ComponentT const>::apply(shared_from_this());
  // return std::shared_ptr<ComponentT const>(shared_from_this());
}

// ----------------------------------------------------------------------------

template <typename T>
std::shared_ptr<T> Component::create_component(const std::string &name)
{
  std::shared_ptr<T> comp = make_component<T>(name);
  add_component(comp);
  return std::shared_ptr<T>(comp);
}

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
