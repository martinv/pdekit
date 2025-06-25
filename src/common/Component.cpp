#include <boost/algorithm/string.hpp>

#include "common/Component.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

Component::Component(const std::string &name) : m_name(name)
{
  m_child_components.resize(0);
}

// ----------------------------------------------------------------------------

Component::~Component()
{
}

// ----------------------------------------------------------------------------

const std::string &Component::name() const
{
  return m_name;
}

// ----------------------------------------------------------------------------

const URI Component::uri() const
{
  // First create shared pointer
  auto sptr = m_parent.lock();

  if (!sptr)
  {
    return URI(std::string("/"), URIType::CPATH);
  }
  return sptr->uri() / URI(name(), URIType::CPATH);
}

// ----------------------------------------------------------------------------

std::shared_ptr<Component> Component::access_component(const URI &path)
{
  // Return self for trivial path or at end of recursion.
  if (path.path() == "." || path.empty())
    return shared_from_this();

  // If the path is absolute, make it relative and pass it to the root
  if (path.is_absolute())
  {
    // String without protocol
    std::string new_path = path.path();

    // Remove any leading /
    boost::algorithm::trim_left_if(new_path, boost::algorithm::is_any_of("/"));

    if (new_path.empty())
    {
      return root();
    }

    // Pass the rest to root
    return root()->access_component(URI(new_path, URIType::CPATH));
  }

  // Relative path
  std::string path_str = path.path();
  // Remove trailing /
  boost::algorithm::trim_right_if(path_str, boost::algorithm::is_any_of("/"));

  // Find the first separator
  const std::size_t first_sep = path_str.find("/");
  const bool has_no_separator = (first_sep == std::string::npos);

  // Current part of the path to parse
  const std::string current_part = has_no_separator ? path_str : path_str.substr(0, first_sep);

  // Remainder of the path to parse, set to "." if there were no more parts
  const std::string next_part =
      has_no_separator ? "." : path_str.substr(first_sep + 1, path_str.size());

  // Dispatch to self
  if (current_part == "." || current_part.empty())
    return access_component(next_part);

  // Dispatch to parent
  if (current_part == "..")
  {
    std::shared_ptr<Component> shared_ptr_parent = m_parent.lock();
    return shared_ptr_parent ? shared_ptr_parent->access_component(next_part)
                             : std::shared_ptr<Component>();
  }

  // Dispatch to child
  std::shared_ptr<Component> child = get_child(current_part);
  if (child)
  {
    return child->access_component(next_part);
  }

  // Return null if not found
  return std::shared_ptr<Component>();
}

// ----------------------------------------------------------------------------

std::shared_ptr<Component> Component::handle()
{
  return shared_from_this();
}

// ----------------------------------------------------------------------------

std::shared_ptr<Component const> Component::handle() const
{
  return std::shared_ptr<Component const>(shared_from_this());
}

// ----------------------------------------------------------------------------

std::shared_ptr<Component> Component::parent()
{
  std::shared_ptr<Component> shared_ptr_parent = m_parent.lock();
  return shared_ptr_parent;
}

// ----------------------------------------------------------------------------

std::shared_ptr<const Component> Component::root() const
{
  std::shared_ptr<Component const> result            = shared_from_this();
  std::shared_ptr<Component const> shared_ptr_parent = m_parent.lock();

  while (shared_ptr_parent)
  {
    result            = shared_ptr_parent;
    shared_ptr_parent = result->m_parent.lock();
  }
  return result;
}

// ----------------------------------------------------------------------------

std::shared_ptr<Component> Component::root()
{
  std::shared_ptr<Component> result            = shared_from_this();
  std::shared_ptr<Component> shared_ptr_parent = m_parent.lock();

  while (shared_ptr_parent)
  {
    result            = shared_ptr_parent;
    shared_ptr_parent = result->m_parent.lock();
  }
  return result;
}

// ----------------------------------------------------------------------------

std::shared_ptr<Component const> Component::get_child(const std::string &name) const
{
  for (Uint i = 0; i < m_child_components.size(); ++i)
  {
    if (m_child_components[i]->name() == name)
    {
      std::shared_ptr<Component const> child_ptr(m_child_components[i]);
      return child_ptr;
    }
  }
  return std::shared_ptr<Component const>();
}

// ----------------------------------------------------------------------------

std::shared_ptr<Component> Component::get_child(const std::string &name)
{
  for (Uint i = 0; i < m_child_components.size(); ++i)
  {
    if (m_child_components[i]->name() == name)
    {
      return m_child_components[i];
    }
  }
  return std::shared_ptr<Component>();
}

// ----------------------------------------------------------------------------

Component &Component::add_component(const std::shared_ptr<Component> &subcomp)
{
  m_child_components.push_back(subcomp); // add to all component list
  subcomp->m_parent = shared_from_this();
  return *subcomp;
}

// ----------------------------------------------------------------------------

std::string Component::tree(Uint depth, Uint recursion_level) const
{
  std::string tree;
  if (recursion_level <= depth || depth == 0)
  {
    for (Uint i = 0; i < recursion_level; i++)
      tree += "  ";
    tree += name();

    tree += " (" + derived_type_name() + ")\n";

    /*
    for(const Component &c : *this)
    {
      tree += c.tree(depth, recursion_level + 1);
    }
    */

    for (Uint i = 0; i < m_child_components.size(); ++i)
    {
      tree += m_child_components[i]->tree(depth, recursion_level + 1);
    }
  }
  return tree;
}

// ----------------------------------------------------------------------------

size_t Component::count_children() const
{
  return m_child_components.size();
}

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit
