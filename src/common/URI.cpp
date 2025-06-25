#include <ios>

#include "boost/algorithm/string.hpp"

#include "common/URI.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

const std::string URITypeInfo::Name[URITypeInfo::NbValues] = {"Invalid", "http", "https", "cpath",
                                                              "file"};

const URIType URITypeInfo::Value[URITypeInfo::NbValues] = {
    URIType::INVALID, URIType::HTTP, URIType::HTTPS, URIType::CPATH, URIType::FILE};

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const URIType uri_type)
{
  const std::underlying_type<URIType>::type uri_type_idx =
      static_cast<std::underlying_type<URIType>::type>(uri_type);
  if (uri_type_idx >= URITypeInfo::NbValues)
  {
    os << URITypeInfo::Name[0];
  }
  else
  {
    os << URITypeInfo::Name[uri_type_idx];
  }
  return os;
}

// ----------------------------------------------------------------------------

URI::URI() : m_path(), m_type(URIType::CPATH)
{
}

// ----------------------------------------------------------------------------

URI::URI(const URI &path)
{
  operator=(path);
}

// ----------------------------------------------------------------------------

URI::URI(const std::string &s) : m_path(s), m_type(URIType::CPATH)
{
  split_path(s, m_type, m_path);
  cleanup();
}

// ----------------------------------------------------------------------------

URI::URI(const char *c) : m_path(c), m_type(URIType::CPATH)
{
  std::string s(c);
  split_path(s, m_type, m_path);
  cleanup();
}

// ----------------------------------------------------------------------------

URI::URI(const std::string &s, URIType p) : m_path(s), m_type(p)
{
  /// @todo check path
  // throw NotImplemented(FromHere(), "Implement this");
  cleanup();
}

// ----------------------------------------------------------------------------

URI &URI::operator=(const URI &other_uri)
{
  m_path = other_uri.m_path;
  m_type = other_uri.m_type;
  cleanup();
  return *this;
}

// ----------------------------------------------------------------------------

bool URI::operator==(const URI &right) const
{
  return m_type == right.m_type && m_path == right.m_path;
}

// ----------------------------------------------------------------------------

bool URI::operator!=(const URI &right) const
{
  return !operator==(right);
}

// ----------------------------------------------------------------------------

URI URI::operator/(const URI &uri_rhs) const
{
  return (!m_path.empty() && !uri_rhs.m_path.empty()) ? URI(m_path + "/" + uri_rhs.m_path)
                                                      : // both not empty
             URI(m_path + uri_rhs.m_path);              // one is empty
}

// ----------------------------------------------------------------------------

bool URI::is_absolute() const
{
  return boost::algorithm::starts_with(m_path, "/");
}

// ----------------------------------------------------------------------------

bool URI::is_relative() const
{
  return !is_absolute();
}

// ----------------------------------------------------------------------------

bool URI::is_complete() const
{
  return !(boost::algorithm::starts_with(m_path, ".") || boost::algorithm::contains(m_path, "./") ||
           boost::algorithm::contains(m_path, "/."));
}

// ----------------------------------------------------------------------------

bool URI::empty() const
{
  return m_path.empty();
}

// ----------------------------------------------------------------------------

const std::string URI::string() const
{
  // if the path is not empty, we prepend the protocol
  if (!m_path.empty())
    return URITypeInfo::Name[static_cast<std::underlying_type<URIType>::type>(m_type)] + ':' +
           m_path;

  return m_path;
}

// ----------------------------------------------------------------------------

const URI URI::base_path() const
{
  using namespace boost::algorithm;

  if (m_path == "/")
    return *this;

  if (!contains(m_path, "/"))
    return URI("./", m_type);
  else
  {
    std::string rpath = m_path;
    rpath.erase(find_last(rpath, "/").begin(), rpath.end());
    if (rpath.empty())
    {
      assert(is_absolute()); // this case should only happen on
                             // first-level absolute paths such as /Model
      rpath = "/";
    }
    return URI(rpath, m_type);
  }
}

// ----------------------------------------------------------------------------

const std::string URI::name() const
{
  using namespace boost::algorithm;
  std::string name = string(); // This calls the METHOD 'string()'
  if (find_last(name, "/").begin() == name.end())
    name.erase(name.begin(), find_last(name, ":").begin() + 1);
  else
    name.erase(name.begin(), find_last(name, "/").begin() + 1);
  return name;
}

// ----------------------------------------------------------------------------

URIType URI::uri_type() const
{
  return m_type;
}

// ----------------------------------------------------------------------------

void URI::set_uri_type(const URIType type)
{
  m_type = type;
}

// ----------------------------------------------------------------------------

const std::string URI::path() const
{
  return m_path;
}

// ----------------------------------------------------------------------------

void URI::set_path(const std::string &path)
{
  m_path = path;
  cleanup();
}

// ----------------------------------------------------------------------------

void URI::split_path(const std::string &path, URIType &protocol, std::string &real_path)
{
  // by default the protocol is CPATH
  protocol  = URIType::CPATH;
  real_path = path;

  size_t colon_pos = path.find_first_of(':');

  // if the colon has been found
  if (colon_pos != std::string::npos)
  {
    // extract the procotol
    std::string protocol_str = path.substr(0, colon_pos);

    // extract the path
    real_path = real_path.substr(colon_pos + 1, path.length() - colon_pos - 1);

    // check that the protocol is valid
    for (Uint i = 0; i < URITypeInfo::NbValues; ++i)
    {
      if (URITypeInfo::Name[i] == protocol_str)
      {
        protocol = URITypeInfo::Value[i];
        break;
      }
    }
    /*
    if(protocol == URIType::INVALID)
      throw ProtocolError(FromHere(), "\'" + protocol_str + "\' is not a
    supported protocol");
    */
  }
}

// ----------------------------------------------------------------------------

void URI::cleanup()
{
  if (m_type == URIType::CPATH && !m_path.empty())
  {
    const Uint path_size = m_path.size();
    std::string cleaned_path;
    cleaned_path.reserve(path_size);
    cleaned_path.push_back(m_path[0]);
    for (Uint i = 1; i != path_size; ++i)
    {
      if (m_path[i] != '/' || m_path[i - 1] != '/')
        cleaned_path.push_back(m_path[i]);
    }
    m_path = cleaned_path;
  }
}

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit
