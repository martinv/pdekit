#ifndef PDEKIT_Common_URI_hpp
#define PDEKIT_Common_URI_hpp

#include <iosfwd>
#include <string>

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

enum class URIType : short unsigned
{
  INVALID = 0,
  HTTP    = 1,
  HTTPS   = 2,
  CPATH   = 3,
  FILE    = 4
};

// ----------------------------------------------------------------------------

struct URITypeInfo
{
  static const Uint NbValues = 5;
  static const std::string Name[NbValues];
  static const URIType Value[NbValues];
};

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const URIType uri_type);

// ----------------------------------------------------------------------------

class URI
{

  public:
  /// Default constructor
  URI();

  /// Copy constructor from other path object
  /// @param path object
  URI(const URI &path);

  /// Constructor from string object
  /// @param s string with path
  URI(const std::string &s);

  /// Constructor from const char*
  /// @param c C string with path
  URI(const char *c);

  /// Constructor from string object and separate protocol
  /// @pre assumes that string does not have a protocol, just the path
  /// @param s string with path
  /// @param p scheme type e.g. (HTTP,CPATH,FILE)
  URI(const std::string &s, URIType type);

  /// Assignment operator
  URI &operator=(const URI &other_uri);

  /// comparison operator
  bool operator==(const URI &right) const;
  bool operator!=(const URI &right) const;

  /// concatenation operator with URI
  URI operator/(const URI &uri_rhs) const;

  // accessors

  /// Check if path is absolute.
  /// Should start with "//"
  /// @returns true if the path is absolute
  bool is_absolute() const;

  /// Check is path is relative.
  /// Should not start with "//"
  /// @returns true if the path is a relative path
  bool is_relative() const;

  /// check this path is complete
  /// @post true if does not contain ".." or "."
  bool is_complete() const;

  /// @return if the path is empty
  bool empty() const;

  /// @return the full URI as a string
  const std::string string() const;

  /// @return the base path (i.e. path to the parent)
  /// If the current uri is cpath:/dir1/dir2/dir3,
  /// then base_path is cpath:/dir1/dir2
  const URI base_path() const;

  /// @return the name of the object, without the path
  /// For cpath:/dir1/dir2/dir3, name() should return 'dir3'
  const std::string name() const;

  /// Gives the protocol (if any).
  /// @return Returns the protocol. May return @c URIType::INVALID if no
  /// protocol has been specified.
  URIType uri_type() const;

  /// Changes the protocol to the supplied scheme
  /// @post scheme() will return the supplied protocol
  void set_uri_type(const URIType type);

  /// Gives the URI path, which is the URI without the scheme (protocol)
  /// @return Returns the URI path
  const std::string path() const;

  /// Changes the URI path
  /// @post path() will return the supplied path
  void set_path(const std::string &path);

  /// Splits a given path into the URI protocol and the real path.
  /// @param path Path to split.
  /// @param protocol Variable where to store the found protocol. The value is
  /// set to @c URI::Protocol::INVALID if no protocol found.
  /// @param real_path Variable where to store the path.
  /// @throw ProtocolError If a an unknown protocol is found.
  static void split_path(const std::string &path, URIType &protocol, std::string &real_path);

  private:
  /// Cleans up the stored string, i.e. remove multiple / in sequence, ...
  void cleanup();

  /// Path represented by this URI
  std::string m_path;

  /// Current URI protocol
  URIType m_type;
};

} // namespace common

} // namespace pdekit

#endif
