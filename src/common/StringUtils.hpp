#ifndef PDEKIT_Common_StringUtils_hpp
#define PDEKIT_Common_StringUtils_hpp

#include <sstream>
#include <string>
#include <vector>

namespace pdekit
{

namespace common
{

class StringUtils
{
  public:
  /// Split a string into a vector of strings given a delimiter
  /// @param input_str ... input string that should be split
  /// @param delim     ... delimiter
  /// @return elemds   ... vector of strings split based on delim
  static void split_string(const std::string &input_str, char delim,
                           std::vector<std::string> &elems);

  // --------------------------------------------------------------------------

  /// Convert an instance of type T (can be integer, for example) to string
  /// @param t ... value that should be converted to string
  /// @return  ... string representation of t
  template <typename T>
  static std::string to_string(const T &t)
  {
    std::ostringstream os;
    os << t;
    return os.str();
  }

  // --------------------------------------------------------------------------

  /// Convert a string to value of given type
  /// For example, we want to convert the string "123"
  /// to integer 123 - we can write
  /// const Uint n = StringUtils::from_string<Uint>("123");
  /// @param s ... string that should be converted to a value of type T
  /// @return  ... a value of type t
  template <typename T>
  static T from_string(const std::string &s)
  {
    std::istringstream is(s);
    T t;
    is >> t;
    return t;
  }
};

// --------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
