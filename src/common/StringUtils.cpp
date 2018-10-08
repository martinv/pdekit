#include "common/StringUtils.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

void StringUtils::split_string(const std::string &input_str, char delim,
                               std::vector<std::string> &elems)
{
  elems.clear();
  std::stringstream ss(input_str);
  std::string item;
  while (std::getline(ss, item, delim))
  {
    elems.push_back(item);
  }
}

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit
