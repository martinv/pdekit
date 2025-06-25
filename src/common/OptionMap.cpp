#include <fstream>
#include <iostream>
#include <tuple>

#include "common/OptionMap.hpp"
#include "common/StringUtils.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

/*
const std::vector<std::string> OptionMap::registered_type_names = {
  "Real", "Float", "Int", "Uint", "String", "StringTuple2", "StringTuple3"
};
*/

// ----------------------------------------------------------------------------

OptionMap::~OptionMap()
{
  clear();
}

// ----------------------------------------------------------------------------

void OptionMap::clear()
{
  for (size_t i = 0; i < m_options.size(); ++i)
  {
    delete m_options[i];
  }
  m_options.resize(0);
}

// ----------------------------------------------------------------------------

void OptionMap::initialize_from_file(const std::string &filename)
{
  clear();
  std::ifstream ifile;
  ifile.open(filename.c_str());

  std::string buffer;
  std::vector<std::string> words;
  std::string option_name;
  std::string option_type;

  while (getline(ifile, buffer))
  {
    // Lines starting with '#' are considered comments
    if (buffer.front() != '#')
    {
      StringUtils::split_string(buffer, ' ', words);

      if (words.size() < 2)
      {
        std::cout << "OptionMap::initialize_from_file: line " << buffer << " is malformatted"
                  << std::endl;
        break;
      }
      option_name = words[0];
      option_type = words[1];

      if (option_type == "Real")
      {
        this->create<Real>(option_name);
      }
      else if (option_type == "Float")
      {
        this->create<Float>(option_name);
      }
      else if (option_type == "Int")
      {
        this->create<Int>(option_name);
      }
      else if (option_type == "Uint")
      {
        this->create<Uint>(option_name);
      }
      else if (option_type == "String")
      {
        this->create<std::string>(option_name);
      }
      else if (option_type == "StringTuple2")
      {
        using value_type = std::tuple<std::string, std::string>;
        this->create<value_type>(option_name);
      }
      else if (option_type == "StringTuple3")
      {
        using value_type = std::tuple<std::string, std::string, std::string>;
        this->create<value_type>(option_name);
      }

      else
      {
        std::cerr << "OptionMap::initalize_from_file: option type \"" << option_type
                  << "\" not recognized" << std::endl;
      }
    }
  }

  ifile.close();
}

// ----------------------------------------------------------------------------

void OptionMap::configure_from_file(const std::string &filename)
{
  std::ifstream ifile;
  ifile.open(filename.c_str());

  std::string buffer;
  std::vector<std::string> words;
  std::string option_name;
  std::string option_type;

  while (getline(ifile, buffer))
  {
    // Lines starting with '#' are considered comments
    if (buffer.front() != '#')
    {
      StringUtils::split_string(buffer, ' ', words);

      if (words.size() < 3)
      {
        std::cout << "OptionMap::configure_from_file: line " << buffer << " is malformatted"
                  << std::endl;
        break;
      }
      option_name = words[0];
      option_type = words[1];

      if (option_type == "Real")
      {
        const Real value = StringUtils::from_string<Real>(words[2]);
        this->set<Real>(option_name, value);
      }
      else if (option_type == "Float")
      {
        const Float value = StringUtils::from_string<Float>(words[2]);
        this->set<Float>(option_name, value);
      }
      else if (option_type == "Int")
      {
        const Int value = StringUtils::from_string<Int>(words[2]);
        this->set<Int>(option_name, value);
      }
      else if (option_type == "Uint")
      {
        const Uint value = StringUtils::from_string<Uint>(words[2]);
        this->set<Uint>(option_name, value);
      }
      else if (option_type == "String")
      {
        this->set<std::string>(option_name, words[2]);
      }
      else if (option_type == "StringTuple2")
      {
        using value_type       = std::tuple<std::string, std::string>;
        const value_type value = value_type(words[2], words[3]);
        this->set<value_type>(option_name, value_type(words[2], words[3]));
      }
      else if (option_type == "StringTuple3")
      {
        using value_type       = std::tuple<std::string, std::string, std::string>;
        const value_type value = value_type(words[2], words[3], words[4]);
        this->set<value_type>(option_name, value_type(words[2], words[3], words[4]));
      }

      else
      {
        std::cerr << "OptionMap::configure_from_file: option type \"" << option_type
                  << "\" not recognized" << std::endl;
      }
    }
  }

  ifile.close();
}

// ----------------------------------------------------------------------------

void OptionMap::print() const
{
  for (detail::OptionArrayBase *base_ptr : m_options)
  {
    base_ptr->print();
    std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------

/*
void OptionMap::test(const std::string &name) const
{
  DataHelper<std::tuple_size<registered_types>::value, 0>::get_pos(name,
registered_type_names);
}
*/

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit
