#ifndef PDEKIT_Common_Options_hpp
#define PDEKIT_Common_Options_hpp

#include <vector>

#include "common/Meta.hpp"
#include "common/PDEKit.hpp"

namespace pdekit
{

namespace common
{

namespace detail
{

// ----------------------------------------------------------------------------

class OptionPairBase
{
  public:
  OptionPairBase()
  {
  }

  virtual ~OptionPairBase()
  {
  }

  virtual void print() const = 0;

  private:
};

// ----------------------------------------------------------------------------

template <typename T>
class OptionPair : public OptionPairBase
{
  public:
  // Default constructor
  OptionPair() : OptionPairBase()
  {
  }

  // Default constructor
  OptionPair(const std::string &name, const T &value) : OptionPairBase()
  {
    m_name = name;
    m_val  = value;
  }

  // Copy constructor
  OptionPair(const OptionPair &rhs)
  {
    m_name = rhs.m_name;
    m_val  = rhs.m_val;
  }

  // Assignment operator
  OptionPair &operator=(const OptionPair &rhs)
  {
    m_name = rhs.m_name;
    m_val  = rhs.m_val;
    return *this;
  }

  // Destructor
  ~OptionPair() override
  {
  }

  // Return option name
  const std::string &name() const
  {
    return m_name;
  }

  // Return option value
  const T &value() const
  {
    return m_val;
  }

  void print() const override
  {
    std::cout << m_name << ": ";
    to_stream(std::cout, m_val);
    std::cout << std::endl;
  }

  private:
  template <typename U>
  static std::ostream &to_stream(std::ostream &os, const U &value)
  {
    os << value;
    return os;
  }

  /// Overload printing of tuples
  template <typename U1, typename U2>
  static std::ostream &to_stream(std::ostream &os, const std::tuple<U1, U2> &tuple)
  {
    os << "[" << std::get<0>(tuple) << "," << std::get<1>(tuple) << "]";
    return os;
  }

  template <typename U1, typename U2, typename U3>
  static std::ostream &to_stream(std::ostream &os, const std::tuple<U1, U2, U3> &tuple)
  {
    os << "[" << std::get<0>(tuple) << "," << std::get<1>(tuple) << "," << std::get<2>(tuple)
       << "]";
    return os;
  }

  // Option name
  std::string m_name;

  // Actual value
  T m_val;
};

// ----------------------------------------------------------------------------

class OptionArrayBase
{
  public:
  OptionArrayBase() = default;

  virtual ~OptionArrayBase() = default;

  virtual void insert_value(const OptionPairBase &value) = 0;

  virtual void print() const = 0;
};

// ----------------------------------------------------------------------------

template <typename T>
class OptionArray : public OptionArrayBase
{
  public:
  // Default constructor
  OptionArray() : OptionArrayBase()
  {
  }

  // Copy constructor
  OptionArray(const OptionArray &rhs)
  {
    m_names.resize(rhs.m_names.size());
    m_names = rhs.m_names;

    m_vals.resize(rhs.m_vals.size());
    m_vals = rhs.m_vals;
  }

  // Assignment operator
  OptionArray &operator=(const OptionArray &rhs)
  {
    m_names.resize(rhs.m_names.size());
    m_names = rhs.m_names;

    m_vals.resize(rhs.m_vals.size());
    m_vals = rhs.m_vals;
    return *this;
  }

  // Destructor
  ~OptionArray() override
  {
  }

  // Return the size of the array
  size_t size() const
  {
    return m_names.size();
  }

  // Return the i-th option name
  const std::string &name(int i) const
  {
    return m_names[i];
  }

  // Return the i-th option value
  const T &value(int i) const
  {
    return m_vals[i];
  }

  // Change value of option with given position
  void change_value(const int i, const T &new_value)
  {
    m_vals[i] = new_value;
  }

  void insert_value(const OptionPairBase &value) override
  {
    const OptionPair<T> *concrete_pair = dynamic_cast<const OptionPair<T> *>(&value);
    if (concrete_pair)
    {
      // std::cout << "Dynamic cast successful." << std::endl;
      m_names.push_back(concrete_pair->name());
      m_vals.push_back(concrete_pair->value());
    }
    else
    {
      // std::cout << "Dynamic cast failed." << std::endl;
    }
  }

  void print() const override
  {
    for (size_t i = 0; i < m_names.size(); ++i)
    {
      std::cout << m_names[i] << ": ";
      to_stream(std::cout, m_vals[i]);
      std::cout << std::endl;
    }
  }

  private:
  template <typename U>
  static std::ostream &to_stream(std::ostream &os, const U &value)
  {
    os << value;
    return os;
  }

  /// Overload printing of tuples
  template <typename U1, typename U2>
  static std::ostream &to_stream(std::ostream &os, const std::tuple<U1, U2> &tuple)
  {
    os << "[" << std::get<0>(tuple) << "," << std::get<1>(tuple) << "]";
    return os;
  }

  template <typename U1, typename U2, typename U3>
  static std::ostream &to_stream(std::ostream &os, const std::tuple<U1, U2, U3> &tuple)
  {
    os << "[" << std::get<0>(tuple) << "," << std::get<1>(tuple) << "," << std::get<2>(tuple)
       << "]";
    return os;
  }

  // Option name
  std::vector<std::string> m_names;

  // Actual value
  std::vector<T> m_vals;
};

// ----------------------------------------------------------------------------

} // namespace detail

// ----------------------------------------------------------------------------

class OptionMap
{
  public:
  /// Default constructor
  OptionMap() = default;

  /// Delete copy constructor
  OptionMap(const OptionMap &other) = delete;

  /// Destructor
  ~OptionMap();

  /// Deleted assignment operator
  OptionMap &operator=(const OptionMap &rhs) = delete;

  /// Create option and set value
  template <typename T>
  void create(const std::string &name, const T &value);

  /// Create option
  template <typename T>
  void create(const std::string &name);

  /// Get option value
  template <typename T>
  T get(const std::string &name) const;

  /// Set option value
  template <typename T>
  void set(const std::string &name, const T &value);

  /// Clear the map - remove all options
  void clear();

  /// Initialize from file - fills the map with options and their
  /// Default values
  void initialize_from_file(const std::string &filename);

  /// Initialize from file - fills the map with options and their
  /// Default values
  void configure_from_file(const std::string &filename);

  /// Print all options
  void print() const;

  // void test(const std::string &name) const;

  private:
  /*
  /// TYPES
  using registered_types =
      std::tuple<Real, Float, Int, Uint, std::string, std::tuple<std::string,
  std::string>, std::tuple<std::string, std::string, std::string>>;

  static const std::vector<std::string> registered_type_names;

  /// METHODS

  template <Uint Size, Uint Pos>
  struct DataHelper
  {
    static void get_pos(const std::string &type_name,
                        const std::vector<std::string> &registered_names)
    {
      if (registered_names[Pos] == type_name)
      {
        std::cout << "Name " << type_name << " found on position " << Pos <<
  std::endl;
      }
      else
      {
        DataHelper<Size, Pos + 1>::get_pos(type_name, registered_names);
      }
    }
  };

  template <Uint Size>
  struct DataHelper<Size, Size - 1>
  {
    static void get_pos(const std::string &type_name,
                        const std::vector<std::string> &registered_names)
    {
      if (registered_names[Size - 1] == type_name)
      {
        std::cout << "Name " << type_name << " found on position " << Size - 1
  << std::endl;
      }
    }
  };
  */

  /// DATA
  using option_storage = std::vector<detail::OptionArrayBase *>;
  option_storage m_options;
};

// ----------------------------------------------------------------------------

template <typename T>
void OptionMap::create(const std::string &name, const T &value)
{
  // First go through existing option arrays and try to find out if there's
  // already an array with existing data type
  for (detail::OptionArrayBase *base_ptr : m_options)
  {
    detail::OptionArray<T> *arr_ptr = dynamic_cast<detail::OptionArray<T> *>(base_ptr);
    if (arr_ptr)
    {
      detail::OptionPair<T> opt_pair(name, value);
      arr_ptr->insert_value(opt_pair);
      return;
    }
  }

  // If option array with suitable type does not exist yet,
  // create it and insert option

  detail::OptionArrayBase *new_opt_vector = new detail::OptionArray<T>;
  detail::OptionPair<T> opt_pair(name, value);
  new_opt_vector->insert_value(opt_pair);
  m_options.push_back(new_opt_vector);
}

// ----------------------------------------------------------------------------

template <typename T>
void OptionMap::create(const std::string &name)
{
  this->create<T>(name, T());
}

// ----------------------------------------------------------------------------

template <typename T>
T OptionMap::get(const std::string &name) const
{
  for (detail::OptionArrayBase *base_ptr : m_options)
  {
    const detail::OptionArray<T> *arr_ptr = dynamic_cast<detail::OptionArray<T> *>(base_ptr);
    if (arr_ptr)
    {
      for (size_t i = 0; i < arr_ptr->size(); ++i)
      {
        if (arr_ptr->name(i) == name)
        {
          return arr_ptr->value(i);
        }
      }
    }
  }

  return T();
}

// ----------------------------------------------------------------------------

template <typename T>
void OptionMap::set(const std::string &name, const T &value)
{
  for (detail::OptionArrayBase *base_ptr : m_options)
  {
    detail::OptionArray<T> *arr_ptr = dynamic_cast<detail::OptionArray<T> *>(base_ptr);
    if (arr_ptr)
    {
      for (size_t i = 0; i < arr_ptr->size(); ++i)
      {
        if (arr_ptr->name(i) == name)
        {
          arr_ptr->change_value(i, value);
          return;
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit

#endif // PDEKIT_Common_Options_hpp
