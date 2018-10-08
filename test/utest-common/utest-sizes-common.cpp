#include <bitset>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>

#include "common/PDEKit.hpp"

template <typename T>
void print_type_size(const std::string &type_name)
{
  /// Output name of the type and its size
  /// Pad the name so that the first column is 15 characters wide
  std::cout << std::setw(15) << type_name << std::setw(10) << sizeof(T) << std::endl;
}

template <typename T>
void print_type_size_and_limits(const std::string &type_name)
{
  /// Output name of the type, its size and min/max limits
  /// Pad the name so that the first column is 15 characters wide
  std::cout << std::setw(15) << type_name << std::setw(10) << sizeof(T) << "\t["
            << std::numeric_limits<T>::min() << "," << std::numeric_limits<T>::max() << "]"
            << std::endl;
}

using namespace pdekit;

int main()
{
  std::cout << "           TYPE        SIZE       LIMITS" << std::endl;
  std::cout << " ==========================================" << std::endl;

  print_type_size<bool>("bool");

  print_type_size<std::bitset<8>>("std::bitset<8>");

  print_type_size_and_limits<short int>("short int");

  print_type_size_and_limits<unsigned short int>("unsigned short int");

  print_type_size_and_limits<SUint>("SUint");

  print_type_size_and_limits<Uint>("Uint");

  print_type_size_and_limits<LUint>("LUint");

  print_type_size_and_limits<LLUint>("LLUint");

  print_type_size<Uint *>("Uint*");

  print_type_size_and_limits<Float>("Float");

  print_type_size_and_limits<Double>("Double");

  print_type_size_and_limits<LDouble>("Long Double");

  print_type_size_and_limits<Real>("Real");

  print_type_size<Real *>("Real*");

  print_type_size<std::ptrdiff_t>("ptrdiff_t");

  print_type_size<std::map<Uint, Uint *>::iterator>("std::map<Uint,Uint*>::iterator");

  return 0;
}
