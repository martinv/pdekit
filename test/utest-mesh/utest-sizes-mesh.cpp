#include <iomanip>
#include <iostream>

#include "common/PtrHandle.hpp"
#include "mesh/MeshEntity.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"
#include "mesh/std_region/StdRegion.hpp"
#include "mesh/std_region/StdRegionEntity.hpp"

template <typename T>
void print_type_size(const std::string &type_name)
{
  /// Output name of the type and its size
  /// Pad the name so that the first column is 35 characters wide
  std::cout << std::setw(35) << type_name << std::setw(10) << sizeof(T) << std::endl;
}

using namespace pdekit;

int main()
{
  std::cout << std::setw(30) << "TYPE" << std::setw(15) << "SIZE" << std::endl;
  std::cout << " ==================================================" << std::endl;

  print_type_size<mesh::StdRegion>("ReferenceElement");
  print_type_size<mesh::MeshEntity>("MeshEntity");
  print_type_size<mesh::StdRegionEntity>("StdRegionEntity");
  print_type_size<mesh::StdRegionEntity *>("StdRegionEntity*");
  print_type_size<common::PtrHandle<mesh::StdRegionEntity const>>("WrappedPtr<RefEntity const>");
  print_type_size<mesh::IncidenceEntry>("IncidencePair");
  print_type_size<mesh::TraceIncidences>("LocalIncidences");

  return 0;
}
