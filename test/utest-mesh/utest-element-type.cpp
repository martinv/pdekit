#include <iostream>

#include "common/PDEKit.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::mesh::sf;

int main()
{

  SFTag etype(ElemShape::Quad, SFunc::Lagrange, P6, ModalBasis::Modal);

  std::string etype_desc;

  etype_desc = etype.as_string();
  std::cout << etype_desc << std::endl;

  etype_desc = "Tetra-BoundaryAdapted-Nodal-P2-Modal";

  etype = SFTag::string_to_tag(etype_desc);

  std::cout << "Element type:" << std::endl;
  std::cout << "[" << etype.as_string() << "]" << std::endl;

  // ==========================================================================

  /*
  std::cout << "The id of Lagrange Triag p1 in 2d = ";
  std::cout << SFTypeToID<Triag,Lagrange,Equidist,P1,Modal>::value <<
  std::endl;

  std::cout << "Decomposing the element id 255786:" << std::endl;
  std::cout << "Elem shape           = " <<
  ElemShape::Names[SFIDToType<255786>::shape] << std::endl;
  std::cout << "Shape function type  = " <<
  SFunc::Names[SFIDToType<255786>::sf_type] << std::endl;
  std::cout << "Reference topology   = " <<
  RefTopology::Names[SFIDToType<255786>::ref_topo] << std::endl;
  std::cout << "Shape function order = " <<
  PolyOrder::Names[SFIDToType<255786>::sf_order] << std::endl;
  std::cout << "Prime basis          = " <<
  PrimeBasis::Names[SFIDToType<255786>::prime_basis] << std::endl;
  */

  SFTag type2(ElemShape::Triag, SFunc::Lagrange, P1, ModalBasis::Modal);
  std::cout << "Shape function type = " << type2.as_string() << std::endl;

  return 0;
}
