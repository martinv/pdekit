/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE reference_topology_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers

#include <iostream>

#include "common/PDEKit.hpp"
#include "mesh/EntityDofRealign.hpp"
#include "mesh/std_region/StdRegion.hpp"
#include "mesh/std_region/StdRegionFactory.hpp"
#include "mesh/std_region/StdRegionWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

BOOST_AUTO_TEST_CASE(reference_topology_utest)
{
  std::cout << "============================================" << std::endl;
  std::cout << "Reference topology test - creation and" << std::endl
            << "decomposition of reference topology" << std::endl
            << "tags:" << std::endl;
  std::cout << "============================================" << std::endl;

  // Compose an reference topology tag
  PointSetTag reft(ElemShape::Quad, P6, PointSetID::Equidist);

  std::string reft_desc;
  reft_desc = reft.as_string();
  // std::cout << reft_desc << std::endl;
  assert(reft_desc == "Quad-P6-Equidist");

  // Create reference topology tag from string
  reft_desc = "Tetra-P7-Equidist";

  reft = PointSetTag::string_to_tag(reft_desc);

  std::cout << "Reference topology (should be Tetra-P7-Equidist):" << std::endl;
  std::cout << "[" << reft.as_string() << "]" << std::endl << std::endl;

  // ----------------------------------------------------------------------------

  PointSetTag p1_tri(ElemShape::Triag, P1, PointSetID::Equidist);

  std::cout << "The id of Equidist Triag p1 = ";
  std::cout << PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist) << std::endl;
  std::cout << "Decomposing the reference topology id 885:" << std::endl;
  std::cout << "Polynomial order = " << p1_tri.poly_order() << std::endl;

  PointSetTag type2(ElemShape::Triag, P1, PointSetID::Equidist);
  const PointSetTag type3(type2);
  std::cout << "Reference topology type = " << type3.as_string() << std::endl;

  StdRegionFactory::instance_type &rt_factory = StdRegionFactory::instance();

  const StdRegionFactory::instance_type::const_product_base_ptr rt =
      rt_factory.create(PointSetTag::string_to_tag("Quad-P1-Equidist"));

  /*
  std::cout << "Coordinates:\n" << rt->coords() << std::endl;
  std::cout << "Entities:" << std::endl;
  for(Uint i = 0; i < rt->nb_entities(_1D); ++i)
  {
    std::cout << rt->local_connectivity(_1D,i) << std::endl;
  }

  std::cout << "Faces: " << std::endl;
  std::cout << rt->local_connectivity(_2D,0) << std::endl;

  std::cout << "Volumes: " << std::endl;
  std::cout << rt->local_connectivity(_3D,0) << std::endl;
  */

  std::cout << std::endl;

  std::cout << "============================================" << std::endl;
  std::cout << "Reference element test:" << std::endl;
  std::cout << "============================================" << std::endl;

  StdRegion re;
  re.change_type(PointSetTag::string_to_tag("Quad-P3-Equidist"));
  // re.change_type(RefTopologyTag::string_to_tag("Tetra-P4-Equidist"));

  std::cout << "Reference element type = " << re.get().pt_set_id().as_string() << std::endl;
  BOOST_CHECK_EQUAL(re.get().pt_set_id().as_string(), "Quad-P3-Equidist");

  std::cout << "Edges:" << std::endl;
  for (Uint i = 0; i < re.get().nb_entities(_1D); ++i)
  {
    std::cout << *re.get().elem_entity(_1D, i) << std::endl;
  }

  std::cout << "Faces:" << std::endl;
  for (Uint i = 0; i < re.get().nb_entities(_2D); ++i)
  {
    std::cout << *re.get().elem_entity(_2D, i) << std::endl;
  }

  if (re.get().topo_dim() > 2)
  {
    std::cout << "Volumes:" << std::endl;
    for (Uint i = 0; i < re.get().nb_entities(_3D); ++i)
    {
      std::cout << *re.get().elem_entity(_3D, i) << std::endl;
    }
  }

  std::cout << "COORDINATES OF TETRA:" << std::endl;
  re.change_type(PointSetTag::string_to_tag("Tetra-P3-Equidist"));
  std::cout << re.get().coordinates() << std::endl;
  std::cout << std::endl;

  std::cout << "============================================" << std::endl;
  std::cout << "Permutation test:" << std::endl;
  std::cout << "============================================" << std::endl;

  EntityDofRealign permutation;

  std::vector<Uint> p3triag(10);

  p3triag[0] = 100;
  p3triag[1] = 101;
  p3triag[2] = 102;
  p3triag[3] = 20;
  p3triag[4] = 21;
  p3triag[5] = 30;
  p3triag[6] = 31;
  p3triag[7] = 40;
  p3triag[8] = 41;
  p3triag[9] = 1;

  EntityRealignCode perm_code = EntityRealignCode::identity(ElemShape::Line);
  perm_code.add_rotation();
  permutation.change_type(PointSetTag::string_to_tag("Line-P3-Equidist"), perm_code);

  std::cout << "Rotation of Line-P3-Equidist = " << std::endl;
  permutation.get().print();
  std::cout << std::endl;

  perm_code.reset();
  perm_code.set_nb_flips(1);
  perm_code.set_nb_rotations(1);
  permutation.change_type(PointSetTag::string_to_tag("Line-P3-Equidist"), perm_code);

  std::cout << "Rotation and flip of Line-P3-Equidist = " << std::endl;
  permutation.get().print();
  std::cout << std::endl;

  perm_code = EntityRealignCode::identity(ElemShape::Triag);
  perm_code.set_nb_flips(0);
  perm_code.set_nb_rotations(3);
  permutation.change_type(PointSetTag::string_to_tag("Triag-P3-Equidist"), perm_code);

  for (Uint i = 0; i < permutation.get().size(); ++i)
  {
    assert(permutation[i] == i);
  }

  std::cout << "Sample p3 triangle before any permutation:" << std::endl;
  for (Uint i = 0; i < p3triag.size(); ++i)
  {
    std::cout << p3triag[i] << " ";
  }
  std::cout << std::endl << std::endl;

  std::cout << "3 x rotation of Triag-P3-Equidist = " << std::endl;
  permutation.get().print();
  std::cout << std::endl;

  perm_code.set_nb_flips(1);
  perm_code.set_nb_rotations(1);
  permutation.change_type(PointSetTag::string_to_tag("Triag-P3-Equidist"), perm_code);

  std::cout << "Flip and rotation of Triag-P3-Equidist = " << std::endl;
  permutation.get().print();
  std::cout << std::endl;

  BOOST_CHECK_EQUAL(p3triag[permutation.get().vertex(0)], 102u);
  BOOST_CHECK_EQUAL(p3triag[permutation.get().vertex(1)], 101u);
  BOOST_CHECK_EQUAL(p3triag[permutation.get().vertex(2)], 100u);
  BOOST_CHECK_EQUAL(p3triag[permutation.get().vertex(3)], 31u);
  BOOST_CHECK_EQUAL(p3triag[permutation.get().vertex(4)], 30u);
  BOOST_CHECK_EQUAL(p3triag[permutation.get().vertex(5)], 21u);
  BOOST_CHECK_EQUAL(p3triag[permutation.get().vertex(6)], 20u);
  BOOST_CHECK_EQUAL(p3triag[permutation.get().vertex(7)], 41u);
  BOOST_CHECK_EQUAL(p3triag[permutation.get().vertex(8)], 40u);
  BOOST_CHECK_EQUAL(p3triag[permutation.get().vertex(9)], 1u);

  std::cout << "Sample p3 triangle after flip and rotation:" << std::endl;
  for (Uint i = 0; i < p3triag.size(); ++i)
  {
    std::cout << p3triag[permutation.get().vertex(i)] << " ";
  }
  std::cout << std::endl;

  std::cout << "============================================" << std::endl;

  std::cout << "Warpblend test" << std::endl;
  StdRegionBuilder::ptr rtw = rt_factory.create(PointSetTag::string_to_tag("Triag-P9-Warpblend"));

  math::DenseDMat<Real> warp_coords;
  rtw->coordinates(warp_coords);

  std::cout << "Coordinates = " << std::endl << warp_coords << std::endl;

  mesh::StdRegionWriter::write_to_vtk(
      mesh::PointSetTag(ElemShape::Tetra, P9, PointSetID::Warpblend), "reference_topology.vtk");
  mesh::StdRegionWriter::write_to_point3d(
      mesh::PointSetTag(ElemShape::Tetra, P9, PointSetID::Warpblend), "reference_topology.3D");
  mesh::StdRegionWriter::write_to_gmsh(
      mesh::PointSetTag(ElemShape::Tetra, P9, PointSetID::Warpblend), "reference_topology.msh");
}
