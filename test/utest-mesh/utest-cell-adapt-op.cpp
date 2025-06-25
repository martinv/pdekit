/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cell_adapt_op_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <fstream>
#include <iostream>

/// PDEKIT headers
#include "common/StringUtils.hpp"
#include "mesh/MeshConstants.hpp"
#include "mesh/adaptation/CellAdaptOp.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

void write_coords_to_file(const math::DenseDMat<Real> &coords, const std::string &filename)
{
  std::ofstream outfile;
  outfile.open(filename.c_str());
  outfile.precision(15);

  for (Uint r = 0; r < coords.rows(); ++r)
  {
    outfile << coords.const_row_transp(r) << std::endl;
  }
  outfile.close();
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cell_adapt_op_utest)
{
  adapt::CellAdaptOp adapt_op;
  const adapt::CellAdaptOpTag adapt_tag(ElemShape::Triag, CellTransform::UNIFORM_REFINE);

  adapt_op.change_type(adapt_tag);

  math::DenseDMat<Real> coord_in, coord_out;

  mesh::StdRegion std_reg;
  std_reg.change_type(PointSetTag(ElemShape::Triag, P9, PointSetID::Equidist));

  const math::DenseDMat<Real> &coords = std_reg.get().coordinates();

  coord_in.resize(coords.rows(), coords.cols());
  coord_in = coords;

  write_coords_to_file(coord_in, "cell_adapt_op_utest_coord_in.dat");

  const std::vector<Uint> adapt_sub_elem_id = {1, 3, 2};

  for (Uint adapt = 0; adapt < adapt_sub_elem_id.size(); ++adapt)
  {
    adapt_op.get().transform_coords(coord_in, coord_in, adapt_sub_elem_id[adapt], coord_out);

    const std::string outfilename =
        "cell_adapt_op_utest_coord_out" + common::StringUtils::to_string(adapt) + ".dat";
    write_coords_to_file(coord_out, outfilename);

    coord_in = coord_out;
  }
}

// ----------------------------------------------------------------------------
