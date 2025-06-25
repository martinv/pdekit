/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE element_adapt_interpolator_utest
#include <boost/test/unit_test.hpp>

#include <iomanip>
#include <iostream>
#include <vector>

#include "common/PDEKit.hpp"
#include "interpolation/ElementAdaptInterpolator.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::interpolation;

// ----------------------------------------------------------------------------

void check_interpolated_values_on_children(const PointSetTag parent_tag,
                                           const CellTransform adapt_case, const bool print)
{
  StdRegion parent_std_region(parent_tag);

  math::DenseDMat<Real> const &parent_coords = parent_std_region.get().coordinates();

  // 2 Fields: one for x^2 + y^3, the other for x^3 + y^2
  math::DenseDMat<Real> field_values(parent_coords.rows(), 2);
  for (Uint i = 0; i < parent_coords.rows(); ++i)
  {
    const math::DenseConstVecView<Real> node_coords = parent_coords.const_row_transp(i);
    field_values(i, 0) = std::pow(node_coords[X0], 2) + std::pow(node_coords[X1], 3);
    field_values(i, 1) = std::pow(node_coords[X0], 3) + std::pow(node_coords[X1], 2);
  }

  ElementAdaptInterpolator interp;
  const adapt::CellAdaptOpTag cell_adapt_op_tag(parent_tag.elem_shape(), adapt_case);
  interp.compute_child_values(parent_tag, cell_adapt_op_tag, field_values);

  if (print)
  {
    std::cout.precision(8);
    std::cout.setf(std::ios::fixed);
    std::cout << "=========================== Parent data "
                 "==========================="
              << std::endl;
    std::cout << "         x0               x1              f[0]           "
                 "  f[1]    "
              << std::endl;

    for (Uint i = 0; i < parent_coords.rows(); ++i)
    {
      const math::DenseConstVecView<Real> node_coords = parent_coords.const_row_transp(i);
      std::cout << std::setw(15) << node_coords[X0] << "  " << std::setw(15) << node_coords[X1]
                << "  " << std::setw(15) << field_values(i, 0) << "  " << std::setw(15)
                << field_values(i, 1) << std::endl;
    }

    std::cout << std::endl;
  }

  // Ask the local adaptation strategy to compute child coordinates
  mesh::adapt::CellAdaptOp adapt_operation(cell_adapt_op_tag);
  std::vector<math::DenseDMat<Real>> child_coords;
  adapt_operation.get().compute_child_coords(parent_tag, child_coords);

  for (Uint c = 0; c < interp.nb_child_value_blocks(); ++c)
  {
    if (print)
    {
      std::cout << "========================== Child data [" << c
                << "] =========================" << std::endl;
      std::cout << "         x0               x1              f[0]             "
                   "f[1]    "
                << std::endl;
    }

    const math::DenseConstMatView<Real> child_values = interp.interpolated_child_values(c);

    for (Uint i = 0; i < parent_coords.rows(); ++i)
    {
      if (print)
      {
        std::cout << std::setw(15) << child_coords[c](i, X0) << "  " << std::setw(15)
                  << child_coords[c](i, X1) << "  " << std::setw(15) << child_values(i, 0) << "  "
                  << std::setw(15) << child_values(i, 1) << std::endl;
      }

      const Real f0_check =
          std::pow(child_coords[c](i, X0), 2) + std::pow(child_coords[c](i, X1), 3);
      const Real f1_check =
          std::pow(child_coords[c](i, X0), 3) + std::pow(child_coords[c](i, X1), 2);

      BOOST_CHECK_LE(std::abs(f0_check - child_values(i, 0)), 1.e-7);
      BOOST_CHECK_LE(std::abs(f1_check - child_values(i, 1)), 1.e-7);
    }
    if (print)
    {
      std::cout << std::endl;
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(element_adapt_interpolator_basic_test)
{
  for (Uint p = 3; p <= 5; ++p)
  {
    PointSetTag parent_tag(ElemShape::Triag, p, PointSetID::Equidist);
    check_interpolated_values_on_children(parent_tag, CellTransform::UNIFORM_REFINE, false);
  }
}

// ----------------------------------------------------------------------------
