/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE local_interpolation_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "math/DenseVecView.hpp"
#include "mesh/CellGeometry.hpp"
#include "mesh/adaptation/LocalInterpolator.hpp"
#include "mesh/std_region/StdRegion.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

Real f_interp0_2D(const Real x0, const Real x1)
{
  return 7.0 * x0 * x0 - 3.0 * x1 * x1 * x1;
}
Real f_interp1_2D(const Real x0, const Real x1)
{
  return 3.14 * x0 * x1 - 2.0 * x0 * x1 * x1;
}

Real f_interp0_3D(const Real x0, const Real x1, const Real x2)
{
  return 3.14 * x0 * x0 * x1 - 2.0 * x0 * x1 * x1 + 7.14524 * x0 * x1 * x2;
}

Real f_interp1_3D(const Real x0, const Real x1, const Real x2)
{
  return 3.14 * x0 * x0 * x1 + 2.0 * x0 * x1 * x1 + 7.14524 * x0 * x2 * x2;
}

Real f_interp2_3D(const Real x0, const Real x1, const Real x2)
{
  return 0.00014 * x0 * x0 * x0 - 8.0 * x1 * x1 * x1 + 7.14333 * x2 * x2 * x2;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(local_interpolator_triag_utest)
{
  adapt::LocalInterpolator interpolator;

  mesh::StdRegion std_region_in;
  mesh::StdRegion std_region_out;

  for (Uint p = 3; p <= 8; ++p)
  {
    const PointSetTag std_reg_tag_in(ElemShape::Triag, p, PointSetID::Warpblend);
    const PointSetTag std_reg_tag_out(ElemShape::Triag, p + 2, PointSetID::Equidist);

    /*
    std::cout << "Interpolating from " << std_reg_tag_in.as_string() << " to
    "
              << std_reg_tag_out.as_string() << std::endl;
    */

    std_region_in.change_type(std_reg_tag_in);
    std_region_out.change_type(std_reg_tag_out);

    const math::DenseDMat<Real> &std_reg_in_coords  = std_region_in.get().coordinates();
    const math::DenseDMat<Real> &std_reg_out_coords = std_region_out.get().coordinates();

    math::DenseDMat<Real> data_in(std_region_in.get().nb_nodes(), 2);

    for (Uint n = 0; n < std_region_in.get().nb_nodes(); ++n)
    {
      data_in(n, 0) = f_interp0_2D(std_reg_in_coords(n, X0), std_reg_in_coords(n, X1));
      data_in(n, 1) = f_interp1_2D(std_reg_in_coords(n, X0), std_reg_in_coords(n, X1));
    }

    math::DenseDMat<Real> data_out_ref(std_region_out.get().nb_nodes(), 2);

    for (Uint n = 0; n < std_region_out.get().nb_nodes(); ++n)
    {
      data_out_ref(n, 0) = f_interp0_2D(std_reg_out_coords(n, X0), std_reg_out_coords(n, X1));
      data_out_ref(n, 1) = f_interp1_2D(std_reg_out_coords(n, X0), std_reg_out_coords(n, X1));
    }

    const math::DenseConstMatView<Real> data_out =
        interpolator.transfer_data(std_reg_tag_in, std_reg_tag_out, data_in);

    /*
    std::cout << "Data out - reference:" << std::endl;
    std::cout << data_out_ref << std::endl;

    std::cout << "Data out - interpolated:" << std::endl;
    std::cout << data_out << std::endl;
    */

    BOOST_CHECK_EQUAL(data_out_ref.rows(), data_out.rows());
    BOOST_CHECK_EQUAL(data_out_ref.cols(), data_out.cols());

    for (Uint r = 0; r < data_out_ref.rows(); ++r)
    {
      for (Uint c = 0; c < data_out_ref.cols(); ++c)
      {
        const Real diff = std::abs(data_out_ref(r, c) - data_out(r, c));
        BOOST_CHECK_LE(diff, 1.e-13);
      }
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(local_interpolator_quad_utest)
{
  adapt::LocalInterpolator interpolator;

  mesh::StdRegion std_region_in;
  mesh::StdRegion std_region_out;

  for (Uint p = 3; p <= 8; ++p)
  {
    const PointSetTag std_reg_tag_in(ElemShape::Quad, p, PointSetID::Warpblend);
    const PointSetTag std_reg_tag_out(ElemShape::Quad, p + 2, PointSetID::Equidist);

    /*
    std::cout << "Interpolating from " << std_reg_tag_in.as_string() << " to
    "
              << std_reg_tag_out.as_string() << std::endl;
    */

    std_region_in.change_type(std_reg_tag_in);
    std_region_out.change_type(std_reg_tag_out);

    const math::DenseDMat<Real> &std_reg_in_coords  = std_region_in.get().coordinates();
    const math::DenseDMat<Real> &std_reg_out_coords = std_region_out.get().coordinates();

    math::DenseDMat<Real> data_in(std_region_in.get().nb_nodes(), 2);

    for (Uint n = 0; n < std_region_in.get().nb_nodes(); ++n)
    {
      data_in(n, 0) = f_interp0_2D(std_reg_in_coords(n, X0), std_reg_in_coords(n, X1));
      data_in(n, 1) = f_interp1_2D(std_reg_in_coords(n, X0), std_reg_in_coords(n, X1));
    }

    math::DenseDMat<Real> data_out_ref(std_region_out.get().nb_nodes(), 2);

    for (Uint n = 0; n < std_region_out.get().nb_nodes(); ++n)
    {
      data_out_ref(n, 0) = f_interp0_2D(std_reg_out_coords(n, X0), std_reg_out_coords(n, X1));
      data_out_ref(n, 1) = f_interp1_2D(std_reg_out_coords(n, X0), std_reg_out_coords(n, X1));
    }

    const math::DenseConstMatView<Real> data_out =
        interpolator.transfer_data(std_reg_tag_in, std_reg_tag_out, data_in);

    /*
    std::cout << "Data out - reference:" << std::endl;
    std::cout << data_out_ref << std::endl;

    std::cout << "Data out - interpolated:" << std::endl;
    std::cout << data_out << std::endl;
    */

    BOOST_CHECK_EQUAL(data_out_ref.rows(), data_out.rows());
    BOOST_CHECK_EQUAL(data_out_ref.cols(), data_out.cols());

    for (Uint r = 0; r < data_out_ref.rows(); ++r)
    {
      for (Uint c = 0; c < data_out_ref.cols(); ++c)
      {
        const Real diff = std::abs(data_out_ref(r, c) - data_out(r, c));
        BOOST_CHECK_LE(diff, 1.e-14);
      }
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(local_interpolator_tetra_utest)
{
  adapt::LocalInterpolator interpolator;

  mesh::StdRegion std_region_in;
  mesh::StdRegion std_region_out;

  for (Uint p = 3; p <= 7; ++p)
  {
    const PointSetTag std_reg_tag_in(ElemShape::Tetra, p, PointSetID::Warpblend);
    const PointSetTag std_reg_tag_out(ElemShape::Tetra, p + 2, PointSetID::Equidist);

    /*
    std::cout << "Interpolating from " << std_reg_tag_in.as_string() << " to
    "
              << std_reg_tag_out.as_string() << std::endl;
    */

    std_region_in.change_type(std_reg_tag_in);
    std_region_out.change_type(std_reg_tag_out);

    const math::DenseDMat<Real> &std_reg_in_coords  = std_region_in.get().coordinates();
    const math::DenseDMat<Real> &std_reg_out_coords = std_region_out.get().coordinates();

    math::DenseDMat<Real> data_in(std_region_in.get().nb_nodes(), 3);

    for (Uint n = 0; n < std_region_in.get().nb_nodes(); ++n)
    {
      data_in(n, 0) = f_interp0_3D(std_reg_in_coords(n, X0), std_reg_in_coords(n, X1),
                                   std_reg_in_coords(n, X2));
      data_in(n, 1) = f_interp1_3D(std_reg_in_coords(n, X0), std_reg_in_coords(n, X1),
                                   std_reg_in_coords(n, X2));
      data_in(n, 2) = f_interp2_3D(std_reg_in_coords(n, X0), std_reg_in_coords(n, X1),
                                   std_reg_in_coords(n, X2));
    }

    math::DenseDMat<Real> data_out_ref(std_region_out.get().nb_nodes(), 3);

    for (Uint n = 0; n < std_region_out.get().nb_nodes(); ++n)
    {
      data_out_ref(n, 0) = f_interp0_3D(std_reg_out_coords(n, X0), std_reg_out_coords(n, X1),
                                        std_reg_out_coords(n, X2));
      data_out_ref(n, 1) = f_interp1_3D(std_reg_out_coords(n, X0), std_reg_out_coords(n, X1),
                                        std_reg_out_coords(n, X2));
      data_out_ref(n, 2) = f_interp2_3D(std_reg_out_coords(n, X0), std_reg_out_coords(n, X1),
                                        std_reg_out_coords(n, X2));
    }

    const math::DenseConstMatView<Real> data_out =
        interpolator.transfer_data(std_reg_tag_in, std_reg_tag_out, data_in);

    /*
    std::cout << "Data out - reference:" << std::endl;
    std::cout << data_out_ref << std::endl;

    std::cout << "Data out - interpolated:" << std::endl;
    std::cout << data_out << std::endl;
    */

    BOOST_CHECK_EQUAL(data_out_ref.rows(), data_out.rows());
    BOOST_CHECK_EQUAL(data_out_ref.cols(), data_out.cols());

    for (Uint r = 0; r < data_out_ref.rows(); ++r)
    {
      for (Uint c = 0; c < data_out_ref.cols(); ++c)
      {
        const Real diff = std::abs(data_out_ref(r, c) - data_out(r, c));
        BOOST_CHECK_LE(diff, 1.e-13);
      }
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(local_interpolator_hexa_utest)
{
  adapt::LocalInterpolator interpolator;

  mesh::StdRegion std_region_in;
  mesh::StdRegion std_region_out;

  for (Uint p = 3; p <= 4; ++p)
  {
    const PointSetTag std_reg_tag_in(ElemShape::Hexa, p, PointSetID::Equidist);
    const PointSetTag std_reg_tag_out(ElemShape::Hexa, p + 2, PointSetID::Equidist);

    /*
    std::cout << "Interpolating from " << std_reg_tag_in.as_string() << " to
    "
              << std_reg_tag_out.as_string() << std::endl;
    */

    std_region_in.change_type(std_reg_tag_in);
    std_region_out.change_type(std_reg_tag_out);

    const math::DenseDMat<Real> &std_reg_in_coords  = std_region_in.get().coordinates();
    const math::DenseDMat<Real> &std_reg_out_coords = std_region_out.get().coordinates();

    math::DenseDMat<Real> data_in(std_region_in.get().nb_nodes(), 3);

    for (Uint n = 0; n < std_region_in.get().nb_nodes(); ++n)
    {
      data_in(n, 0) = f_interp0_3D(std_reg_in_coords(n, X0), std_reg_in_coords(n, X1),
                                   std_reg_in_coords(n, X2));
      data_in(n, 1) = f_interp1_3D(std_reg_in_coords(n, X0), std_reg_in_coords(n, X1),
                                   std_reg_in_coords(n, X2));
      data_in(n, 2) = f_interp2_3D(std_reg_in_coords(n, X0), std_reg_in_coords(n, X1),
                                   std_reg_in_coords(n, X2));
    }

    math::DenseDMat<Real> data_out_ref(std_region_out.get().nb_nodes(), 3);

    for (Uint n = 0; n < std_region_out.get().nb_nodes(); ++n)
    {
      data_out_ref(n, 0) = f_interp0_3D(std_reg_out_coords(n, X0), std_reg_out_coords(n, X1),
                                        std_reg_out_coords(n, X2));
      data_out_ref(n, 1) = f_interp1_3D(std_reg_out_coords(n, X0), std_reg_out_coords(n, X1),
                                        std_reg_out_coords(n, X2));
      data_out_ref(n, 2) = f_interp2_3D(std_reg_out_coords(n, X0), std_reg_out_coords(n, X1),
                                        std_reg_out_coords(n, X2));
    }

    const math::DenseConstMatView<Real> data_out =
        interpolator.transfer_data(std_reg_tag_in, std_reg_tag_out, data_in);

    /*
    std::cout << "Data out - reference:" << std::endl;
    std::cout << data_out_ref << std::endl;

    std::cout << "Data out - interpolated:" << std::endl;
    std::cout << data_out << std::endl;
    */

    BOOST_CHECK_EQUAL(data_out_ref.rows(), data_out.rows());
    BOOST_CHECK_EQUAL(data_out_ref.cols(), data_out.cols());

    for (Uint r = 0; r < data_out_ref.rows(); ++r)
    {
      for (Uint c = 0; c < data_out_ref.cols(); ++c)
      {
        const Real diff = std::abs(data_out_ref(r, c) - data_out(r, c));
        BOOST_CHECK_LE(diff, 1.e-13);
      }
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(local_interpolator_triag_coordinates_utest)
{
  adapt::LocalInterpolator interpolator;

  mesh::StdRegion std_region_in;
  mesh::StdRegion std_region_out;

  for (Uint p = 3; p <= 8; ++p)
  {
    const PointSetTag std_reg_tag_in(ElemShape::Triag, p, PointSetID::Warpblend);
    const PointSetTag std_reg_tag_out(ElemShape::Triag, p + 2, PointSetID::Equidist);

    /*
    std::cout << "Interpolating from " << std_reg_tag_in.as_string() << " to
    "
              << std_reg_tag_out.as_string() << std::endl;
    */

    std_region_in.change_type(std_reg_tag_in);
    std_region_out.change_type(std_reg_tag_out);

    const math::DenseDMat<Real> &std_reg_in_coords  = std_region_in.get().coordinates();
    const math::DenseDMat<Real> &std_reg_out_coords = std_region_out.get().coordinates();

    std::vector<Uint> raw_nodes_in((p + 1) * (p + 2) / 2);
    std::iota(raw_nodes_in.begin(), raw_nodes_in.end(), 0);

    const common::ArrayView<const Uint, _1D, Uint> nodes_view_in(raw_nodes_in.data(),
                                                                 raw_nodes_in.size());

    const mesh::MeshEntity entity_in(nodes_view_in, 0, std_reg_tag_in);

    std::vector<Real> raw_coords(std_reg_in_coords.rows() * std_reg_in_coords.cols());

    for (Uint r = 0; r < std_reg_in_coords.rows(); ++r)
    {
      for (Uint c = 0; c < std_reg_in_coords.cols(); ++c)
      {
        raw_coords[r * std_reg_in_coords.cols() + c] = std_reg_in_coords(r, c);
      }
    }

    const math::DenseVecView<const Real> coord_view_in(raw_coords.data(), raw_coords.size());
    mesh::CellGeometry<_2D> triag_geo_in(coord_view_in, entity_in);

    const math::DenseConstMatView<Real> data_out =
        interpolator.transfer_coords(std_reg_tag_in, std_reg_tag_out, triag_geo_in);

    /*
    std::cout << "Data out - reference:" << std::endl;
    std::cout << std_reg_out_coords << std::endl;

    std::cout << "Data out - interpolated:" << std::endl;
    std::cout << data_out << std::endl;
    */

    BOOST_CHECK_EQUAL(std_reg_out_coords.rows(), data_out.rows());
    BOOST_CHECK_EQUAL(std_reg_out_coords.cols(), data_out.cols());

    for (Uint r = 0; r < std_reg_out_coords.rows(); ++r)
    {
      for (Uint c = 0; c < std_reg_out_coords.cols(); ++c)
      {
        const Real diff = std::abs(std_reg_out_coords(r, c) - data_out(r, c));
        BOOST_CHECK_LE(diff, 1.e-13);
      }
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(local_interpolator_facet_codim1_utest)
{
  adapt::LocalInterpolator interpolator;

  mesh::StdRegion std_region_in;
  mesh::StdRegion std_region_out;

  for (Uint p = 3; p <= 8; ++p)
  {
    const PointSetTag std_reg_tag_in(ElemShape::Line, p, PointSetID::Warpblend);
    const PointSetTag std_reg_tag_out(ElemShape::Line, p + 2, PointSetID::Equidist);

    /*
    std::cout << "Interpolating from " << std_reg_tag_in.as_string() << " to
    "
              << std_reg_tag_out.as_string() << std::endl;
    */

    std_region_in.change_type(std_reg_tag_in);
    std_region_out.change_type(std_reg_tag_out);

    const math::DenseDMat<Real> &std_reg_in_coords  = std_region_in.get().coordinates();
    const math::DenseDMat<Real> &std_reg_out_coords = std_region_out.get().coordinates();

    std::vector<Uint> raw_nodes_in(p + 1);
    std::iota(raw_nodes_in.begin(), raw_nodes_in.end(), 0);

    const common::ArrayView<const Uint, _1D, Uint> nodes_view_in(raw_nodes_in.data(),
                                                                 raw_nodes_in.size());

    const mesh::MeshEntity entity_in(nodes_view_in, 0, std_reg_tag_in);

    const Uint nb_nodes_in = std_reg_in_coords.rows();
    const Uint geo_dim_in  = std_reg_in_coords.cols() + 1;

    std::vector<Real> raw_coords(nb_nodes_in * geo_dim_in);

    for (Uint r = 0; r < nb_nodes_in; ++r)
    {
      for (Uint c = 0; c < (geo_dim_in - 1); ++c)
      {
        raw_coords[r * geo_dim_in + c] = std_reg_in_coords(r, c) + 0.1 * (c + 1);
      }
      raw_coords[r * geo_dim_in + geo_dim_in - 1] = 2.0;
    }

    const math::DenseVecView<const Real> coord_view_in(raw_coords.data(), raw_coords.size());
    mesh::CellGeometry<_2D> line_geo_in(coord_view_in, entity_in);

    const math::DenseConstMatView<Real> data_out =
        interpolator.transfer_coords(std_reg_tag_in, std_reg_tag_out, line_geo_in);

    /*
    std::cout << "Data out - reference:" << std::endl;
    std::cout << std_reg_out_coords << std::endl;

    std::cout << "Data out - interpolated:" << std::endl;
    std::cout << data_out << std::endl;
    */

    const Uint nb_nodes_out = std_reg_out_coords.rows();
    const Uint geo_dim_out  = std_reg_out_coords.cols() + 1;

    math::DenseDMat<Real> std_reg_out_coords_check;
    std_reg_out_coords_check.resize(nb_nodes_out, geo_dim_out);

    for (Uint r = 0; r < nb_nodes_out; ++r)
    {
      for (Uint c = 0; c < (geo_dim_out - 1); ++c)
      {
        std_reg_out_coords_check(r, c) = std_reg_out_coords(r, c) + 0.1 * (c + 1);
      }
      std_reg_out_coords_check(r, geo_dim_out - 1) = 2.0;
    }

    /*
    std::cout << "Geometry in:" << std::endl;
    std::cout << line_geo_in << std::endl;
    std::cout << "Computed:" << std::endl;
    std::cout << data_out << std::endl;
    std::cout << "Reference:" << std::endl;
    std::cout << std_reg_out_coords_check << std::endl;
    */

    BOOST_CHECK_EQUAL(std_reg_out_coords_check.rows(), data_out.rows());
    BOOST_CHECK_EQUAL(std_reg_out_coords_check.cols(), data_out.cols());

    for (Uint r = 0; r < std_reg_out_coords_check.rows(); ++r)
    {
      for (Uint c = 0; c < std_reg_out_coords_check.cols(); ++c)
      {
        const Real diff = std::abs(std_reg_out_coords_check(r, c) - data_out(r, c));
        BOOST_CHECK_LE(diff, 1.e-13);
      }
    }
  }
}

// ----------------------------------------------------------------------------
