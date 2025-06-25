#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/// TYPES
struct Coordinate
{
  double x;
  double y;
};

/// FUNCTIONS
void ComputeYCoord(const double tparam, double &cy)
{
  cy = 0.12 / 0.2 *
       (0.2969 * std::sqrt(tparam) - 0.126 * tparam - 0.3516 * tparam * tparam +
        0.2843 * tparam * tparam * tparam - 0.1015 * tparam * tparam * tparam * tparam);
}

int main(int argc, char *argv[])
{

  std::ofstream outfile;
  std::stringstream ss;
  std::string outfile_name;
  //   std::string input_arg;

  if (argc != 2)
  {
    std::cerr << "Wrong number of arguments.\n";
    std::cerr << "Usage: " << argv[0] << " Filename where Filename = output *.geo file.\n";
    std::cerr << "Exiting.\n";
    exit(1);
  }
  else
  {
    outfile_name = std::string(argv[1]);
  }

  /// ============================================================================
  /// Constants: change the following to
  /// geometry of the output geo file

  // Characteristic lengths of points on the airfoil
  const double lc = 1.5;

  // Half of the width of the domain
  const double W = 2.0;

  // Number of intervals, in which the choord is divided
  // when computing the contour of the airfoil
  const unsigned int nb_x_intervals = 100;

  /// ============================================================================
  /// Extra variables needed for generating the aifoil

  std::vector<Coordinate> coordinates;
  double coord_x, coord_y;
  double dx, dx_ref;
  double weight;
  const double weight_threshold = 0.01;
  unsigned int nb_nodes, pt_idx;
  unsigned int le_idx, te_idx;
  unsigned int temp_idx, div_pt1, div_pt2;

  /// Start writing the output geo file

  outfile.precision(3);
  outfile.setf(std::ios_base::fixed);
  outfile.open(outfile_name.c_str());

  outfile << "// =============================================================" << std::endl;
  outfile << "// USER-DEFINED PARAMETERS" << std::endl;
  outfile << "// =============================================================" << std::endl;

  outfile << "// Characteristic length at airfoil points" << std::endl;
  outfile << "lc = " << lc << ";\n\n";

  outfile << "// Half-width of the domain" << std::endl;
  outfile << "W = " << W << ";\n\n";

  dx_ref = 1.0 / nb_x_intervals;

  coord_x  = 1.009;
  nb_nodes = 1;

  weight = 1.0 - 4 * (coord_x - 0.5) * (coord_x - 0.5);
  weight = std::max(weight, weight_threshold);
  dx     = dx_ref * weight;

  while (coord_x - dx > 0)
  {
    coord_x -= dx;
    nb_nodes++;

    weight = 1.0 - 4 * (coord_x - 0.5) * (coord_x - 0.5);
    weight = std::max(weight, weight_threshold);

    dx = dx_ref * weight;
    if (coord_x < 0.1)
      dx *= 0.5;
    if (coord_x < 0.05)
      dx *= 0.5;
  }

  le_idx = nb_nodes; // C++ indexing: first point has index 0
  nb_nodes *= 2;     // Include points on pressure side of the airfoil
  te_idx = 0;        // C++ indexing: first point has index 0

  coordinates.resize(nb_nodes);

  //====================================================
  //         STORE AIRFOIL POINT COORDINATES
  //====================================================

  // Points at trailing edge
  coord_x = 1.0089;
  coord_y = 0.0;
  pt_idx  = 0;

  coordinates[0].x = coord_x;
  coordinates[0].y = coord_y;

  weight = 1.0 - 4 * (coord_x - 0.5) * (coord_x - 0.5);
  weight = std::max(weight, weight_threshold);

  dx = dx_ref * weight;

  while (coord_x - dx > 0)
  {

    coord_x -= dx;
    ComputeYCoord(coord_x, coord_y);

    pt_idx++;
    coordinates[pt_idx].x = coord_x;
    coordinates[pt_idx].y = coord_y;

    coordinates[nb_nodes - pt_idx].x = coord_x;
    coordinates[nb_nodes - pt_idx].y = -coord_y;

    weight = 1.0 - 4 * (coord_x - 0.5) * (coord_x - 0.5);
    weight = std::max(weight, weight_threshold);

    dx = dx_ref * weight;
    if (coord_x < 0.1)
      dx *= 0.5;
    if (coord_x < 0.05)
      dx *= 0.5;

  } // while

  pt_idx++;
  coordinates[pt_idx].x = 0.0;
  coordinates[pt_idx].y = 0.0;

  // Write the coordinates to file
  outfile.setf(std::ios_base::scientific);
  outfile.precision(10);
  for (unsigned int i = 0; i < coordinates.size(); ++i)
  {
    outfile << "Point(" << i + 1 << ") = {" << coordinates[i].x << ", " << coordinates[i].y
            << ", -W, lc};\n";
  }
  outfile << std::endl;

  // The contour of the airfoil consists of 6 splines (3 on top, 3 on bottom).
  // Each side (bottom/top) is therefore divided by 2 points (div_pt1 and
  // div_pt2) into three segments. The INDEX of this points is around chord
  // length 0.049 (i.e. close to the leading edge) and 0.59 ( ~ in the
  // 'middle' of the airfoil )
  temp_idx = 0;
  while (coordinates[temp_idx].x > 0.049)
    ++temp_idx;
  div_pt1 = temp_idx - 1;

  temp_idx = 0;
  while (coordinates[temp_idx].x > 0.59)
    ++temp_idx;
  div_pt2 = temp_idx - 1;

  outfile << "trans_pt1 = " << le_idx - div_pt1 << ";\n";
  outfile << "trans_pt2 = " << div_pt1 - div_pt2 << ";\n";
  outfile << "le_idx = " << le_idx + 1 << ";\n";
  outfile << "nb_airfoil_nodes = " << nb_nodes << ";\n\n";
  // outfile << "Spline(" << airfoil_contour_idx[0] << ") =
  // {le_idx:le_idx-trans_pt1:-1};\n"; outfile << "Spline(" <<
  // airfoil_contour_idx[1] << ") =
  // {le_idx-trans_pt1:le_idx-trans_pt1-trans_pt2:-1};\n";
  // outfile << "Spline(" << airfoil_contour_idx[2] << ") =
  // {le_idx-trans_pt1-trans_pt2:1:-1};\n"; outfile << "Spline(" <<
  // airfoil_contour_idx[3] << ") = {1,
  // nb_airfoil_nodes:le_idx+trans_pt1+trans_pt2:-1};\n";
  // outfile << "Spline(" << airfoil_contour_idx[4] << ") =
  // {le_idx+trans_pt1+trans_pt2:le_idx+trans_pt1:-1};\n";
  // outfile << "Spline(" << airfoil_contour_idx[5] << ") =
  // {le_idx+trans_pt1:le_idx:-1};\n";

  outfile << "// Lines defining the airfoil contour:" << std::endl;
  outfile << "Spline(1) = {le_idx:le_idx-trans_pt1:-1};\n";
  outfile << "Spline(2) = {le_idx-trans_pt1:le_idx-trans_pt1-trans_pt2:-1};\n";
  outfile << "Spline(3) = {le_idx-trans_pt1-trans_pt2:1:-1};\n";
  outfile << "Spline(4) = {1, nb_airfoil_nodes:le_idx+trans_pt1+trans_pt2:-1};\n";
  outfile << "Spline(5) = {le_idx+trans_pt1+trans_pt2:le_idx+trans_pt1:-1};\n";
  outfile << "Spline(6) = {le_idx+trans_pt1:le_idx:-1};\n";

  outfile << std::endl;

  outfile.close();

  return 0;
}
