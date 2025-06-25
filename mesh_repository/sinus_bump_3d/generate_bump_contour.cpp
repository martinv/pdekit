#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

int main()
{

  const double PI = 3.14159265358979323846;
  std::ofstream outfile;
  const int nb_pts = 200;

  const double dx = 2.0 / (nb_pts - 1);

  outfile.precision(10);
  outfile.setf(std::ios_base::scientific);
  outfile.open("bump_contour.geo");

  outfile << "//Characteristic length at points" << std::endl;
  outfile << "lc = 0.2;\n\n";
  outfile << "// Half-width of the domain" << std::endl;
  outfile << "W = 0.5;\n\n";

  outfile << "Point(1) = {0.0,  0.0, -W, lc};\n";

  for (int i = 0; i < nb_pts; ++i)
  {
    const double x = 1.0 + i * dx;
    const double y = 0.1 * std::sin(PI * (x - 1.5)) + 0.1;

    outfile << "Point(" << 2 + i << ") = {" << x << ", " << y << ", -W, lc};\n";
  }

  outfile << "Point(" << nb_pts + 2 << ") = {4.0, 0.0, -W, lc};\n";
  outfile << "Point(" << nb_pts + 3 << ") = {4.0, 1.0, -W, lc};\n";
  outfile << "Point(" << nb_pts + 4 << ") = {0.0, 1.0, -W, lc};\n";

  outfile << "\n";

  outfile << "Line(1) = {1,2};\n";
  outfile << "Spline(2) = {2:" << nb_pts + 1 << "};\n";
  outfile << "Line(3) = {" << nb_pts + 1 << "," << nb_pts + 2 << "};\n";
  outfile << "Line(4) = {" << nb_pts + 2 << "," << nb_pts + 3 << "};\n";
  outfile << "Line(5) = {" << nb_pts + 3 << "," << nb_pts + 4 << "};\n";
  outfile << "Line(6) = {" << nb_pts + 4 << ",1};\n";

  outfile << "\n";

  /*
  outfile << "Physical Line( \"LowerWall\" ) = { 1 };\n";
  outfile << "Physical Line( \"SubOutlet\" ) = { 2 };\n";
  outfile << "Physical Line( \"UpperWall\" ) = { 3 };\n";
  outfile << "Physical Line( \"SubInlet\" ) = { 4 };\n";

  outfile << "\n";

  outfile << "Line Loop(1) = {1,2,3,4};\n";

  outfile << "\n";

  outfile << "Plane Surface(1) = {1};\n";
  outfile << "Physical Surface( \"InnerCells\" ) = {1};\n";

  outfile << "\n";

  outfile << "pts_hor = 15; //number of points on horizontal wall\n";
  outfile << "pts_vert = 4; //number of points on vertical wall\n\n";


  outfile << "Transfinite Line{1} = pts_hor;\n";
  outfile << "Transfinite Line{2} = pts_vert;\n";
  outfile << "Transfinite Line{3} = pts_hor Using Bump 0.8;\n\n";
  outfile << "Transfinite Line{4} = pts_vert;\n";

  outfile << "Transfinite Surface { 1 } = {  } Right;\n";
  */

  outfile.close();

  return 0;
}
