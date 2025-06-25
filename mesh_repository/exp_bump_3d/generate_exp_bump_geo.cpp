#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

void insert_point_to_file(const unsigned int pt_idx, const double x, const double y, const double z,
                          const std::string &lc_name, std::ofstream &fs)
{
  fs << "Point(" << pt_idx << ") = {" << x << ", " << y << ", " << z << ", " << lc_name << "};"
     << std::endl;
}

int main()
{
  const unsigned int N = 200;

  const double x_min = -1.5;
  const double x_max = 1.5;
  const double y_min = 0.0;
  const double y_max = 0.8;
  const double z_min = -0.4;
  const double z_max = 0.4;

  const double dx = (x_max - x_min) / (N - 1);

  // Ids of important vertices
  const unsigned int BottomLeftIn   = 1;
  const unsigned int BottomLeftOut  = N;
  const unsigned int BottomRightIn  = N + 1;
  const unsigned int BottomRightOut = 2 * N;

  const unsigned int TopLeftIn   = 2 * N + 1;
  const unsigned int TopLeftOut  = 2 * N + 2;
  const unsigned int TopRightIn  = 2 * N + 3;
  const unsigned int TopRightOut = 2 * N + 4;

  std::ofstream outfile;

  outfile.setf(std::ios::fixed);
  outfile.precision(15);

  outfile.open("exp_bump.geo");

  outfile << "lc = 0.2;\n" << std::endl;

  for (unsigned int i = 0; i < N; ++i)
  {
    const double x = x_min + i * dx;
    const double y = 0.0625 * std::exp(-25. * x * x);

    insert_point_to_file(i + 1, x, y, z_min, "lc", outfile);
  }

  for (unsigned int i = 0; i < N; ++i)
  {
    const double x = x_min + i * dx;
    const double y = 0.0625 * std::exp(-25. * x * x);

    insert_point_to_file(N + i + 1, x, y, z_max, "lc", outfile);
  }

  insert_point_to_file(TopLeftIn, x_min, y_max, z_min, "lc", outfile);
  insert_point_to_file(TopRightIn, x_min, y_max, z_max, "lc", outfile);
  insert_point_to_file(TopLeftOut, x_max, y_max, z_min, "lc", outfile);
  insert_point_to_file(TopRightOut, x_max, y_max, z_max, "lc", outfile);

  outfile << std::endl;
  outfile << "Spline(1) = {" << BottomLeftIn << ":" << BottomLeftOut << "};" << std::endl;
  outfile << "Spline(2) = {" << BottomRightIn << ":" << BottomRightOut << "};" << std::endl;
  outfile << "Line(3)   = {" << BottomLeftIn << "," << BottomRightIn << "};" << std::endl;
  outfile << "Line(4)   = {" << BottomRightOut << "," << BottomLeftOut << "};" << std::endl;

  outfile << "Line(5)   = {" << TopLeftIn << "," << TopLeftOut << "};" << std::endl;
  outfile << "Line(6)   = {" << TopRightIn << "," << TopRightOut << "};" << std::endl;
  outfile << "Line(7)   = {" << TopLeftIn << "," << TopRightIn << "};" << std::endl;
  outfile << "Line(8)   = {" << TopRightOut << "," << TopLeftOut << "};" << std::endl;

  outfile << "Line(9)   = {" << BottomRightIn << "," << TopRightIn << "};" << std::endl;
  outfile << "Line(10)  = {" << BottomRightOut << "," << TopRightOut << "};" << std::endl;
  outfile << "Line(11)  = {" << BottomLeftOut << "," << TopLeftOut << "};" << std::endl;
  outfile << "Line(12)  = {" << BottomLeftIn << "," << TopLeftIn << "};" << std::endl;

  outfile << std::endl;

  outfile << "Line Loop(1) = {3,2,4,-1};" << std::endl;
  outfile << "Line Loop(2) = {6,8,-5,7};" << std::endl;
  outfile << "Line Loop(3) = {3,9,-7,-12};" << std::endl;
  outfile << "Line Loop(4) = {2,10,-6,-9};" << std::endl;
  outfile << "Line Loop(5) = {4,11,-8,-10};" << std::endl;
  outfile << "Line Loop(6) = {12,5,-11,-1};" << std::endl;

  outfile << std::endl;
  outfile << "Ruled Surface(1) = {1};" << std::endl;

  for (unsigned int i = 1; i < 6; ++i)
  {
    outfile << "Plane Surface(" << i + 1 << ") = {" << i + 1 << "};" << std::endl;
  }

  outfile << "Surface Loop(1) = {1:6};" << std::endl;
  outfile << std::endl;
  outfile << "Volume(1) = {1};";
  outfile << std::endl;

  outfile << std::endl;
  outfile << "Physical Surface(\"Inlet\")         = {3};" << std::endl;
  outfile << "Physical Surface(\"Outlet\")        = {5};" << std::endl;
  outfile << "Physical Surface(\"SymmetryLeft\")  = {6};" << std::endl;
  outfile << "Physical Surface(\"SymmetryRight\") = {4};" << std::endl;
  outfile << "Physical Surface(\"Top\")           = {2};" << std::endl;
  outfile << "Physical Surface(\"Bottom\")        = {1};" << std::endl;

  outfile << std::endl;
  outfile << "Physical Volume(\"Fluid\") = {1};" << std::endl;

  outfile.close();

  return 0;
}
