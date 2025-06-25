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

  const double dx = (x_max - x_min) / (N - 1);

  // Ids of important vertices
  const unsigned int BottomIn  = 1;
  const unsigned int BottomOut = N;

  const unsigned int TopOut = N + 1;
  const unsigned int TopIn  = N + 2;

  std::ofstream outfile;

  outfile.setf(std::ios::fixed);
  outfile.precision(15);

  outfile.open("exp_bump_2d.geo");

  outfile << "lc = 0.2;\n" << std::endl;

  for (unsigned int i = 0; i < N; ++i)
  {
    const double x = x_min + i * dx;
    const double y = 0.0625 * std::exp(-25. * x * x);

    insert_point_to_file(i + 1, x, y, 0.0, "lc", outfile);
  }

  insert_point_to_file(TopOut, x_max, y_max, 0.0, "lc", outfile);
  insert_point_to_file(TopIn, x_min, y_max, 0.0, "lc", outfile);

  outfile << std::endl;
  outfile << "Spline(1) = {" << BottomIn << ":" << BottomOut << "};" << std::endl;
  outfile << "Line(2)   = {" << BottomOut << "," << TopOut << "};" << std::endl;
  outfile << "Line(3)   = {" << TopOut << "," << TopIn << "};" << std::endl;
  outfile << "Line(4)   = {" << TopIn << "," << BottomIn << "};" << std::endl;

  outfile << std::endl;

  outfile << "Line Loop(1) = {1,2,3,4};" << std::endl;

  outfile << std::endl;
  outfile << "Plane Surface(1) = {1};" << std::endl;

  outfile << std::endl;
  outfile << "Physical Line(\"inlet\")  = {4};" << std::endl;
  outfile << "Physical Line(\"bottom_wall\") = {1};" << std::endl;
  outfile << "Physical Line(\"outlet\") = {2};" << std::endl;
  outfile << "Physical Line(\"top_wall\")    = {3};" << std::endl;

  outfile << std::endl;
  outfile << "Physical Surface(\"domain\") = {1};" << std::endl;

  outfile.close();

  return 0;
}
