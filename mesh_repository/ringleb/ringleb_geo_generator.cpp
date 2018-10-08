#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

void place_wall_points(const double k, const size_t n_points, std::vector<double> &x,
                       std::vector<double> &y)
{
  const double gamma = 1.4;
  const double V0    = 0.43;
  const double dV    = (k - V0) / (n_points - 1);

  x.resize(2 * n_points - 1);
  y.resize(2 * n_points - 1);

  for (int i = 0; i < n_points; ++i)
  {
    const double V   = V0 + dV * i;
    const double a   = sqrt(1.0 - 0.5 * (gamma - 1) * V * V);
    const double rho = std::pow(a, (2. / (gamma - 1)));
    const double p   = 1. / gamma * std::pow(a, 2 * gamma / (gamma - 1));
    const double J   = 1. / a + 1. / (3. * a * a * a) + 1. / (5. * a * a * a * a * a) -
                     0.5 * std::log((1. + a) / (1. - a));

    x[i] = 1. / (2. * rho) * (1. / (V * V) - 2. / (k * k)) + J / 2.;
    y[i] = 1. / (k * rho * V) * std::sqrt(1. - (V / k) * (V / k));
  }

  for (int i = 0; i < (n_points - 1); ++i)
  {
    x[n_points + i] = x[n_points - 2 - i];
    y[n_points + i] = -y[n_points - 2 - i];
  }
}

void place_inlet_outlet_points(const double V, const size_t n_points, const double sign,
                               std::vector<double> &x, std::vector<double> &y)
{
  const double gamma = 1.4;
  const double k_min = 0.6;
  const double k_max = 0.98;
  const double dk    = (k_max - k_min) / (n_points - 1);

  x.resize(n_points);
  y.resize(n_points);

  for (int i = 0; i < n_points; ++i)
  {
    const double k   = k_min + dk * i;
    const double a   = sqrt(1.0 - 0.5 * (gamma - 1) * V * V);
    const double rho = std::pow(a, (2. / (gamma - 1)));
    const double p   = 1. / gamma * std::pow(a, 2 * gamma / (gamma - 1));
    const double J   = 1. / a + 1. / (3. * a * a * a) + 1. / (5 * a * a * a * a * a) -
                     0.5 * std::log((1. + a) / (1. - a));

    x[i] = 1. / (2. * rho) * (1. / (V * V) - 2. / (k * k)) + J / 2.;
    y[i] = sign * 1. / (k * rho * V) * std::sqrt(1. - (V / k) * (V / k));
  }
}

int main()
{
  const double k_min = 0.6;
  const double k_max = 0.98;

  const double V0 = 0.43;

  std::vector<double> x, y;

  // Generate gmsh geometry file

  // Number of generating points
  const size_t n_gen_points = 100;

  std::ofstream outfile;
  outfile.setf(std::ios::fixed);
  outfile.precision(15);

  outfile.open("ringleb.geo");

  outfile << "lc = 0.3;" << std::endl;

  size_t point_counter = 1;

  // Inner wall
  place_wall_points(k_min, n_gen_points, x, y);
  const size_t num_inner_wall_pts = x.size();

  for (int i = 0; i < x.size(); ++i)
  {
    outfile << "Point(" << point_counter++ << ") = {" << x[i] << "," << y[i] << ", 0.0, lc};"
            << std::endl;
  }

  // Outlet
  const double outlet_sign = -1.0;
  place_inlet_outlet_points(V0, n_gen_points, outlet_sign, x, y);

  // Remove first point - first shift points on position 1,2, ... to the left
  // This will overwrite x[0] and y[0]
  for (int i = 0; i < x.size(); ++i)
  {
    x[i] = x[i + 1];
    y[i] = y[i + 1];
  }

  // Now shrink the vectors by 1
  x.pop_back();
  y.pop_back();

  const size_t num_outlet_pts = x.size();

  for (int i = 0; i < num_outlet_pts; ++i)
  {
    outfile << "Point(" << point_counter++ << ") = {" << x[i] << "," << y[i] << ", 0.0, lc};"
            << std::endl;
  }

  // Outer wall
  place_wall_points(k_max, n_gen_points, x, y);
  // Skip the last point - already written as inner wall point
  x.pop_back();
  y.pop_back();

  // Reverse the coordinates
  std::reverse(x.begin(), x.end());
  std::reverse(y.begin(), y.end());

  const size_t num_outer_wall_pts = x.size();

  for (int i = 0; i < x.size(); ++i)
  {
    outfile << "Point(" << point_counter++ << ") = {" << x[i] << "," << y[i] << ", 0.0, lc};"
            << std::endl;
  }

  // Inlet
  const double inlet_sign = 1.0;
  place_inlet_outlet_points(V0, n_gen_points, inlet_sign, x, y);

  // Reverse the coordinates
  std::reverse(x.begin(), x.end());
  std::reverse(y.begin(), y.end());

  // Remove first point - first shift points on position 1,2, ... to the left
  // This will overwrite x[0] and y[0]
  for (int i = 0; i < x.size(); ++i)
  {
    x[i] = x[i + 1];
    y[i] = y[i + 1];
  }

  // Now shrink the vectors by 1
  x.pop_back();
  y.pop_back();

  // In addition, remove the last point, because
  // that's the first point of the inner wall - we're closing the loop here
  x.pop_back();
  y.pop_back();

  const size_t num_inlet_pts = x.size();

  for (int i = 0; i < x.size(); ++i)
  {
    outfile << "Point(" << point_counter++ << ") = {" << x[i] << "," << y[i] << ", 0.0, lc};"
            << std::endl;
  }

  // Write down lines
  outfile << "Spline(1) = {1:" << num_inner_wall_pts << "};" << std::endl;
  outfile << "Spline(2) = {" << num_inner_wall_pts << ":" << num_inner_wall_pts + num_outlet_pts
          << "};" << std::endl;
  outfile << "Spline(3) = {" << num_inner_wall_pts + num_outlet_pts << ":"
          << num_inner_wall_pts + num_outlet_pts + num_outer_wall_pts << "};" << std::endl;
  outfile << "Spline(4) = {" << num_inner_wall_pts + num_outlet_pts + num_outer_wall_pts << ":"
          << num_inner_wall_pts + num_outlet_pts + num_outer_wall_pts + num_inlet_pts << ",1};"
          << std::endl;

  outfile << "Line Loop(1) = {1,2,3,4};" << std::endl;
  outfile << "Plane Surface(1) = {1};" << std::endl;

  outfile << "Physical Line(\"outer_wall\") = {1};" << std::endl;
  outfile << "Physical Line(\"outlet\") = {2};" << std::endl;
  outfile << "Physical Line(\"inner_wall\") = {3};" << std::endl;
  outfile << "Physical Line(\"inlet\") = {4};" << std::endl;

  outfile << "Physical Surface(\"interior\") = {1};" << std::endl;
  outfile << std::endl;

  outfile.close();

  return 0;
}
