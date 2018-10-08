#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

// ----------------------------------------------------------------------------

double F_eval(const double x, const double y, const double a)
{
  const double gamma = 1.4;
  // Note that the domain of log((1+a)/(1-a)) is the interval -1 < a < 1
  // a is speed of sound -> negative values don't make physical sense -> 0 < a
  // <
  // 1
  const double J = 1. / a + 1. / (3. * a * a * a) + 1. / (5. * a * a * a * a * a) -
                   0.5 * std::log((1. + a) / (1. - a));
  const double rho = std::pow(a, (2. / (gamma - 1)));
  const double V2  = 2. * (1. - a * a) / (gamma - 1);

  const double result = (x - 0.5 * J) * (x - 0.5 * J) + y * y - 1. / (4. * rho * rho * V2 * V2);

  return result;
}

// ----------------------------------------------------------------------------

double Newton_isotach_solve(const double x, const double y)
{
  // First find a good initial solve
  double a_test        = 0.1;
  const double da_test = 0.01;

  double err_min = std::abs(F_eval(x, y, a_test));

  double a   = a_test;
  double err = err_min;

  while (a_test < 1.0) // See remark in function F_eval about admissible values of a_test
  {
    const double err_test = std::abs(F_eval(x, y, a_test));
    if (err_test < err)
    {
      a   = a_test;
      err = err_test;
    }
    a_test += da_test;
  }

  // std::cout << "Initial a = " << a << std::endl;
  // std::cout << "Initial error = " << err << std::endl;

  const double tol      = 1.e-14;
  const size_t max_iter = 20;
  err                   = 2. * tol;

  // double a = 0.5;
  // double err = 2. * tol;
  size_t iter = 0;

  while ((err > tol) && (iter < max_iter))
  {
    const double da   = std::max(1.e-10, 1.e-7 * std::abs(a));
    const double dFda = (F_eval(x, y, a + da) - F_eval(x, y, a - da)) / (2. * da);

    // std::cout << "a[" << iter << "] = " << a << ", err = " << err <<
    // std::endl;;

    a   = a - F_eval(x, y, a) / dFda;
    err = std::abs(F_eval(x, y, a));
    iter++;
  }
  // std::cout << "a[final] = " << a << ", err = " << err << std::endl;;

  return a;
}

// ----------------------------------------------------------------------------

void write_node_data_header(std::ostream &file, const std::string &field_name)
{
  file << "$NodeData" << std::endl;
  file << "1" << std::endl;
  file << "\"" << field_name << "\"" << std::endl;
  file << "1" << std::endl;
  file << "0" << std::endl;
  file << "3" << std::endl;
  file << "0" << std::endl;
  file << "1" << std::endl;
}

// ----------------------------------------------------------------------------

int main(int argc, const char *argv[])
{
  if (argc != 2)
  {
    std::cout << "  Use as " << argv[0] << " file.msh \n" << std::endl;
    exit(2);
  }

  std::string mesh_infilename = std::string(argv[1]);

  std::ifstream infile;
  std::ofstream outfile;
  std::string buffer = "";

  size_t n_nodes;
  int node_id;
  double x, y, z;

  std::vector<double> density;
  std::vector<std::tuple<double, double>> vel;
  std::vector<double> pressure;
  std::vector<double> Ma;

  infile.open(mesh_infilename.c_str());

  while (buffer != "$Nodes")
  {
    getline(infile, buffer);
  }

  infile >> n_nodes;

  density.resize(n_nodes);
  vel.resize(n_nodes);
  pressure.resize(n_nodes);
  Ma.resize(n_nodes);

  const double gamma = 1.4;

  for (int i = 0; i < n_nodes; ++i)
  {
    infile >> node_id;
    infile >> x;
    infile >> y;
    infile >> z;

    const double a   = Newton_isotach_solve(x, y);
    const double rho = std::pow(a, 2. / (gamma - 1.0));

    const double J = 1. / a + 1. / (3. * a * a * a) + 1. / (5. * a * a * a * a * a) -
                     0.5 * std::log((1. + a) / (1. - a));
    const double V2 = 2. * (1. - a * a) / (gamma - 1);
    const double V  = std::sqrt(V2);

    // Use expression for x(,V) to solve for k:
    const double k2 = 2. / (1. / V2 - 2. * rho * x + J * rho);

    // We assume only positive k
    const double k = std::sqrt(k2);

    // We know k = 1/psi = V/sin(theta), where theta ... hodograph
    // transformation angle Due to roundoff errors, it can happen that
    // sin_theta is slightly bigger than 1
    const double sin_theta = std::min(1.0, V / k);
    const double cos_theta = std::sqrt(1. - sin_theta * sin_theta);

    const double u_sign = (y >= 0.0) ? -1.0 : 1.0;

    const double u = u_sign * V * cos_theta;
    const double v = -V * sin_theta;

    density[node_id - 1]          = rho;
    std::get<0>(vel[node_id - 1]) = u;
    std::get<1>(vel[node_id - 1]) = v;

    pressure[node_id - 1] = 1. / gamma * std::pow(a, 2. * gamma / (gamma - 1));
    Ma[node_id - 1]       = V / a;
  }
  infile.close();

  outfile.open("ringleb_exact.msh");

  // Write density
  write_node_data_header(outfile, "rho");
  outfile << density.size() << std::endl;
  for (int i = 0; i < density.size(); ++i)
  {
    outfile << i + 1 << " " << density[i] << std::endl;
  }
  outfile << "$EndNodeData" << std::endl;

  // Write u-velocity component
  write_node_data_header(outfile, "u");
  outfile << vel.size() << std::endl;
  for (int i = 0; i < vel.size(); ++i)
  {
    outfile << i + 1 << " " << std::get<0>(vel[i]) << std::endl;
  }
  outfile << "$EndNodeData" << std::endl;

  // Write v-velocity component
  write_node_data_header(outfile, "v");
  outfile << vel.size() << std::endl;
  for (int i = 0; i < vel.size(); ++i)
  {
    outfile << i + 1 << " " << std::get<1>(vel[i]) << std::endl;
  }
  outfile << "$EndNodeData" << std::endl;

  // Write pressure
  write_node_data_header(outfile, "pressure");
  outfile << pressure.size() << std::endl;
  for (int i = 0; i < pressure.size(); ++i)
  {
    outfile << i + 1 << " " << pressure[i] << std::endl;
  }
  outfile << "$EndNodeData" << std::endl;

  // Write Mach number
  write_node_data_header(outfile, "Ma");
  outfile << Ma.size() << std::endl;
  for (int i = 0; i < Ma.size(); ++i)
  {
    outfile << i + 1 << " " << Ma[i] << std::endl;
  }
  outfile << "$EndNodeData" << std::endl;

  outfile.close();

  return 0;
}
