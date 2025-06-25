#include <iomanip>

#include "solver/SolverIO.hpp"

namespace pdekit
{

namespace solver
{

// ----------------------------------------------------------------------------

void SolverIO::print_iter_and_res_norm(const Uint iter, const Real CFL,
                                       const math::DenseDVec<Real> &norm)
{
  std::cout << "Iter = " << std::setw(5) << iter << " , CFL = " << std::setw(5)
            << std::setprecision(3) << CFL << " , ";

  std::cout.precision(10);
  std::cout << "res =";
  for (Uint i = 0; i < norm.size(); ++i)
  {
    std::cout << " " << std::setw(15) << norm[i];
  }

  std::cout << " , log(res) =";
  for (Uint i = 0; i < norm.size(); ++i)
  {
    std::cout << " " << std::setw(15) << std::log(norm[i]);
  }
  std::cout << std::endl;
}

// ----------------------------------------------------------------------------

void SolverIO::print_iter_and_res_norm_w_timing(const Uint iter, const Real CFL,
                                                const math::DenseDVec<Real> &norm,
                                                const double cpu_duration)
{
  std::cout << "Iter = " << std::setw(5) << iter << " , CFL = " << std::setw(5)
            << std::setprecision(3) << CFL << " , ";

  std::cout.precision(10);
  std::cout << "res =";
  for (Uint i = 0; i < norm.size(); ++i)
  {
    std::cout << " " << std::setw(15) << norm[i];
  }

  std::cout << " , log(res) =";
  for (Uint i = 0; i < norm.size(); ++i)
  {
    std::cout << " " << std::setw(15) << std::log(norm[i]);
  }
  std::cout << "  cpu_time = " << cpu_duration;
  std::cout << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace solver

} // namespace pdekit
