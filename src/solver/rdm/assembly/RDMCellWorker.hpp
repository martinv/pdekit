#ifndef PDEKIT_Solver_RDM_Assembly_RDM_Cell_Worker_hpp
#define PDEKIT_Solver_RDM_Assembly_RDM_Cell_Worker_hpp

#include "interpolation/FunctionSpace.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

namespace detail
{

template <typename MeshConfig, typename Physics, typename SchemeTraits>
class RDMCellWorker
{
  public:
  using f_space   = interpolation::FunctionSpace<MeshConfig>;
  RDMCellWorker() = default;

  virtual ~RDMCellWorker();

  private:
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
RDMCellWorker<MeshConfig, Physics, SchemeTraits>::~RDMCellWorker()
{
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit
