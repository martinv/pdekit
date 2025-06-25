#ifndef PDEKIT_RDM_Lax_Friedrichs_Residual_Limiter_hpp
#define PDEKIT_RDM_Lax_Friedrichs_Residual_Limiter_hpp

#include "solver/rdm/RDMethodScratchData.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

namespace internal
{
template <typename Physics, Uint NEQ>
struct LFResidualLimiter
{
  inline static void limit_cell_residuals(typename Physics::FluxV const &cell_residual,
                                          typename Physics::FluxV &sum_xp,
                                          RDMethodScratchData<Physics> &MD)
  {
  }

  inline static void limit_facet_residuals(typename Physics::FluxV const &facet_residual,
                                           typename Physics::FluxV &sum_xp,
                                           RDMethodScratchData<Physics, Physics::DIM - 1> &MDLeft,
                                           RDMethodScratchData<Physics, Physics::DIM - 1> &MDRight)
  {
  }
};

// --------------------------------------------------------------------------

template <typename Physics>
struct LFResidualLimiter<Physics, 1u>
{
  inline static void limit_cell_residuals(typename Physics::FluxV const &cell_residual,
                                          typename Physics::FluxV &sum_xp,
                                          RDMethodScratchData<Physics> &MD)
  {
    // ---------------------------
    // Limiting for HO scheme
    // ---------------------------
    // Compute the sum of x+
    sum_xp.fill(0.0);

#if 0
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      // If cell_residual[eq] is NOT ZERO, update the MD.m_elem_node_res[eq].
      // Otherwise (i.e. cell_residual[eq] == 0), don't update anything, because
      // m_sump_xp[eq] must be a NaN (we computed it as
      //    std::max(MD.m_elem_node_res[n][eq] / cell_residual[eq], 0.0);
      // which means we divided by zero!
      if (std::abs(cell_residual[eq]) > 1.e-14)
      {
        for (Uint n = 0; n < MD.m_nb_nodes; ++n)
        {
          // Compute x+ and store it (temporarily) in MD.m_elem_node_res[n][eq]
          MD.m_elem_node_res[n][eq] = std::max(MD.m_elem_node_res[n][eq] / cell_residual[eq], 0.0);
          sum_xp[eq] += MD.m_elem_node_res[n][eq];
        }

        for (Uint n = 0; n < MD.m_nb_nodes; ++n)
        {
          MD.m_elem_node_res[n][eq] = MD.m_elem_node_res[n][eq] / sum_xp[eq] * cell_residual[eq];
        }
      }
      /*
      else
      {
        for (Uint n = 0; n < MD.m_nb_nodes; ++n)
        {
          MD.m_elem_node_res[n][eq] = 0.0;
        }
      }
      */
    }

#endif

    const Real eps = 1.e-12;

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      // If cell_residual[eq] is NOT ZERO, update the
      // MD.m_elem_node_res[eq]. Otherwise (i.e. cell_residual[eq] == 0),
      // don't update anything, because m_sump_xp[eq] must be a NaN (we
      // computed it as
      //    std::max(MD.m_elem_node_res[n][eq] / cell_residual[eq], 0.0);
      // which means we divided by zero!
      if (std::abs(cell_residual[eq]) > 1.e-14)
      {
        for (Uint n = 0; n < MD.m_nb_nodes; ++n)
        {
          // Compute x+ and store it (temporarily) in
          // MD.m_elem_node_res[n][eq]
          MD.m_elem_node_res[n][eq] =
              std::max(MD.m_elem_node_res[n][eq] / cell_residual[eq], 0.0) + eps;
          sum_xp[eq] += MD.m_elem_node_res[n][eq];
        }

        for (Uint n = 0; n < MD.m_nb_nodes; ++n)
        {
          MD.m_elem_node_res[n][eq] = MD.m_elem_node_res[n][eq] / sum_xp[eq] * cell_residual[eq];
        }
      }
    }
  }

  // --------------------------------------------------------------------------

  inline static void limit_facet_residuals(typename Physics::FluxV const &facet_residual,
                                           typename Physics::FluxV &sum_xp,
                                           RDMethodScratchData<Physics, Physics::DIM - 1> &MDLeft,
                                           RDMethodScratchData<Physics, Physics::DIM - 1> &MDRight)
  {
    // ---------------------------
    // Limiting for HO scheme
    // ---------------------------

    sum_xp.fill(0.0);
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      // If cell_residual[eq] is NOT ZERO, update the
      // MD.m_elem_node_res[eq]. Otherwise (i.e. cell_residual[eq] == 0),
      // don't update anything, because m_sump_xp[eq] must be a NaN (we
      // computed it as
      //    std::max(MD.m_elem_node_res[n][eq] / cell_residual[eq], 0.0);
      // which means we divided by zero!
      if (std::abs(facet_residual[eq]) > 1.e-14)
      {
        for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
        {
          // Compute x+ and store it (temporarily) in
          // MD.m_elem_node_res[n][eq]
          MDLeft.m_elem_node_res[n][eq] =
              std::max(MDLeft.m_elem_node_res[n][eq] / facet_residual[eq], 0.0);
          sum_xp[eq] += MDLeft.m_elem_node_res[n][eq];
        }
        for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
        {
          // Compute x+ and store it (temporarily) in
          // MDRight.m_elem_node_res[n][eq]
          MDRight.m_elem_node_res[n][eq] =
              std::max(MDRight.m_elem_node_res[n][eq] / facet_residual[eq], 0.0);
          sum_xp[eq] += MDRight.m_elem_node_res[n][eq];
        }

        for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
        {
          MDLeft.m_elem_node_res[n][eq] =
              MDLeft.m_elem_node_res[n][eq] / sum_xp[eq] * facet_residual[eq];
        }

        for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
        {
          MDRight.m_elem_node_res[n][eq] =
              MDRight.m_elem_node_res[n][eq] / sum_xp[eq] * facet_residual[eq];
        }
      }
      /*
      else
      {
        for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
        {
          MDLeft.m_elem_node_res[n][eq] = 0.0;
        }
        for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
        {
          MDRight.m_elem_node_res[n][eq] = 0.0;
        }
      }
      */
    }
  }
};

// --------------------------------------------------------------------------

} // namespace internal

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
