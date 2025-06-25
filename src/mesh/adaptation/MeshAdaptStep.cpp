#include "mesh/adaptation/MeshAdaptStep.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

MeshAdaptStep::MeshAdaptStep()
    : m_adapt_type(AdaptationType::h), m_h_strategy(h_AdaptStrategy::w_hanging_nodes)
{
  m_cell_transform.resize(0);
  m_cell_p_order.resize(0);
  m_nb_coarsening_ops = 0;
  m_nb_do_nothing_ops = 0;
}

// ----------------------------------------------------------------------------

MeshAdaptStep::~MeshAdaptStep()
{
}

// ----------------------------------------------------------------------------

AdaptationType MeshAdaptStep::adapt_type() const
{
  return m_adapt_type;
}

// ----------------------------------------------------------------------------

h_AdaptStrategy MeshAdaptStep::h_adapt_strategy() const
{
  return m_h_strategy;
}

// ----------------------------------------------------------------------------

void MeshAdaptStep::define_p_adapt_ops(const std::vector<Uint> &p_order)
{
  m_adapt_type = AdaptationType::p;

  m_cell_p_order.resize(p_order.size());
  std::copy(p_order.cbegin(), p_order.cend(), m_cell_p_order.begin());

  // Set the h-adaptation vector size to 0 - only p adaptation is considered
  m_cell_transform.resize(0);

  m_nb_coarsening_ops = 0;
  m_nb_do_nothing_ops = 0;
}

// ----------------------------------------------------------------------------

void MeshAdaptStep::define_p_adapt_ops(const std::vector<Uint> &old_p_order,
                                       const std::vector<Uint> &proposed_p_order,
                                       const std::vector<Real> &quality_measure,
                                       const Real refinement_ratio)
{
  m_adapt_type = AdaptationType::p;

  // First entry: value of quality measure
  // Second entry: active cell index corresponding to that quality measure
  std::vector<std::pair<Real, Uint>> quality_measure_sort;

  quality_measure_sort.resize(old_p_order.size());

  for (Uint i = 0; i < old_p_order.size(); ++i)
  {
    quality_measure_sort[i].first  = quality_measure[i];
    quality_measure_sort[i].second = i;
  }

  std::sort(quality_measure_sort.begin(), quality_measure_sort.end());

  const Uint last_cell_to_refine = static_cast<Uint>(refinement_ratio * old_p_order.size());

  m_cell_p_order.resize(old_p_order.size());

  for (Uint i = 0; i < quality_measure_sort.size(); ++i)
  {
    const Uint active_cell_idx = quality_measure_sort[i].second;
    if (i < last_cell_to_refine)
    {
      m_cell_p_order[active_cell_idx] = proposed_p_order[active_cell_idx];
    }
    else
    {
      m_cell_p_order[active_cell_idx] = old_p_order[active_cell_idx];
    }
  }

  // Set the h-adaptation vector size to 0 - only p adaptation is considered
  m_cell_transform.resize(0);

  m_nb_coarsening_ops = 0;
  m_nb_do_nothing_ops = 0;
}

// ----------------------------------------------------------------------------

const std::vector<CellTransform> &MeshAdaptStep::adapt_ops() const
{
  return m_cell_transform;
}

// ----------------------------------------------------------------------------

const std::vector<Uint> &MeshAdaptStep::cell_poly_orders() const
{
  return m_cell_p_order;
}

// ----------------------------------------------------------------------------

bool MeshAdaptStep::is_pure_coarsening() const
{
  if (m_adapt_type != AdaptationType::h)
  {
    return false;
  }

  return ((m_cell_transform.size() == (m_nb_coarsening_ops + m_nb_do_nothing_ops)));
}
// ----------------------------------------------------------------------------

Uint MeshAdaptStep::nb_coarsening_ops() const
{
  return m_nb_coarsening_ops;
}

// ----------------------------------------------------------------------------

void MeshAdaptStep::print_ops() const
{
  if (m_adapt_type == AdaptationType::h)
  {
    Uint nb_coarsening_ops       = 0;
    Uint nb_do_nothing_ops       = 0;
    Uint nb_refinement_ops       = 0;
    Uint nb_aniso_refinement_ops = 0;

    for (Uint i = 0; i < m_cell_transform.size(); ++i)
    {
      if (m_cell_transform[i] == CellTransform::COARSEN)
      {
        nb_coarsening_ops++;
      }
      if (m_cell_transform[i] == CellTransform::NO_TRANS)
      {
        nb_do_nothing_ops++;
      }
      if (CellTransformTraits::is_refinement(m_cell_transform[i]))
      {
        nb_refinement_ops++;
      }

      if (CellTransformTraits::is_aniso_refinement(m_cell_transform[i]))
      {
        nb_aniso_refinement_ops++;
      }
    }

    std::cout << "MeshAdaptStep: operation statistics:" << std::endl;
    std::cout << "   > " << nb_coarsening_ops << " coarsening operations" << std::endl;
    std::cout << "   > " << nb_do_nothing_ops << " do-nothing operations" << std::endl;
    std::cout << "   > " << nb_refinement_ops << " refinement operations ("
              << nb_aniso_refinement_ops << " anisotropic, "
              << nb_refinement_ops - nb_aniso_refinement_ops << " isotropic)" << std::endl;
    std::cout << "   > TOTAL: " << m_cell_transform.size() << " h-adaptation operations"
              << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
