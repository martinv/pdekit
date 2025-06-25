#ifndef PDEKIT_Mesh_Adaptation_Mesh_Adapt_Sequence_hpp
#define PDEKIT_Mesh_Adaptation_Mesh_Adapt_Sequence_hpp

#include "mesh/adaptation/MeshAdaptStep.hpp"

namespace pdekit
{

namespace mesh
{

template <typename MeshConfig>
class Mesh;

template <typename MeshConfig>
class MeshAdaptSequence
{
  public:
  /// Default constructor
  MeshAdaptSequence();

  /// Default destructor
  ~MeshAdaptSequence();

  /// Return adaptation type
  AdaptationType adapt_type() const;

  /// Return h_adaptation strategy
  h_AdaptStrategy h_adapt_strategy() const;

  /// Define what h-type adaptation
  /// operation should be applied to each cell in mesh
  void define_h_adapt_ops(Tria<MeshConfig> const &mesh, std::vector<CellTransform> const &adapt_op,
                          const h_AdaptStrategy strategy);

  /// Define what h-type adaptation
  /// operation should be applied to each cell in mesh
  void define_h_adapt_ops(Tria<MeshConfig> const &mesh, std::vector<CellTransform> const &adapt_op,
                          std::vector<Real> const &indicator, const Real ratio,
                          const h_AdaptStrategy strategy);

  /// Define what p-type adaptation
  /// operation should be applied to each cell in mesh
  void define_p_adapt_ops(const std::vector<Uint> &p_order);

  /// Define what p-type adaptation
  /// operation should be applied to each cell in mesh
  void define_p_adapt_ops(const std::vector<Uint> &old_p_order,
                          const std::vector<Uint> &proposed_p_order,
                          const std::vector<Real> &quality_measure, const Real refinement_ratio);

  /// Return a vector of adaptation operations
  const std::vector<CellTransform> &adapt_ops() const;

  /// Return a vector of polynomial orders
  const std::vector<Uint> &cell_poly_orders() const;

  /// Return a value indicating whether this is a pure coarsening operation
  bool is_pure_coarsening() const;

  /// Return the number of cells to which coarsening will be applied
  Uint nb_coarsening_ops() const;

  /// Print stored adaptation operations on the screen
  void print_ops() const;

  private:
  MeshAdaptStep m_adapt_step;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
MeshAdaptSequence<MeshConfig>::MeshAdaptSequence()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
MeshAdaptSequence<MeshConfig>::~MeshAdaptSequence()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
AdaptationType MeshAdaptSequence<MeshConfig>::adapt_type() const
{
  return m_adapt_step.adapt_type();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
h_AdaptStrategy MeshAdaptSequence<MeshConfig>::h_adapt_strategy() const
{
  return m_adapt_step.h_adapt_strategy();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshAdaptSequence<MeshConfig>::define_h_adapt_ops(Tria<MeshConfig> const &mesh,
                                                       std::vector<CellTransform> const &adapt_op,
                                                       const h_AdaptStrategy strategy)
{
  m_adapt_step.define_h_adapt_ops(mesh, adapt_op, strategy);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshAdaptSequence<MeshConfig>::define_h_adapt_ops(Tria<MeshConfig> const &mesh,
                                                       std::vector<CellTransform> const &adapt_op,
                                                       std::vector<Real> const &indicator,
                                                       const Real ratio,
                                                       const h_AdaptStrategy strategy)
{
  m_adapt_step.define_h_adapt_ops(mesh, adapt_op, indicator, ratio, strategy);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshAdaptSequence<MeshConfig>::define_p_adapt_ops(const std::vector<Uint> &p_order)
{
  m_adapt_step.define_p_adapt_ops(p_order);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshAdaptSequence<MeshConfig>::define_p_adapt_ops(const std::vector<Uint> &old_p_order,
                                                       const std::vector<Uint> &proposed_p_order,
                                                       const std::vector<Real> &quality_measure,
                                                       const Real refinement_ratio)
{
  m_adapt_step.define_p_adapt_ops(old_p_order, proposed_p_order, quality_measure, refinement_ratio);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const std::vector<CellTransform> &MeshAdaptSequence<MeshConfig>::adapt_ops() const
{
  return m_adapt_step.adapt_ops();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const std::vector<Uint> &MeshAdaptSequence<MeshConfig>::cell_poly_orders() const
{
  return m_adapt_step.cell_poly_orders();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
bool MeshAdaptSequence<MeshConfig>::is_pure_coarsening() const
{
  return m_adapt_step.is_pure_coarsening();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint MeshAdaptSequence<MeshConfig>::nb_coarsening_ops() const
{
  return m_adapt_step.nb_coarsening_ops();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshAdaptSequence<MeshConfig>::print_ops() const
{
  m_adapt_step.print_ops();
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
