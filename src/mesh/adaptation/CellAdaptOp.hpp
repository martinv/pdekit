#ifndef PDEKIT_Mesh_Adaptation_Cell_Split_Strategy_hpp
#define PDEKIT_Mesh_Adaptation_Cell_Split_Strategy_hpp

#include <array>
#include <functional>
#include <vector>

#include "common/ArrayView.hpp"
#include "common/Flyweight.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseVecView.hpp"
#include "mesh/adaptation/CellAdaptOpTag.hpp"
#include "mesh/local_topology/CellSubdomainTag.hpp"
#include "mesh/local_topology/TraceEntityTuple.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"
#include "mesh/std_region/StdRegion.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

class CellAdaptOpInstance
{
  public:
  /// STATIC VARIABLES

  static const CellAdaptOpTag undefined;

  /// Static method to fill a permutation
  static void construct(const CellAdaptOpTag key, CellAdaptOpInstance &concrete_op);

  /// Default constructor
  CellAdaptOpInstance();

  /// Constructor permutation given interpolation point set id and permutation
  /// code
  CellAdaptOpInstance(CellAdaptOpTag &key);

  /// Copy constructor
  CellAdaptOpInstance(const CellAdaptOpInstance &other_split_strategy);

  /// Assignment operator
  CellAdaptOpInstance &operator=(const CellAdaptOpInstance &other_split_strategy);

  /// Destructor
  ~CellAdaptOpInstance();

  /// Return the interpolation point set type that
  /// this reference element represents
  CellAdaptOpTag cell_adapt_op_tag() const;

  /// Return the number of parent facets
  Uint nb_parent_facets() const;

  /// Return the number of child elements (obtained after splitting the
  /// parent)
  Uint nb_child_elems() const;

  /// Return the adaptation operation id that will be applied to facet f
  CellTransform parent_facet_adapt_op_id(const Uint f) const;

  /// Return the shape of parent facet
  ElemShape parent_facet_shape(const Uint f) const;

  /// Return the vector of child shapes
  ElemShape child_elem_shape(const Uint c) const;

  /// Return the number of child facets embedded in parent face
  Uint nb_trace_child_facets(const Uint f) const;

  /// Return the number of internal child facets
  Uint nb_internal_child_facets() const;

  /// Return one internal child facet
  const TraceIncidences internal_child_facet(const Uint f) const;

  /// Return the vector of parent-child incidences on facet
  const TraceIncidences parent_child_incidences_on_facet(const Uint f) const;

  /// Get the local id of parent face which contains given facet of child
  /// @return a tuple holding the local id of the parent face (1st entry)
  ///         and also the position of child facet in the parent facet (2nd
  ///         entry)
  const std::tuple<SUint, SUint> containing_parent_facet_id(const Uint child_id,
                                                            const Uint child_facet_id) const;

  /// Get the trace incidence block which contains given facet of child
  /// @return a tuple holding the local id of the parent face (1st entry)
  ///         and also the position of child facet in the parent facet (2nd
  ///         entry)
  Uint containing_parent_interior_facet(const Uint child_id, const Uint child_facet_id) const;

  /// Compute nodal coordinates of child elements in reference space
  /// @param parent_type ... type of parent cell. This decides polynomial
  /// order and also
  ///                       type of nodal distribution (Equidist or Warpblend)
  /// @param child_coords ... vector of matrices with coordinates of child
  /// elements
  ///                         The length of the vector is equal to the number
  ///                         of child elements
  void compute_child_coords(const PointSetTag parent_type,
                            std::vector<math::DenseDMat<Real>> &child_coords) const;

  /// Perform linear transformation of coordinates so that they are all mapped
  /// inside one of parent's children
  /// @param coords_in ... matrix that stores (row-wise) the coordinates of
  /// each input
  ///                      node before transformation
  /// @param child_id   ... index of child whose transformer should be used
  ///                       A transformer applies a linear transform that maps
  ///                       coordinates in parent reference space inside the
  ///                       child region
  /// @param coords_out ... vector of matrices with coordinates of child
  /// elements
  ///                         The length of the vector is equal to the number
  ///                         of child elements
  void transform_coords(const math::DenseDMat<Real> &parent_ref_coords,
                        const math::DenseDMat<Real> &coords_in, const Uint child_id,
                        math::DenseDMat<Real> &coords_out) const;

  /// Print the internal data
  void print() const;

  private:
  /// TYPEDEFS
  typedef std::function<void(math::DenseDMat<Real> const &ref_coords,
                             math::DenseConstVecView<Real> const &coord_in,
                             math::DenseVecView<Real> &coord_out)>
      coord_transform_function;

  /// Type of reference topology
  CellAdaptOpTag m_cell_adapt_op_tag;

  /// Number of parent facets
  Uint m_nb_parent_facets;

  /// Number of child elements
  Uint m_nb_child_elems;

  /// Reference coordinates of parent (stored row-wise)
  math::DenseDMat<Real> m_parent_ref_coords;

  /// Type of adaptation operation applied to each facet
  std::vector<CellTransform> m_facet_adapt_op_id;

  /// Cell shapes after the split is performed
  std::vector<ElemShape> m_child_elem_shapes;

  /// Vector of ALL local child incidences within parent element
  std::vector<IncidenceEntry> m_local_internal_child_incidences;

  /// Vector of all permutations of all local child facets which do not lie
  /// on the facets of the parent
  std::vector<EntityDofRealign> m_local_internal_child_permutations;

  /// Vector of offsets marking multiple incidences of each internal face
  std::vector<Uint> m_local_facet_offsets;

  /// Vector of local child incidences on the facets of parent element
  std::vector<std::vector<IncidenceEntry>> m_parent_child_incidences_on_facets;

  /// Vector of local facet permutations on the facets of parent element
  std::vector<std::vector<EntityDofRealign>> m_parent_child_permutations_on_facets;

  /// Vector of functions that transform the node coordinates of parent
  /// element into node coordinates of child element
  std::vector<coord_transform_function> m_coord_transformers;
};

// ----------------------------------------------------------------------------

inline CellAdaptOpTag CellAdaptOpInstance::cell_adapt_op_tag() const
{
  return m_cell_adapt_op_tag;
}

// ----------------------------------------------------------------------------

inline Uint CellAdaptOpInstance::nb_parent_facets() const
{
  return m_nb_parent_facets;
}

// ----------------------------------------------------------------------------

inline Uint CellAdaptOpInstance::nb_child_elems() const
{
  return m_nb_child_elems;
}

// ----------------------------------------------------------------------------

inline CellTransform CellAdaptOpInstance::parent_facet_adapt_op_id(const Uint f) const
{
  return m_facet_adapt_op_id[f];
}

// ----------------------------------------------------------------------------

inline ElemShape CellAdaptOpInstance::parent_facet_shape(const Uint f) const
{
  const PointSetTag parent_tag(m_cell_adapt_op_tag.elem_shape(), P1, PointSetID::Equidist);
  const StdRegion parent_std_reg(parent_tag);
  const Uint parent_dim       = parent_std_reg.get().topo_dim();
  const PointSetTag facet_tag = (*parent_std_reg.get().elem_entity(parent_dim - 1, f)).pt_set_id();
  return facet_tag.elem_shape();
}

// ----------------------------------------------------------------------------

inline ElemShape CellAdaptOpInstance::child_elem_shape(const Uint c) const
{
  return m_child_elem_shapes[c];
}

// ----------------------------------------------------------------------------

inline Uint CellAdaptOpInstance::nb_trace_child_facets(const Uint f) const
{
  return m_parent_child_incidences_on_facets[f].size();
}

// ----------------------------------------------------------------------------

inline Uint CellAdaptOpInstance::nb_internal_child_facets() const
{
  return (m_local_facet_offsets.size() <= 1 ? 0 : m_local_facet_offsets.size() - 1);
}

// ----------------------------------------------------------------------------

inline const TraceIncidences CellAdaptOpInstance::internal_child_facet(const Uint f) const
{
  const Uint pos  = m_local_facet_offsets[f];
  const Uint size = m_local_facet_offsets[f + 1] - m_local_facet_offsets[f];

  common::ArrayView<const IncidenceEntry, _1D, Uint> incidences(
      &m_local_internal_child_incidences[pos], size);
  common::ArrayView<const EntityDofRealign, _1D, Uint> permutations(
      &m_local_internal_child_permutations[pos], size);

  const TraceIncidences child_facet(incidences, permutations, f);
  return child_facet;
}

// ----------------------------------------------------------------------------

inline const TraceIncidences CellAdaptOpInstance::parent_child_incidences_on_facet(
    const Uint f) const
{
  common::ArrayView<const IncidenceEntry, _1D, Uint> incidences(
      m_parent_child_incidences_on_facets[f].data(), m_parent_child_incidences_on_facets[f].size());
  common::ArrayView<const EntityDofRealign, _1D, Uint> permutations(
      m_parent_child_permutations_on_facets[f].data(),
      m_parent_child_incidences_on_facets[f].size());

  const TraceIncidences child_facet_incidences(incidences, permutations, f);
  return child_facet_incidences;
}

// ----------------------------------------------------------------------------
// Flyweight policy and flyweight definition of CellAdaptOp
// ----------------------------------------------------------------------------

struct CellAdaptOpFlyweightPolicy
{
  typedef CellAdaptOpTag key_type;
};

// ----------------------------------------------------------------------------

typedef common::Flyweight<CellAdaptOpInstance, CellAdaptOpFlyweightPolicy> CellAdaptOp;

// ----------------------------------------------------------------------------

} // namespace adapt

} // namespace mesh

} // namespace pdekit

#endif
