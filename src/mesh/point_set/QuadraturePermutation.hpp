#ifndef PDEKIT_Mesh_Quadrature_Permutation_hpp
#define PDEKIT_Mesh_Quadrature_Permutation_hpp

#include <iostream>
#include <tuple>
#include <vector>

#include "common/Flyweight.hpp"
#include "mesh/EntityRealignCode.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

class QuadraturePermutationInstance
{
  public:
  /// STATIC VARIABLES

  static const std::tuple<mesh::PointSetTag, Uint, mesh::EntityRealignCode> undefined;

  /// Static method to fill a permutation
  static void construct(const std::tuple<mesh::PointSetTag, Uint, mesh::EntityRealignCode> &key,
                        QuadraturePermutationInstance &quad_permutation);

  /// Default constructor
  QuadraturePermutationInstance();

  /// Constructor permutation given interpolation point set id and permutation
  /// code
  QuadraturePermutationInstance(
      const std::tuple<mesh::PointSetTag, Uint, mesh::EntityRealignCode> &key);

  /// Copy constructor
  QuadraturePermutationInstance(const QuadraturePermutationInstance &other_permutation);

  /// Assignment operator
  QuadraturePermutationInstance &operator=(const QuadraturePermutationInstance &other_permutation);

  /// Destructor
  ~QuadraturePermutationInstance();

  /// Return the interpolation point set type that
  /// this reference element represents
  inline mesh::PointSetTag type_id() const
  {
    return m_quad_tag;
  }

  /// Return the local id in case this is quadrature on faces
  /// of a reference element
  /// Local id is then the number of one of the faces of the parent
  /// element
  inline Uint local_id() const
  {
    return m_local_id;
  }

  /// Return the size of the permutation vector
  inline Uint size() const
  {
    return m_permutation.size();
  }

  /// Return the i-th value
  inline const Uint &vertex(const Uint i) const
  {
    return m_permutation[i];
  }

  /// Return the code of this permutation

  inline const mesh::EntityRealignCode &code() const
  {
    return m_permutation_code;
  }

  /// Print the internal data
  void print() const;

  private:
  /// Type of reference topology
  mesh::PointSetTag m_quad_tag;
  Uint m_local_id;
  mesh::EntityRealignCode m_permutation_code;
  std::vector<Uint> m_permutation;
};

// ----------------------------------------------------------------------------

struct QuadraturePermutationFlyweightPolicy
{
  typedef std::tuple<mesh::PointSetTag, Uint, mesh::EntityRealignCode> key_type;
};

// ----------------------------------------------------------------------------

typedef common::Flyweight<QuadraturePermutationInstance, QuadraturePermutationFlyweightPolicy>
    QuadraturePermutation;

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
