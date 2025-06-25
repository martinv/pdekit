#ifndef PDEKIT_Mesh_Point_Set_Base_hpp
#define PDEKIT_Mesh_Point_Set_Base_hpp

#include <memory>
#include <string>

#include "math/DenseConstMatView.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "mesh/EntityRealignCode.hpp"

namespace pdekit
{

namespace mesh
{

class StdPointSetBase
{
  public:
  /// TYPEDEFS:
  typedef std::shared_ptr<StdPointSetBase> ptr;
  typedef std::shared_ptr<StdPointSetBase const> const_ptr;

  /// Default constructor
  StdPointSetBase();

  /// Deleted copy constructor
  StdPointSetBase(const StdPointSetBase &other_base) = delete;

  /// Destructor
  virtual ~StdPointSetBase() = 0;

  /// Deleted assignment operator
  StdPointSetBase &operator=(const StdPointSetBase &rhs) = delete;

  static std::string type_name()
  {
    return "StdPointSetBase";
  }

  virtual std::string name() const = 0;

  /// Order of polynomial which this quadrature integrates exactly
  virtual Uint order() const = 0;

  /// Topological dimension of element for which this quadrature can be
  /// applied
  virtual Uint dim() const = 0;

  /// Topological codimension of quadrature. For example, if the quadrature
  /// has points on faces of 3D element, i.e. on 2D boundary of 3D object,
  /// then the codimension is codim = 3 - 2 = 1
  virtual Uint codim() const = 0;

  /// Return the number of local entities on which this quadrature
  /// has points defined
  virtual Uint nb_local_entities() const = 0;

  /// Return the number of quadrature points
  virtual Uint size(const Uint local_idx = 0) const = 0;

  /// Return a matrix containing all the reference coordinates
  virtual void reference_coords(math::DenseDMat<Real> &coords, const Uint local_idx = 0) const = 0;

  /// Return a vector containing the reference coordinates
  virtual void weights(math::DenseDVec<Real> &wgt, const Uint local_idx = 0) const = 0;

  /// Fill a vector which represents a permutation of the nodes of a point set
  virtual void permutation(const Uint local_id, const mesh::EntityRealignCode &permutation_code,
                           std::vector<Uint> &permutation_vec) = 0;

  protected:
};

} // namespace mesh

} // namespace pdekit

#endif
