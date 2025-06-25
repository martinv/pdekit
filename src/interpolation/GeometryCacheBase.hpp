#ifndef PDEKIT_Interpolation_Geometry_Cache_Base_hpp
#define PDEKIT_Interpolation_Geometry_Cache_Base_hpp

#include <unordered_map>

#include "common/DataMap.hpp"
#include "common/Meta.hpp"
#include "interpolation/FEValues.hpp"
#include "math/DenseConstMatView.hpp"
#include "mesh/CellGeometry.hpp"
#include "mesh/DiscreteElemKey.hpp"
#include "mesh/DofCoordinates.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

struct CacheInsertManual
{
};

struct CacheInsertAutomatic
{
};

// ----------------------------------------------------------------------------

template <typename Derived, typename GeoCacheMode>
class GeometryCacheBase;

// ----------------------------------------------------------------------------

template <typename Derived>
class GeometryCacheBase<Derived, CacheInsertManual>
{
  public:
  /// Default constructor
  GeometryCacheBase() = default;

  /// Default constructor
  ~GeometryCacheBase() = default;

  /// Set the size of buffer
  template <typename DiscreteElemKeyIterator>
  void allocate(const DiscreteElemKeyIterator keys_begin, const DiscreteElemKeyIterator keys_end,
                Uint const nb_blocks)
  {
    this->as_derived().allocate_manual(keys_begin, keys_end, nb_blocks);
  }

  template <Uint GDIM>
  void push_back_to_buffer(const mesh::DofCoordinates<GDIM> &cell_coords,
                           mesh::DiscreteElemKey const key)
  {
    this->as_derived().push_back_to_buffer_manual(cell_coords, key);
  }

  /// Insert generic data to buffer
  template <Uint GDIM>
  void push_back_to_buffer(const mesh::CellGeometry<GDIM> &cell_coords,
                           mesh::DiscreteElemKey const key)
  {
    this->as_derived().push_back_to_buffer_manual(cell_coords, key);
  }

  template <Uint GDIM>
  void push_back_to_buffer_and_interpolate(const mesh::CellGeometry<GDIM> &cell_coords,
                                           const mesh::DiscreteElemKey key_coords,
                                           const mesh::DiscreteElemKey key_interp)
  {
    this->as_derived().push_back_to_buffer_manual_and_interpolate(cell_coords, key_coords,
                                                                  key_interp);
  }

  void push_back_to_buffer(const math::DenseConstMatView<Real> &cell_coords,
                           mesh::DiscreteElemKey const key)
  {
    this->as_derived().push_back_to_buffer_manual(cell_coords, key);
  }

  private:
  inline Derived &as_derived()
  {
    return static_cast<Derived &>(*this);
  }

  inline const Derived &as_derived() const
  {
    return static_cast<Derived const &>(*this);
  }
};

// ----------------------------------------------------------------------------

template <typename Derived>
class GeometryCacheBase<Derived, CacheInsertAutomatic>
{
  public:
  /// Default constructor
  GeometryCacheBase() = default;

  /// Default constructor
  ~GeometryCacheBase() = default;

  template <typename DiscreteElemKeyInserter>
  void record_insertion_pattern(DiscreteElemKeyInserter &inserter, Uint const nb_blocks)
  {
    this->as_derived().record_insertion_pattern_automatic(inserter, nb_blocks);
  }

  private:
  inline Derived &as_derived()
  {
    return static_cast<Derived &>(*this);
  }

  inline const Derived &as_derived() const
  {
    return static_cast<Derived const &>(*this);
  }
};

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
