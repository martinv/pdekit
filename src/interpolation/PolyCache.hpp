#ifndef PDEKIT_Interpolation_Polynomial_Cache_hpp
#define PDEKIT_Interpolation_Polynomial_Cache_hpp

#include "common/DataMap.hpp"
#include "interpolation/FEValues.hpp"
#include "interpolation/PolyCacheCellView.hpp"
#include "interpolation/PolyCacheDescription.hpp"

namespace pdekit
{

namespace interpolation
{

class PolyCache
{
  private:
  /// DATA TYPES

  struct PolyData
  {
    /// Constructor
    PolyData();

    /// Copy constructor
    PolyData(const PolyData &rhs);

    /// Destructor
    ~PolyData();

    /// Type of standard region
    mesh::DiscreteElemKey pdata_type;

    /// Shape function values - Vandermonde matrices
    interpolation::FEValues ref_fe_values;

    /// Number of entries of this type
    Uint nb_cached_values;
  };

  public:
  using cellwise_view = PolyCacheCellView<PolyData>;

  /// Default constructor
  PolyCache();

  /// Default constructor
  ~PolyCache() = default;

  /// Set the maximum number of entries
  void allocate(const PolyCacheDescription &cache_desc);

  /// Put new value in the cache
  void push_back_to_buffer(mesh::DiscreteElemKey key);

  /// Empty the data in the buffer
  void flush();

  /// Clear all internal data
  void clear();

  /// Return nb. values that have been pushed to buffer so far
  Uint nb_values_in_buffer() const;

  /// Return the max number of values that the buffer can take
  Uint capacity() const;

  /// Return the values of i-th polynomial stored
  const interpolation::FEValues &fe_values(const Uint i) const;

  /// Get one block of values
  cellwise_view const cellwise_values(const Uint idx) const;

  /// Get key for the i-th value
  const mesh::DiscreteElemKey key(const Uint idx) const;

  /// Position of i-th cell in its data block
  Uint position_in_block(const Uint idx) const;

  /// Print for debugging
  void print() const;

  private:
  enum
  {
    PDATA        = 0,
    POS_IN_PDATA = 1
  };

  /// METHODS
  void create_new_poly_data(mesh::DiscreteElemKey const &block_type);

  /// DATA

  /// All available poly data (each only once)
  std::vector<PolyData> m_poly_data;

  /// Vector of all cached values:
  /// The i-th entry holds the index in m_poly_data
  /// to know which type of poly_data this entry corresponds
  std::vector<Uint> m_cached_values;

  /// Linear index of metric data, consists of a vector of pairs
  /// [PointSetTag,Uint]
  /// Value m_mdata_index[i] is a pair [FIRST,SECOND] and corresponds to
  /// buffer data which can be found in m_poly_data[FIRST],
  /// cell block on position [SECOND] in this particular PolyData
  std::vector<std::tuple<Uint, Uint>> m_data_index;
};

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
