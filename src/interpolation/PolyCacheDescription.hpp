#ifndef PDEKIT_Interpolation_Poly_Cache_Description_hpp
#define PDEKIT_Interpolation_Poly_Cache_Description_hpp

#include <tuple>
#include <vector>

#include "common/DataMap.hpp"
#include "mesh/DiscreteElemKey.hpp"
#include "mesh/std_region/PointSetTagExt.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------
// PolyCacheDesription keeps track of how many blocks (where by 'block' we
// mean a block of data corresponding to one DiscreteElemKey) and how
// many entries in each block are present
// ----------------------------------------------------------------------------

class PolyCacheDescription
{
  public:
  /// Constructor
  PolyCacheDescription();

  /// Copy constructor
  PolyCacheDescription(const PolyCacheDescription &rhs);

  /// Destructor
  ~PolyCacheDescription();

  /// ASsignement operator
  const PolyCacheDescription &operator=(const PolyCacheDescription &rhs);

  /// Increase the number of counted for given block type
  void add_to_block_count(const mesh::DiscreteElemKey &type, const Uint nb_values = 1);

  /// Increase the number of counted for given block type
  void add_to_block_count(const mesh::PointSetTagExt &std_reg_key, const mesh::sf::SFTag sf_type,
                          const mesh::PointSetTagExt quad_tag, const Uint nb_values = 1);

  /// Get count for one block type
  Uint block_size(const mesh::DiscreteElemKey &type) const;

  /// Return all available cache block types and their sizes
  std::vector<std::tuple<mesh::DiscreteElemKey, Uint>> all_block_sizes() const;

  /// Get block type of i-th block
  const mesh::DiscreteElemKey &block_type(const Uint block_idx) const;

  /// Get the number of all blocks
  Uint nb_blocks() const;

  /// Get the number of all cached entries
  Uint nb_cached_entries() const;

  /// Reset everything - clear internal data
  void clear();

  /// Print contents
  void print() const;

  private:
  /// For each block, store how many entries are in the block
  std::vector<std::tuple<mesh::DiscreteElemKey, Uint>> m_blocks;
};

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
