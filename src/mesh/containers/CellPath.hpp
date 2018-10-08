#ifndef PDEKIT_Mesh_Cell_Path_hpp
#define PDEKIT_Mesh_Cell_Path_hpp

#include <memory>
#include <vector>

#include "common/PDEKit.hpp"
#include "mesh/MeshIndex.hpp"

namespace pdekit
{

namespace mesh
{

namespace detail
{

/// ---------------------------------------------------------------------------
/// @brief
/// Class representing one part (segment) of cell path. Compresses
/// the information into a single integer
/// ---------------------------------------------------------------------------

class PathSegment
{
  public:
  enum
  {
    max_num_entries = 8,
    max_entry_value = 15
  };

  PathSegment();

  ~PathSegment();

  void set_entry(const Uint pos, const Uint value);

  Uint entry(const Uint pos) const;

  private:
  uint32_t m_value;
};

} // namespace detail

std::ostream &operator<<(std::ostream &os, const detail::PathSegment &segment);

/// ---------------------------------------------------------------------------
/// @brief
/// A path of a cell in mesh as seen from zero level cell. Cells obtained
/// by mesh refinement can be uniquely identified by a sequence of their child
/// indices.
/// Example: suppose that cell nr. 34 is split into 4 children (0,1,2,3) and we
/// are interested in child with index 2. This child is again split into 4
/// children and this time, we're interested in child with index 1. This cell
/// is split into 3 children, and we descend into cell 0. This location of this
/// cell can be recorded as a sequence of numbers { 34 | 2-1-0 }: we start
/// with base-level cell 34 and repeatedly pass through its children
/// and children of children 2,1,0
/// ---------------------------------------------------------------------------

class CellPath;

std::ostream &operator<<(std::ostream &os, const CellPath &path);

class CellPath
{
  private:
  enum
  {
    max_len = detail::PathSegment::max_entry_value * detail::PathSegment::max_num_entries
  };

  public:
  /// Default constructor
  CellPath();

  /// Copy constructor
  CellPath(const CellPath &rhs);

  /// Construct from zero-level cell index and path entries
  CellPath(const FlatIdx zero_level_cell_id, const std::vector<Uint> &entries);

  /// Destructor
  ~CellPath();

  /// Assignment operator
  CellPath &operator=(const CellPath &rhs);

  /// Build cell path given zero-level cell and path entries
  void build(const FlatIdx zero_level_cell_id, const std::vector<Uint> &entries);

  /// Return the size (length) of the path
  Uint size() const;

  /// Return the capacity (max number of entries that the path can store)
  Uint capacity() const;

  /// Reserve sufficient capacity
  void reserve(const Uint capacity);

  /// Set the base cell index
  void set_base_cell_id(const FlatIdx cell_id);

  /// Get the cell id of base cell
  FlatIdx base_cell_id() const;

  /// Push back value into path
  void push_back(const Uint value);

  /// Remove the last entry from the path
  void pop();

  /// Get one entry in the path
  Uint path_entry(const Uint i) const;

  /// Print - debugging format
  void print() const;

  /// Print the path
  friend std::ostream &operator<<(std::ostream &os, const CellPath &path);

  private:
  void reallocate_storage_blocks(const Uint num_blocks);

  /// Index of cell in base mesh
  FlatIdx m_zero_level_cell;

  /// Number of entries in path
  Uint m_num_entries;

  /// The actual path: a series of segments
  // Unique ptr to an array is std::unique_ptr<T[]>,
  // not std::unique_ptr<T> !
  std::unique_ptr<detail::PathSegment[]> m_path;
};

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
