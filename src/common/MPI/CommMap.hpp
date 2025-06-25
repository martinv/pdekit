#ifndef PDEKIT_Common_MPI_Comm_Map_hpp
#define PDEKIT_Common_MPI_Comm_Map_hpp

#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace common
{

namespace mpi
{

class CommMap
{
  public:
  /// TYPEDEFS

  using idx_type = Int;

  /// STATIC VARIABLES

  static const idx_type INVALID_ID;

  /// Constructor
  CommMap(const std::string &name);

  /// Disabled copy constructor
  CommMap(const CommMap &other) = delete;

  /// Default destructor
  ~CommMap();

  /// Disabled assigment operator
  CommMap &operator=(const CommMap &rhs) = delete;

  /// Get the global dof id corresponding to local dof id
  inline idx_type global_id(const idx_type local_id)
  {
    return m_local_to_global_id_map[local_id];
  }

  /// Get the local id based corresponding to the global id
  inline idx_type local_id(const idx_type global_id) const
  {
    const std::unordered_map<idx_type, idx_type>::const_iterator it =
        m_global_to_local_id_map.find(global_id);
    if (it != m_global_to_local_id_map.end())
    {
      return it->second;
    }
    return INVALID_ID;
  }

  /// Get the name of this CommMap
  const std::string &name() const;

  /// Number of entries in the map
  Uint size() const;

  /// Clear all the internal data
  void clear();

  /// Fill the dof map - we need to know how to construct the relation
  /// local_id
  /// => global_id
  /// and how many degrees of freedom per mesh entity are stored
  void fill(std::vector<idx_type> &global_id);

  private:
  /// Name of this dof map
  std::string m_name;

  /// Direct mapping of dofs - the number of i-th
  /// dof is stored in the i-th row m_direct_dof_id
  std::vector<idx_type> m_local_to_global_id_map;

  /// The inverse map to the map above
  std::unordered_map<idx_type, idx_type> m_global_to_local_id_map;
};

} // namespace mpi

} // namespace common

} // namespace pdekit

#endif
