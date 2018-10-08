#ifndef PDEKIT_Mesh_Cell_Marker_hpp
#define PDEKIT_Mesh_Cell_Marker_hpp

#include <map>
#include <string>
#include <vector>

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace mesh
{

class CellMarker
{
  public:
  /// Constructor
  CellMarker();

  /// Copy constructor, intentionally not implemented
  CellMarker(const CellMarker &other_marker) = delete;

  /// Assignment operator, intentionally not implemented
  CellMarker &operator=(const CellMarker &other_marker) = delete;

  /// Destructor
  ~CellMarker();

  /// Static copy function to copy one marker to another
  static void copy(const CellMarker &source, CellMarker &target);

  /// Set the name which describes what this marker is used for
  void set_type_name(const std::string &type_name);

  /// Resize the data with maker values
  void fill(const std::vector<std::pair<Uint, std::string>> &marker_names,
            std::vector<Uint> &marker_values);

  /// Resize the data with maker values
  void fill(const std::map<Uint, std::string> &marker_names, std::vector<Uint> &marker_values);

  /// Push back some more values
  void push_back_values(const std::vector<Uint> &new_values);

  /// Return the string representing the marker type
  const std::string type_name() const;

  /// Return the name of marker for given cell
  const std::string name(const Uint c) const;

  /// Return the value for given cell
  Uint value(const Uint c) const;

  /// Number of values of markers
  Uint nb_values() const;

  /// Number of types of markers
  Uint nb_types() const;

  /// Get all marker names
  void all_marker_names(std::map<Uint, std::string> &marker_names) const;

  private:
  /// String representing marker type (e.g. 'material', 'subdomain')
  std::string m_marker_type;

  /// Map marker => name corresponding to marker
  std::map<Uint, std::string> m_names;

  /// Value of marker for each cell
  std::vector<Uint> m_values;
};

// ----------------------------------------------------------------------------

inline Uint CellMarker::value(const Uint c) const
{
  return m_values[c];
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
