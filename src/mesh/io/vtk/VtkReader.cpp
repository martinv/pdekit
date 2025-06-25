#include "mesh/io/vtk/VtkReader.hpp"
#include "common/StringUtils.hpp"

namespace pdekit
{

namespace mesh
{

namespace vtk
{

// ----------------------------------------------------------------------------

VtkReader::VtkReader()
{
}

// ----------------------------------------------------------------------------

VtkReader::~VtkReader()
{
}

// ----------------------------------------------------------------------------
//                              PRIVATE METHODS
// ----------------------------------------------------------------------------

void VtkReader::reset()
{
  m_nb_nodes    = 0;
  m_nb_elements = 0;
  m_nb_nodes_per_dim.resize(0);
  m_elem_topo_dim.resize(0);
  m_element_type.resize(0);
  m_coords.resize(0);
  m_connectivity.resize(0);
}

void VtkReader::read_coordinates(const Uint geo_dim)
{
  std::cout << "[VTKReader] Reading node coordinates ..." << std::endl;

  m_coords.resize(geo_dim * m_nb_nodes);
  // VTK format stores 3 coordinate component per node, even for 2D meshes,
  // so 'one_vector_node' has to have length 3
  math::DenseDVec<Real> one_node_coords(_3D);

  for (Uint i_node = 0; i_node < m_nb_nodes; i_node++)
  {

    m_infile >> one_node_coords[X];
    m_infile >> one_node_coords[Y];
    m_infile >> one_node_coords[Z];

    for (Uint component = 0; component < geo_dim; component++)
    {

      m_coords[i_node * geo_dim + component] = one_node_coords[component];
    }
  }

  std::cout << "[VTKReader] ... finished reading coords" << std::endl;
}

// ----------------------------------------------------------------------------

void VtkReader::read_celltypes_unstructured_grid(const Uint topo_dim)
{
  std::cout << "[VTKReader] Reading cell types ..." << std::endl;

  PDEKitToVTK pdekit_to_vtk;

  m_nb_nodes_per_dim.resize(topo_dim + 1);
  m_nb_nodes_per_dim.assign(topo_dim + 1, 0);

  m_elem_topo_dim.resize(m_nb_elements);

  m_element_type.resize(topo_dim + 1);

  for (Uint d = 0; d < (topo_dim + 1); ++d)
  {
    m_element_type[d] = std::unique_ptr<std::vector<StdRegion>>(new std::vector<StdRegion>());
  }

  Uint vtk_elem_type;
  StdRegion std_region;

  for (Uint i_elem = 0; i_elem < m_nb_elements; i_elem++)
  {
    // m_infile >> m_element_type[i_elem];
    m_infile >> vtk_elem_type;

    const ElemShape shape = pdekit_to_vtk.element_shape(vtk_elem_type);
    const Uint order      = pdekit_to_vtk.elem_order(vtk_elem_type);
    const PointSetTag std_region_tag(shape, order, PointSetID::Equidist);

    std_region.change_type(std_region_tag);

    const Uint dim      = std_region.get().topo_dim();
    const Uint nb_nodes = std_region.get().nb_nodes();

    m_elem_topo_dim[i_elem] = dim;

    m_element_type[dim]->push_back(std_region);
    m_nb_nodes_per_dim[dim] += nb_nodes;
  }

  /*
  for (Uint dim = 0; dim < m_element_type.size(); ++dim)
  {
    for (Uint i = 0; i < m_element_type[dim]->size(); ++i)
    {
      std::cout << (*m_element_type[dim])[i].get().type_id().as_string() <<
  std::endl;
    }
  }
  */

  std::cout << "[VTKReader] ... finished reading cell types" << std::endl;
}

// ----------------------------------------------------------------------------

} // Namespace vtk

} // Namespace mesh

} // Namespace pdekit
