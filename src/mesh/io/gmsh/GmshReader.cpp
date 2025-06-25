#include "mesh/io/gmsh/GmshReader.hpp"
#include "math/DenseDVec.hpp"

namespace pdekit
{

namespace mesh
{

namespace gmsh
{

// ----------------------------------------------------------------------------

GmshReader::GmshReader()
{
}

// ----------------------------------------------------------------------------

GmshReader::~GmshReader()
{
}

// ----------------------------------------------------------------------------
//                              PRIVATE METHODS
// ----------------------------------------------------------------------------

void GmshReader::reset()
{
  m_max_topo_dim = 0;

  m_nb_nodes = 0;

  m_nb_elements = 0;
  m_coords.resize(0);

  m_inner_domains.clear();
  m_inner_domains.resize(0);

  m_boundary_domains.clear();
  m_boundary_domains.resize(0);

  for (Uint d = 0; d < 4; ++d)
  {
    m_tag_to_elem_count[d].clear();
  }

  for (Uint i = 0; i < m_nb_elem_of_type.size(); ++i)
  {
    m_nb_elem_of_type[i] = 0;
  }
}

// ----------------------------------------------------------------------------

// Read the gmsh file and find out where are the sections containing info about
// physical
// tags, node coordinates and cell connectivity
void GmshReader::get_record_positions()
{
  std::string subdomain_names("$PhysicalNames");
  std::string nodes("$Nodes");
  std::string elements("$Elements");
  // std::string element_data("$ElementData");
  // std::string node_data("$NodeData");
  // std::string element_node_data("$ElementNodeData");

  int p;
  while (!m_infile.eof())
  {
    p = m_infile.tellg();
    getline(m_infile, tempstr);
    if (tempstr.find(subdomain_names) != std::string::npos)
    {
      m_subdomain_names_filepos = p;

      Uint nb_subdomains;
      m_infile >> nb_subdomains;

      std::string tempstr;

      std::vector<SubDomainHeaderData> temp_domains(nb_subdomains);

      for (Uint is = 0; is < temp_domains.size(); ++is)
      {
        m_infile >> temp_domains[is].topo_dim;
        m_infile >> temp_domains[is].idx;
        // The original name of the subdomain in the mesh file has
        // quotes, we want to strip them off
        m_infile >> tempstr;
        temp_domains[is].name = tempstr.substr(1, tempstr.length() - 2);
      }

      for (Uint is = 0; is < temp_domains.size(); ++is)
      {
        std::cout << " Subdomain " << temp_domains[is].name << " , " << temp_domains[is].idx
                  << " , " << temp_domains[is].topo_dim << std::endl;
      }

      // Sort the domains into inner and boundary domains
      Uint max_domain_dim = temp_domains[0].topo_dim;
      for (Uint is = 1; is < temp_domains.size(); ++is)
      {
        max_domain_dim = std::max(max_domain_dim, temp_domains[is].topo_dim);
      }

      for (Uint is = 0; is < temp_domains.size(); ++is)
      {
        if (temp_domains[is].topo_dim == max_domain_dim)
        {
          m_inner_domains.push_back(temp_domains[is]);
        }
        else
        {
          m_boundary_domains.push_back(temp_domains[is]);
        }
      }

    } // case subdomain names found

    else if (tempstr.find(nodes) != std::string::npos)
    {
      m_node_coord_filepos = p;
      m_infile >> m_nb_nodes;
      std::cout << "[GmshReader] The total number of nodes is " << m_nb_nodes << std::endl;
    }

    else if (tempstr.find(elements) != std::string::npos)
    {
      m_elements_filepos = p;

      // Reset the number of all element types to 0
      m_nb_elem_of_type.fill(0);

      m_infile >> m_nb_elements;
      std::cout << "[GmshReader] The total number of elements is " << m_nb_elements << std::endl;

      Uint elem_idx, elem_type, nb_tags, phys_tag;

      /// Let's count:
      /// 1) how many elements are in each domain (how many elements are
      /// associated with
      ///    each physical tag)
      /// 2) how many elements of each type are present

      // Clearing of variables
      for (Uint i = 0; i < m_tag_to_elem_count.size(); ++i)
      {
        m_tag_to_elem_count[i].clear();
      }

      for (Uint et = 0; et < m_nb_elem_of_type.size(); ++et)
      {
        m_nb_elem_of_type[et] = 0;
      }

      PDEKitToGmsh pdekit_to_gmsh;

      // Read the gmsh connectivity the first time
      for (Uint ie = 0; ie < m_nb_elements; ++ie)
      {
        m_infile >> elem_idx;
        m_infile >> elem_type;
        m_infile >> nb_tags;
        m_infile >> phys_tag;
        getline(m_infile, tempstr);

        // Count the number of elements in each domain
        std::map<Uint, Uint> &tag_to_elem_count =
            m_tag_to_elem_count[pdekit_to_gmsh.elem_dim(elem_type)];
        tag_to_elem_count[phys_tag]++;

        // Count elements of each type
        m_nb_elem_of_type[elem_type]++;
      }

      /// Find out what is the maximum topological dimension in this mesh
      m_max_topo_dim = _0D;

      for (Uint et = 0; et < pdekit_to_gmsh.nb_elem_types(); ++et)
      {
        if (m_nb_elem_of_type[et] > 0)
        {
          m_max_topo_dim = std::max(pdekit_to_gmsh.elem_dim(et), m_max_topo_dim);
        }
      }

      std::cout << "The maximum topological dimension of this mesh is " << m_max_topo_dim
                << std::endl;

    } // case reading connectivity

  } // while

  m_infile.clear();
}

// ----------------------------------------------------------------------------

void GmshReader::read_coordinates(const Uint geo_dim)
{
  std::cout << "[GmshReader] Reading node coordinates ..." << std::endl;
  tempstr = "";

  m_infile.seekg(m_node_coord_filepos, std::ios::beg);
  getline(m_infile, tempstr);
  getline(m_infile, tempstr);

  // For the moment, we suppose that the coordinates are numbered
  // starting from 1
  m_coords.resize(m_nb_nodes * geo_dim);
  // Gmsh stores 3 coordinate component per node, even for 2D meshes,
  // so 'one_vector_node' has to have length 3
  math::DenseDVec<Real> one_node_coords(_3D);

  Uint node_id;

  for (Uint node_count = 0; node_count < m_nb_nodes; ++node_count)
  {
    m_infile >> node_id;
    node_id--;

    m_infile >> one_node_coords[X];
    m_infile >> one_node_coords[Y];
    m_infile >> one_node_coords[Z];

    for (Uint component = 0; component < geo_dim; ++component)
    {
      m_coords[node_id * geo_dim + component] = one_node_coords[component];
    }
  }
  //   std::cout << geo_coords << std::endl;
  std::cout << "[GmshReader] ... finished reading coords" << std::endl;
}

// ----------------------------------------------------------------------------

} // Namespace gmsh

} // Namespace mesh

} // Namespace pdekit
