#include <ctime>
#include <iostream>
#include <memory>

#include "math/MathConstants.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/MeshManipulator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "solver/rdm/RDTimeUpdate.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::solver;

typedef Cart2D MeshConfig;
typedef Tria<MeshConfig> MeshType;

// ----------------------------------------------------------------------------

int main()
{
  // ------------------------------------------------------
  // READ GEOMETRY MESH, PREPARE SOLUTION MESH
  // ------------------------------------------------------

  MeshType::shared_ptr mesh = std::make_shared<MeshType>("mesh");
  MeshCreator::make_unit_quad(*mesh, "dofs", 20, true);

  const result_of::dof_map_t<MeshConfig> &dofs = *(mesh->dof_storage("dofs"));

  rdm::RDTimeUpdate time_step;

  time_step.setup<MeshConfig>(*mesh, dofs);

  return 0;
}
