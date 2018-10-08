#ifndef PDEKIT_Mesh_Shape_Function_Modal_Basis_Factory_hpp
#define PDEKIT_Mesh_Shape_Function_Modal_Basis_Factory_hpp

#include "common/FactoryT.hpp"
#include "common/Singleton.hpp"
#include "mesh/shape_function/ModalBasisTag.hpp"
#include "mesh/shape_function/ModalExpansion.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

namespace detail
{

class ModalBasisFactorySetup
{
  public:
  static void setup_instance(common::FactoryT<ModalExpansion, mesh::ModalBasisTag> &instance);
};

} // namespace detail

} // namespace sf

} // namespace mesh

using ModalBasisFactory =
    common::Singleton<common::FactoryT<mesh::sf::ModalExpansion, mesh::ModalBasisTag>,
                      mesh::sf::detail::ModalBasisFactorySetup>;

} // namespace pdekit

#endif
