#include "mesh/shape_function/ModalBasisFactory.hpp"
#include "common/FactoryPool.hpp"
#include "mesh/shape_function/CarnevaliExpansionLine.hpp"
#include "mesh/shape_function/CarnevaliExpansionTriag.hpp"
#include "mesh/shape_function/DubinerExpansionTetra.hpp"
#include "mesh/shape_function/DubinerExpansionTriag.hpp"
#include "mesh/shape_function/ModalExpansionHexa.hpp"
#include "mesh/shape_function/ModalExpansionLine.hpp"
#include "mesh/shape_function/ModalExpansionQuad.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

namespace detail
{

void ModalBasisFactorySetup::setup_instance(
    common::FactoryT<sf::ModalExpansion, mesh::ModalBasisTag> &instance)
{
  // Register itself in FactoryPool:
  common::FactoryPool::instance_type &fpool = common::FactoryPool::instance();
  std::shared_ptr<common::AbstractFactoryBase> prime_basis_factory =
      ModalBasisFactory::shared_ptr_to_instance();
  fpool.register_factory(sf::ModalExpansion::type_name(), prime_basis_factory);

  // Fill the factory

  using namespace mesh::sf;

  instance.register_builder<ModalExpansionLine>(ModalBasisTag(ModalBasis::Modal, ElemShape::Line));
  instance.register_builder<CarnevaliExpansionLine>(
      ModalBasisTag(ModalBasis::Carnevali, ElemShape::Line));
  instance.register_builder<DubinerExpansionTriag>(
      ModalBasisTag(ModalBasis::Modal, ElemShape::Triag));
  instance.register_builder<CarnevaliExpansionTriag>(
      ModalBasisTag(ModalBasis::Carnevali, ElemShape::Triag));
  instance.register_builder<ModalExpansionQuad>(ModalBasisTag(ModalBasis::Modal, ElemShape::Quad));
  instance.register_builder<DubinerExpansionTetra>(
      ModalBasisTag(ModalBasis::Modal, ElemShape::Tetra));
  instance.register_builder<ModalExpansionHexa>(ModalBasisTag(ModalBasis::Modal, ElemShape::Hexa));
}

} // namespace detail

} // namespace sf

} // namespace mesh

// Force instantantiation of certain classes
template class common::FactoryT<mesh::sf::ModalExpansion, mesh::ModalBasisTag>;
template class common::Singleton<common::FactoryT<mesh::sf::ModalExpansion, mesh::ModalBasisTag>,
                                 mesh::sf::detail::ModalBasisFactorySetup>;

} // namespace pdekit
