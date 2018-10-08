#include "mesh/adaptation/CellAdaptOpFactory.hpp"
#include "common/FactoryPool.hpp"
#include "mesh/CellTransform.hpp"
#include "mesh/adaptation/CellAdaptOpQuad.hpp"
#include "mesh/adaptation/CellAdaptOpTriag.hpp"
#include "mesh/adaptation/DoNothingCellAdaptOp.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

namespace internal
{

void CellAdaptOpFactorySetup::setup_instance(
    common::FactoryT<CellAdaptOpBase, CellAdaptOpTag> &instance)
{
  //  std::cout << "SF setup" << std::endl;

  // Register itself in FactoryPool:
  common::FactoryPool::instance_type &fpool = common::FactoryPool::instance();
  std::shared_ptr<common::AbstractFactoryBase> facet_incidence_factory =
      CellAdaptOpFactory::shared_ptr_to_instance();
  fpool.register_factory(CellAdaptOpBase::type_name(), facet_incidence_factory);

  // ------------------------------------------
  // Default: do nothing (do not split the cell
  // ------------------------------------------

  instance.register_builder<DoNothingCellAdaptOpLine>(
      CellAdaptOpTag(ElemShape::Line, CellTransform::NO_TRANS));

  instance.register_builder<DoNothingCellAdaptOpTriag>(
      CellAdaptOpTag(ElemShape::Triag, CellTransform::NO_TRANS));

  instance.register_builder<DoNothingCellAdaptOpQuad>(
      CellAdaptOpTag(ElemShape::Quad, CellTransform::NO_TRANS));

  instance.register_builder<DoNothingCellAdaptOpTetra>(
      CellAdaptOpTag(ElemShape::Tetra, CellTransform::NO_TRANS));

  instance.register_builder<DoNothingCellAdaptOpHexa>(
      CellAdaptOpTag(ElemShape::Hexa, CellTransform::NO_TRANS));

  instance.register_builder<DoNothingCellAdaptOpPrism>(
      CellAdaptOpTag(ElemShape::Prism, CellTransform::NO_TRANS));

  instance.register_builder<DoNothingCellAdaptOpPyramid>(
      CellAdaptOpTag(ElemShape::Pyramid, CellTransform::NO_TRANS));

  // -----------------------------------
  // Strategies to split a triangle cell
  // -----------------------------------

  instance.register_builder<CellAdaptOpTriagUniformRefine>(
      CellAdaptOpTag(ElemShape::Triag, CellTransform::UNIFORM_REFINE));

  instance.register_builder<CellAdaptOpTriagAnisoRefineOrthoFace0>(
      CellAdaptOpTag(ElemShape::Triag, CellTransform::ANISO_REFINE_ORTHO_FACE_0));

  instance.register_builder<CellAdaptOpTriagAnisoRefineOrthoFace1>(
      CellAdaptOpTag(ElemShape::Triag, CellTransform::ANISO_REFINE_ORTHO_FACE_1));

  instance.register_builder<CellAdaptOpTriagAnisoRefineOrthoFace2>(
      CellAdaptOpTag(ElemShape::Triag, CellTransform::ANISO_REFINE_ORTHO_FACE_2));

  // ----------------------------------------
  // Strategies to split a quadrilateral cell
  // ----------------------------------------

  instance.register_builder<CellAdaptOpQuadUniformRefine>(
      CellAdaptOpTag(ElemShape::Quad, CellTransform::UNIFORM_REFINE));
}

} // namespace internal

} // namespace adapt

} // namespace mesh

// Force instantantiation of certain classes
// This has to be done in a namespace that is enclosing common:: and mesh:: (?)
template class common::FactoryT<mesh::adapt::internal::CellAdaptOpBase,
                                mesh::adapt::CellAdaptOpTag>;
template class common::Singleton<
    common::FactoryT<mesh::adapt::internal::CellAdaptOpBase, mesh::adapt::CellAdaptOpTag>,
    mesh::adapt::internal::CellAdaptOpFactorySetup>;

} // namespace pdekit
