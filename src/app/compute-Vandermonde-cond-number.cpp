#include <iomanip>
#include <iostream>

#include "common/PDEKit.hpp"
#include "math/DenseDMat.hpp"
#include "math/unary_ops/MatrixCondNumber.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"
#include "mesh/std_region/PointSetTag.hpp"
#include "mesh/std_region/StdRegion.hpp"

using namespace pdekit;

int main()
{

  const ElemShape elem_shape    = ElemShape::Triag;
  const PointSetID ref_topology = PointSetID::Equidist;
  const SFunc sf_type_id        = SFunc::Modal;

  const Uint p_min = 1;
  const Uint p_max = 15;

  mesh::StdRegion ref_elem;
  mesh::sf::ShapeFunction sf;
  math::DenseDMat<Real> Vandermonde;

  std::cout.precision(10);
  std::cout.setf(std::ios::fixed);

  for (Uint p = p_min; p <= p_max; ++p)
  {
    const mesh::PointSetTag ref_elem_tag = mesh::PointSetTag(elem_shape, p, ref_topology);
    ref_elem.change_type(ref_elem_tag);

    const mesh::sf::SFTag sf_tag = mesh::sf::SFTag(elem_shape, sf_type_id, p, ModalBasis::Modal);
    sf.change_type(ref_elem_tag, sf_tag);

    Vandermonde.resize(ref_elem.get().nb_nodes(), sf.get().nb_dof());

    sf.get().compute_ref_values(ref_elem.get().coordinates(), Vandermonde);

    // std::cout << "Nb. of dof in ref. elem = " <<
    // ref_elem.get().nb_nodes() << std::endl; std::cout << "Nb. of shape
    // functions = " << sf.get().nb_dof()
    // << std::endl; std::cout << "Vandermonde = " << std::endl <<
    // Vandermonde
    // << std::endl;

    std::cout << p << " " << std::setfill(' ') << std::setw(25)
              << math::cond_number_1_norm(Vandermonde) << std::endl;
  }

  return 0;
}
