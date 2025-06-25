#include "mesh/shape_function/ShapeFunction.hpp"
#include "mesh/shape_function/ModalBasisFactory.hpp"
#include "mesh/std_region/StdRegion.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

namespace detail
{

const std::tuple<PointSetTag, SFTag> ShapeFunctionInstance::undefined =
    std::tuple<PointSetTag, SFTag>(
        PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined),
        SFTag(ElemShape::Undefined, SFunc::Undefined, P0, ModalBasis::Undefined));

// ----------------------------------------------------------------------------

void ShapeFunctionInstance::construct(const std::tuple<PointSetTag, SFTag> sf_key,
                                      ShapeFunctionInstance &sf_instance)
{
  const PointSetTag tmp_std_reg_tag = std::get<0>(sf_key);
  const SFTag sf_tag                = std::get<1>(sf_key);
  sf_instance.m_tag                 = sf_tag;

  // Get the correct prime basis
  ModalBasisFactory::instance_type &prime_basis_factory = ModalBasisFactory::instance();
  const ModalBasisTag prime_basis_tag = ModalBasisTag(sf_tag.prime_basis(), sf_tag.elem_shape());
  // const Uint prime_basis_uid =
  // PrimeBasisTag::keys_to_uid(Modal,tag.key<ElemShape>());
  sf_instance.m_prime_basis = prime_basis_factory.create(prime_basis_tag);

  // THIS IS IMPORTANT: IT SETS THE POLYNOMIAL ORDER OF THE PRIME EXPANSION!!!
  sf_instance.m_prime_basis->set_polynomial_order(sf_tag.poly_order());

  // THE VANDERMONDE MATRIX IS EVALUATED IN THE DOF OF GIVEN REFERENCE
  // ELEMENT: Here we fill the coordinates
  const PointSetTag rt_tag(sf_tag.elem_shape(), sf_tag.poly_order(),
                           tmp_std_reg_tag.ref_topology());
  mesh::StdRegion const std_region        = mesh::StdRegion(rt_tag);
  math::DenseDMat<Real> const &ref_coords = std_region.get().coordinates();
  sf_instance.m_prime_coords.resize(ref_coords.rows(), ref_coords.cols());
  sf_instance.m_prime_coords = ref_coords;

  // And here we compute the Vandermonde matrix
  sf_instance.m_prime_basis->Vandermonde_matrix(sf_instance.m_prime_coords, sf_instance.m_V);
  // sf_instance.m_V.transpose_in_place();

  sf_instance.m_invV.resize(sf_instance.m_V.rows(), sf_instance.m_V.cols());
  sf_instance.m_V.inv(sf_instance.m_invV);

  if (sf_tag.shape_function() == SFunc::Modal)
  {
    const ModalBasisFactory::instance_type::const_product_base_ptr modal_expansion =
        prime_basis_factory.create(ModalBasisTag(ModalBasis::Modal, sf_tag.elem_shape()));
    // THIS IS IMPORTANT: IT SETS THE POLYNOMIAL ORDER OF THE PRIME
    // EXPANSION!!!
    modal_expansion->set_polynomial_order(sf_tag.poly_order());
    modal_expansion->is_leading_expansion_term(sf_instance.m_is_leading_expansion_term);

    // Fill the polynomial degree of each mode
    modal_expansion->mode_poly_deg(sf_instance.m_mode_poly_deg);
  }
  else if (sf_tag.shape_function() == SFunc::Lagrange)
  {
    sf_instance.m_is_leading_expansion_term.resize(sf_instance.m_prime_basis->nb_modes());
    sf_instance.m_is_leading_expansion_term.fill(false);

    sf_instance.m_mode_poly_deg.resize(sf_instance.m_prime_basis->nb_modes());
    sf_instance.m_mode_poly_deg.fill(sf_tag.poly_order());
  }
}

// ----------------------------------------------------------------------------

ShapeFunctionInstance::ShapeFunctionInstance()
{
  m_tag.set_elem_shape(ElemShape::Undefined);
  m_tag.set_shape_function(SFunc::Undefined);
  m_tag.set_poly_order(0);
  m_tag.set_prime_basis(ModalBasis::Undefined);

  m_V.resize(0, 0);
  m_invV.resize(0, 0);
  //  m_V_pt_set.resize(0,0);
  //  m_dV_pt_set.resize(0);

  m_prime_coords.resize(0, 0);
}

// ----------------------------------------------------------------------------

ShapeFunctionInstance::ShapeFunctionInstance(const std::tuple<PointSetTag, SFTag> sf_key,
                                             const SFTag sf_tag)
{
  construct(sf_key, *this);
}

// ----------------------------------------------------------------------------

ShapeFunctionInstance::~ShapeFunctionInstance()
{
}

// ----------------------------------------------------------------------------

const SFTag &ShapeFunctionInstance::tag() const
{
  return m_tag;
}

// ----------------------------------------------------------------------------

Uint ShapeFunctionInstance::nb_dof() const
{
  return m_prime_coords.rows();
}

// ----------------------------------------------------------------------------

Uint ShapeFunctionInstance::topo_dim() const
{
  return m_prime_coords.cols();
}

// ----------------------------------------------------------------------------

const ShapeFunctionInstance::coord_t &ShapeFunctionInstance::ref_coords() const
{
  return m_prime_coords;
}

// ----------------------------------------------------------------------------

void ShapeFunctionInstance::ref_coords(coord_t &coord) const
{
  coord.resize(m_prime_coords.rows(), m_prime_coords.cols());

  coord = m_prime_coords;
}

// ----------------------------------------------------------------------------

void ShapeFunctionInstance::compute_ref_values(const math::DenseDMat<Real> &ref_coords,
                                               math::DenseDMat<Real> &values) const
{
  // For nodal shape functions: first build prime basis, then transform to
  // to desired nodal basis
  if (m_tag.shape_function() == SFunc::Lagrange)
  {
    math::DenseDMat<Real> V_pt_set;
    m_prime_basis->Vandermonde_matrix(ref_coords, V_pt_set);
    // V_pt_set.transpose_in_place();
    values.resize(V_pt_set.rows(), m_invV.cols());
    values = V_pt_set * m_invV;

    // values.transpose_in_place();

    /// CLIP OFF CLOSE-TO-ZERO VALUES AND MAKE SURE THAT THE SUM OF ALL
    /// SF IN EVERY POINT IS EQUAL TO ONE
    Real row_sum = 0.0;
    for (Uint row = 0; row < values.rows(); ++row)
    {
      row_sum = 0.0;
      for (Uint col = 0; col < values.cols(); ++col)
      {
        if (std::abs(values(row, col)) < 1.e-14)
        {
          values(row, col) = 0.0;
        }
        row_sum += values(row, col);
      }

      if (row_sum != 1.0)
      {
        const Real factor = 1.0 / row_sum;
        for (Uint col = 0; col < values.cols(); ++col)
        {
          values(row, col) *= factor;
        }
      } // If sum of sf values in one row is not = 1

    } // Loop over all rows of values
  }

  // -----------------------

  // Modal shape functions: directly compute values
  else if (sf_is_modal(m_tag.shape_function()))
  {
    // Get the correct prime basis
    ModalBasisFactory::instance_type &prime_basis_factory = ModalBasisFactory::instance();
    const ModalBasisTag prime_basis_tag =
        ModalBasisTag(cast_modal_sf_to_modal_basis(m_tag.shape_function()), m_tag.elem_shape());

    const ModalBasisFactory::instance_type::const_product_base_ptr modal_basis =
        prime_basis_factory.create(prime_basis_tag);
    modal_basis->set_polynomial_order(m_tag.poly_order());
    modal_basis->Vandermonde_matrix(ref_coords, values);
  }

  // -----------------------
}

// ----------------------------------------------------------------------------

void ShapeFunctionInstance::compute_ref_derivatives(
    const math::DenseDMat<Real> &ref_coords, std::vector<math::DenseDMat<Real>> &values) const
{
  if (m_tag.shape_function() == SFunc::Lagrange)
  {
    std::vector<math::DenseDMat<Real>> dV_pt_set;
    m_prime_basis->Vandermonde_matrix_derivatives(ref_coords, dV_pt_set);

    values.resize(dV_pt_set.size());
    for (Uint d = 0; d < values.size(); ++d)
    {
      values[d].resize(dV_pt_set[0].rows(), m_invV.cols());
      values[d] = dV_pt_set[d] * m_invV;
    }
  }
  else if (sf_is_modal(m_tag.shape_function()))
  {
    // Get the correct prime basis
    ModalBasisFactory::instance_type &prime_basis_factory = ModalBasisFactory::instance();
    const ModalBasisTag prime_basis_tag =
        ModalBasisTag(cast_modal_sf_to_modal_basis(m_tag.shape_function()), m_tag.elem_shape());

    const ModalBasisFactory::instance_type::const_product_base_ptr modal_basis =
        prime_basis_factory.create(prime_basis_tag);
    modal_basis->set_polynomial_order(m_tag.poly_order());
    modal_basis->Vandermonde_matrix_derivatives(ref_coords, values);
  }
}

// ----------------------------------------------------------------------------

void ShapeFunctionInstance::compute_ref_values(const math::DenseDMat<Real> &ref_coords,
                                               const Uint dim, const Uint local_entity_id,
                                               math::DenseDMat<Real> &values) const
{
  this->compute_ref_values(ref_coords, values);
}

// ----------------------------------------------------------------------------

void ShapeFunctionInstance::compute_ref_derivatives(
    const math::DenseDMat<Real> &ref_coords, const Uint dim, const Uint local_entity_id,
    std::vector<math::DenseDMat<Real>> &values) const
{
  this->compute_ref_derivatives(ref_coords, values);
}

// ----------------------------------------------------------------------------

const math::DenseConstVecView<bool> ShapeFunctionInstance::is_leading_expansion_term() const
{
  const math::DenseConstVecView<bool> let(&m_is_leading_expansion_term[0],
                                          m_is_leading_expansion_term.size(), 1u);
  return let;
}

// ----------------------------------------------------------------------------

const math::DenseConstVecView<Uint> ShapeFunctionInstance::mode_poly_deg() const
{
  const math::DenseConstVecView<Uint> poly_deg(&m_mode_poly_deg[0], m_mode_poly_deg.size(), 1u);
  return poly_deg;
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace sf

} // namespace mesh

} // namespace pdekit
