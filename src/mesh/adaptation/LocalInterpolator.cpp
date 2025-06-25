#include "mesh/adaptation/LocalInterpolator.hpp"
#include "mesh/std_region/StdRegion.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

// ----------------------------------------------------------------------------

const math::DenseConstMatView<Real> LocalInterpolator::transfer_data(
    PointSetTag const std_reg_tag_in, PointSetTag const std_reg_tag_out,
    const math::DenseDMat<Real> &data_in)
{
  interpolator_key search_key(std_reg_tag_in, std_reg_tag_out);

  map_iterator interp_mat_it = m_interp_matrices.find(search_key);

  if (interp_mat_it == m_interp_matrices.end())
  {
    interp_mat_it = build_interpolation_work_matrix(search_key);
  } // If interpolation matrices not found and had to be created

  const math::DenseDMat<Real> &I = interp_mat_it->second;
  m_work_data_out.resize(I.rows() * data_in.cols());
  math::DenseMatView<Real> data_out(m_work_data_out.data(), data_in.cols(), I.rows(),
                                    data_in.cols());

  data_out = I * data_in;

  math::DenseConstMatView<Real> data_out_const(m_work_data_out.data(), data_in.cols(), I.rows(),
                                               data_in.cols());
  return data_out_const;
}

// ----------------------------------------------------------------------------

const math::DenseConstMatView<Real> LocalInterpolator::transfer_data(
    PointSetTag const std_reg_tag_in, PointSetTag const std_reg_tag_out,
    const math::DenseConstMatView<Real> &data_in)
{
  interpolator_key search_key(std_reg_tag_in, std_reg_tag_out);

  map_iterator interp_mat_it = m_interp_matrices.find(search_key);

  if (interp_mat_it == m_interp_matrices.end())
  {
    interp_mat_it = build_interpolation_work_matrix(search_key);
  } // If interpolation matrices not found and had to be created

  const math::DenseDMat<Real> &I = interp_mat_it->second;
  m_work_data_out.resize(I.rows() * data_in.cols());
  math::DenseMatView<Real> data_out(m_work_data_out.data(), data_in.cols(), I.rows(),
                                    data_in.cols());

  data_out = I * data_in;

  math::DenseConstMatView<Real> data_out_const(m_work_data_out.data(), data_in.cols(), I.rows(),
                                               data_in.cols());
  return data_out_const;
}

// ----------------------------------------------------------------------------

typename LocalInterpolator::map_iterator LocalInterpolator::build_interpolation_work_matrix(
    interpolator_key const &search_key)
{
  const PointSetTag std_reg_tag_in  = search_key.first;
  const PointSetTag std_reg_tag_out = search_key.second;

  std::pair<map_iterator, bool> out = m_interp_matrices.insert(
      std::pair<interpolator_key, math::DenseDMat<Real>>(search_key, math::DenseDMat<Real>()));

  StdRegion std_reg_in(std_reg_tag_in);
  StdRegion std_reg_out(std_reg_tag_out);

  map_iterator interp_mat_it = out.first;

  math::DenseDMat<Real> &Imat = interp_mat_it->second;
  Imat.resize(std_reg_out.get().nb_nodes(), std_reg_in.get().nb_nodes());

  if (std_reg_tag_in == std_reg_tag_out)
  {
    Imat.fill(0.0);
    for (Uint i = 0; i < Imat.rows(); ++i)
    {
      Imat(i, i) = 1.0;
    }
  }
  else
  {
    sf::ShapeFunction shape_func;
    const sf::SFTag sf_tag_in(std_reg_tag_in.elem_shape(), SFunc::Lagrange,
                              std_reg_tag_in.poly_order(), ModalBasis::Modal);
    shape_func.change_type(std_reg_tag_in, sf_tag_in);
    shape_func.get().compute_ref_values(std_reg_out.get().coordinates(), Imat);
  }

  return interp_mat_it;
}

// ----------------------------------------------------------------------------

} // namespace adapt

} // namespace mesh

} // namespace pdekit
