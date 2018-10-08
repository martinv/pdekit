#include <cmath>
#include <limits>

#include "interpolation/FEValues.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"
#include "mesh/std_region/StdRegion.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

void FEValues::copy(FEValues const &fe_source, FEValues &fe_target)
{
  fe_target.m_std_region_tag = fe_source.m_std_region_tag;
  fe_target.m_sf_type        = fe_source.m_sf_type;

  fe_target.m_topo_dim = fe_source.m_topo_dim;

  // Copy reference coordinates
  fe_target.m_ref_coord.resize(fe_source.m_ref_coord.rows(), fe_source.m_ref_coord.cols());

  for (Uint r = 0; r < fe_target.m_ref_coord.rows(); ++r)
  {
    for (Uint c = 0; c < fe_target.m_ref_coord.cols(); ++c)
    {
      fe_target.m_ref_coord(r, c) = fe_source.m_ref_coord(r, c);
    }
  }

  // Copy Vandermonde matrix of shape functions
  fe_target.m_V.resize(fe_source.m_V.rows(), fe_source.m_V.cols());
  for (Uint r = 0; r < fe_target.m_V.rows(); ++r)
  {
    for (Uint c = 0; c < fe_target.m_V.cols(); ++c)
    {
      fe_target.m_V(r, c) = fe_source.m_V(r, c);
    }
  }

  // Copy Vandermonde matrices of shape function derivatives
  fe_target.m_dV.resize(fe_source.m_dV.size());
  for (Uint d = 0; d < fe_target.m_dV.size(); ++d)
  {
    fe_target.m_dV[d].resize(fe_source.m_dV[d].rows(), fe_source.m_dV[d].cols());
    for (Uint r = 0; r < fe_target.m_dV[d].rows(); ++r)
    {
      for (Uint c = 0; c < fe_target.m_dV[d].cols(); ++c)
      {
        fe_target.m_dV[d](r, c) = fe_source.m_dV[d](r, c);
      }
    }
  }

  // Copy point set
  fe_target.m_point_set.resize(fe_source.m_point_set.rows(), fe_source.m_point_set.cols());
  for (Uint r = 0; r < fe_target.m_point_set.rows(); ++r)
  {
    for (Uint c = 0; c < fe_target.m_point_set.cols(); ++c)
    {
      fe_target.m_point_set(r, c) = fe_source.m_point_set(r, c);
    }
  }

  // Copy point set
  fe_target.m_pt_weights.resize(fe_source.m_pt_weights.size());
  for (Uint i = 0; i < fe_target.m_pt_weights.size(); ++i)
  {
    fe_target.m_pt_weights[i] = fe_source.m_pt_weights[i];
  }
}

// ----------------------------------------------------------------------------

FEValues::FEValues()
{
}

// ----------------------------------------------------------------------------

FEValues::FEValues(const mesh::PointSetTag std_region_tag, const mesh::sf::SFTag &sf_type)
{
  m_std_region_tag = std_region_tag;
  m_sf_type        = sf_type;
  setup_type();
}

// ----------------------------------------------------------------------------

FEValues::~FEValues()
{
}

// ----------------------------------------------------------------------------

void FEValues::configure(const mesh::PointSetTag std_region_tag, const mesh::sf::SFTag &sf_type)
{
  m_std_region_tag = std_region_tag;
  m_sf_type        = sf_type;
  setup_type();
}

// ----------------------------------------------------------------------------

mesh::PointSetTag FEValues::std_region_id() const
{
  return m_std_region_tag;
}

// ----------------------------------------------------------------------------

mesh::sf::SFTag FEValues::sf_type()
{
  return m_sf_type;
}

// ----------------------------------------------------------------------------

mesh::sf::SFTag FEValues::sf_type() const
{
  return m_sf_type;
}

// ----------------------------------------------------------------------------

Uint FEValues::topo_dim() const
{
  return m_topo_dim;
}

// ----------------------------------------------------------------------------

Uint FEValues::nb_nodes() const
{
  return m_ref_coord.rows();
}

// ----------------------------------------------------------------------------

Uint FEValues::nb_qd_pts() const
{
  return m_pt_weights.size();
}

// ----------------------------------------------------------------------------

const FEValues::coord_t &FEValues::ref_coord() const
{
  return m_ref_coord;
}

// ----------------------------------------------------------------------------

void FEValues::fill_Vandermonde(const math::DenseDMat<Real> &point_set,
                                const math::DenseDVec<Real> &weights, const bool apply_filter)
{
  // I) perform standard filling of Vandermonde matrices

  mesh::sf::ShapeFunction shape_func;
  // shape_func.change_type(m_sf_type.elem_shape(),
  // m_sf_type.shape_function(), Equidist,
  //                       m_sf_type.poly_order(), Modal);

  shape_func.change_type(m_std_region_tag, m_sf_type);

  const Uint nb_nodes = m_ref_coord.rows();

  m_V.resize(point_set.rows(), nb_nodes);
  for (Uint dim = 0; dim < m_topo_dim; ++dim)
  {
    m_dV[dim].resize(point_set.rows(), nb_nodes);
  }

  shape_func.get().compute_ref_values(point_set, m_V);
  shape_func.get().compute_ref_derivatives(point_set, m_dV);

  m_point_set.resize(point_set.rows(), point_set.cols());
  for (Uint r = 0; r < point_set.rows(); ++r)
  {
    for (Uint c = 0; c < point_set.cols(); ++c)
    {
      m_point_set(r, c) = point_set(r, c);
    }
  }

  m_pt_weights.resize(weights.size());

  for (Uint i = 0; i < weights.size(); ++i)
  {
    m_pt_weights[i] = weights[i];
  }

  // If filtering not needed, we're done -> exit here
  if (!apply_filter)
  {
    return;
  }

  // II) Apply filter

  // For the moment filtering is implemented only on quads (other modal bases
  // do not know how to evaluate mode_poly_deg() ....
  if ((m_sf_type.shape_function() == SFunc::Lagrange) &&
      (m_std_region_tag.elem_shape() == ElemShape::Quad))
  {
    mesh::sf::SFTag modal_sf_type = m_sf_type;
    modal_sf_type.set_shape_function(SFunc::Modal);

    shape_func.change_type(m_std_region_tag, modal_sf_type);

    mesh::StdRegion std_region;
    std_region.change_type(m_std_region_tag);

    math::DenseDMat<Real> V_modal, V_modal_inv, filter;
    const Uint nb_modes = shape_func.get().nb_dof();
    V_modal.resize(nb_modes, nb_modes);
    V_modal_inv.resize(nb_modes, nb_modes);
    filter.resize(nb_modes, nb_modes);
    filter.fill(0.0);

    // Compute 'modal' Vandermonde matrix in degrees of freedom...
    shape_func.get().compute_ref_values(std_region.get().coordinates(), V_modal);
    // ... and its inverse
    V_modal.inv(V_modal_inv);

    const math::DenseConstVecView<Uint> mode_p_deg = shape_func.get().mode_poly_deg();

    Uint p_max = 0;
    for (Uint i = 0; i < mode_p_deg.size(); ++i)
    {
      p_max = std::max(p_max, mode_p_deg[i]);
    }

    // This means that only the modes with polynomial order p = p_max - 2
    // or higher will be considered for filtering
    const Uint cutoff = 2;

    const Real eps   = std::numeric_limits<Real>::epsilon();
    const Real alpha = -std::log(eps);

    for (Uint i = 0; i < mode_p_deg.size(); ++i)
    {
      if ((mode_p_deg[i] + cutoff) >= p_max)
      {
        filter(i, i) = std::exp(-alpha * std::pow(1. * mode_p_deg[i] / p_max, 18));
      }
      else
      {
        filter(i, i) = 1.0;
      }
    }

    /*
    std::cout << "V_mod = " << V_modal << std::endl;
    std::cout << "V_mod_inv = " << V_modal_inv << std::endl;
    std::cout << "Filter = " << filter << std::endl;
    */

    math::DenseDMat<Real> V_temp(m_V.rows(), m_V.cols());
    V_temp = m_V;

    /*
    std::cout << "Size of V_mod = (" << V_modal.rows() << " x " <<
    V_modal.cols() << ")" << std::endl; std::cout << "Size of V_mod_inv = ("
    << V_modal_inv.rows() << " x " << V_modal_inv.cols() << ")"
    << std::endl;
    std::cout << "Size of filter= (" << filter.rows() << " x " <<
    filter.cols()
    << ")" << std::endl; std::cout << "Size of Vandermonde matrix = (" <<
    m_V.rows() << " x " << m_V.cols() << ")" << std::endl;
    */

    m_V = V_temp * (V_modal * filter * V_modal_inv);

  } // If this is a Lagrange function
}

// ----------------------------------------------------------------------------

math::DenseDMat<Real> const &FEValues::Vandermonde() const
{
  return m_V;
}

// ----------------------------------------------------------------------------

math::DenseDMat<Real> const &FEValues::deriv_Vandermonde(const Uint variable) const
{
  return m_dV[variable];
}

// ----------------------------------------------------------------------------

const math::DenseDVec<Real> &FEValues::qw() const
{
  return m_pt_weights;
}

// ----------------------------------------------------------------------------

const math::DenseDMat<Real> &FEValues::qp() const
{
  return m_point_set;
}

// ----------------------------------------------------------------------------

void FEValues::print() const
{
  std::cout << " === FE VALUES === " << std::endl;
  std::cout << "SF type: " << m_sf_type.as_string() << std::endl;
  std::cout << "Reference coordinates: " << std::endl << m_ref_coord << std::endl;
  std::cout << "Vandermonde matrix:\n";
  std::cout << "size = [" << m_V.rows() << " x " << m_V.cols() << "] " << std::endl;
  std::cout << m_V << std::endl;

  std::cout << "Derivatives of Vandermonde matrix:\n";
  for (Uint d = 0; d < m_dV.size(); ++d)
  {
    std::cout << "dim: " << d << ", size = [" << m_dV[d].rows() << " x " << m_dV[d].cols() << "] "
              << std::endl;
    std::cout << m_dV[d] << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Reference coordinates in which the Vandermonde matrix was computed:" << std::endl
            << m_point_set << std::endl;

  std::cout << " ============== " << std::endl;
}

// ----------------------------------------------------------------------------

void FEValues::setup_type()
{
  mesh::StdRegion ref_elem;
  ref_elem.change_type(m_std_region_tag);

  m_topo_dim = ref_elem.get().topo_dim();

  m_dV.resize(m_topo_dim);

  mesh::sf::ShapeFunction shape_func;
  // shape_func.change_type(m_type.key<ElemShape>(),m_type.key<SFunc>(),Equidist,m_type.key<PolyOrder>(),Modal);
  shape_func.change_type(m_std_region_tag, m_sf_type);
  shape_func.get().ref_coords(m_ref_coord);
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit
