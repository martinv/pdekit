#ifndef PDEKIT_RD_Method_Scratch_Data_hpp
#define PDEKIT_RD_Method_Scratch_Data_hpp

#include <array>

#include "common/ArrayView.hpp"
#include "interpolation/FEValues.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

// ----------------------------------------------------------------------------
//       CONTAINER FOR THE DATA IN REFERENCE SPACE
//       CONCERNING THE SOLUTION IN THE ELEMENT
// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM = PhysModel::DIM>
class RDMethodScratchData
{
  public:
  /// TYPEDEFS

  typedef PhysModel physics;

  enum
  {
    NEQ = PhysModel::NEQ
  };

  /// Constructor
  RDMethodScratchData();

  /// Destructor
  ~RDMethodScratchData();

  /// This method resizes all member variables
  /// Depending on what's the finite element passed to it
  void resize_variables(interpolation::FEValues const &fe);

  /// This method resizes all member variables
  void resize_variables(mesh::StdRegionEntity const &sub_entity, interpolation::FEValues const &fe);

  Uint nb_nodes()
  {
    return m_nb_nodes;
  }

  math::DenseConstVecView<Uint> sub_idx()
  {
    return math::DenseConstVecView<Uint>(m_sub_idx.data(), m_sub_idx.size());
  }

  math::DenseVecView<Real> elem_node_res()
  {
    return m_elem_node_res.block(0, m_elem_node_res.size());
  }

  math::DenseVecView<Real> elem_node_res_block(const Uint first_idx, const Uint size)
  {
    return m_elem_node_res.block(first_idx, size);
  }

  math::DenseVecView<Real> n_node_res()
  {
    return m_n_node_res.block(0, m_n_node_res.size());
  }

  math::DenseVecView<Real> n_node_res_block(const Uint first_idx, const Uint size)
  {
    return m_n_node_res.block(first_idx, size);
  }

  math::DenseVecView<Real> lda_node_res()
  {
    return m_lda_node_res.block(0, m_lda_node_res.size());
  }

  math::DenseVecView<Real> lda_node_res_block(const Uint first_idx, const Uint size)
  {
    return m_lda_node_res.block(first_idx, size);
  }

  math::DenseMatView<Real> elem_jacobian()
  {
    return m_elem_jacobian.block(0, 0, m_elem_jacobian.rows(), m_elem_jacobian.cols());
  }

  math::DenseMatView<Real> elem_jacobian_stab()
  {
    return m_elem_jacobian_stab.block(0, 0, m_elem_jacobian_stab.rows(),
                                      m_elem_jacobian_stab.cols());
  }

  inline typename PhysModel::JM &Rv()
  {
    return m_spectral_mats[0];
  }

  inline typename PhysModel::JM &Lv()
  {
    return m_spectral_mats[1];
  }

  inline typename PhysModel::JM &Dvp()
  {
    return m_spectral_mats[2];
  }

  inline typename PhysModel::JM &Dvm()
  {
    return m_spectral_mats[3];
  }

  inline typename PhysModel::JM &sum_Kp()
  {
    return m_spectral_mats[4];
  }

  inline typename PhysModel::JM &inv_sum_Kp()
  {
    return m_spectral_mats[5];
  }

  common::ArrayView<typename PhysModel::JM, _1D, Uint> Kp()
  {
    common::ArrayView<typename PhysModel::JM, _1D, Uint> result(m_Kp.data(), m_Kp.size());
    return result;
  }

  common::ArrayView<typename PhysModel::JM, _1D, Uint> Km()
  {
    common::ArrayView<typename PhysModel::JM, _1D, Uint> result(m_Km.data(), m_Km.size());
    return result;
  }

  private:
  /// Number of nodes for given element type
  Uint m_nb_nodes;

  /// Vector of local indices of facet nodes as seen
  /// from the master (volume) reference element
  /// Example: first facet of a P2 triangle will have
  /// sub-indexes {0,1,3}, second facet will have {1,2,4},
  /// third facet will have {2,0,5}
  std::vector<Uint> m_sub_idx;

  /// Nodal residuals
  math::DenseDVec<Real> m_elem_node_res;

  /// TODO: this should be scheme specific
  /// Residuals of the LDA scheme
  math::DenseDVec<Real> m_n_node_res;

  /// TODO: this should be scheme specific
  /// Residuals of the N scheme
  math::DenseDVec<Real> m_lda_node_res;

  /// Element Jacobi matrix
  math::DenseDMat<Real> m_elem_jacobian;

  /// Element Jacobi matrix - stabilization part
  math::DenseDMat<Real> m_elem_jacobian_stab;

  /*
  /// Values of flux Jacobian, it's inverse and the diagonal
  /// matrix of Jacobian eigenvalues
  typename PhysModel::JM Rv;
  typename PhysModel::JM Lv;
  /// Matrix of positive eigenvalues
  typename PhysModel::JM Dvp;
  /// Matrix of negative eigenvalues
  typename PhysModel::JM Dvm;
  /// Matrix to store sum(K+)
  typename PhysModel::JM sum_Kp;
  /// Matrix to store inv(sum(K+))
  typename PhysModel::JM inv_sum_Kp;
  */

  /// Stores the following matrices (in this order):
  /// R(right eigenvectors), L( = inv(R)), D+, D-, sum(K+), inv(sum(K+))
  std::array<typename PhysModel::JM, 6> m_spectral_mats;

  /// The matrix K+ in every quadrature node
  std::vector<typename PhysModel::JM> m_Kp;

  /// The matrix K- in every quadrature node
  std::vector<typename PhysModel::JM> m_Km;

  public:
  /// Integral of flux in the element
  typename PhysModel::FluxV m_flux_integral_in_elem;

  /// Update coefficients for the vertices of the element
  /// Rows:    element nodes [nb_nodes]
  math::DenseDVec<Real> m_elem_wave_speed;

  /// Values of solution gradient in one point
  typename PhysModel::SolGradM m_grad_u_at_point;

  /// Flux Jacobians in one point
  std::array<typename PhysModel::JM, PhysModel::DIM> m_dFdu;

  /// Residual vector at one point
  typename PhysModel::FluxV m_res_at_point;

  /// Jacobi block contribution at one point
  typename PhysModel::JM m_jacobi_at_point;

  /// Fluxes computed at one point
  typename PhysModel::FluxV m_norm_flux_at_point;

  /// Vandermonde matrix of solution shape function VALUES
  /// in reference space
  math::DenseDMat<Real> m_V_u;

  /// Vandermonde matrices of solution shape function DERIVATIVES
  /// in reference space
  std::array<math::DenseDMat<Real>, PhysModel::DIM> m_dV_u;

  /// Values of shape function (vector-valued) at one point
  typename PhysModel::JM m_sf_at_pt;

  /// Gradient of solution with respect to physical coordinates
  std::array<math::DenseConstMatView<Real>, PhysModel::DIM> m_grad_u;

  /// Gradient of fluxes with respect to physical coordinates
  std::array<math::DenseConstMatView<Real>, PhysModel::DIM> m_grad_F;

  /// Gradient of one shape function in one point (in reference space)
  math::DenseSVec<Real, PhysModel::DIM> m_grad_sf_at_pt_ref;

  /// Gradient of one shape function in one point (in physical space)
  math::DenseSVec<Real, PhysModel::DIM> m_grad_sf_at_pt_phys;

  /// Physical properties at one point
  typename PhysModel::Properties m_props;

  /// Maximum eigenvalues of the decomposed jacobian \grad(F) \cdot \nabla
  /// \varphi where F is the tensor of inviscid fluxes
  Real m_max_eigenvalue;

  /// Values of blending coefficient for nonlinear schemes
  math::DenseDVec<Real> m_blending_coeff;

  /// Artificial viscosity
  math::DenseDVec<Real> m_art_visc;

  /// Determine whether externally computed blending
  /// coefficient should be used
  bool m_use_external_theta;

  private:
};

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
RDMethodScratchData<PhysModel, TopoDIM>::RDMethodScratchData()
    : m_nb_nodes(0u), m_max_eigenvalue(0.0)
{
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
RDMethodScratchData<PhysModel, TopoDIM>::~RDMethodScratchData()
{
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
void RDMethodScratchData<PhysModel, TopoDIM>::resize_variables(interpolation::FEValues const &fe)

{
  m_nb_nodes = fe.nb_nodes();

  // This is a case where FEValues 'fe' represents Finite Element values
  // in a 'volume' element and no node sub-indexing is needed.
  // m_sub_idx is resized here just for safety and filled with 'identity'
  // values
  m_sub_idx.resize(m_nb_nodes);
  for (Uint v = 0; v < m_nb_nodes; ++v)
  {
    m_sub_idx[v] = v;
  }

  m_elem_node_res.resize(NEQ * m_nb_nodes);

  m_n_node_res.resize(NEQ * m_nb_nodes);

  m_lda_node_res.resize(NEQ * m_nb_nodes);

  m_elem_jacobian.resize(NEQ * m_nb_nodes, NEQ * m_nb_nodes);
  m_elem_jacobian_stab.resize(NEQ * m_nb_nodes, NEQ * m_nb_nodes);

  /*
  Rv.fill(0.0);
  Lv.fill(0.0);
  Dvp.fill(0.0);
  Dvm.fill(0.0);
  sum_Kp.fill(0.0);
  inv_sum_Kp.fill(0.0);
  */

  for (Uint i = 0; i < m_spectral_mats.size(); ++i)
  {
    m_spectral_mats[i].fill(0.0);
  }

  m_Kp.resize(m_nb_nodes);
  m_Km.resize(m_nb_nodes);

  for (Uint n = 0; n < m_nb_nodes; ++n)
  {
    m_Kp[n].fill(0.0);
    m_Km[n].fill(0.0);
  }

  const math::DenseDMat<Real> &V = fe.Vandermonde();
  m_V_u.resize(V.rows(), V.cols());
  m_V_u = V;

  m_elem_wave_speed.resize(m_nb_nodes);

  /// NOTE THAT m_dV IS AN ARRAY OF LENGTH PhysModel::DIM, BUT THE FEValues
  /// OBJECT CAN HAVE LOWER TOPOLOGICAL DIMENSION!
  // for (Uint d = 0; d < PhysModel::DIM; ++d)
  for (Uint d = 0; d < TopoDIM; ++d)
  {
    const math::DenseDMat<Real> &dV = fe.deriv_Vandermonde(d);
    m_dV_u[d].resize(dV.rows(), dV.cols());
    m_dV_u[d] = dV;
  }

  m_blending_coeff.resize(m_nb_nodes);
  m_blending_coeff.fill(0.0);

  m_art_visc.resize(fe.nb_qd_pts());
  m_art_visc.fill(0.0);

  m_use_external_theta = false;
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
void RDMethodScratchData<PhysModel, TopoDIM>::resize_variables(
    mesh::StdRegionEntity const &sub_entity, interpolation::FEValues const &fe)

{
  /*
  m_nb_nodes = sub_entity.nb_vert();

  m_sub_idx.resize(sub_entity.nb_vert());
  for (Uint v = 0; v < sub_entity.nb_vert(); ++v)
  {
    m_sub_idx[v] = sub_entity.vertex(v);
  }

  m_elem_node_res.resize(NEQ*m_nb_nodes);

  m_n_node_res.resize(NEQ*m_nb_nodes);

  m_lda_node_res.resize(NEQ*m_nb_nodes);

  m_elem_update_coeff.resize(m_nb_nodes);

  Rv.fill(0.0);
  Lv.fill(0.0);
  Dvp.fill(0.0);
  Dvm.fill(0.0);

  m_Kp.resize(m_nb_nodes);
  m_Km.resize(m_nb_nodes);

  for (Uint n = 0; n < m_nb_nodes; ++n)
  {
    m_Kp[n].fill(0.0);
    m_Km[n].fill(0.0);
  }

  sum_Kp.fill(0.0);
  inv_sum_Kp.fill(0.0);

  const math::DynamicMatrix<Real> &V = fe.Vandermonde();
  m_V.resize(V.rows(), V.cols());
  m_V = V;

  /// NOTE THAT m_dV IS AN ARRAY OF LENGTH PhysModel::DIM, BUT THE FEValues
  /// OBJECT CAN HAVE LOWER TOPOLOGICAL DIMENSION!
  // for (Uint d = 0; d < PhysModel::DIM; ++d)
  for (Uint d = 0; d < TopoDIM; ++d)
  {
    const math::DynamicMatrix<Real> &dV = fe.deriv_Vandermonde(d);
    m_dV[d].resize(dV.rows(), dV.cols());
    m_dV[d] = dV;
  }

  m_blending_coeff.resize(m_nb_nodes);

  m_max_edge_length = 0.0;

  m_art_visc = 0.0;

  */
  resize_variables(fe);
  m_sub_idx.resize(sub_entity.nb_vert());
  for (Uint v = 0; v < sub_entity.nb_vert(); ++v)
  {
    m_sub_idx[v] = sub_entity.vertex(v);
  }
}

// ----------------------------------------------------------------------------
// Special variant for the B scheme
// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM = PhysModel::DIM>
class PGRDBMethodData : public RDMethodScratchData<PhysModel, TopoDIM>
{
  public:
  /// Constructor
  PGRDBMethodData();

  /// Destructor
  ~PGRDBMethodData();

  void resize_variables(interpolation::FEValues const &fe);

  void resize_variables(mesh::StdRegionEntity const &sub_entity, interpolation::FEValues const &fe);

  typename PhysModel::SolV theta_cellwise;
  Real abs_phi;
  Real abs_phiN;

  private:
};

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
PGRDBMethodData<PhysModel, TopoDIM>::PGRDBMethodData() : RDMethodScratchData<PhysModel, TopoDIM>()
{
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
PGRDBMethodData<PhysModel, TopoDIM>::~PGRDBMethodData()
{
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
void PGRDBMethodData<PhysModel, TopoDIM>::resize_variables(interpolation::FEValues const &fe)
{
  RDMethodScratchData<PhysModel, TopoDIM>::resize_variables(fe);
  theta_cellwise.fill(0.0);
  abs_phi  = 0.0;
  abs_phiN = 0.0;
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
void PGRDBMethodData<PhysModel, TopoDIM>::resize_variables(mesh::StdRegionEntity const &sub_entity,
                                                           interpolation::FEValues const &fe)
{
  RDMethodScratchData<PhysModel, TopoDIM>::resize_variables(sub_entity, fe);
  theta_cellwise.fill(0.0);
  abs_phi  = 0.0;
  abs_phiN = 0.0;
}

// ----------------------------------------------------------------------------
// Special variant for the Lax-Friedrichs scheme
// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM = PhysModel::DIM>
class PGRDLFMethodData : public RDMethodScratchData<PhysModel, TopoDIM>
{
  public:
  /// Constructor
  PGRDLFMethodData();

  /// Destructor
  ~PGRDLFMethodData();

  void resize_variables(interpolation::FEValues const &fe);

  void resize_variables(mesh::StdRegionEntity const &sub_entity, interpolation::FEValues const &fe);

  /// Total cell residual in one element - needed to compute limited residuals
  /// in LLF scheme
  typename PhysModel::FluxV m_total_cell_res;

  /// Average state
  typename PhysModel::FluxV m_u_avg;

  /// Sum of x+
  typename PhysModel::FluxV m_sum_x_p;

  /// Vector of Jacobian flux eigenvalues
  typename PhysModel::FluxV m_eigenvalues;

  /// The i-th entry of this vector is 1/DIM * int{ \nabla \varphi_i }
  std::vector<typename PhysModel::CoordV> m_LF_proj_vec;

  /// The factor k in LF: k_i = \int_K { max_eigenvalue(dF/dU \cdot \nabla
  /// \varphi_i) }
  math::DenseDVec<Real> m_k_LF;

  /// Values of upwind stabilization coefficient
  math::DenseDVec<Real> m_LF_upwind_stab;

  /// Jacobians projected on (some) direction
  typename PhysModel::JM m_J_proj;

  /// Coordinates of barycenter
  typename PhysModel::CoordV m_barycenter;

  /// Cell volume
  Real m_cell_volume;

  /// Stabilisation factor for Lax-Friedrichs scheme
  Real m_LF_alpha;

  private:
};

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
PGRDLFMethodData<PhysModel, TopoDIM>::PGRDLFMethodData()
    : RDMethodScratchData<PhysModel, TopoDIM>(), m_cell_volume(0.0), m_LF_alpha(0.0)
{
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
PGRDLFMethodData<PhysModel, TopoDIM>::~PGRDLFMethodData()
{
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
void PGRDLFMethodData<PhysModel, TopoDIM>::resize_variables(interpolation::FEValues const &fe)
{
  RDMethodScratchData<PhysModel, TopoDIM>::resize_variables(fe);
  m_total_cell_res.fill(0.0);
  m_sum_x_p.fill(0.0);
  m_LF_proj_vec.resize(RDMethodScratchData<PhysModel, TopoDIM>::m_nb_nodes);
  m_k_LF.resize(RDMethodScratchData<PhysModel, TopoDIM>::m_nb_nodes);
  m_LF_upwind_stab.resize(fe.nb_qd_pts());
  m_LF_alpha = 0.0;
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
void PGRDLFMethodData<PhysModel, TopoDIM>::resize_variables(mesh::StdRegionEntity const &sub_entity,
                                                            interpolation::FEValues const &fe)
{
  RDMethodScratchData<PhysModel, TopoDIM>::resize_variables(sub_entity, fe);
  m_total_cell_res.fill(0.0);
  m_sum_x_p.fill(0.0);
  m_LF_proj_vec.resize(RDMethodScratchData<PhysModel, TopoDIM>::m_nb_nodes);
  m_k_LF.resize(RDMethodScratchData<PhysModel, TopoDIM>::m_nb_nodes);
  m_LF_upwind_stab.resize(fe.nb_qd_pts());
  m_LF_alpha = 0.0;
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
