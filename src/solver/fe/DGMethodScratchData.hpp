#ifndef PDEKIT_DG_Method_Scratch_Data_hpp
#define PDEKIT_DG_Method_Scratch_Data_hpp

#include <array>

#include "interpolation/FEValues.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseDVec.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

// ----------------------------------------------------------------------------
//       CONTAINER FOR THE DATA IN REFERENCE SPACE
//       CONCERNING THE SOLUTION IN THE ELEMENT
// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM = PhysModel::DIM>
class DGMethodScratchData
{
  public:
  /// TYPEDEFS

  using physics = PhysModel;

  enum
  {
    NEQ = PhysModel::NEQ
  };

  /// Constructor
  DGMethodScratchData();

  /// Destructor
  ~DGMethodScratchData();

  /// This method resizes all member variables
  /// Depending on what's the finite element passed to it
  void resize_variables(interpolation::FEValues const &fe);

  /// This method resizes all member variables
  void resize_variables(mesh::StdRegionEntity const &sub_entity, interpolation::FEValues const &fe);

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

  /// Integral of flux in the element
  typename PhysModel::FluxV m_flux_integral_in_elem;

  /// Update coefficients for the vertices of the element
  /// Rows:    element nodes [nb_nodes]
  math::DenseDVec<Real> m_elem_wave_speed;

  /// Values of flux Jacobian, it's inverse and the diagonal
  /// matrix of Jacobian eigenvalues
  typename PhysModel::JM Rv;
  typename PhysModel::JM Lv;
  /// Matrix of positive eigenvalues
  typename PhysModel::JM Dvp;
  /// Matrix of negative eigenvalues
  typename PhysModel::JM Dvm;

  /// The matrix K+ in every quadrature node
  std::vector<typename PhysModel::JM> m_Kp;

  /// The matrix K- in every quadrature node
  std::vector<typename PhysModel::JM> m_Km;

  typename PhysModel::JM sum_Kp;
  typename PhysModel::JM inv_sum_Kp;

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
  math::DenseDMat<Real> m_V;

  /// Vandermonde matrices of solution shape function DERIVATIVES
  /// in reference space
  std::array<math::DenseDMat<Real>, PhysModel::DIM> m_dV;

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
};

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
DGMethodScratchData<PhysModel, TopoDIM>::DGMethodScratchData()
{
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
DGMethodScratchData<PhysModel, TopoDIM>::~DGMethodScratchData()
{
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
void DGMethodScratchData<PhysModel, TopoDIM>::resize_variables(interpolation::FEValues const &fe)

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

  m_elem_wave_speed.resize(m_nb_nodes);

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

  const math::DenseDMat<Real> &V = fe.Vandermonde();
  m_V.resize(V.rows(), V.cols());
  m_V = V;

  /// NOTE THAT m_dV IS AN ARRAY OF LENGTH PhysModel::DIM, BUT THE FEValues
  /// OBJECT CAN HAVE LOWER TOPOLOGICAL DIMENSION!
  // for (Uint d = 0; d < PhysModel::DIM; ++d)
  for (Uint d = 0; d < TopoDIM; ++d)
  {
    const math::DenseDMat<Real> &dV = fe.deriv_Vandermonde(d);
    m_dV[d].resize(dV.rows(), dV.cols());
    m_dV[d] = dV;
  }

  m_blending_coeff.resize(m_nb_nodes);
  m_blending_coeff.fill(0.0);

  m_art_visc.resize(fe.nb_qd_pts());
  m_art_visc.fill(0.0);

  m_use_external_theta = false;
}

// ----------------------------------------------------------------------------

template <typename PhysModel, Uint TopoDIM>
void DGMethodScratchData<PhysModel, TopoDIM>::resize_variables(
    mesh::StdRegionEntity const &sub_entity, interpolation::FEValues const &fe)

{
  resize_variables(fe);
  m_sub_idx.resize(sub_entity.nb_vert());
  for (Uint v = 0; v < sub_entity.nb_vert(); ++v)
  {
    m_sub_idx[v] = sub_entity.vertex(v);
  }
}

// ----------------------------------------------------------------------------

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
