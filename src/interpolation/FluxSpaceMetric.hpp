#ifndef PDEKIT_Interpolation_Flux_Space_Metric_hpp
#define PDEKIT_Interpolation_Flux_Space_Metric_hpp

#include <set>

#include "common/Meta.hpp"
#include "interpolation/CellFluxMetric.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "interpolation/MetricFlags.hpp"
#include "interpolation/SolutionSpaceMetric.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM = MeshConfig::TDIM>
class FluxSpaceMetric
{
  private:
  enum
  {
    evaluate_boundary_normals = MeshConfig::GDIM > DIM ? 1 : 0
  };

  struct MetricData
  {
    enum
    {
      GDIM = MeshConfig::GDIM
    };
    enum
    {
      TDIM = DIM
    };

    /// Type of this element
    mesh::PointSetTagExt std_region_type;

    /// Solution shape function values in reference space
    math::DenseDMat<Real> Vs;

    /// Solution shape function derivatives in reference space
    std::array<math::DenseDMat<Real>, DIM> dVs;

    /// Reference element values for flux shape functions
    math::DenseDMat<Real> Vf;

    /// Flux shape function derivatives in reference space
    std::array<math::DenseDMat<Real>, DIM> dVf;

    /// Reference element values for conversion from geometry
    /// to flux nodes - this matrix is used to interpolate
    /// geometry coordinates to flux DOF
    math::DenseDMat<Real> Vgf;

    /// Derivatives of reference element values for conversion from geometry
    /// to flux nodes - this matrix is used to compute jacobians of
    /// transformation reference space -> physical space in the nodes
    /// of the flux element
    std::array<math::DenseDMat<Real>, DIM> dVgf;

    /// Reference element values for conversion from solution
    /// to flux nodes - this matrix is used to interpolate
    /// solution to flux DOF
    math::DenseDMat<Real> Vsf;

    /// Derivatives of shape functions for conversion from solution
    /// to flux nodes in reference space
    std::array<math::DenseDMat<Real>, DIM> dVsf_ref;

    /// Derivatives of shape functions for conversion from solution
    /// to flux nodes in physical space
    std::array<math::DenseDMat<Real>, DIM> dVsf_phys;

    /// Coordinates of flux nodes (of ONE cell)
    math::DenseDMat<Real> Gf;

    /// Derivatives of coordinate transformation
    std::array<math::DenseDMat<Real>, DIM> dGf;

    /// Solution values in flux nodes (of ONE cell)
    math::DenseDMat<Real> Sf;

    /// Derivativs of solution values in flux nodes (of ONE cell)
    std::array<math::DenseDMat<Real>, DIM> dSf;

    /// How many elements have been filled
    Uint nb_elem_filled;

    /// Number of dof per elem
    Uint nb_dof_per_elem;

    /// Number of flux dof per elem
    Uint nb_flux_dof_per_elem;

    /// Number of points in which we evaluate (quadrature points)
    Uint nb_qd_pts_per_elem;

    /// Matrix of flux values in nodes of a BLOCK of elements
    /// This is an array of matrices: the first matrix is for flux in the
    /// x-direction, the second matrix for flux in the y-direction etc.
    /// Each matrix should be of size [ (nb. flux DOF) x (nb_blocks * NEQ) ]
    std::array<math::DenseDMat<Real>, Physics::DIM> m_flux_in_flux_nodes;

    /// Vector of matrices of flux values in nodes of a BLOCK of elements
    /// m_flux_deriv_in_flux_nodes[X0] is a matrix of size
    /// [ (nb. flux DOF) x (nb_blocks * NEQ) ]
    /// and it holds all values of dF0/dX0 for a block of cells
    std::array<math::DenseDMat<Real>, Physics::DIM> m_flux_deriv_in_flux_nodes;

    /// Matrix of flux values in quadrature points of a BLOCK of elements
    /// m_flux_in_qd_pts[X0] is a matrix of size
    /// [ (nb. quadrature points) x (nb_blocks * NEQ) ]
    /// and it holds all values of F0(u) for a block of cells
    /// Similarly m_flux_in_qd_pts[X1] has values of F1(u) etc.
    std::array<math::DenseDMat<Real>, Physics::DIM> m_flux_in_qd_pts;

    /// Vector of matrices of flux values in quadrature points of a BLOCK of
    /// elements m_flux_deriv_in_qd_pts[X0] is a matrix of size [ (nb.
    /// quadrature points) x (nb_blocks * NEQ) ] and it holds all values of
    /// dF0(u)/dX0 for a block of cells Similarly m_flux_deriv_in_qd_pts[X1]
    /// has values of dF1(u)/dX1 etc.
    std::array<math::DenseDMat<Real>, Physics::DIM> m_flux_deriv_in_qd_pts;

    /// Temporary extra storage for flux derivatives
    std::array<math::DenseDMat<Real>, Physics::DIM> m_tmp_flux_deriv;
  };

  public:
  // typedef CellFluxMetric<MetricData, Physics> cellwise_metric;

  using cellwise_metric =
      typename common::SelectType<evaluate_boundary_normals, FacetFluxMetric<MetricData, Physics>,
                                  CellFluxMetric<MetricData, Physics>>::type;

  /// Default constructor
  FluxSpaceMetric();

  /// Default destructor
  ~FluxSpaceMetric();

  /// Return the number of elements that fit in buffer
  Uint max_nb_blocks_in_buffer() const;

  /// Return the number of fields per dof
  Uint nb_fields() const;

  /// Set the size of buffer
  template <typename FeValsIterator>
  void allocate_buffer(const SFunc geo_sf_type, const Uint geo_sf_order,
                       const FeValsIterator fe_vals_begin, const FeValsIterator fe_vals_end,
                       Uint const nb_blocks);

  /// Empty the data in the buffer
  void empty_buffer();

  /// Clear all internal data
  void clear();

  /// Return nb. values that have been pushed to buffer so far
  Uint nb_values_in_buffer() const;

  /// Collocation - interpolate the values and their derivatives
  void evaluate(GeometryCache<MeshConfig::GDIM> const &gcache,
                GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &geo_met,
                SolutionCache const &scache, SolutionSpaceMetric<MeshConfig, DIM> const &sol_met,
                const RebuildMetricIndex rebuild_index);

  /// Get one block of values
  cellwise_metric const cellwise_values(const Uint idx) const;

  private:
  /// TYPES

  /// This class actually knows how to evaluate the derivatives, jacobians ...
  /// It is a helper class visible only to FunctionSpaceMetric. Client code
  /// of FunctionSpaceMetric will never need to use it directly.

  template <Uint GEODIM, Uint TDIM, Uint Dummy>
  class MetricComputer
  {
  };

  /// DATA

  /// Maximum number of cell entries the buffer can take
  Uint m_nb_blocks;

  /// A map of buffers: for each element type we have one buffer
  common::DataMap<mesh::PointSetTagExt, MetricData> m_metric_data_map;

  /// Linear index of metric data, consists of a vector of pairs
  /// [PointSetTag,Uint]
  /// Value m_mdata_index[i] is a pair [FIRST,SECOND] and corresponds to
  /// metric data which can be found in
  /// m_metric_data_map.std_region_data(FIRST), cell block on position
  /// [SECOND] in this particular MetricData
  std::vector<std::pair<common::PtrHandle<MetricData>, Uint>> m_mdata_index;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
FluxSpaceMetric<MeshConfig, Physics, DIM>::FluxSpaceMetric()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
FluxSpaceMetric<MeshConfig, Physics, DIM>::~FluxSpaceMetric()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
Uint FluxSpaceMetric<MeshConfig, Physics, DIM>::max_nb_blocks_in_buffer() const
{
  return m_nb_blocks;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
Uint FluxSpaceMetric<MeshConfig, Physics, DIM>::nb_fields() const
{
  return Physics::NEQ;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
template <typename FeValsIterator>
void FluxSpaceMetric<MeshConfig, Physics, DIM>::allocate_buffer(const SFunc geo_sf_type,
                                                                const Uint geo_sf_order,
                                                                const FeValsIterator fe_vals_begin,
                                                                const FeValsIterator fe_vals_end,
                                                                Uint const nb_blocks)
{
  /// NOTE: for discontinuous method, there is no need to allocate
  /// m_coord_buffer and copy data into it, but instead m_coord_buffer could
  /// be just a MatrixBlock<Real>. To make this class work with both
  /// continuous and discontinuous fields, each metric_data should have an
  /// array of 'raw' values and m_coord_buffer should become a
  /// MatrixBlock<Real>. Then the procedure would be as follows: I)  For
  /// continous method, 'raw' values would allocate sufficient memory
  ///     and let m_coord_buffer be a proxy to this block
  /// II) For discontinous method, no 'raw' values would be allocated, but
  /// only
  ///     m_coord_buffer would serve as proxy to wherever the data is stored.
  ///     Since the data is already discontinous, it does not need to be
  ///     copied, but the work can be directly done on it by means of the
  ///     proxy.

  m_metric_data_map.clear();

  m_nb_blocks = nb_blocks;

  FeValsIterator it = fe_vals_begin;
  Uint n_elem_types = 0;

  for (; it != fe_vals_end; ++it)
  {
    n_elem_types++;
  }

  m_mdata_index.reserve(n_elem_types * nb_blocks);
  m_mdata_index.resize(0);

  mesh::StdRegion flux_ref_element;
  mesh::sf::ShapeFunction shape_function;

  std::vector<math::DenseDMat<Real>> tmp_derivatives(DIM);

  it = fe_vals_begin;

  for (; it != fe_vals_end; ++it)
  {
    common::PtrHandle<MetricData> met_data = m_metric_data_map.create(it.key_value());

    (*met_data).std_region_type = it.key_value();

    // Copy solution shape function values
    math::DenseDMat<Real> const &V = (*it.data_ptr()).Vandermonde();
    (*met_data).Vs.resize(V.rows(), V.cols());
    (*met_data).Vs = V;

    // ------------------------------------------------------------
    // DERIVATIVES are used and allocated only in case GDIM == TDIM
    // ------------------------------------------------------------

    if (MeshConfig::GDIM == DIM)
    {
      for (Uint d = 0; d < DIM; ++d)
      {
        math::DenseDMat<Real> const &dV = (*it.data_ptr()).deriv_Vandermonde(d);

        (*met_data).dVs[d].resize(dV.rows(), dV.cols());
        (*met_data).dVs[d] = dV;
      }
    }

    // Shape function for solution
    const mesh::sf::SFTag solution_sf_tag = (*it.data_ptr()).sf_type();

    // Shape function for geometry
    mesh::sf::SFTag geometry_sf_tag = solution_sf_tag;
    geometry_sf_tag.set_shape_function(geo_sf_type);
    geometry_sf_tag.set_poly_order(geo_sf_order);

    // Flux shape function will be of the same type as solution shape
    // function, but one polynomial order higher
    mesh::sf::SFTag flux_sf_tag = solution_sf_tag;
    flux_sf_tag.set_poly_order(solution_sf_tag.poly_order() + P1);
    // NOTE! This will make sure that Lagrange basis is used to compute flux
    // values and flux derivatives in quadrature points even if the solution
    // interpolant is a modal basis
    flux_sf_tag.set_shape_function(SFunc::Lagrange);

    const mesh::PointSetTag tmp_tag = (*met_data).std_region_type.std_region_tag();
    const mesh::PointSetTag flux_std_reg_tag(tmp_tag.elem_shape(), tmp_tag.poly_order() + P1,
                                             tmp_tag.ref_topology());

    shape_function.change_type(flux_std_reg_tag, flux_sf_tag);

    // Vf should have dimension [ (nb. quad. pts) x (nb. DOF per flux elem)
    // ]
    shape_function.get().compute_ref_values((*it.data_ptr()).qp(), (*met_data).Vf);

    // ------------------------------------------------------------
    // DERIVATIVES are used and allocated only in case GDIM == TDIM
    // ------------------------------------------------------------

    if (MeshConfig::GDIM == DIM)
    {
      shape_function.get().compute_ref_derivatives((*it.data_ptr()).qp(), tmp_derivatives);

      for (Uint d = 0; d < DIM; ++d)
      {
        math::DenseDMat<Real> const &dV = tmp_derivatives[d];
        (*met_data).dVf[d].resize(dV.rows(), dV.cols());

        // dVf[d] should have dimension [ (nb. quad. pts) x (nb. DOF per
        // flux elem) ]
        (*met_data).dVf[d] = dV;
      }
    }

    flux_ref_element.change_type(flux_std_reg_tag);

    shape_function.change_type(flux_std_reg_tag, geometry_sf_tag);
    shape_function.get().compute_ref_values(flux_ref_element.get().coordinates(), (*met_data).Vgf);

    // ------------------------------------------------------------
    // DERIVATIVES are used and allocated only in case GDIM == TDIM
    // ------------------------------------------------------------

    if (MeshConfig::GDIM == DIM)
    {
      shape_function.get().compute_ref_derivatives(flux_ref_element.get().coordinates(),
                                                   tmp_derivatives);

      for (Uint d = 0; d < DIM; ++d)
      {
        math::DenseDMat<Real> const &dV = tmp_derivatives[d];
        (*met_data).dVgf[d].resize(dV.rows(), dV.cols());
        (*met_data).dVgf[d] = dV;
      }
    }

    shape_function.change_type(flux_std_reg_tag, solution_sf_tag);
    shape_function.get().compute_ref_values(flux_ref_element.get().coordinates(), (*met_data).Vsf);

    // ------------------------------------------------------------
    // DERIVATIVES are used and allocated only in case GDIM == TDIM
    // ------------------------------------------------------------

    if (MeshConfig::GDIM == DIM)
    {
      shape_function.get().compute_ref_derivatives(flux_ref_element.get().coordinates(),
                                                   tmp_derivatives);

      for (Uint d = 0; d < DIM; ++d)
      {
        math::DenseDMat<Real> const &dV = tmp_derivatives[d];
        (*met_data).dVsf_ref[d].resize(dV.rows(), dV.cols());

        // dVsf[d] should have dimension [ (nb. DOF per flux elem) x
        // (nb. DOF per solution elem) ]
        (*met_data).dVsf_ref[d] = dV;
      }

      for (Uint d = 0; d < DIM; ++d)
      {
        (*met_data).dVsf_phys[d].resize((*met_data).dVsf_ref[d].rows(),
                                        (*met_data).dVsf_ref[d].cols());
      }
    }

    const Uint nb_nodes_in_flux_elem = flux_ref_element.get().nb_nodes();

    // Coordinates of flux nodes in one cell
    (*met_data).Gf.resize(nb_nodes_in_flux_elem, MeshConfig::GDIM);
    // Solution values in flux nodes of one cell
    (*met_data).Sf.resize(nb_nodes_in_flux_elem, Physics::NEQ);

    // ------------------------------------------------------------
    // DERIVATIVES are used and allocated only in case GDIM == TDIM
    // ------------------------------------------------------------

    if (MeshConfig::GDIM == DIM)
    {
      for (Uint d = 0; d < DIM; ++d)
      {
        (*met_data).dGf[d].resize(nb_nodes_in_flux_elem, MeshConfig::GDIM);
        (*met_data).dSf[d].resize(nb_nodes_in_flux_elem, Physics::NEQ);
      }
    }

    // Vandermonde matrix:
    // number of rows corresponds to the number of Gauss points
    // number of columns corresponds to the number of shape functions

    (*met_data).nb_dof_per_elem      = (*it.data_ptr()).Vandermonde().cols();
    (*met_data).nb_flux_dof_per_elem = flux_ref_element.get().nb_nodes();
    (*met_data).nb_qd_pts_per_elem   = (*it.data_ptr()).Vandermonde().rows();

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      (*met_data).m_flux_in_flux_nodes[d].resize((*met_data).nb_flux_dof_per_elem,
                                                 nb_blocks * Physics::NEQ);
      (*met_data).m_flux_in_qd_pts[d].resize((*met_data).nb_qd_pts_per_elem,
                                             nb_blocks * Physics::NEQ);

      // DERIVATIVES are used and allocated only in case GDIM == TDIM
      if (MeshConfig::GDIM == DIM)
      {
        (*met_data).m_flux_deriv_in_flux_nodes[d].resize((*met_data).nb_flux_dof_per_elem,
                                                         nb_blocks * Physics::NEQ);

        (*met_data).m_flux_deriv_in_qd_pts[d].resize((*met_data).nb_qd_pts_per_elem,
                                                     nb_blocks * Physics::NEQ);

        (*met_data).m_tmp_flux_deriv[d].resize((*met_data).nb_qd_pts_per_elem,
                                               nb_blocks * Physics::NEQ);
      }
    }

    (*met_data).nb_elem_filled = 0;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
void FluxSpaceMetric<MeshConfig, Physics, DIM>::empty_buffer()
{
  for (typename common::DataMap<mesh::PointSetTagExt, MetricData>::iterator it =
           m_metric_data_map.begin();
       it != m_metric_data_map.end(); ++it)
  {
    (*it.data_ptr()).nb_elem_filled = 0;
  }
  m_mdata_index.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
void FluxSpaceMetric<MeshConfig, Physics, DIM>::clear()
{
  m_nb_blocks = 0;

  m_metric_data_map.clear();

  m_mdata_index.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
Uint FluxSpaceMetric<MeshConfig, Physics, DIM>::nb_values_in_buffer() const
{
  return m_mdata_index.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
void FluxSpaceMetric<MeshConfig, Physics, DIM>::evaluate(
    GeometryCache<MeshConfig::GDIM> const &gcache,
    GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &geo_met,
    SolutionCache const &scache, SolutionSpaceMetric<MeshConfig, DIM> const &sol_met,
    const RebuildMetricIndex rebuild_index)
{
  // 1) Set the number of filled elements to zero for all cell blocks
  for (typename common::DataMap<mesh::PointSetTagExt, MetricData>::iterator it =
           m_metric_data_map.begin();
       it != m_metric_data_map.end(); ++it)
  {
    (*it.data_ptr()).nb_elem_filled = 0;
  }

  // 2) Rebuild index if necessary
  if (rebuild_index)
  {
    // Rebuild the data index for fast access to metric data
    m_mdata_index.resize(0);

    for (Uint i = 0; i < scache.nb_values_in_buffer(); ++i)
    {
      common::PtrHandle<MetricData> md =
          m_metric_data_map.std_region_data(scache.std_region_type(i));
      (*md).nb_elem_filled++;
      m_mdata_index.push_back(
          std::pair<common::PtrHandle<MetricData>, Uint>(md, scache.position_in_block(i)));
    }
  }

  // 3) Compute flux derivatives in quadrature points
  MetricComputer<MeshConfig::GDIM, DIM, 0>::evaluate_metric_terms(gcache, geo_met, scache, sol_met,
                                                                  *this);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
typename FluxSpaceMetric<MeshConfig, Physics, DIM>::cellwise_metric const FluxSpaceMetric<
    MeshConfig, Physics, DIM>::cellwise_values(const Uint idx) const
{
  cellwise_metric cm(m_mdata_index[idx].first, m_mdata_index[idx].second);
  return cm;
}

// ----------------------------------------------------------------------------
// Specializations of FluxSpaceMetric::MetricComputer for some combinations
// of dimensions  - here 2D volume metric (GDIM = 2, TDIM = 2)
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
template <Uint Dummy>
class FluxSpaceMetric<MeshConfig, Physics, DIM>::MetricComputer<_2D, _2D, Dummy>
{
  public:
  static void evaluate_metric_terms(
      GeometryCache<MeshConfig::GDIM> const &gcache,
      GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &geo_met,
      SolutionCache const &scache, SolutionSpaceMetric<MeshConfig, DIM> const &sol_met,
      FluxSpaceMetric<MeshConfig, Physics, DIM> &flux_met)
  {
    using geo_cell_metric_type =
        typename GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM>::cellwise_metric;

    using flux_data_type = typename FluxSpaceMetric<MeshConfig, Physics, DIM>::MetricData;

    math::DenseSMat<Real, _2D, _2D> J;
    // Derivatives of shape function in one point
    math::DenseSVec<Real, _2D> d_phi_ref;

    typename Physics::Properties phys_properties;
    math::DenseSMat<Real, Physics::NEQ, Physics::DIM> tmp_flux;
    math::DenseSMat<Real, Physics::NEQ, Physics::DIM> solution_grad;

    // ============================================================
    // 1) LOOP OVER ALL CELLS AND COMPUTE FLUX VALUES IN FLUX NODES
    // ============================================================

    for (Uint c = 0; c < flux_met.m_mdata_index.size(); ++c)
    {
      const math::DenseConstMatView<Real> geo_dof_coords = gcache.cell_values(c);
      const math::DenseConstMatView<Real> sol_dof_values = scache.cell_values(c);

      flux_data_type &flux_data = (*flux_met.m_mdata_index[c].first);

      // Compute physical coordinates in the nodes of flux element
      // [ (nb. flux DOF) x GDIM] = [ (nb. flux DOF) x (nb. sol. DOF)] *
      // [(nb. sol. DOF) x GDIM ]
      flux_data.Gf = flux_data.Vgf * geo_dof_coords;

      // Compute derivatives of transformation physical - reference space
      // in the nodes of flux element
      for (Uint d = 0; d < DIM; ++d)
      {
        flux_data.dGf[d] = flux_data.dVgf[d] * geo_dof_coords;
      }

      // Compute solution values in the nodes of flux element
      // [ (nb. flux DOF) x NEQ]  =  [ (nb. flux DOF) x (nb. sol. DOF)] *
      // [(nb. sol. DOF) x NEQ ]
      flux_data.Sf = flux_data.Vsf * sol_dof_values;

      /// DEBUG: THIS SEEMS TO BE CORRECT
      /*
      std::cout << "\n\n-----------------------------------------" <<
      std::endl; std::cout << "Solution in flux nodes:" << std::endl;
      std::cout << flux_data.Sf << std::endl;
      */

      for (Uint nf = 0; nf < flux_data.nb_flux_dof_per_elem; ++nf)
      {
        J(XI0, X0) = flux_data.dGf[XI0](nf, X0); // dX/dXi
        J(XI0, X1) = flux_data.dGf[XI0](nf, X1); // dY/dXi

        J(XI1, X0) = flux_data.dGf[XI1](nf, X0); // dX/dEta
        J(XI1, X1) = flux_data.dGf[XI1](nf, X1); // dY/dEta

        const Real inv_det = 1.0 / (J(XI0, X0) * J(XI1, X1) - J(XI0, X1) * J(XI1, X0));

        for (Uint ns = 0; ns < flux_data.nb_dof_per_elem; ++ns)
        {
          d_phi_ref[XI0] = flux_data.dVsf_ref[XI0](nf, ns);
          d_phi_ref[XI1] = flux_data.dVsf_ref[XI1](nf, ns);

          flux_data.dVsf_phys[XI0](nf, ns) =
              inv_det * (J(XI1, X1) * d_phi_ref[X0] - J(XI0, X1) * d_phi_ref[X1]);

          flux_data.dVsf_phys[XI1](nf, ns) =
              inv_det * (-J(XI1, X0) * d_phi_ref[X0] + J(XI0, X0) * d_phi_ref[X1]);
        }
      }

      // Compute derivatives of solution (in the physical space)
      // in the nodes of flux element and store them in dSf
      for (Uint d = 0; d < DIM; ++d)
      {
        // [(nb. flux dof) x  NEQ] = [(nb. flux dof) x (nb. sol. dof)] *
        // [(nb. sol. dof) x NEQ]
        flux_data.dSf[d] = flux_data.dVsf_phys[d] * sol_dof_values;
      }

      // Offset determines where in m_flux_in_flux_nodes are stored values
      // for cell with 'global' buffer index c, which in given block of
      // cells is on position m_mdata_index[c].second
      const Uint offset = flux_met.m_mdata_index[c].second * Physics::NEQ;

      // Compute fluxes in flux dofs of one cell
      for (Uint n = 0; n < flux_data.nb_flux_dof_per_elem; ++n)
      {
        solution_grad.insert_col(X0, flux_data.dSf[X0].const_row(n));
        solution_grad.insert_col(X1, flux_data.dSf[X1].const_row(n));

        Physics::compute_properties(flux_data.Gf.const_row_transp(n),
                                    flux_data.Sf.const_row_transp(n), solution_grad,
                                    phys_properties);

        Physics::flux(phys_properties, tmp_flux);

        // m_flux_in_flux_nodes is of size [nb. qd_pts x (nb_blocks NEQ)]

        for (Uint eq = 0; eq < Physics::NEQ; ++eq)
        {
          flux_data.m_flux_in_flux_nodes[X0](n, offset + eq) = tmp_flux(eq, X0);
          flux_data.m_flux_in_flux_nodes[X1](n, offset + eq) = tmp_flux(eq, X1);
        }

      } // Loop over nodes of one flux element

    } // Loop over all cells in buffer

    // ============================================================
    // 2) Compute fluxes in quadrature points
    // ============================================================

    for (typename common::DataMap<mesh::PointSetTagExt, flux_data_type>::iterator it =
             flux_met.m_metric_data_map.begin();
         it != flux_met.m_metric_data_map.end(); ++it)
    {
      flux_data_type &fdata = (*it.data_ptr());

      if (fdata.nb_elem_filled > 0)
      {
        // [(nb. qd. pts) x NEQ] = [ (nb. qd. pts) x (nb. flux DOF) ]  *
        // [(nb. flux DOF) x NEQ]
        fdata.m_flux_in_qd_pts[X0] = fdata.Vf * fdata.m_flux_in_flux_nodes[X0];
        fdata.m_flux_in_qd_pts[X1] = fdata.Vf * fdata.m_flux_in_flux_nodes[X1];
      }
    }

    // ============================================================
    // 3) Compute derivatives of fluxes.
    //    For each coordinate component icomp, compute the derivative
    //    of the i-th flux component with respect to i-th coordinate
    //    variable.
    //    In other words, compute dF0/dx, dF1/dy, dF2/dz (in 3D)
    // ============================================================

    math::DenseSVec<Real, _2D> d_flux_ref, d_flux_phys;

    for (Uint icomp = 0; icomp < Physics::DIM; ++icomp)
    {
      for (typename common::DataMap<mesh::PointSetTagExt, flux_data_type>::iterator it =
               flux_met.m_metric_data_map.begin();
           it != flux_met.m_metric_data_map.end(); ++it)
      {
        flux_data_type &fdata = (*it.data_ptr());

        if (fdata.nb_elem_filled > 0)
        {
          // Flux derivatives
          fdata.m_tmp_flux_deriv[X0] = fdata.dVf[X0] * fdata.m_flux_in_flux_nodes[icomp];
          fdata.m_tmp_flux_deriv[X1] = fdata.dVf[X1] * fdata.m_flux_in_flux_nodes[icomp];
        }
      }

      // The jacobians can be provided by geometry metric
      for (Uint c = 0; c < flux_met.m_mdata_index.size(); ++c)
      {
        const geo_cell_metric_type cell_geo_values = geo_met.cellwise_values(c);

        flux_data_type &flux_data = (*flux_met.m_mdata_index[c].first);

        const Uint offset = flux_met.m_mdata_index[c].second * Physics::NEQ;

        for (Uint q = 0; q < flux_data.nb_qd_pts_per_elem; ++q)
        {
          const math::DenseConstMatView<Real> J_inv = cell_geo_values.inv_jacobi(q);

          for (Uint eq = 0; eq < Physics::NEQ; ++eq)
          {
            d_flux_ref[XI0] = flux_data.m_tmp_flux_deriv[XI0](q, offset + eq);
            d_flux_ref[XI1] = flux_data.m_tmp_flux_deriv[XI1](q, offset + eq);

            // d_flux_phys[X0] = J_inv(0,0) * d_flux_ref[X0] +
            // J_inv(0,1) * d_flux_ref[X1]; d_flux_phys[X1] =
            // J_inv(1,0) * d_flux_ref[X0] + J_inv(1,1) *
            // d_flux_ref[X1];

            d_flux_phys = J_inv * d_flux_ref;

            flux_data.m_flux_deriv_in_qd_pts[icomp](q, offset + eq) = d_flux_phys[icomp];
            // flux_data.m_flux_deriv_in_qd_pts[XI1](q, offset + eq)
            // = d_flux_phys[X1];
          }
        } // Loop over quadrature points

      } // Loop over all cells stored in m_mdata_index

    } // Loop over coordinate components

  } // Static method 'evaluate_metric_terms'
};

// ----------------------------------------------------------------------------
// Specializations of FluxSpaceMetric::MetricComputer for some combinations
// of dimensions  - here 3D volume metric (GDIM = 3, TDIM = 3)
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
template <Uint Dummy>
class FluxSpaceMetric<MeshConfig, Physics, DIM>::MetricComputer<_3D, _3D, Dummy>
{
  public:
  static void evaluate_metric_terms(
      GeometryCache<MeshConfig::GDIM> const &gcache,
      GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &geo_met,
      SolutionCache const &scache, SolutionSpaceMetric<MeshConfig, DIM> const &sol_met,
      FluxSpaceMetric<MeshConfig, Physics, DIM> &flux_met)
  {
    using geo_cell_metric_type =
        typename GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM>::cellwise_metric;

    // typedef typename SolutionSpaceMetric<MeshConfig,
    // DIM>::cellwise_metric sol_cell_metric_type;

    using flux_data_type = typename FluxSpaceMetric<MeshConfig, Physics, DIM>::MetricData;

    math::DenseSMat<Real, _3D, _3D> J, inv_J;
    math::DenseSVec<Real, _3D> d_phi_ref; // Derivatives of shape function
                                          // in one point

    typename Physics::Properties phys_properties;

    math::DenseSMat<Real, Physics::NEQ, Physics::DIM> tmp_flux;

    math::DenseSMat<Real, Physics::NEQ, Physics::DIM> solution_grad;

    // ============================================================
    // 1) LOOP OVER ALL CELLS AND COMPUTE FLUX VALUES IN FLUX NODES
    // ============================================================

    for (Uint c = 0; c < flux_met.m_mdata_index.size(); ++c)
    {
      const math::DenseConstMatView<Real> geo_dof_coords = gcache.cell_values(c);
      const math::DenseConstMatView<Real> sol_dof_values = scache.cell_values(c);

      flux_data_type &flux_data = (*flux_met.m_mdata_index[c].first);

      // Compute physical coordinates in the nodes of flux element
      // [ (nb. flux DOF) x GDIM] = [ (nb. flux DOF) x (nb. sol. DOF)] *
      // [(nb. sol. DOF) x GDIM ]
      flux_data.Gf = flux_data.Vgf * geo_dof_coords;

      // Compute derivatives of transformation physical - reference space
      // in the nodes of flux element
      for (Uint d = 0; d < DIM; ++d)
      {
        flux_data.dGf[d] = flux_data.dVgf[d] * geo_dof_coords;
      }

      // Compute solution values in the nodes of flux element
      // [ (nb. flux DOF) x NEQ]  =  [ (nb. flux DOF) x (nb. sol. DOF)] *
      // [(nb. sol. DOF) x NEQ ]
      flux_data.Sf = flux_data.Vsf * sol_dof_values;

      for (Uint nf = 0; nf < flux_data.nb_flux_dof_per_elem; ++nf)
      {
        J(XI0, X0) = flux_data.dGf[XI0](nf, X0); // dX/dXi
        J(XI0, X1) = flux_data.dGf[XI0](nf, X1); // dY/dXi
        J(XI0, X2) = flux_data.dGf[XI0](nf, X2); // dZ/dXi

        J(XI1, X0) = flux_data.dGf[XI1](nf, X0); // dX/dEta
        J(XI1, X1) = flux_data.dGf[XI1](nf, X1); // dY/dEta
        J(XI1, X2) = flux_data.dGf[XI1](nf, X2); // dZ/dEta

        J(XI2, X0) = flux_data.dGf[XI2](nf, X0); // dX/dZeta
        J(XI2, X1) = flux_data.dGf[XI2](nf, X1); // dY/dZeta
        J(XI2, X2) = flux_data.dGf[XI2](nf, X2); // dZ/dZeta

        J.inv(inv_J);

        for (Uint ns = 0; ns < flux_data.nb_dof_per_elem; ++ns)
        {
          d_phi_ref[XI0] = flux_data.dVsf_ref[XI0](nf, ns);
          d_phi_ref[XI1] = flux_data.dVsf_ref[XI1](nf, ns);
          d_phi_ref[XI2] = flux_data.dVsf_ref[XI2](nf, ns);

          flux_data.dVsf_phys[XI0](nf, ns) = inv_J(XI0, X0) * d_phi_ref[X0] +
                                             inv_J(XI0, X1) * d_phi_ref[X1] +
                                             inv_J(XI0, X2) * d_phi_ref[X2];

          flux_data.dVsf_phys[XI1](nf, ns) = inv_J(XI1, X0) * d_phi_ref[X0] +
                                             inv_J(XI1, X1) * d_phi_ref[X1] +
                                             inv_J(XI1, X2) * d_phi_ref[X2];

          flux_data.dVsf_phys[XI2](nf, ns) = inv_J(XI2, X0) * d_phi_ref[X0] +
                                             inv_J(XI2, X1) * d_phi_ref[X1] +
                                             inv_J(XI2, X2) * d_phi_ref[X2];
        }
      }

      // Compute derivatives of solution (in the physical space)
      // in the nodes of flux element and store them in dSf
      for (Uint d = 0; d < DIM; ++d)
      {
        // [(nb. flux dof) x  NEQ] = [(nb. flux dof) x (nb. sol. dof)] *
        // [(nb. sol. dof) x NEQ]
        flux_data.dSf[d] = flux_data.dVsf_phys[d] * sol_dof_values;
      }

      // Offset determines where in m_flux_in_flux_nodes are stored values
      // for cell with 'global' buffer index c, which in given block of
      // cells is on position m_mdata_index[c].second
      const Uint offset = flux_met.m_mdata_index[c].second * Physics::NEQ;

      // Compute fluxes in flux dofs of one cell
      for (Uint n = 0; n < flux_data.nb_flux_dof_per_elem; ++n)
      {
        solution_grad.insert_col(X0, flux_data.dSf[X0].const_row(n));
        solution_grad.insert_col(X1, flux_data.dSf[X1].const_row(n));
        solution_grad.insert_col(X1, flux_data.dSf[X2].const_row(n));

        Physics::compute_properties(flux_data.Gf.const_row_transp(n),
                                    flux_data.Sf.const_row_transp(n), solution_grad,
                                    phys_properties);

        Physics::flux(phys_properties, tmp_flux);

        // m_flux_in_flux_nodes is of size [nb. qd_pts x (nb_blocks
        // NEQ)]

        for (Uint eq = 0; eq < Physics::NEQ; ++eq)
        {
          flux_data.m_flux_in_flux_nodes[X0](n, offset + eq) = tmp_flux(eq, X0);
          flux_data.m_flux_in_flux_nodes[X1](n, offset + eq) = tmp_flux(eq, X1);
          flux_data.m_flux_in_flux_nodes[X2](n, offset + eq) = tmp_flux(eq, X2);
        }

      } // Loop over nodes of one flux element

    } // Loop over all cells in buffer

    // ============================================================
    // 2) Compute fluxes in quadrature points
    // ============================================================

    for (typename common::DataMap<mesh::PointSetTagExt, flux_data_type>::iterator it =
             flux_met.m_metric_data_map.begin();
         it != flux_met.m_metric_data_map.end(); ++it)
    {
      flux_data_type &fdata = (*it.data_ptr());

      if (fdata.nb_elem_filled > 0)
      {
        // [(nb. qd. pts) x NEQ] = [ (nb. qd. pts) x (nb. flux DOF) ]  *
        // [(nb. flux DOF) x NEQ]
        fdata.m_flux_in_qd_pts[X0] = fdata.Vf * fdata.m_flux_in_flux_nodes[X0];
        fdata.m_flux_in_qd_pts[X1] = fdata.Vf * fdata.m_flux_in_flux_nodes[X1];
        fdata.m_flux_in_qd_pts[X2] = fdata.Vf * fdata.m_flux_in_flux_nodes[X2];
      }
    }

    // ============================================================
    // 3) Compute derivatives of fluxes.
    //    For each coordinate component icomp, compute the derivative
    //    of the i-th flux component with respect to i-th coordinate
    //    variable.
    //    In other words, compute dF0/dx, dF1/dy, dF2/dz (in 3D)
    // ============================================================

    math::DenseSVec<Real, _3D> d_flux_ref, d_flux_phys;

    for (Uint icomp = 0; icomp < Physics::DIM; ++icomp)
    {

      for (typename common::DataMap<mesh::PointSetTagExt, flux_data_type>::iterator it =
               flux_met.m_metric_data_map.begin();
           it != flux_met.m_metric_data_map.end(); ++it)
      {
        flux_data_type &fdata = (*it.data_ptr());

        if (fdata.nb_elem_filled > 0)
        {
          // Flux derivatives
          fdata.m_tmp_flux_deriv[X0] = fdata.dVf[X0] * fdata.m_flux_in_flux_nodes[icomp];
          fdata.m_tmp_flux_deriv[X1] = fdata.dVf[X1] * fdata.m_flux_in_flux_nodes[icomp];
          fdata.m_tmp_flux_deriv[X2] = fdata.dVf[X2] * fdata.m_flux_in_flux_nodes[icomp];
        }
      }

      // The jacobians can be provided by geometry metric
      for (Uint c = 0; c < flux_met.m_mdata_index.size(); ++c)
      {
        const geo_cell_metric_type cell_geo_values = geo_met.cellwise_values(c);

        flux_data_type &flux_data = (*flux_met.m_mdata_index[c].first);

        const Uint offset = flux_met.m_mdata_index[c].second * Physics::NEQ;

        for (Uint q = 0; q < flux_data.nb_qd_pts_per_elem; ++q)
        {
          const math::DenseConstMatView<Real> J_inv = cell_geo_values.inv_jacobi(q);

          for (Uint eq = 0; eq < Physics::NEQ; ++eq)
          {
            d_flux_ref[XI0] = flux_data.m_tmp_flux_deriv[XI0](q, offset + eq);
            d_flux_ref[XI1] = flux_data.m_tmp_flux_deriv[XI1](q, offset + eq);
            d_flux_ref[XI2] = flux_data.m_tmp_flux_deriv[XI2](q, offset + eq);

            // d_flux_phys[X0] = J_inv(0,0) * d_flux_ref[X0] +
            // J_inv(0,1) * d_flux_ref[X1]; d_flux_phys[X1] =
            // J_inv(1,0) * d_flux_ref[X0] + J_inv(1,1) *
            // d_flux_ref[X1];

            d_flux_phys = J_inv * d_flux_ref;

            flux_data.m_flux_deriv_in_qd_pts[icomp](q, offset + eq) = d_flux_phys[icomp];
            // flux_data.m_flux_deriv_in_qd_pts[XI1](q, offset + eq)
            // = d_flux_phys[X1];
          }
        } // Loop over quadrature points

      } // Loop over all cells stored in m_mdata_index

    } // Loop over coordinate components

  } // Static method 'evaluate_metric_terms'
};

// ----------------------------------------------------------------------------
// Specializations of FluxSpaceMetric::MetricComputer for some combinations
// of dimensions - here edge metric in 2D space (GDIM = 2, TDIM = 1)
// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
template <Uint Dummy>
class FluxSpaceMetric<MeshConfig, Physics, DIM>::MetricComputer<_2D, _1D, Dummy>
{
  public:
  static void evaluate_metric_terms(
      GeometryCache<MeshConfig::GDIM> const &gcache,
      GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &geo_met,
      SolutionCache const &scache, SolutionSpaceMetric<MeshConfig, DIM> const &sol_met,
      FluxSpaceMetric<MeshConfig, Physics, DIM> &flux_met)
  {
    using flux_data_type = typename FluxSpaceMetric<MeshConfig, Physics, DIM>::MetricData;

    typename Physics::Properties phys_properties;

    math::DenseSMat<Real, Physics::NEQ, Physics::DIM> tmp_flux;

    math::DenseSMat<Real, Physics::NEQ, Physics::DIM> solution_grad;
    solution_grad.fill(0.0);

    // ============================================================
    // 1) LOOP OVER ALL CELLS AND COMPUTE FLUX VALUES IN FLUX NODES
    // ============================================================

    for (Uint c = 0; c < flux_met.m_mdata_index.size(); ++c)
    {
      const math::DenseConstMatView<Real> geo_dof_coords = gcache.cell_values(c);
      const math::DenseConstMatView<Real> sol_dof_values = scache.cell_values(c);

      flux_data_type &flux_data = (*flux_met.m_mdata_index[c].first);

      // Compute physical coordinates in the nodes of flux element
      // [ (nb. flux DOF) x GDIM] = [ (nb. flux DOF) x (nb. sol. DOF)] *
      // [(nb. sol. DOF) x GDIM ]
      flux_data.Gf = flux_data.Vgf * geo_dof_coords;

      // Compute solution values in the nodes of flux element
      // [ (nb. flux DOF) x NEQ]  =  [ (nb. flux DOF) x (nb. sol. DOF)] *
      // [(nb. sol. DOF) x NEQ ]
      flux_data.Sf = flux_data.Vsf * sol_dof_values;

      // Offset determines where in m_flux_in_flux_nodes are stored values
      // for cell with 'global' buffer index c, which in given block of
      // cells is on position m_mdata_index[c].second
      const Uint offset = flux_met.m_mdata_index[c].second * Physics::NEQ;

      // Compute fluxes in flux dofs of one cell
      for (Uint n = 0; n < flux_data.nb_flux_dof_per_elem; ++n)
      {
        Physics::compute_properties(flux_data.Gf.const_row_transp(n),
                                    flux_data.Sf.const_row_transp(n), solution_grad,
                                    phys_properties);

        Physics::flux(phys_properties, tmp_flux);

        // m_flux_in_flux_nodes is of size [nb. qd_pts x (nb_blocks
        // NEQ)]

        for (Uint eq = 0; eq < Physics::NEQ; ++eq)
        {
          flux_data.m_flux_in_flux_nodes[X0](n, offset + eq) = tmp_flux(eq, X0);
          flux_data.m_flux_in_flux_nodes[X1](n, offset + eq) = tmp_flux(eq, X1);
        }

      } // Loop over nodes of one flux element

    } // Loop over all cells in buffer

    // ============================================================
    // 2) Compute fluxes in quadrature points
    // ============================================================

    for (typename common::DataMap<mesh::PointSetTagExt, flux_data_type>::iterator it =
             flux_met.m_metric_data_map.begin();
         it != flux_met.m_metric_data_map.end(); ++it)
    {
      flux_data_type &fdata = (*it.data_ptr());

      if (fdata.nb_elem_filled > 0)
      {
        // [(nb. qd. pts) x NEQ] = [ (nb. qd. pts) x (nb. flux DOF) ]  *
        // [(nb. flux DOF) x NEQ]
        fdata.m_flux_in_qd_pts[X0] = fdata.Vf * fdata.m_flux_in_flux_nodes[X0];
        fdata.m_flux_in_qd_pts[X1] = fdata.Vf * fdata.m_flux_in_flux_nodes[X1];
      }
    }

  } // Static method 'evaluate_metric_terms'
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint DIM>
template <Uint Dummy>
class FluxSpaceMetric<MeshConfig, Physics, DIM>::MetricComputer<_3D, _2D, Dummy>
{
  public:
  static void evaluate_metric_terms(
      GeometryCache<MeshConfig::GDIM> const &gcache,
      GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> const &geo_met,
      SolutionCache const &scache, SolutionSpaceMetric<MeshConfig, DIM> const &sol_met,
      FluxSpaceMetric<MeshConfig, Physics, DIM> &flux_met)
  {
    using flux_data_type = typename FluxSpaceMetric<MeshConfig, Physics, DIM>::MetricData;

    typename Physics::Properties phys_properties;

    math::DenseSMat<Real, Physics::NEQ, Physics::DIM> tmp_flux;

    math::DenseSMat<Real, Physics::NEQ, Physics::DIM> solution_grad;
    solution_grad.fill(0.0);

    // ============================================================
    // 1) LOOP OVER ALL CELLS AND COMPUTE FLUX VALUES IN FLUX NODES
    // ============================================================

    for (Uint c = 0; c < flux_met.m_mdata_index.size(); ++c)
    {
      const math::DenseConstMatView<Real> geo_dof_coords = gcache.cell_values(c);
      const math::DenseConstMatView<Real> sol_dof_values = scache.cell_values(c);

      flux_data_type &flux_data = (*flux_met.m_mdata_index[c].first);

      // Compute physical coordinates in the nodes of flux element
      // [ (nb. flux DOF) x GDIM] = [ (nb. flux DOF) x (nb. sol. DOF)] *
      // [(nb. sol. DOF) x GDIM ]
      flux_data.Gf = flux_data.Vgf * geo_dof_coords;

      // Compute solution values in the nodes of flux element
      // [ (nb. flux DOF) x NEQ]  =  [ (nb. flux DOF) x (nb. sol. DOF)] *
      // [(nb. sol. DOF) x NEQ ]
      flux_data.Sf = flux_data.Vsf * sol_dof_values;

      // Offset determines where in m_flux_in_flux_nodes are stored values
      // for cell with 'global' buffer index c, which in given block of
      // cells is on position m_mdata_index[c].second
      const Uint offset = flux_met.m_mdata_index[c].second * Physics::NEQ;

      // Compute fluxes in flux dofs of one cell
      for (Uint n = 0; n < flux_data.nb_flux_dof_per_elem; ++n)
      {
        Physics::compute_properties(flux_data.Gf.const_row_transp(n),
                                    flux_data.Sf.const_row_transp(n), solution_grad,
                                    phys_properties);

        Physics::flux(phys_properties, tmp_flux);

        // m_flux_in_flux_nodes is of size [nb. qd_pts x (nb_blocks
        // NEQ)]

        for (Uint eq = 0; eq < Physics::NEQ; ++eq)
        {
          flux_data.m_flux_in_flux_nodes[X0](n, offset + eq) = tmp_flux(eq, X0);
          flux_data.m_flux_in_flux_nodes[X1](n, offset + eq) = tmp_flux(eq, X1);
          flux_data.m_flux_in_flux_nodes[X2](n, offset + eq) = tmp_flux(eq, X2);
        }

      } // Loop over nodes of one flux element

    } // Loop over all cells in buffer

    // ============================================================
    // 2) Compute fluxes in quadrature points
    // ============================================================

    for (typename common::DataMap<mesh::PointSetTagExt, flux_data_type>::iterator it =
             flux_met.m_metric_data_map.begin();
         it != flux_met.m_metric_data_map.end(); ++it)
    {
      flux_data_type &fdata = (*it.data_ptr());

      if (fdata.nb_elem_filled > 0)
      {
        // [(nb. qd. pts) x NEQ] = [ (nb. qd. pts) x (nb. flux DOF) ]  *
        // [(nb. flux DOF) x NEQ]
        fdata.m_flux_in_qd_pts[X0] = fdata.Vf * fdata.m_flux_in_flux_nodes[X0];
        fdata.m_flux_in_qd_pts[X1] = fdata.Vf * fdata.m_flux_in_flux_nodes[X1];
        fdata.m_flux_in_qd_pts[X2] = fdata.Vf * fdata.m_flux_in_flux_nodes[X2];
      }
    }

  } // Static method 'evaluate_metric_terms'
};

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
