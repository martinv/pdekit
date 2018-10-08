#ifndef PDEKIT_PG_RDM_Picard_Jac_Cell_Assembly_Worker_hpp
#define PDEKIT_PG_RDM_Picard_Jac_Cell_Assembly_Worker_hpp

#include "solver/rdm/assembly/PGRDMImplicitCellAssemblyWorker.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

namespace detail
{

#define DO_TIMINGS 0

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
class PGRDMPicardJacCellAssemblyWorker : public PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>
{
  public:
  /// TYPEDEFS
  using tria_t    = typename result_of::tria_t<MeshConfig>;
  using dof_map_t = typename result_of::dof_map_t<MeshConfig>;
  using method_data_type =
      typename internal::CellSchemeSelector<Physics, SchemeTraits>::type::method_data;
  using f_space_cells = interpolation::FunctionSpace<MeshConfig>;

  PGRDMPicardJacCellAssemblyWorker(const Uint worker_id);

  ~PGRDMPicardJacCellAssemblyWorker();

  void configure_cell_spaces(tria_t const &tria, dof_map_t const &sol_dofs,
                             const Uint first_cell_idx, const Uint last_cell_idx,
                             const Uint nb_blocks, const SFunc sf_type, const PointSetID quad_type,
                             const Uint quad_order) override;

  void assemble_mat_and_rhs_part(const tria_t &tria, const dof_map_t &sol_dofs,
                                 const interpolation::VectorMeshFunction<Real> &solution,
                                 RDTimeUpdate &time_update,
                                 const std::vector<bool> &is_Dirichlet_node,
                                 ls::TpetraCrsMatrix<Real> &mat,
                                 ls::TpetraMultiVector<Real> &rhs) override;

  void assemble_rhs_part(const tria_t &tria, const dof_map_t &sol_dofs,
                         const interpolation::VectorMeshFunction<Real> &solution,
                         RDTimeUpdate &time_update, const std::vector<bool> &is_Dirichlet_node,
                         ls::TpetraMultiVector<Real> &rhs) override;

  private:
  using base = PGRDMImplicitCellAssemblyWorker<MeshConfig, Physics>;

  enum
  {
    NEQ = Physics::NEQ
  };

  using geo_metric_type  = typename base::geo_metric_type;
  using sol_metric_type  = typename base::sol_metric_type;
  using flux_metric_type = typename base::flux_metric_type;

  /// A map which associates to each element type in solution space a
  /// corresponding reference element for fluxes
  common::DataMap<mesh::PointSetTagExt, method_data_type> m_rdm_method_data;

  /// Splitter to be used
  typename internal::CellSchemeSelector<Physics, SchemeTraits>::type m_scheme;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
PGRDMPicardJacCellAssemblyWorker<
    MeshConfig, Physics, SchemeTraits>::PGRDMPicardJacCellAssemblyWorker(const Uint worker_id)
    : base(worker_id)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
PGRDMPicardJacCellAssemblyWorker<MeshConfig, Physics,
                                 SchemeTraits>::~PGRDMPicardJacCellAssemblyWorker()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMPicardJacCellAssemblyWorker<MeshConfig, Physics, SchemeTraits>::configure_cell_spaces(
    tria_t const &tria, dof_map_t const &sol_dofs, const Uint first_cell_idx,
    const Uint last_cell_idx, const Uint nb_blocks, const SFunc sf_type, const PointSetID quad_type,
    const Uint quad_order)
{
  base::configure_cell_spaces(tria, sol_dofs, first_cell_idx, last_cell_idx, nb_blocks, sf_type,
                              quad_type, quad_order);

  common::DataMap<mesh::PointSetTagExt, interpolation::FEValues> const &elem_type_map =
      base::m_sol_cell_space.reference_elements();

  m_rdm_method_data.clear();

  for (common::DataMap<mesh::PointSetTagExt, interpolation::FEValues>::const_iterator it =
           elem_type_map.cbegin();
       it != elem_type_map.cend(); ++it)
  {
    common::PtrHandle<method_data_type> method_data = m_rdm_method_data.create(it.key_value());

    (*method_data).resize_variables(*it.data_ptr());
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMPicardJacCellAssemblyWorker<MeshConfig, Physics, SchemeTraits>::assemble_mat_and_rhs_part(
    const tria_t &tria, const dof_map_t &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &solution, RDTimeUpdate &time_update,
    const std::vector<bool> &is_Dirichlet_node, ls::TpetraCrsMatrix<Real> &mat,
    ls::TpetraMultiVector<Real> &rhs)
{
#if DO_TIMINGS
  std::ofstream outfile;
  std::stringstream ss;
  ss << "cell_worker" << std::setfill('0') << std::setw(3) << m_worker_id << "_lhs_rhs_timing.dat";
  outfile.open(ss.str(), std::ios::app);
#endif

  /*
  auto sf_generator = [this](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, base::m_sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [this](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, base::m_quad_order, base::m_quad_type);
  };
  */

  const bool has_source_term     = base::m_sources ? true : false;
  const bool has_blending_coeff  = base::m_blending_coeff ? true : false;
  const bool has_artif_viscosity = base::m_art_viscosity ? true : false;

  const interpolation::ScalarMeshFunction<Real> &nodal_dual_volume =
      time_update.nodal_dual_volume();

  for (const typename dof_map_t::const_dof_range_typed &dof_group :
       sol_dofs.all_active_dof_groups())
  {
    typename dof_map_t::const_dof_iterator_typed sol_dof_block_begin = dof_group.begin();
    typename dof_map_t::const_dof_iterator_typed sol_dof_block_end   = dof_group.begin();

    const Uint nb_dof_per_elem = sol_dof_block_begin->mesh_entity().nb_vert();

    const Uint required_mat_cache_capacity =
        base::m_nb_blocks * (nb_dof_per_elem * NEQ) * (nb_dof_per_elem * NEQ);

    const Uint required_vec_cache_capacity = base::m_nb_blocks * (nb_dof_per_elem * NEQ);

    if (base::m_mat_buffer.capacity() < required_mat_cache_capacity)
    {
      base::m_mat_buffer.reserve(required_mat_cache_capacity);
    }
    if (base::m_rhs_buffer.capacity() < required_vec_cache_capacity)
    {
      base::m_rhs_buffer.reserve(required_vec_cache_capacity);
    }
    if (base::m_time_update_buffer.capacity() < required_vec_cache_capacity)
    {
      base::m_time_update_buffer.reserve(required_vec_cache_capacity);
    }

    bool process_cells_in_group = true;

    while (process_cells_in_group)
    {
      base::m_mat_buffer.resize(0);
      base::m_rhs_buffer.resize(0);
      base::m_time_update_buffer.resize(0);

      // **************************************************************************
      // 1a) FILL GEOMETRY CACHE
      // 1b) FILL SOLUTION CACHE AND OPTIONALLY SOURCE AND STABILIZATION
      // CACHE
      // **************************************************************************

      base::m_geo_cache.flush();
      base::m_sol_cache.flush();
      base::m_res_cache.flush();

      if (has_source_term)
      {
        base::m_source_cache.flush();
      }
      if (has_artif_viscosity)
      {
        base::m_stab_coeff_cache.flush();
      }

      Uint cell_idx_in_metric = 0;

      while ((sol_dof_block_end != dof_group.end()) && (cell_idx_in_metric < base::m_nb_blocks))
      {
        const mesh::CellTopologyView<MeshConfig> tcell_view = sol_dof_block_end->tcell();
        const mesh::MeshEntity sol_cell                     = sol_dof_block_end->mesh_entity();

        if ((base::m_first_cell_idx <= sol_cell.idx()) && (sol_cell.idx() <= base::m_last_cell_idx))
        {
          const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell_view.coordinates();

#if 0
          const mesh::PointSetTag geo_pt_set_tag = geo_cell.std_region_id();
          const mesh::PointSetTagExt active_geo_cell_tag_ext(geo_pt_set_tag, P0,
                                                             mesh::CellTransform::NO_TRANS, 0u);
          const mesh::PointSetTagExt active_sol_cell_tag_ext(sol_cell.std_region_id(), P0,
                                                             mesh::CellTransform::NO_TRANS, 0u);

          const mesh::sf::SFTag geo_sf_tag =
              sf_generator(geo_pt_set_tag.elem_shape(), geo_pt_set_tag.poly_order());
          const mesh::PointSetTag geo_quad_tag =
              quad_generator(geo_pt_set_tag.elem_shape(), geo_pt_set_tag.poly_order());
          const mesh::PointSetTagExt geo_quad_tag_ext(geo_quad_tag, P0,
                                                      mesh::CellTransform::NO_TRANS, 0);
          const mesh::DiscreteElemKey geo_key(active_geo_cell_tag_ext, geo_sf_tag,
                                              geo_quad_tag_ext);
#endif

          const mesh::DiscreteElemKey geo_key = base::m_geo_key_cache.key(cell_idx_in_metric);
          const mesh::PointSetTagExt active_sol_cell_tag_ext(sol_cell.pt_set_id(), P0,
                                                             mesh::CellTransform::NO_TRANS, 0u);

          base::m_geo_cache.push_back_to_buffer(geo_cell_coords, geo_key);

          // std::cout << "Solution cell [" << cell_iter->idx() << "]
          // " << *cell_iter << std::endl;
          base::m_sol_cache.push_back_to_buffer(sol_cell, solution, active_sol_cell_tag_ext);

          if (has_source_term)
          {
            base::m_source_cache.push_back_to_buffer(sol_cell, *base::m_sources,
                                                     active_sol_cell_tag_ext);
          }
          if (has_artif_viscosity)
          {
            base::m_stab_coeff_cache.push_back_to_buffer(sol_cell, *base::m_art_viscosity,
                                                         active_sol_cell_tag_ext);
          }
          cell_idx_in_metric++;
        } // If solution cell index is within range

        ++sol_dof_block_end;
      }

      // **************************************************************************
      // 2) EVALUATE METRIC TERMS
      // **************************************************************************

      base::m_geo_metric.empty_buffer();
      base::m_sol_metric.empty_buffer();
      base::m_flux_metric.empty_buffer();

      base::m_geo_metric.evaluate(base::m_geo_cache, interpolation::RebuildMetricIndex{true});
      base::m_sol_metric.evaluate(base::m_geo_metric, base::m_sol_cache,
                                  interpolation::ComputeMetricDerivs{true},
                                  interpolation::RebuildMetricIndex{true});
      base::m_flux_metric.evaluate(base::m_geo_cache, base::m_geo_metric, base::m_sol_cache,
                                   base::m_sol_metric, interpolation::RebuildMetricIndex{true});

      if (has_source_term)
      {
        base::m_source_metric.empty_buffer();
        base::m_source_metric.evaluate(base::m_geo_metric, base::m_source_cache,
                                       interpolation::ComputeMetricDerivs{false},
                                       interpolation::RebuildMetricIndex{true});
      }

      if (has_artif_viscosity)
      {
        base::m_stab_coeff_metric.empty_buffer();
        base::m_stab_coeff_metric.evaluate(base::m_geo_metric, base::m_stab_coeff_cache,
                                           interpolation::ComputeMetricDerivs{false},
                                           interpolation::RebuildMetricIndex{true});
      }

      cell_idx_in_metric = 0;

      // **************************************************************************
      // 3) COMPUTE THE RESIDUALS AND CELL JACOBIANS AND STORE THEM
      // **************************************************************************

      // Since we are still looping over element groups, we can get the
      // method data for a particular element type before the actual loop
      // over cells
      const mesh::PointSetTag cell_type_id = dof_group.begin()->pt_set_id();

      common::PtrHandle<method_data_type> method_data = m_rdm_method_data.std_region_data(
          mesh::PointSetTagExt(cell_type_id, P0, mesh::CellTransform::NO_TRANS, 0u));

      Uint rhs_fill = 0;

      for (typename dof_map_t::const_dof_iterator_typed cell_iter = sol_dof_block_begin;
           cell_iter != sol_dof_block_end; ++cell_iter)
      {
        const mesh::MeshEntity solution_elem = cell_iter->mesh_entity();

        if ((base::m_first_cell_idx <= solution_elem.idx()) &&
            (solution_elem.idx() <= base::m_last_cell_idx))
        {
          const math::DenseConstMatView<Real> sol_nodal_values =
              base::m_sol_cache.cell_values(cell_idx_in_metric);

          /*
          const typename geo_metric_type::cellwise_metric geo_cell_met
          = base::m_geo_metric.cellwise_values(cell_idx_in_metric);

          const typename sol_metric_type::cellwise_metric sol_cell_met
          = base::m_sol_metric.cellwise_values(cell_idx_in_metric);

          const typename flux_metric_type::cellwise_metric
          flux_cell_met =
              base::m_flux_metric.cellwise_values(cell_idx_in_metric);
          */

          base::m_const_method_data.CGM = base::m_geo_metric.cellwise_values(cell_idx_in_metric);
          base::m_const_method_data.CSM = base::m_sol_metric.cellwise_values(cell_idx_in_metric);
          base::m_const_method_data.CFM = base::m_flux_metric.cellwise_values(cell_idx_in_metric);

          (*method_data).m_use_external_theta = false;
          if (has_blending_coeff)
          {
            (*method_data).m_use_external_theta = true;
            //(*method_data).base::m_blending_coeff.fill((*base::m_blending_coeff)[solution_elem.idx()]);

            for (Uint v = 0; v < solution_elem.nb_vert(); ++v)
            {
              (*method_data).m_blending_coeff[v] =
                  (*base::m_blending_coeff)[solution_elem.vertex(v)];
            }
          }

          if (has_artif_viscosity)
          {
            /*
            const Real art_visc =
            (*base::m_art_viscosity)[solution_elem.idx()];
            (*method_data).m_art_visc.fill(art_visc);
            */

            const typename sol_metric_type::cellwise_metric stab_coeff_cell_met =
                base::m_stab_coeff_metric.cellwise_values(cell_idx_in_metric);

            const math::DenseConstMatView<Real> stab_coeff_at_qd_pt =
                stab_coeff_cell_met.field_values();

            (*method_data).m_art_visc = stab_coeff_at_qd_pt.col(0);
          }

          if (has_source_term)
          {
            base::m_const_method_data.CSrcM =
                base::m_source_metric.cellwise_values(cell_idx_in_metric);
            m_scheme.compute_adv_reaction_residuals(sol_nodal_values, base::m_const_method_data,
                                                    *method_data);
          }
          else
          {
            m_scheme.compute_adv_res_and_jacobian(sol_nodal_values, base::m_const_method_data,
                                                  *method_data);
          }

          math::DenseMatView<Real> const elem_jacobian   = (*method_data).elem_jacobian();
          math::DenseDVec<Real> const &elem_update_coeff = (*method_data).m_elem_wave_speed;

          for (Uint m = 0; m < (*method_data).nb_nodes(); ++m)
          {
            // Get a slice of all nodal residuals computer in this
            // element. The 'slice' is only the portion
            // corresponding to node n
            math::DenseVecView<Real> nodal_res = (*method_data).elem_node_res_block(m * NEQ, NEQ);

            const Real inv_nodal_volume = 1. / nodal_dual_volume[solution_elem.vertex(m)];

            for (Uint eq_m = 0; eq_m < NEQ; ++eq_m)
            {
              nodal_res[eq_m] *= inv_nodal_volume;

              const Uint dof_idx = solution_elem.vertex(m) * NEQ + eq_m;

              if (!is_Dirichlet_node[dof_idx])
              {
                // rhs(dof_idx, 0) -= elem_res[n][eq];
                base::m_rhs_buffer.push_back(std::tuple<Uint, Real>(dof_idx, -nodal_res[eq_m]));
                rhs_fill++;
              }
            }

            // -------------------------------------
            // Put cell Jacobian entries into buffer
            // -------------------------------------

            for (Uint n = 0; n < (*method_data).nb_nodes(); ++n)
            {
              for (Uint eq_m = 0; eq_m < NEQ; ++eq_m)
              {
                const Uint global_row_idx = NEQ * solution_elem.vertex(m) + eq_m;

                if (!is_Dirichlet_node[global_row_idx])
                {
                  for (Uint eq_n = 0; eq_n < NEQ; ++eq_n)
                  {
                    const Uint global_col_idx = NEQ * solution_elem.vertex(n) + eq_n;

                    const Real jac_entry =
                        inv_nodal_volume * elem_jacobian(m * NEQ + eq_m, n * NEQ + eq_n);
                    base::m_mat_buffer.push_back(
                        std::tuple<Uint, Uint, Real>(global_row_idx, global_col_idx, jac_entry));
                  }
                } // If this is not Dirichlet DOF
              }
            } // Loop over local columns

            base::m_time_update_buffer.push_back(
                std::tuple<Uint, Real>(solution_elem.vertex(m), elem_update_coeff[m]));
          } // Loop over m local nodes

          const math::DenseVecView<Real> elem_node_res = (*method_data).elem_node_res();
          base::m_res_cache.push_vec_to_buffer(solution_elem, elem_node_res,
                                               mesh::PointSetTagExt(solution_elem.pt_set_id(), P0,
                                                                    mesh::CellTransform::NO_TRANS,
                                                                    0u));

          cell_idx_in_metric++;

        } // If the index of solution cell is within range

      } // Loop over all cells of one block

      // -----------------------------------------------
      // Accumulate the buffer values into system matrix
      // -----------------------------------------------
      mat.add_values(base::m_mat_buffer);

// --------------------------------------------
// Accumulate the buffer values into RHS vector
// --------------------------------------------
#if DO_TIMINGS
      std::clock_t c_start = std::clock();
      auto t_start         = std::chrono::high_resolution_clock::now();
#endif

      rhs.add_values(0, base::m_rhs_buffer);

#if DO_TIMINGS
      std::clock_t c_end  = std::clock();
      auto t_end          = std::chrono::high_resolution_clock::now();
      double cpu_duration = (c_end - c_start) / (double)CLOCKS_PER_SEC;
      auto wclock_duration =
          0.001 * std::chrono::duration<double, std::milli>(t_end - t_start).count();
      outfile << "rhs accumulation cpu/wall  " << std::setw(10) << std::setprecision(10)
              << cpu_duration << " " << std::setw(10) << std::setprecision(10) << wclock_duration
              << std::endl;
#endif

      // --------------------------------------------
      // Accumulate the time update values
      // --------------------------------------------
      time_update.accumulate_wave_speeds(base::m_time_update_buffer);

      if (sol_dof_block_end == dof_group.end())
      {
        process_cells_in_group = false;
      }
      else
      {
        sol_dof_block_begin = sol_dof_block_end;
      }

    } // while(process_cells_in_group)

  } // Lop over all cell groups

#if DO_TIMINGS
  outfile.close();
#endif
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMPicardJacCellAssemblyWorker<MeshConfig, Physics, SchemeTraits>::assemble_rhs_part(
    const tria_t &tria, const dof_map_t &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &solution, RDTimeUpdate &time_update,
    const std::vector<bool> &is_Dirichlet_node, ls::TpetraMultiVector<Real> &rhs)
{
#if DO_TIMINGS
  std::ofstream outfile;
  std::stringstream ss;
  ss << "cell_worker" << std::setfill('0') << std::setw(3) << m_worker_id << "_rhs_timing.dat";
  outfile.open(ss.str(), std::ios::app);
#endif

  /*
  auto sf_generator = [this](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, base::m_sf_type, order, ModalBasis::Modal);
  };

  auto quad_generator = [this](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, base::m_quad_order, base::m_quad_type);
  };
  */

  const bool has_source_term     = base::m_sources ? true : false;
  const bool has_blending_coeff  = base::m_blending_coeff ? true : false;
  const bool has_artif_viscosity = base::m_art_viscosity ? true : false;

  const interpolation::ScalarMeshFunction<Real> &nodal_dual_volume =
      time_update.nodal_dual_volume();

  for (const typename dof_map_t::const_dof_range_typed &dof_group :
       sol_dofs.all_active_dof_groups())
  {
    typename dof_map_t::const_dof_iterator_typed sol_dof_block_begin = dof_group.begin();
    typename dof_map_t::const_dof_iterator_typed sol_dof_block_end   = dof_group.begin();

    const Uint nb_dof_per_elem = sol_dof_block_begin->mesh_entity().nb_vert();

    const Uint required_vec_cache_capacity = base::m_nb_blocks * (nb_dof_per_elem * NEQ);

    if (base::m_rhs_buffer.capacity() < required_vec_cache_capacity)
    {
      base::m_rhs_buffer.reserve(required_vec_cache_capacity);
    }
    if (base::m_time_update_buffer.capacity() < required_vec_cache_capacity)
    {
      base::m_time_update_buffer.reserve(required_vec_cache_capacity);
    }

    bool process_cells_in_group = true;

    while (process_cells_in_group)
    {
      base::m_rhs_buffer.resize(0);
      base::m_time_update_buffer.resize(0);

      // **************************************************************************
      // 1a) FILL GEOMETRY CACHE
      // 1b) FILL SOLUTION CACHE AND OPTIONALLY SOURCE CACHE
      // **************************************************************************

      base::m_geo_cache.flush();
      base::m_sol_cache.flush();

      if (has_source_term)
      {
        base::m_source_cache.flush();
      }
      if (has_artif_viscosity)
      {
        base::m_stab_coeff_cache.flush();
      }

      Uint cell_idx_in_metric = 0;

      while ((sol_dof_block_end != dof_group.end()) && (cell_idx_in_metric < base::m_nb_blocks))
      {
        const mesh::CellTopologyView<MeshConfig> tcell = sol_dof_block_end->tcell();
        const mesh::MeshEntity sol_cell                = sol_dof_block_end->mesh_entity();

        if ((base::m_first_cell_idx <= sol_cell.idx()) && (sol_cell.idx() <= base::m_last_cell_idx))
        {
          const mesh::CellGeometry<MeshConfig::GDIM> geo_cell_coords = tcell.coordinates();

#if 0
          const mesh::PointSetTag geo_pt_set_tag = geo_cell.std_region_id();
          const mesh::PointSetTagExt active_geo_cell_tag_ext(geo_pt_set_tag, P0,
                                                             mesh::CellTransform::NO_TRANS, 0u);
          const mesh::PointSetTagExt active_sol_cell_tag_ext(sol_cell.std_region_id(), P0,
                                                             mesh::CellTransform::NO_TRANS, 0u);

          const mesh::sf::SFTag geo_sf_tag =
              sf_generator(geo_pt_set_tag.elem_shape(), geo_pt_set_tag.poly_order());
          const mesh::PointSetTag geo_quad_tag =
              quad_generator(geo_pt_set_tag.elem_shape(), geo_pt_set_tag.poly_order());
          const mesh::PointSetTagExt geo_quad_tag_ext(geo_quad_tag, P0,
                                                      mesh::CellTransform::NO_TRANS, 0);
          const mesh::DiscreteElemKey geo_key(active_geo_cell_tag_ext, geo_sf_tag,
                                              geo_quad_tag_ext);
#endif

          const mesh::DiscreteElemKey geo_key = base::m_geo_key_cache.key(cell_idx_in_metric);
          const mesh::PointSetTagExt active_sol_cell_tag_ext(sol_cell.pt_set_id(), P0,
                                                             mesh::CellTransform::NO_TRANS, 0u);

          base::m_geo_cache.push_back_to_buffer(geo_cell_coords, geo_key);

          // std::cout << "Solution cell [" << cell_iter->idx() << "]
          // " << *cell_iter << std::endl;
          base::m_sol_cache.push_back_to_buffer(sol_cell, solution, active_sol_cell_tag_ext);

          if (has_source_term)
          {
            base::m_source_cache.push_back_to_buffer(sol_cell, *base::m_sources,
                                                     active_sol_cell_tag_ext);
          }
          if (has_artif_viscosity)
          {
            base::m_stab_coeff_cache.push_back_to_buffer(sol_cell, *base::m_art_viscosity,
                                                         active_sol_cell_tag_ext);
          }

          cell_idx_in_metric++;
        } // If solution cell index is within range

        ++sol_dof_block_end;
      }

      // **************************************************************************
      // 2) EVALUATE METRIC TERMS
      // **************************************************************************

      base::m_geo_metric.empty_buffer();
      base::m_sol_metric.empty_buffer();

      base::m_flux_metric.empty_buffer();

      base::m_geo_metric.evaluate(base::m_geo_cache, interpolation::RebuildMetricIndex{true});
      base::m_sol_metric.evaluate(base::m_geo_metric, base::m_sol_cache,
                                  interpolation::ComputeMetricDerivs{true},
                                  interpolation::RebuildMetricIndex{true});
      base::m_flux_metric.evaluate(base::m_geo_cache, base::m_geo_metric, base::m_sol_cache,
                                   base::m_sol_metric, interpolation::RebuildMetricIndex{true});

      if (has_source_term)
      {
        base::m_source_metric.empty_buffer();
        base::m_source_metric.evaluate(base::m_geo_metric, base::m_source_cache,
                                       interpolation::ComputeMetricDerivs{false},
                                       interpolation::RebuildMetricIndex{true});
      }

      if (has_artif_viscosity)
      {
        base::m_stab_coeff_metric.empty_buffer();
        base::m_stab_coeff_metric.evaluate(base::m_geo_metric, base::m_stab_coeff_cache,
                                           interpolation::ComputeMetricDerivs{false},
                                           interpolation::RebuildMetricIndex{true});
      }

      cell_idx_in_metric = 0;

      // **************************************************************************
      // 3) COMPUTE THE RESIDUALS AND STORE THEM
      // **************************************************************************

      const mesh::PointSetTag cell_type_id = dof_group.begin()->pt_set_id();

      common::PtrHandle<method_data_type> method_data = m_rdm_method_data.std_region_data(
          mesh::PointSetTagExt(cell_type_id, P0, mesh::CellTransform::NO_TRANS, 0u));

      Uint rhs_fill = 0;

      for (typename dof_map_t::const_dof_iterator_typed cell_iter = sol_dof_block_begin;
           cell_iter != sol_dof_block_end; ++cell_iter)
      {
        const mesh::MeshEntity solution_elem = cell_iter->mesh_entity();

        if ((base::m_first_cell_idx <= solution_elem.idx()) &&
            (solution_elem.idx() <= base::m_last_cell_idx))
        {
          const math::DenseConstMatView<Real> sol_nodal_values =
              base::m_sol_cache.cell_values(cell_idx_in_metric);

          /*
          const typename geo_metric_type::cellwise_metric geo_cell_met
          = base::m_geo_metric.cellwise_values(cell_idx_in_metric);

          const typename sol_metric_type::cellwise_metric sol_cell_met
          = base::m_sol_metric.cellwise_values(cell_idx_in_metric);

          const typename flux_metric_type::cellwise_metric
          flux_cell_met =
              base::m_flux_metric.cellwise_values(cell_idx_in_metric);
          */

          base::m_const_method_data.CGM = base::m_geo_metric.cellwise_values(cell_idx_in_metric);
          base::m_const_method_data.CSM = base::m_sol_metric.cellwise_values(cell_idx_in_metric);
          base::m_const_method_data.CFM = base::m_flux_metric.cellwise_values(cell_idx_in_metric);

          (*method_data).m_use_external_theta = false;
          if (has_blending_coeff)
          {
            (*method_data).m_use_external_theta = true;
            //(*method_data).base::m_blending_coeff.fill((*base::m_blending_coeff)[solution_elem.idx()]);

            for (Uint v = 0; v < solution_elem.nb_vert(); ++v)
            {
              (*method_data).m_blending_coeff[v] =
                  (*base::m_blending_coeff)[solution_elem.vertex(v)];
            }
          }

          if (has_artif_viscosity)
          {
            /*
            const Real art_visc =
            (*base::m_art_viscosity)[solution_elem.idx()];
            (*method_data).m_art_visc.fill(art_visc);
            */

            const typename sol_metric_type::cellwise_metric stab_coeff_cell_met =
                base::m_stab_coeff_metric.cellwise_values(cell_idx_in_metric);

            const math::DenseConstMatView<Real> stab_coeff_at_qd_pt =
                stab_coeff_cell_met.field_values();

            (*method_data).m_art_visc = stab_coeff_at_qd_pt.col(0);
          }

          if (has_source_term)
          {
            base::m_const_method_data.CSrcM =
                base::m_source_metric.cellwise_values(cell_idx_in_metric);
            m_scheme.compute_adv_reaction_residuals(sol_nodal_values, base::m_const_method_data,
                                                    *method_data);
          }
          else
          {
            m_scheme.compute_adv_residuals(sol_nodal_values, base::m_const_method_data,
                                           *method_data);
          }

          math::DenseDVec<Real> const &elem_update_coeff = (*method_data).m_elem_wave_speed;

          for (Uint n = 0; n < (*method_data).nb_nodes(); ++n)
          {
            // Get a slice of all nodal residuals computer in this
            // element. The 'slice' is only the portion
            // corresponding to node n
            math::DenseVecView<Real> nodal_res = (*method_data).elem_node_res_block(n * NEQ, NEQ);

            const Real inv_nodal_volume = 1. / nodal_dual_volume[solution_elem.vertex(n)];

            for (Uint eq = 0; eq < NEQ; ++eq)
            {
              nodal_res[eq] *= inv_nodal_volume;

              const Uint dof_idx = solution_elem.vertex(n) * NEQ + eq;

              if (!is_Dirichlet_node[dof_idx])
              {
                // rhs(dof_idx, 0) -= elem_res[n][eq];
                base::m_rhs_buffer.push_back(std::tuple<Uint, Real>(dof_idx, -nodal_res[eq]));
                rhs_fill++;
              }
            }

            // time_update.accumulate_nodal_wave_speed(solution_elem.vertex(n),
            // elem_update_coeff[n]);
            base::m_time_update_buffer.push_back(
                std::tuple<Uint, Real>(solution_elem.vertex(n), elem_update_coeff[n]));
          }

          cell_idx_in_metric++;

        } // If the index of solution cell is within range

      } // Loop over all cells of one block

// --------------------------------------------
// Accumulate the buffer values into RHS vector
// --------------------------------------------
#if DO_TIMINGS
      const std::clock_t c_start = std::clock();
      auto t_start               = std::chrono::high_resolution_clock::now();
#endif

      rhs.add_values(0, base::m_rhs_buffer);

#if DO_TIMINGS
      std::clock_t c_end        = std::clock();
      auto t_end                = std::chrono::high_resolution_clock::now();
      const double cpu_duration = (c_end - c_start) / (double)CLOCKS_PER_SEC;
      auto wclock_duration =
          0.001 * std::chrono::duration<double, std::milli>(t_end - t_start).count();
      outfile << "rhs  accumulation cpu/wall " << std::setw(10) << std::setprecision(10)
              << cpu_duration << " " << std::setw(10) << std::setprecision(10) << wclock_duration
              << std::endl;
#endif

      // --------------------------------------------
      // Accumulate nodal wave speeds
      // --------------------------------------------
      time_update.accumulate_wave_speeds(base::m_time_update_buffer);

      if (sol_dof_block_end == dof_group.end())
      {
        process_cells_in_group = false;
      }
      else
      {
        sol_dof_block_begin = sol_dof_block_end;
      }

    } // while (process_cells_in_group)

  } // Lop over all cell groups

#if DO_TIMINGS
  outfile.close();
#endif
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
