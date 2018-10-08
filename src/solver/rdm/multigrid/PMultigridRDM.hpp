#ifndef PDEKIT_Solver_RDM_Multigrid_P_Multigrid_RDM_hpp
#define PDEKIT_Solver_RDM_Multigrid_P_Multigrid_RDM_hpp

#include <cmath>
#include <memory>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

#include "interpolation/FunctionSpace.hpp"
#include "interpolation/L2ProjectionGlobal.hpp"
#include "interpolation/L2ProjectionLocal.hpp"
//#include "interpolation/OperatorResidual.hpp"
#include "interpolation/ResidualRestrictionLocal.hpp"
#include "interpolation/mesh_function/function_ops/MeshFunctionNorm.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "solver/InitialCondition.hpp"
#include "solver/rdm/PGRDMethodImplicit.hpp"
#include "solver/rdm/bc/RDMBCBase.hpp"
#include "solver/rdm/cellsplitters/CellSchemeSelector.hpp"
#include "solver/time/ExplicitTimeStepper.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
class PGRDMethodExplicit;

// ----------------------------------------------------------------------------

namespace multigrid
{

template <typename MeshConfig, typename Physics, typename SchemeTraits>
class PMultigridRDM
{
  private:
  typedef typename internal::CellSchemeSelector<Physics, SchemeTraits>::type cell_scheme_type;

  public:
  /// TYPEDEFS
  typedef mesh::Tria<MeshConfig> mesh_type;
  typedef Physics phys_model;
  typedef RDMBCBase<MeshConfig, phys_model, phys_model::DIM - 1> bc_base_type;

  /// Default constructor
  PMultigridRDM();

  /// Default destructor
  ~PMultigridRDM();

  /// Return the number of levels
  Uint nb_levels() const;

  /// Build a sequence of nested meshes with varying polynomial
  /// order of elements
  void build_mesh_p_sequence(typename mesh::Tria<MeshConfig>::shared_ptr &geo_mesh,
                             const std::string &geo_dof_handler_name,
                             const std::string &cycle_description,
                             const PointSetID topology_type = PointSetID::Equidist);

  /// Define how many smoothing sweeps should be performed at each level
  /// The parameter is a vector of integers of the same length as
  /// 'cycle_description' in the previous method and it defines number of
  /// sweeps on each level
  void define_nb_smoothing_iters_per_level(const std::vector<Uint> &nb_iter_per_level);

  /// Apply initial condition to top-level mesh
  void apply_initial_condition(
      const std::string &domain_name,
      Real (*expr_ptr)(const math::DenseConstVecView<Real> &,
                       const interpolation::VectorMeshFunction<Real>::const_entry_type &,
                       const Uint));

  /// Apply initial condition to top-level mesh
  void apply_initial_condition(const std::string &domain_name,
                               const math::DenseDVec<Real> &state_vector);

  /// Add a new boundary condition to all levels
  /// @param bc_type_name ... type of this boundary condition (e.g.
  /// WeakSubInlet)
  /// @param boundary_name  ... name of the mesh boundary segment where the bc
  /// should be applied
  std::shared_ptr<bc_base_type> add_boundary_condition(const std::string &bc_type_name,
                                                       const std::string &condition_name,
                                                       const std::string &boundary_name,
                                                       const Uint level);

  /// Set the expression that defines the boundary condition values
  /// Set an analytical expression if needed for the bc
  /// @param bc_type_name ... type of this boundary condition (e.g.
  /// WeakDirichletBC)
  /// @param expr_ptr ... expression which defines the values of this boundary
  /// condition
  /*
  void set_expression(const std::string &bc_type_name,
                      const Real (*expr_ptr)(const math::VectorBlock<Real> &,
                                             const typename
  interpolation::FunctionSpace< MeshConfig>::function_t::column_type &, const
  Uint));
  */

  /// Perform the initial iterations to prepare the full multigrid algorithm
  void mg_startup(const Uint nb_iter = 10);

  /// Make one multigrid cycle
  void make_cycle();

  /// Compute the residual norm on given level
  void compute_residual_norm(const Uint level, math::DenseDVec<Real> &norm) const;

  /// Write all the underlying meshes to files
  void write_sequence_to_files() const;

  void smooth_single_level(const Uint level, const Uint nb_iter);

  private:
  /// TYPEDEFS
  typedef std::shared_ptr<mesh::Tria<MeshConfig>> mesh_shared_ptr;
  typedef interpolation::FunctionSpace<MeshConfig> function_space;
  typedef interpolation::ScalarMeshFunction<Real> scal_mesh_function;
  typedef interpolation::VectorMeshFunction<Real> vect_mesh_function;

  typedef typename vect_mesh_function::entry_type node_value_type;
  typedef typename vect_mesh_function::const_entry_type const_node_value_type;

  typedef interpolation::L2ProjectionGlobal sol_restriction_op_type;
  typedef interpolation::L2ProjectionGlobal defect_prolongation_op_type;
  typedef interpolation::L2ProjectionGlobal residual_restriction_op_type;
  // typedef interpolation::ResidualRestrictionLocal
  // residual_restriction_op_type;

  /// TYPES

  struct MGLevelData
  {
    MGLevelData()
    {
    }

    ~MGLevelData()
    {
    }

    typename mesh_type::shared_ptr mesh_ptr;
    common::PtrHandle<typename mesh_type::dof_storage_type> geo_dofs;
    common::PtrHandle<typename mesh_type::dof_storage_type> sol_dofs;

    Uint poly_order;
    typename vect_mesh_function::ptr u;
    typename vect_mesh_function::ptr u0;
    typename vect_mesh_function::ptr rdm_res;
    typename vect_mesh_function::ptr res;
    typename vect_mesh_function::ptr tmp_work;
    // PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits> method;
    std::unique_ptr<RDMethod<MeshConfig, Physics>> method;
    // typename interpolation::OperatorResidual<MeshConfig, Physics>
    // op_residual;

    // ExplicitEuler m_time_stepper;
    time::RK3TVD m_time_stepper;
    // solver::L2Projection m_sol_restriction;
    // solver::L2Projection m_defect_prolongation;
  };

  struct PolyRule
  {
    // Return: a tuple consisting of
    // 1) Shape function tag of source shape function (which we're
    // projecting FROM) 2) Shape function tag of target polynomial (i.e
    // space we're projecting TO)
    std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> sf_tags(const mesh::PointSetTag &geo_tag,
                                                         const mesh::PointSetTag &sol_tag_src,
                                                         const mesh::PointSetTag &sol_tag_tgt) const
    {
      const mesh::sf::SFTag src_poly_space_tag(sol_tag_src.elem_shape(), SFunc::Lagrange,
                                               sol_tag_src.poly_order(), ModalBasis::Modal);

      const mesh::sf::SFTag tgt_poly_space_tag(sol_tag_tgt.elem_shape(), SFunc::Lagrange,
                                               sol_tag_tgt.poly_order(), ModalBasis::Modal);

      const std::tuple<mesh::sf::SFTag, mesh::sf::SFTag> result =
          std::make_tuple(src_poly_space_tag, tgt_poly_space_tag);
      return result;
    }

    // Return: a tuple consisting of
    // 1) Quadrature order to be used in order to construct LHS matrix
    // operator 2) Quadrature order to be used in order to construct RHS
    std::tuple<Uint, Uint> quadrature_orders(const mesh::PointSetTag &geo_tag,
                                             const mesh::PointSetTag &sol_tag_src,
                                             const mesh::PointSetTag &sol_tag_tgt) const
    {
      const Uint src_quad_order =
          std::max(sol_tag_src.poly_order() + sol_tag_tgt.poly_order(), geo_tag.poly_order());
      const Uint tgt_quad_order = std::max(2 * sol_tag_tgt.poly_order(), geo_tag.poly_order());

      const std::tuple<Uint, Uint> result = std::make_tuple(src_quad_order, tgt_quad_order);
      return result;
    }
  };

  typedef std::pair<Uint, Uint> transfer_op_map_key;
  struct UnsignedIntPairKeyHasher
  {
    inline std::size_t operator()(const transfer_op_map_key &key) const
    {
      return key.first ^ key.second;
    }
  };

  typedef std::unordered_map<transfer_op_map_key, std::unique_ptr<sol_restriction_op_type>,
                             UnsignedIntPairKeyHasher>
      sol_restrict_op_map;

  typedef std::unordered_map<transfer_op_map_key, std::unique_ptr<residual_restriction_op_type>,
                             UnsignedIntPairKeyHasher>
      res_restrict_op_map;

  typedef std::unordered_map<transfer_op_map_key, std::unique_ptr<defect_prolongation_op_type>,
                             UnsignedIntPairKeyHasher>
      defect_prolongation_op_map;

  /// METHODS
  /// Prepare the function spaces - multigrid setup
  void build_solvers();

  /// Prepare transfer operators
  void configure_transfer_operators();

  /// Clear member data
  void clear_mesh_sequence();

  /// Run the smoother
  void smooth(MGLevelData &mg_data, const Uint nb_iter);

  /// Move one level down in the multigrid cycle
  /// (i.e. go to coarser level)
  void move_down(MGLevelData &mg_data_fine, MGLevelData &mg_data_coarse,
                 sol_restriction_op_type &sol_restriction_operator,
                 residual_restriction_op_type &res_restriction_operator);

  /// Move one level up in the multigrid cycle
  /// (i.e. go to finer level)
  void move_up(MGLevelData &mg_data_coarse, MGLevelData &mg_data_fine,
               defect_prolongation_op_type &defect_prolongation_operator);

  /// DATA
  /// Description of the multigrid cycle
  /// This vector holds values of polynomial orders
  /// as they go in one multigrid cycle
  std::vector<Uint> m_cycle_description;

  /// Number of smoothing iterations on each level
  std::vector<Uint> m_nb_smoothing_iters_per_level;

  /// Cache the minimum and maximum polynomial order
  Uint min_p;
  Uint max_p;

  /// Solution restriction operators
  sol_restrict_op_map m_sol_restriction_ops;

  res_restrict_op_map m_res_restriction_ops;

  /// Defect prolongation operators
  defect_prolongation_op_map m_defect_prolongation_ops;

  /// Sequence of meshes for multigrid
  std::unordered_map<Uint, std::unique_ptr<MGLevelData>> m_mesh_sequence;

  /// Say if the output should be verbose
  bool m_verbose;

  std::string m_geo_dof_handler_name;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
PMultigridRDM<MeshConfig, Physics, SchemeTraits>::PMultigridRDM()
    : min_p(0), max_p(0), m_verbose(false), m_geo_dof_handler_name("")
{
  m_cycle_description.resize(0);
  m_nb_smoothing_iters_per_level.resize(0);
  m_mesh_sequence.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
PMultigridRDM<MeshConfig, Physics, SchemeTraits>::~PMultigridRDM()
{
  clear_mesh_sequence();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
Uint PMultigridRDM<MeshConfig, Physics, SchemeTraits>::nb_levels() const
{
  return m_mesh_sequence.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::build_mesh_p_sequence(
    typename mesh::Tria<MeshConfig>::shared_ptr &geo_mesh, const std::string &geo_dof_handler_name,
    const std::string &cycle_description, const PointSetID topology_type)
{
  m_geo_dof_handler_name = geo_dof_handler_name;

  std::vector<std::string> strs;
  boost::split(strs, cycle_description, boost::is_any_of("-"));

  m_cycle_description.resize(0);

  for (Uint s = 0; s < strs.size(); ++s)
  {
    for (Uint p = 0; p < PolyOrder::nb_instances(); ++p)
    {
      if (PolyOrder::name(p) == strs[s])
      {
        m_cycle_description.push_back(PolyOrder::value(p));
        break;
      }
    }
  }

  std::cout << "Cycle configured with the following polynomial orders: " << std::endl;
  for (Uint i = 0; i < m_cycle_description.size(); ++i)
  {
    std::cout << m_cycle_description[i] << " ";
  }
  std::cout << std::endl;

  min_p = m_cycle_description[0];
  max_p = m_cycle_description[0];

  for (Uint i = 1; i < m_cycle_description.size(); ++i)
  {
    min_p = std::min(min_p, m_cycle_description[i]);
    max_p = std::max(max_p, m_cycle_description[i]);
  }

  std::cout << "Minimum p: " << min_p << std::endl;
  std::cout << "Maximum p: " << max_p << std::endl;

  // ----------------
  // Configure meshes
  // ----------------
  m_mesh_sequence.clear();

  common::PtrHandle<typename mesh_type::dof_storage_type> geo_dofs =
      geo_mesh->dof_storage(geo_dof_handler_name);

  typedef typename mesh_type::dof_storage_type dof_storage_type;

  for (Uint l = 0; l < m_cycle_description.size(); ++l)
  {
    typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator data_it =
        m_mesh_sequence.find(m_cycle_description[l]);

    // If mesh data object has not been found, insert a new one
    if (data_it == m_mesh_sequence.end())
    {
      // In case the required order is P0, we set the mesh to use linear
      // elements, hence std::max(1, m_cycle_description[l])
      common::PtrHandle<typename mesh_type::dof_storage_type> sol_dofs =
          geo_mesh->create_dof_storage("sol_dofs_p" +
                                       common::StringUtils::to_string(m_cycle_description[l]));
      dof_storage_type::clone_continuous(*geo_mesh, *geo_dofs, *sol_dofs,
                                         std::max(1u, m_cycle_description[l]), topology_type);

      // This is C++14
      // std::unique_ptr<MGLevelData> mg_data =
      // std::make_unique<MGLevelData>();
      std::unique_ptr<MGLevelData> mg_data(new MGLevelData());

      mg_data->mesh_ptr   = geo_mesh;
      mg_data->geo_dofs   = geo_dofs;
      mg_data->sol_dofs   = sol_dofs;
      mg_data->poly_order = m_cycle_description[l];

      m_mesh_sequence.insert(std::pair<Uint, std::unique_ptr<MGLevelData>>(m_cycle_description[l],
                                                                           std::move(mg_data)));
    }
  }

  /*
  std::cout << "Data inserted to mesh sequence:" << std::endl;
  for (typename std::unordered_map<Uint,
  std::unique_ptr<MGLevelData>>::const_iterator data_iter =
           m_mesh_sequence.cbegin();
       data_iter != m_mesh_sequence.cend(); ++data_iter)
  {
    std::cout << "[" << data_iter->first << "]" << std::endl;
  }
  */

  // Configure function spaces and data for the solver
  build_solvers();

  // Prepare transfer operators (solution restriction, residual restriction
  // and defect prolongation
  configure_transfer_operators();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::define_nb_smoothing_iters_per_level(
    const std::vector<Uint> &nb_iter_per_level)
{
  m_nb_smoothing_iters_per_level.resize(nb_iter_per_level.size());
  m_nb_smoothing_iters_per_level = nb_iter_per_level;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::apply_initial_condition(
    const std::string &domain_name,
    Real (*expr_ptr)(const math::DenseConstVecView<Real> &,
                     const interpolation::VectorMeshFunction<Real>::const_entry_type &, const Uint))
{
  for (typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator data_iter =
           m_mesh_sequence.cbegin();
       data_iter != m_mesh_sequence.cend(); ++data_iter)
  {
    MGLevelData &mg_data = (*data_iter->second);

    solver::InitialCondition<MeshConfig> ic("L0_MG_initial_condition", mg_data.mesh_ptr);

    ic.set_domain(MeshConfig::TDIM, domain_name);
    ic.set_expression(expr_ptr);
    ic.apply_function((*mg_data.sol_dofs).name(), *mg_data.u);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::apply_initial_condition(
    const std::string &domain_name, const math::DenseDVec<Real> &state_vector)
{
  for (typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator data_iter =
           m_mesh_sequence.cbegin();
       data_iter != m_mesh_sequence.cend(); ++data_iter)
  {
    MGLevelData &mg_data = (*data_iter->second);

    solver::InitialCondition<MeshConfig> ic_sol("L0_MG_initial_condition", mg_data.mesh_ptr);

    ic_sol.set_domain(MeshConfig::TDIM, domain_name);
    ic_sol.apply_values((*mg_data.sol_dofs).name(), state_vector, *mg_data.u);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
std::shared_ptr<typename PMultigridRDM<MeshConfig, Physics, SchemeTraits>::bc_base_type> PMultigridRDM<
    MeshConfig, Physics, SchemeTraits>::add_boundary_condition(const std::string &bc_type_name,
                                                               const std::string &condition_name,
                                                               const std::string &boundary_name,
                                                               const Uint level)
{
  typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator data_it =
      m_mesh_sequence.find(level);

  MGLevelData &mg_data = *(data_it->second);

  std::shared_ptr<bc_base_type> new_bc =
      mg_data.method->add_boundary_condition(bc_type_name, condition_name, boundary_name);

  return new_bc;
}

// ----------------------------------------------------------------------------

/*
template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PGRDMultigrid<MeshConfig, RDScheme>::set_expression(
    const std::string &bc_type_name,
    const Real (*expr_ptr)(
        const math::VectorBlock<Real> &,
        const typename
interpolation::FunctionSpace<MeshConfig>::function_t::column_type &, const
Uint))
{
}
*/

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::mg_startup(const Uint nb_iter)
{

  // Use a set to store all polynomial orders
  // They should be stored
  // a) uniquely (no p is repeated twice)
  // b) in increasing order
  std::set<Uint> cycle_poly_orders;

  for (Uint i = 0; i < m_cycle_description.size(); ++i)
  {
    cycle_poly_orders.insert(m_cycle_description[i]);
  }

  // Now we put them in a vector
  std::vector<Uint> cycle_poly_order_vec(0);
  for (std::set<Uint>::const_iterator p_iter = cycle_poly_orders.cbegin();
       p_iter != cycle_poly_orders.cend(); ++p_iter)
  {
    cycle_poly_order_vec.push_back(*p_iter);
  }

  typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator
      mg_data_iter_coarse;
  typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator mg_data_iter_fine;

  for (Uint i = 0; i < (cycle_poly_order_vec.size() - 1); ++i)
  {
    const Uint P_in  = cycle_poly_order_vec[i];
    const Uint P_out = cycle_poly_order_vec[i + 1];

    mg_data_iter_coarse = m_mesh_sequence.find(P_in);
    mg_data_iter_fine   = m_mesh_sequence.find(P_out);

    MGLevelData &mg_data_coarse = (*mg_data_iter_coarse->second);
    MGLevelData &mg_data_fine   = (*mg_data_iter_fine->second);

    smooth(mg_data_coarse, nb_iter);

    typename mesh_type::shared_ptr fine_sol_mesh   = mg_data_fine.sol_mesh;
    typename mesh_type::shared_ptr coarse_sol_mesh = mg_data_coarse.sol_mesh;

    const transfer_op_map_key key = transfer_op_map_key(P_in, P_out);

    const typename defect_prolongation_op_map::const_iterator defect_prolongation_it =
        m_defect_prolongation_ops.find(key);
    defect_prolongation_op_type &defect_prolongation = *(defect_prolongation_it->second);
    defect_prolongation.project(*coarse_sol_mesh, *fine_sol_mesh, *mg_data_coarse.u,
                                *mg_data_fine.u, m_verbose);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::make_cycle()
{
  for (Uint i = 0; i < (m_cycle_description.size() - 1); ++i)
  {
    const Uint P_in  = m_cycle_description[i];
    const Uint P_out = m_cycle_description[i + 1];

    const typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator
        data_iter_in = m_mesh_sequence.find(P_in);

    const typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator
        data_iter_out = m_mesh_sequence.find(P_out);

    MGLevelData &data_in  = *(data_iter_in->second);
    MGLevelData &data_out = *(data_iter_out->second);

    const transfer_op_map_key key = transfer_op_map_key(P_in, P_out);

    // First perform smoothing
    smooth(data_in, m_nb_smoothing_iters_per_level[i]);

    // If the initial polynomial order is higher than the final, we're
    // moving down in the mesh hierarchy

    if (P_in > P_out)
    {
      // Number of iterations on fine level before solution
      // and residual restriction
      // const Uint pre_smooth_fine_iters = 35; // 100

      // Pre-smooth data on fine level
      // smooth(data_in, pre_smooth_fine_iters);

      // Get solution restriction operator
      const typename sol_restrict_op_map::const_iterator sol_restrict_it =
          m_sol_restriction_ops.find(key);
      sol_restriction_op_type &sol_restriction = *(sol_restrict_it->second);

      // Get residual restriction operator
      const typename res_restrict_op_map::const_iterator res_restrict_it =
          m_res_restriction_ops.find(key);
      residual_restriction_op_type &res_restriction = (*res_restrict_it->second);

      // Perform the restriction
      move_down(data_in, data_out, sol_restriction, res_restriction);
    }
    else
    {
      // Number of smoother iterations on coarse level
      // const Uint smooth_coarse_iters = 35; // 100

      // Perform iterations on the coarse level
      // smooth(data_in, smooth_coarse_iters);

      // Prolongate the defect
      const typename defect_prolongation_op_map::const_iterator defect_prolongation_it =
          m_defect_prolongation_ops.find(key);
      defect_prolongation_op_type &defect_prolongation = *(defect_prolongation_it->second);
      move_up(data_in, data_out, defect_prolongation);

      // Number of post-smooth iterations on fine level
      // const Uint post_smooth_fine_iters = 35; // 100;

      // Post-smooth on the fine level. The source term here should be
      // ZERO! smooth(data_out, post_smooth_fine_iters);
    }
  } // Loop over cycle description
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::compute_residual_norm(
    const Uint level, math::DenseDVec<Real> &norm) const
{
  const typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator data_iter =
      m_mesh_sequence.find(level);

  if (data_iter == m_mesh_sequence.end())
  {
    std::cerr << "PGRDMultigrid::compute_residual_norm: mesh level " << level << " does not exist."
              << std::endl;
    return;
  }

  const MGLevelData &mg_data = *(data_iter->second);
  interpolation::norm_L2(*mg_data.rdm_res, norm);
  // interpolation::norm_L2((*mg_data.res) - (*mg_data.rdm_res), norm);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::write_sequence_to_files() const
{
  mesh::gmsh::GmshWriter meshwriter;

  for (typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator data_it =
           m_mesh_sequence.cbegin();
       data_it != m_mesh_sequence.cend(); ++data_it)
  {
    MGLevelData &mg_data = (*data_it->second);

    mesh_type &mg_mesh = *mg_data.mesh_ptr;

    const std::string level_name = common::StringUtils::to_string(data_it->first);

    const std::string geo_outfilename = mg_mesh.name() + "_geo_L" + level_name + ".msh";
    const std::string sol_outfilename = mg_mesh.name() + "_sol_L" + level_name + ".msh";

    meshwriter.write_mesh_to_file(mg_mesh, m_geo_dof_handler_name, geo_outfilename);
    meshwriter.write_mesh_to_file(mg_mesh, (*mg_data.sol_dofs).name(), sol_outfilename);

    meshwriter.append_nodal_function_to_file(mg_mesh, sol_outfilename, *mg_data.u, "solution");
    meshwriter.append_nodal_function_to_file(mg_mesh, sol_outfilename, *mg_data.rdm_res, "rdm_res");
    meshwriter.append_nodal_function_to_file(mg_mesh, sol_outfilename, *mg_data.res, "mg_res");
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::smooth_single_level(const Uint level,
                                                                           const Uint nb_iter)
{
  const typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator data_iter =
      m_mesh_sequence.find(level);

  if (data_iter == m_mesh_sequence.end())
  {
    std::cerr << "PGRDMultigrid::smooth_single_level: mesh level " << level << " does not exist."
              << std::endl;
    return;
  }

  MGLevelData &mg_data = *(data_iter->second);

  smooth(mg_data, nb_iter);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::build_solvers()
{
  for (typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator data_it =
           m_mesh_sequence.cbegin();
       data_it != m_mesh_sequence.cend(); ++data_it)
  {
    // If the required order is P0, then we set the approximation space to
    // P1, and we will use N scheme If the required order is different from
    // P0, keep it as it is
    const Uint element_poly_order = (data_it->first == P0) ? P1 : data_it->first;

    // This is 'scheme' order: P
    const Uint scheme_order = data_it->first;

    MGLevelData &mg_data = (*data_it->second);

    // const PolyOrderID quadrature_order = P4;
    const Uint quadrature_order = element_poly_order * 2;

    if (scheme_order == P0)
    {
      PGRDMethodExplicit<MeshConfig, Physics, PGN> *method_ptr =
          new PGRDMethodExplicit<MeshConfig, Physics, PGN>();
      mg_data.method =
          std::move(std::unique_ptr<PGRDMethodExplicit<MeshConfig, Physics, PGN>>(method_ptr));

      /*
      PGRDMethodImplicit<MeshConfig, Physics, PGN> *method_ptr =
          new PGRDMethodImplicit<MeshConfig, Physics, PGN>();
      mg_data.method =
          std::move(std::unique_ptr<PGRDMethodImplicit<MeshConfig,
      Physics, PGN>>(method_ptr));
      */
    }
    else
    {
      PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits> *method_ptr =
          new PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>();
      mg_data.method = std::move(
          std::unique_ptr<PGRDMethodExplicit<MeshConfig, Physics, SchemeTraits>>(method_ptr));

      /*
      PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits> *method_ptr =
          new PGRDMethodImplicit<MeshConfig, Physics, SchemeTraits>();
      mg_data.method = std::move(
          std::unique_ptr<PGRDMethodImplicit<MeshConfig, Physics,
      SchemeTraits>>(method_ptr));
      */
    }

    mg_data.method->configure_mesh_data(mg_data.mesh_ptr, mg_data.geo_dofs, mg_data.sol_dofs);

    mg_data.method->initialize_work_data(SFunc::Lagrange, PointSetID::Gauss, quadrature_order, 1,
                                         100);

    // mg_data.op_residual.setup(*mg_data.geo_dofs, *mg_data.sol_dofs);

    mg_data.u        = std::make_shared<vect_mesh_function>("", "u");
    mg_data.u0       = std::make_shared<vect_mesh_function>("", "u0");
    mg_data.rdm_res  = std::make_shared<vect_mesh_function>("", "rdm_res");
    mg_data.res      = std::make_shared<vect_mesh_function>("", "rhs");
    mg_data.tmp_work = std::make_shared<vect_mesh_function>("", "tmp_work");

    const Uint nb_nodes = (*mg_data.sol_dofs).nb_nodes();
    mg_data.u->resize(phys_model::NEQ, nb_nodes);
    mg_data.u0->resize(phys_model::NEQ, nb_nodes);
    mg_data.rdm_res->resize(phys_model::NEQ, nb_nodes);
    mg_data.res->resize(phys_model::NEQ, nb_nodes);
    mg_data.tmp_work->resize(phys_model::NEQ, nb_nodes);

    mg_data.u->fill(0.0);
    mg_data.u0->fill(0.0);
    mg_data.rdm_res->fill(0.0);
    mg_data.res->fill(0.0);
    mg_data.tmp_work->fill(0.0);

    /*
    mg_data.u = mg_data.scheme.get_solution();
    mg_data.rdm_nodal_res = mg_data.scheme.get_residuals();
    mg_data.update_coeff = mg_data.scheme.get_update_coeffs();
    */

    mg_data.method->set_vec_function(SolverVecFn::solution, mg_data.u);
    mg_data.method->set_vec_function(SolverVecFn::residuals, mg_data.rdm_res);
    if (element_poly_order != max_p)
    {
      mg_data.method->set_vec_function(SolverVecFn::sources, mg_data.res);
    }

    if (mg_data.method->is_explicit())
    {
      mg_data.m_time_stepper.setup(phys_model::NEQ, nb_nodes);
    }

  } // Loop over all levels in hierarchy
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::configure_transfer_operators()
{
  // Configure transfer operators
  for (Uint i = 0; i < (m_cycle_description.size() - 1); ++i)
  {
    const Uint P_in  = m_cycle_description[i];
    const Uint P_out = m_cycle_description[i + 1];

    const typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator
        data_iter_in = m_mesh_sequence.find(P_in);
    common::PtrHandle<typename mesh_type::dof_storage_type> sol_dofs_in =
        data_iter_in->second->sol_dofs;

    const typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::const_iterator
        data_iter_out = m_mesh_sequence.find(P_out);
    common::PtrHandle<typename mesh_type::dof_storage_type> sol_dofs_out =
        data_iter_out->second->sol_dofs;

    const transfer_op_map_key key = transfer_op_map_key(P_in, P_out);

    // If the 'starting' polynomial order is bigger than the 'final'
    // polyorder, we need to define solution restriction
    if (P_in > P_out)
    {
      typename sol_restrict_op_map::const_iterator it = m_sol_restriction_ops.find(key);
      if (it == m_sol_restriction_ops.end())
      {
        const typename mesh_type::shared_ptr mesh_in = data_iter_in->second->mesh_ptr;

        common::PtrHandle<typename mesh_type::dof_storage_type> geo_dofs_in =
            data_iter_in->second->geo_dofs;

        std::unique_ptr<sol_restriction_op_type> sol_restriction(new sol_restriction_op_type());
        sol_restriction->build_projection_operator<MeshConfig, PolyRule>(*mesh_in, *sol_dofs_in,
                                                                         *sol_dofs_out, PolyRule{});

        m_sol_restriction_ops.insert(
            std::pair<transfer_op_map_key, std::unique_ptr<sol_restriction_op_type>>(
                key, std::move(sol_restriction)));

        std::unique_ptr<residual_restriction_op_type> res_restriction(
            new residual_restriction_op_type());
        res_restriction->build_projection_operator<MeshConfig, PolyRule>(*mesh_in, *sol_dofs_in,
                                                                         *sol_dofs_out, PolyRule{});

        m_res_restriction_ops.insert(
            std::pair<transfer_op_map_key, std::unique_ptr<residual_restriction_op_type>>(
                key, std::move(res_restriction)));
      }
    }
    else
    {
      typename defect_prolongation_op_map::const_iterator it = m_defect_prolongation_ops.find(key);
      if (it == m_defect_prolongation_ops.end())
      {
        const typename mesh_type::shared_ptr mesh_out = data_iter_out->second->mesh_ptr;

        common::PtrHandle<typename mesh_type::dof_storage_type> geo_dofs_out =
            data_iter_out->second->geo_dofs;

        std::unique_ptr<defect_prolongation_op_type> defect_prolongation(
            new defect_prolongation_op_type());
        defect_prolongation->build_projection_operator<MeshConfig, PolyRule>(
            *mesh_out, *sol_dofs_in, *sol_dofs_out, PolyRule{});

        m_defect_prolongation_ops.insert(
            std::pair<transfer_op_map_key, std::unique_ptr<defect_prolongation_op_type>>(
                key, std::move(defect_prolongation)));
      }
    }
  } // Loop over the vector describing one multigrid cycle

  /*
  std::cout << "Inserted solution restriction operators:" << std::endl;
  for (typename transfer_op_map::const_iterator it =
  m_sol_restriction_ops.begin(); it != m_sol_restriction_ops.end(); ++it)
  {
    std::cout << "P" << it->first.first << " -> P" << it->first.second <<
  std::endl;
  }

  std::cout << "Inserted defect prolongation operators:" << std::endl;
  for (typename transfer_op_map::const_iterator it =
  m_defect_prolongation_ops.begin(); it != m_defect_prolongation_ops.end();
  ++it)
  {
    std::cout << "P" << it->first.first << " -> P" << it->first.second <<
  std::endl;
  }
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::clear_mesh_sequence()
{
  for (typename std::unordered_map<Uint, std::unique_ptr<MGLevelData>>::iterator data_it =
           m_mesh_sequence.begin();
       data_it != m_mesh_sequence.end(); ++data_it)
  {
    data_it->second.release();
  }
  m_mesh_sequence.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::smooth(MGLevelData &mg_data,
                                                              const Uint nb_iter)
{
  // If the method is explicit, use RK-type time stepper
  if (mg_data.method->is_explicit())
  {
    for (Uint iter = 0; iter < nb_iter; ++iter)
    {
      mg_data.m_time_stepper.advance_in_time(*mg_data.method, 0.5);
    }
  }
  else
  {
    // Implicit iteration
    mg_data.method->assemble_lhs_and_rhs(5.0);
    mg_data.method->solve();
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::move_down(
    MGLevelData &mg_data_fine, MGLevelData &mg_data_coarse,
    sol_restriction_op_type &sol_restriction_operator,
    residual_restriction_op_type &res_restriction_operator)
{
  const typename mesh_type::shared_ptr mesh_fine                          = mg_data_fine.mesh_ptr;
  common::PtrHandle<typename mesh_type::dof_storage_type> fine_geo_dofs   = mg_data_fine.geo_dofs;
  common::PtrHandle<typename mesh_type::dof_storage_type> fine_sol_dofs   = mg_data_fine.sol_dofs;
  common::PtrHandle<typename mesh_type::dof_storage_type> coarse_sol_dofs = mg_data_coarse.sol_dofs;

  // Data on fine level
  vect_mesh_function &res_fine      = *mg_data_fine.res;
  vect_mesh_function &rdm_res_fine  = *mg_data_fine.rdm_res;
  vect_mesh_function &u_fine        = *mg_data_fine.u;
  vect_mesh_function &tmp_work_fine = *mg_data_fine.tmp_work;

  // Data on coarse level
  vect_mesh_function &res_coarse      = *mg_data_coarse.res;
  vect_mesh_function &rdm_res_coarse  = *mg_data_coarse.rdm_res;
  vect_mesh_function &u_coarse        = *mg_data_coarse.u;
  vect_mesh_function &u0_coarse       = *mg_data_coarse.u0;
  vect_mesh_function &tmp_work_coarse = *mg_data_coarse.tmp_work;

  // ------------------------------------
  // Restrict solution
  // ------------------------------------
  if (m_verbose)
  {
    std::cout << "  PGRDMultigrid:: restricting solution ... " << std::endl;
  }

  u_coarse.fill(0.0);
  sol_restriction_operator.project<MeshConfig, PolyRule>(
      *mesh_fine, *fine_sol_dofs, *coarse_sol_dofs, PolyRule{}, u_fine, u_coarse, m_verbose);

  // ------------------------------------
  // Restrict residuals
  // ------------------------------------
  if (m_verbose)
  {
    std::cout << "  PGRDMultigrid:: restricting residuals ... " << std::endl;
  }

  tmp_work_fine = res_fine - rdm_res_fine;
  tmp_work_coarse.fill(0.0);

  ///@todo consider using residual restriction operator which is different
  /// from
  ///      solution restriction operator
  // tmp_work_coarse now contains restricted fine rhs residual rhs^{k} =>
  // rhs^{k-1}
  res_restriction_operator.project<MeshConfig, PolyRule>(*mesh_fine, *fine_sol_dofs,
                                                         *coarse_sol_dofs, PolyRule{},
                                                         tmp_work_fine, tmp_work_coarse, m_verbose);

  // ------------------------------------
  // 3) Compute the error (defect)
  // ------------------------------------
  if (m_verbose)
  {
    std::cout << "  PGRDMultigrid:: defect evaluation on coarse level ... " << std::endl;
  }

  // The first component in the source term on coarse level is the
  // rdm_nodal_res obtained by performing one iteration on the restricted
  // solution from fine mesh
  u0_coarse = u_coarse;

  // smooth(mg_data_coarse, 1);

  // Run the scheme

  // This step stores A^{k-1} (v^{k-1}) in rdm_res_coarse
  mg_data_coarse.method->assemble_rhs();

  // The new rhs is r^{k-1} + A^{k-1} (v^{k-1}_0)
  res_coarse = rdm_res_coarse + tmp_work_coarse;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, typename SchemeTraits>
void PMultigridRDM<MeshConfig, Physics, SchemeTraits>::move_up(
    MGLevelData &mg_data_coarse, MGLevelData &mg_data_fine,
    defect_prolongation_op_type &defect_prolongation_operator)
{
  const typename mesh_type::shared_ptr mesh_fine                          = mg_data_fine.mesh_ptr;
  common::PtrHandle<typename mesh_type::dof_storage_type> coarse_sol_dofs = mg_data_coarse.sol_dofs;
  common::PtrHandle<typename mesh_type::dof_storage_type> fine_geo_dofs   = mg_data_fine.geo_dofs;
  common::PtrHandle<typename mesh_type::dof_storage_type> fine_sol_dofs   = mg_data_fine.sol_dofs;

  // Data on coarse level
  vect_mesh_function &tmp_work_coarse = *mg_data_coarse.tmp_work;
  // vect_mesh_function &rdm_res_coarse = *mg_data_coarse.rdm_res;
  vect_mesh_function &u_coarse  = *mg_data_coarse.u;
  vect_mesh_function &u0_coarse = *mg_data_coarse.u0;

  // Data on fine level
  vect_mesh_function &tmp_work_fine = *mg_data_fine.tmp_work;
  // vect_mesh_function &rdm_res_fine = *mg_data_fine.rdm_res;
  vect_mesh_function &u_fine = *mg_data_fine.u;
  // vect_mesh_function &u0_fine = *mg_data_fine.u0;

  // ---------------------------------------
  // This is a second part of step 3
  // ---------------------------------------

  // Evaluate the defect after smoothing on coarse level:
  tmp_work_coarse = u_coarse - u0_coarse;

  // ---------------------------------------
  // 4) Prolongate the defect onto fine mesh
  // ---------------------------------------
  if (m_verbose)
  {
    std::cout << "  PGRDMultigrid:: defect prolongation ... " << std::endl;
  }

  defect_prolongation_operator.project<MeshConfig, PolyRule>(
      *mesh_fine, *coarse_sol_dofs, *fine_sol_dofs, PolyRule{}, tmp_work_coarse, tmp_work_fine);

  /*
  // Write prolongated solution and defect to file
  mesh::gmsh::GmshWriter meshwriter;
  const std::string level_name =
  common::StringUtils::to_string(mg_data_fine.poly_order);

  const std::string outfilename = "correction_L_" + level_name + ".msh";

  meshwriter.write_mesh_to_file(*fine_sol_mesh, outfilename);
  meshwriter.append_nodal_function_to_file(*fine_sol_mesh, outfilename,
  u_fine, "u_before");
  meshwriter.append_nodal_function_to_file(*fine_sol_mesh, outfilename,
  tmp_work_fine, "correction");
  */

  u_fine += tmp_work_fine;

  /*
  meshwriter.append_nodal_function_to_file(*fine_sol_mesh, outfilename,
  u_fine, "u_after");
  */

  // ---------------------------------------
  // 5) Smooth data on fine mesh
  // ---------------------------------------
  if (m_verbose)
  {
    std::cout << "  PGRDMultigrid:: post-smoothing on fine level ... " << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace multigrid

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
