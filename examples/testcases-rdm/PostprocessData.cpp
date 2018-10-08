#include <ctime>
#include <iostream>

#include "common/DataMap.hpp"
#include "common/PDEKit.hpp"
#include "mesh/io/MeshManipulator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "physics/euler/Euler2DCons.hpp"
#include "physics/euler/Euler3DCons.hpp"
#include "solver/FEMetric.hpp"
#include "solver/FEMetricData.hpp"
#include "solver/InitialCondition.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::interpolation;
using namespace pdekit::solver;

typedef Cart3D MeshConfig;
typedef Tria<MeshConfig> MeshType;
typedef physics::Euler3DCons PhysicsType;

template <Uint GEODIM, Uint TDIM>
class SurfaceJacobian
{
};

template <>
class SurfaceJacobian<_2D, _1D>
{
  public:
  // typedef std::shared_ptr<FEMetricTerm<_2D,_1D> > geo_data_ptr;

  static void compute(const solver::FEMetricData<_2D, _1D> &geo_data,
                      math::DenseDMat<Real> &normals, math::DenseDVec<Real> &jacobians)
  {
    const math::DenseDMat<Real> &dXqd = geo_data.m_dxqd[XI0];

    for (Uint q = 0; q < geo_data.nb_qd_pts; ++q)
    {
      jacobians[q] = std::sqrt(dXqd(q, X) * dXqd(q, X) + dXqd(q, Y) * dXqd(q, Y));

      normals(q, X) = dXqd(q, Y) / jacobians[q];
      normals(q, Y) = -dXqd(q, X) / jacobians[q];
    }
  }
};

template <>
class SurfaceJacobian<_3D, _2D>
{
  public:
  // typedef std::shared_ptr<FEMetricTerm<_2D,_1D> > geo_data_ptr;

  static void compute(const solver::FEMetricData<_3D, _2D> &geo_data,
                      math::DenseDMat<Real> &normals, math::DenseDVec<Real> &jacobians)
  {
    // Derivatives with respect to xi0
    const math::DenseDMat<Real> &dXI0qd = geo_data.m_dxqd[XI0];
    // Derivatives with respect to xi1
    const math::DenseDMat<Real> &dXI1qd = geo_data.m_dxqd[XI1];

    for (Uint q = 0; q < geo_data.nb_qd_pts; ++q)
    {
      normals(q, X) = dXI0qd(q, Y) * dXI1qd(q, Z) - dXI0qd(q, Z) * dXI1qd(q, Y);
      normals(q, Y) = dXI0qd(q, Z) * dXI1qd(q, X) - dXI0qd(q, X) * dXI1qd(q, Z);
      normals(q, Z) = dXI0qd(q, X) * dXI1qd(q, Y) - dXI0qd(q, Y) * dXI1qd(q, X);

      jacobians[q] = std::sqrt(normals(q, X) * normals(q, X) + normals(q, Y) * normals(q, Y) +
                               normals(q, Z) * normals(q, Z));

      normals(q, X) /= jacobians[q];
      normals(q, Y) /= jacobians[q];
      normals(q, Z) /= jacobians[q];
    }
  }
};

int main()
{
#if 0
  MeshType::shared_ptr mesh = MeshType::shared_ptr(new MeshType("mesh"));

  const std::string infilename = "M6_unstr_euler_c2_p2_sample_result.msh";
  const std::string outfilename = "processed_onera_m6.msh";

  gmsh::GmshReader meshreader;
  gmsh::GmshWriter meshwriter;

  meshreader.read_mesh_from_file(infilename, *mesh);
  MeshManipulator::scale_mesh(*mesh, 1.1963 / 15.1);

  interpolation::FunctionSpace<MeshConfig>::ptr fs_solution =
      std::make_shared<interpolation::FunctionSpace<MeshConfig>>(mesh);

  const PolyOrderID quadrature_order = P3;

  fs_solution->set_reference_fe_values(mesh->topology().dof_storage(), Lagrange, quadrature_order,
                                       Gauss);
  interpolation::FunctionSpace<MeshConfig>::function_t::ptr field_u =
      fs_solution->create_function("u");
  interpolation::FunctionSpace<MeshConfig>::function_t::ptr field_Ma =
      fs_solution->create_function("field_Ma");
  interpolation::FunctionSpace<MeshConfig>::function_t::ptr field_cp =
      fs_solution->create_function("field_cp");


  // This corresponds to angle of attack = 3.06 and Mach number on inlet Ma_in =
  // 0.84
  const Real alpha = 3.06;
  const Real alpha_rad = alpha * PI / 180.0;
  // const Real tan_alpha = 0.05345791105765804881; // tan(3.06*PI/180)
  const Real Ma_in = 0.8395;
  const Real kappa = 1.4;

  const Real rho_in = 1.4;
  const Real u_in = 1.0; // 1.7432/rho;
  const Real v_in = 0.0;
  const Real w_in = u_in * std::tan(PI * alpha / 180.0); // u*tan_alpha;
  const Real speed_in2 = u_in * u_in + v_in * v_in + w_in * w_in;

  const Real e_in =
      speed_in2 / (Ma_in * Ma_in) * rho_in / (kappa * (kappa - 1.)) + 0.5 * rho_in * speed_in2;
  const Real p_in = (kappa - 1.) * (e_in - 0.5 * rho_in * speed_in2);

  math::StaticVector<Real, PhysicsType::DIM> freestream_tangent, freestream_normal;
  freestream_tangent[X] = std::cos(alpha_rad);
  freestream_tangent[Y] = 0.0;
  freestream_tangent[Z] = std::sin(alpha_rad);

  freestream_normal[X] = -std::sin(alpha_rad);
  freestream_normal[Y] = 0.0;
  freestream_normal[Z] = std::cos(alpha_rad);

  meshreader.read_nodal_function_from_file(
      *mesh, infilename, *field_u, { "\"u_0\"", "\"u_1\"", "\"u_2\"", "\"u_3\"", "\"u_4\"" });
  meshreader.read_nodal_function_from_file(*mesh, infilename, *field_Ma, { "\"Ma_0\"" });

  field_cp->resize(1, field_Ma->nb_dof());
  field_cp->fill(0.0);

  typedef interpolation::MeshFunction<MeshConfig, Real>::block_type
      node_value_type;
  typedef interpolation::MeshFunction<MeshConfig, Real>::const_block_type
      const_node_value_type;

  for (Uint n = 0; n < field_Ma->nb_dof(); ++n)
  {
    const_node_value_type u_node = field_u->const_value(n);
    node_value_type cp_node = field_cp->value(n);

    const Real rho = u_node[0];
    const Real speed2 =
        (u_node[1] * u_node[1] + u_node[2] * u_node[2] + u_node[3] * u_node[3]) / (rho * rho);
    const Real p = (kappa - 1.) * (u_node[4] - 0.5 * rho * speed2);

    cp_node[0] = (p - p_in) / (0.5 * rho_in * speed_in2);
  }

  meshwriter.write_mesh_to_file(*mesh, outfilename);
  meshwriter.append_nodal_function_to_file(*mesh, outfilename, *field_u, "u");
  meshwriter.append_nodal_function_to_file(*mesh, outfilename, *field_Ma, "Ma");

  meshwriter.save_mesh_part(*mesh, "Part_Ma.msh", *field_Ma, "Ma", { "Wing", "Symmetry" });
  meshwriter.save_mesh_part(*mesh, "Part_Cp.msh", *field_cp, "cp", { "Wing", "Symmetry" });


  // Integrate over the surface to obtain aerodynamic coefficients

  std::vector<std::string> surface_names = { "Wing" };

  const result_of::geometry<MeshConfig>::type &coordinates = mesh->geometry();

  solver::FEMetric<_3D, _2D> fe_metric;
  solver::FEMetricData<_3D, _2D> metric_data;
  quadrature::Quadrature quad;
  interpolation::FEValues fe_values;
  Uint eshape, poly_order, ref_topology;

  math::DynamicMatrix<Real> coord_at_qd_pt_phys;
  math::DynamicVector<Real> jacobian_in_qd_pt;
  math::DynamicMatrix<Real> element_nodal_coordinates;

  math::DynamicVector<Real> pressure_in_quad_pts;

  math::DynamicMatrix<Real> solution_in_elem_nodes;
  math::DynamicMatrix<Real> solution_in_qd_pts;

  math::DynamicMatrix<Real> normal_in_qd_pt;

  math::StaticVector<Real, PhysicsType::DIM> surface_force;
  surface_force.fill(0.0);

  for (Uint isurf = 0; isurf < surface_names.size(); ++isurf)
  {
    const std::shared_ptr<BoundaryFacets<MeshConfig>> bdry_domain =
        mesh->topology().all_boundaries().domain(surface_names[isurf]);

    for(
        const typename result_of::dof_handler<MeshConfig>::type::const_cell_range_typed &cell_group :
        bdry_domain->all_cell_groups())
    {
      const MeshEntity first_cell = *(cell_group.begin());
      // std::cout << "First cell in group = " << first_cell << std::endl;

      // Detect the element shape, polynomial order and type ('reference
      // topology') by decomposing the cell id
      PointSetTag::decompose_into_fields(first_cell.std_region_id(), eshape, poly_order,
                                          ref_topology);

      quad.change_type(eshape, poly_order, Gauss);

      math::DynamicVector<Real> const &weight_in_qd_pt = quad.get().weights();

      const Uint nb_qd_pts = quad.get().nb_quad_pts();

      fe_values.configure(first_cell.std_region_id(),
                          mesh::sf::SFTag(eshape, Lagrange, Equidist, poly_order, Modal));
      fe_values.fill_Vandermonde(quad.get().coordinates(), quad.get().weights());

      metric_data.resize_variables(fe_values);

      coord_at_qd_pt_phys.resize(nb_qd_pts, first_cell.nb_vert());

      jacobian_in_qd_pt.resize(nb_qd_pts);

      // Resize matrix holding coordinates of one element for given element type
      element_nodal_coordinates.resize(first_cell.nb_vert(), MeshConfig::TDIM);

      solution_in_qd_pts.resize(nb_qd_pts, PhysicsType::NEQ);
      solution_in_elem_nodes.resize(first_cell.nb_vert(), PhysicsType::NEQ);
      pressure_in_quad_pts.resize(nb_qd_pts);

      normal_in_qd_pt.resize(nb_qd_pts, PhysicsType::DIM);

      for (typename result_of::dof_handler<MeshConfig>::type::const_cell_iterator_typed cell_iter =
               cell_group.begin();
           cell_iter != cell_group.end(); ++cell_iter)
      {
        const MeshEntity cell = *cell_iter;

        // Store element coordinates in one matrix
        for (Uint v = 0; v < cell.nb_vert(); ++v)
        {
          element_nodal_coordinates.insert_row(v, coordinates.node(cell.vertex(v)));
        }

        for (Uint v = 0; v < cell.nb_vert(); ++v)
        {
          solution_in_elem_nodes.insert_row(v, field_u->const_value(cell.vertex(v)));
        }

        // Compute metric terms in element
        fe_metric.interpolate(metric_data.V, metric_data.dV, element_nodal_coordinates,
                              coord_at_qd_pt_phys, metric_data.m_dxqd);

        SurfaceJacobian<PhysicsType::DIM, PhysicsType::DIM - 1>::compute(
            metric_data, normal_in_qd_pt, jacobian_in_qd_pt);

        solution_in_qd_pts = metric_data.V * solution_in_elem_nodes;

        for (Uint q = 0; q < nb_qd_pts; ++q)
        {
          const math::ConstVectorBlock<Real> uq = solution_in_qd_pts.const_row_transp(q);

          const Real speed2 = (uq[1] * uq[1] + uq[2] * uq[2] + uq[3] * uq[3]) / (uq[0] * uq[0]);
          const Real p = (kappa - 1.) * (uq[4] - 0.5 * uq[0] * speed2);

          const Real coeff = p * jacobian_in_qd_pt[q] * weight_in_qd_pt[q];

          // '-' sign is here because the normal should be inward-pointing
          surface_force -= coeff * normal_in_qd_pt.const_row_transp(q);
          // surface_force[X] -= coeff * normal_in_qd_pt(q,X);
          // surface_force[Y] -= coeff * normal_in_qd_pt(q,Y);
          // surface_force[Z] -= coeff * normal_in_qd_pt(q,Z);
        }

      } // Loop over all cells in one group

    } // Loop over all cell groups of one surface

  } // Loop over surfaces


  std::cout << "Surface force = " << surface_force << std::endl;

  const Real L =
      (freestream_normal[X] * surface_force[X] + freestream_normal[Y] * surface_force[Y] +
       freestream_normal[Z] * surface_force[Z]);

  const Real D =
      (freestream_tangent[X] * surface_force[X] + freestream_tangent[Y] * surface_force[Y] +
       freestream_tangent[Z] * surface_force[Z]);

  const Real S_ref = 0.7532;

  const Real cL = L / (0.5 * rho_in * speed_in2 * S_ref);
  std::cout << "Lift coefficient = " << cL << std::endl;

  const Real cD = D / (0.5 * rho_in * speed_in2 * S_ref);
  std::cout << "Drag coefficient = " << cD << std::endl;
#endif

  return 0;
}
