#ifndef PDEKIT_Interpolation_Vector_Mesh_Function_hpp
#define PDEKIT_Interpolation_Vector_Mesh_Function_hpp

#include <array>
#include <memory>
#include <ostream>

#include "common/PDEKit.hpp"
#include "interpolation/mesh_function/ScalarMeshFunctionBase.hpp"
#include "interpolation/mesh_function/VectorMeshFunctionOps.hpp"
#include "interpolation/mesh_function/function_ops/VectorMeshFunctionAssign.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseVecView.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/adaptation/MeshAdaptSequence.hpp"
#include "mesh/containers/DofMap.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------
//             Implementation of class VectorMeshFunction
// ----------------------------------------------------------------------------

/// STORAGE LAYOUT:
/// Each entry is a scalar, vector or (rank 2) tensor
/// Entry shape for scalars is simply '0' or should be ignored whatsoever
/// Entry shape for vectors is 'M', where M is the length of each vector
/// Entry shape for tensors is M x N
/// The values of each entry in MeshFunction are stored sequentially and
/// followed by values of the next entry. This means that if TensRank is tensor
/// for example and we store 3x3 tensors, each entry has 9 values which are
/// followed by 9 values of the next entry etc. For total K entries, the
/// MeshFunction has to hold data of size K x ( 3 x 3 ) = 9 x K values

template <typename T = Real>
class VectorMeshFunction : public VectorMeshFunctionBase<VectorMeshFunction<T>>
{

  public:
  /// TYPEDEFS

  typedef std::shared_ptr<VectorMeshFunction<T>> ptr;
  typedef std::shared_ptr<VectorMeshFunction<T> const> const_ptr;
  typedef typename math::DenseConstVecView<T, math::ColumnVector> const_entry_type;
  typedef typename math::DenseVecView<T, math::ColumnVector> entry_type;
  typedef T value_type;

  /// Constructor - takes a function space as an argument. The function is
  /// defined in this function space
  VectorMeshFunction(const std::string &space_name, const std::string &name);

  /// Copy constructor
  VectorMeshFunction(const VectorMeshFunction &f_rhs);

  /// Destructor
  virtual ~VectorMeshFunction();

  /// Get the name of this class
  static std::string type_name()
  {
    return "VectorMeshFunction";
  }

  /// Get the name of the function space this function belongs to
  const std::string &space_name() const;

  /// Get the name of this function
  const std::string &name() const;

  /// Allocate some data
  /// @param number of rows and columns in the table
  void resize(const Uint nb_fields, const Uint nb_entries);

  /// Return one column of the data array
  /// @param index of the row
  /// @return Slice vector representing one row
  entry_type value(const Uint i);

  /// Return one column of the data array
  /// @param index of the row
  /// @return Slice vector representing one row
  const const_entry_type const_value(const Uint i) const;

  void insert_value(const Uint i, const_entry_type const &value);

  /// Get number of rows
  Uint nb_fields() const;

  /// Get number of columns
  Uint nb_entries() const;

  /// Assign one scalar value to all entries in this function
  void fill(const T val);

  /// Adapt the mesh function
  template <typename MeshConfig>
  void adapt(const typename result_of::dof_map_t<MeshConfig> &cell_dofs,
             const typename mesh::MeshAdaptSequence<MeshConfig> &schedule);

  /// Apply dof reordering
  void apply_reordering(const std::vector<Int> &reordering);

  // --------------------------------------------------------------------------
  //             Definition of VectorMeshFunction operators
  // --------------------------------------------------------------------------

  /// Assignment operator between two scalar mesh functions
  VectorMeshFunction &operator=(const VectorMeshFunction &f_rhs);

  /// Assignment operator from a scalar mesh function expression
  /// @param rhs - a scalar expression (for example, a sum of some other
  ///              two scalar mesh functions)
  template <typename VF>
  inline VectorMeshFunction<T> &operator=(const VectorMeshFunctionBase<VF> &rhs);

  /// Accumulation from scalar mesh function expression
  template <typename VF>
  inline VectorMeshFunction<T> &operator+=(const VectorMeshFunctionBase<VF> &rhs);

  /// Subtraction of a scalar mesh function expression
  template <typename VF>
  inline VectorMeshFunction<T> &operator-=(const VectorMeshFunctionBase<VF> &rhs);

  /// Scaling by a SCALAR mesh function expression
  template <typename SF>
  inline VectorMeshFunction<T> &operator*=(const ScalarMeshFunctionBase<SF> &rhs);

  private:
  template <typename MeshConfig>
  void p_adapt(const typename result_of::dof_map_t<MeshConfig> &cell_dofs,
               const std::vector<Uint> &cell_p_orders);

  /// Reference to the underlying function space
  std::string m_space_name;

  /// Name of this function
  std::string m_name;

  /// Number of rows (fields in the data)
  Uint m_nb_fields;

  /// Number of columns (entries in the data -
  /// typically, this should be the number of
  /// nodes in the mesh)
  Uint m_nb_entries;

  /// The actual data
  T *m_data;
};

// ----------------------------------------------------------------------------

template <typename T>
VectorMeshFunction<T>::VectorMeshFunction(const std::string &space_name, const std::string &name)
    : m_space_name(space_name), m_name(name), m_nb_fields(0u), m_nb_entries(0u)
{
  m_data = nullptr;
}

// ----------------------------------------------------------------------------

template <typename T>
VectorMeshFunction<T>::VectorMeshFunction(const VectorMeshFunction<T> &f_rhs)
{
  m_space_name = f_rhs.m_space_name;
  m_name       = f_rhs.m_name;

  m_nb_entries = f_rhs.m_nb_entries;
  m_nb_fields  = f_rhs.m_nb_fields;
  m_data       = new T[m_nb_entries * m_nb_fields];

  for (Uint i = 0; i < m_nb_entries * m_nb_fields; ++i)
  {
    m_data[i] = f_rhs.m_data[i];
  }
}

// ----------------------------------------------------------------------------

template <typename T>
VectorMeshFunction<T>::~VectorMeshFunction()
{
  delete[] m_data;
}

// ----------------------------------------------------------------------------

template <typename T>
const std::string &VectorMeshFunction<T>::space_name() const
{
  return m_space_name;
}

// ----------------------------------------------------------------------------

template <typename T>
const std::string &VectorMeshFunction<T>::name() const
{
  return m_name;
}

// ----------------------------------------------------------------------------

template <typename T>
void VectorMeshFunction<T>::resize(const Uint nb_fields, const Uint nb_entries)
{
  if (nb_fields * nb_entries == 0)
  {
    delete[] m_data;
    m_data       = nullptr;
    m_nb_fields  = 0;
    m_nb_entries = 0;
    return;
  }

  if (m_nb_fields * m_nb_entries != nb_fields * nb_entries)
  {
    delete[] m_data;
    m_data = new T[nb_fields * nb_entries];
  }
  m_nb_fields  = nb_fields;
  m_nb_entries = nb_entries;

  fill(T());
}

// ----------------------------------------------------------------------------

template <typename T>
typename VectorMeshFunction<T>::entry_type inline VectorMeshFunction<T>::value(const Uint i)
{
  // return entry_type(m_data + i, m_nb_fields, m_nb_dof);
  return entry_type(m_data + i * m_nb_fields, m_nb_fields);
}

// ----------------------------------------------------------------------------

template <typename T>
const typename VectorMeshFunction<T>::const_entry_type inline VectorMeshFunction<T>::const_value(
    const Uint i) const
{
  // return const_entry_type(m_data + i, m_nb_fields, m_nb_dof);
  return const_entry_type(m_data + i * m_nb_fields, m_nb_fields);
}

// ----------------------------------------------------------------------------

template <typename T>
inline void VectorMeshFunction<T>::insert_value(const Uint i, const_entry_type const &value)
{
  T *data_ptr = m_data + i * m_nb_fields;

  /*
  for (Uint r = 0; r < m_nb_fields; ++r)
  {
    (*data_ptr) = value[r];
    data_ptr += m_nb_dof;
  }
  */

  for (Uint r = 0; r < m_nb_fields; ++r)
  {
    (*data_ptr) = value[r];
    data_ptr++;
  }
}

// ----------------------------------------------------------------------------

template <typename T>
Uint VectorMeshFunction<T>::nb_fields() const
{
  return m_nb_fields;
}

// ----------------------------------------------------------------------------

template <typename T>
Uint VectorMeshFunction<T>::nb_entries() const
{
  return m_nb_entries;
}

// ----------------------------------------------------------------------------

template <typename T>
void VectorMeshFunction<T>::fill(const T val)
{
  for (Uint i = 0; i < m_nb_fields * m_nb_entries; ++i)
  {
    m_data[i] = val;
  }
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig>
void VectorMeshFunction<T>::adapt(const typename result_of::dof_map_t<MeshConfig> &cell_dofs,
                                  const typename mesh::MeshAdaptSequence<MeshConfig> &schedule)
{
  if (schedule.adapt_type() == mesh::AdaptationType::p)
  {
    p_adapt<MeshConfig>(cell_dofs, schedule.cell_poly_orders());
  }
  else if (schedule.adapt_type() == mesh::AdaptationType::h)
  {
    std::cerr << "VectorMeshFunction::adapt: h adaptation is not implemented." << std::endl;
  }
  else if (schedule.adapt_type() == mesh::AdaptationType::hp)
  {
    std::cerr << "VectorMeshFunction::adapt: hp adaptation is not implemented." << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <typename T>
void VectorMeshFunction<T>::apply_reordering(const std::vector<Int> &reordering)
{
  T *reordered_data = new T[m_nb_entries * m_nb_fields];

  for (Uint i = 0; i < m_nb_entries; ++i)
  {
    for (Uint j = 0; j < m_nb_fields; ++j)
    {
      reordered_data[reordering[i] * m_nb_fields + j] = m_data[i * m_nb_fields + j];
    }
  }

  std::swap(reordered_data, m_data);

  delete[] reordered_data;
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig>
void VectorMeshFunction<T>::p_adapt(const typename result_of::dof_map_t<MeshConfig> &cell_dofs,
                                    const std::vector<Uint> &cell_p_orders)
{
  if (m_nb_entries == 0)
  {
    std::cerr << "VectorMeshFunction::p_adapt::nothing to adapt. Aborting." << std::endl;
    return;
  }

  if (cell_dofs.nb_active_cells() != cell_p_orders.size())
  {
    std::cerr << "VectorMeshFunction::p_adapt::number of active cells and "
                 "number of new polynomial\n"
                 "order specified differ. I need to have new polynomial order "
                 "for each cell. Abortig."
              << std::endl;
    return;
  }

  std::cout << "VectorMeshFunction::adapt::starting p-adaptation of function [" << m_name << "]"
            << std::endl;

  Uint tot_nb_dofs = 0;
  mesh::StdRegion old_cell_type;
  std::vector<mesh::StdRegion> new_cell_type(cell_p_orders.size());
  mesh::sf::ShapeFunction sf;

  // ----------------------------------------------------------------------------
  // PASS 1:
  // a) Estimate the needed storage
  //    Only the number of nodes in the dof handler can potentially change,
  //    the number of cells will remain
  // b) Prepare matrices for local interpolation of vertex coordinates
  //    within each element
  // ----------------------------------------------------------------------------

  typedef std::map<std::pair<Uint, Uint>, std::unique_ptr<math::DenseDMat<Real>>>
      interp_mat_map_type;
  interp_mat_map_type dof_value_interp_map;

  // Cell coordinates before change of polynomial order
  std::map<Uint, std::unique_ptr<math::DenseDMat<Real>>> dof_values_in_map;
  // Cell coordinates after change of polynomial order
  std::map<Uint, std::unique_ptr<math::DenseDMat<Real>>> dof_values_out_map;

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const Uint p_new = cell_p_orders[c];

    const mesh::PointSetTag old_cell_type_tag =
        cell_dofs.active_cell_std_region_id(mesh::ActiveIdx(c));
    const mesh::PointSetTag new_cell_type_tag(old_cell_type_tag.elem_shape(), p_new,
                                              old_cell_type_tag.ref_topology());

    new_cell_type[c].change_type(new_cell_type_tag);
    tot_nb_dofs += new_cell_type[c].get().nb_nodes();

    // Check if new interpolation matrix between the coordinates of new and
    // old cell type (which should only differ by polynomial order) needs to
    // be created. If that is the case, create and fill it.
    if (old_cell_type_tag != new_cell_type_tag)
    {
      const std::pair<Uint, Uint> transfer_matrix_key =
          std::pair<Uint, Uint>(old_cell_type_tag.store_value(), new_cell_type_tag.store_value());

      interp_mat_map_type::const_iterator interp_mat_iter =
          dof_value_interp_map.find(transfer_matrix_key);

      if (interp_mat_iter == dof_value_interp_map.end())
      {
        std::unique_ptr<math::DenseDMat<Real>> interp_mat(new math::DenseDMat<Real>());

        old_cell_type.change_type(old_cell_type_tag);
        interp_mat->resize(new_cell_type[c].get().nb_nodes(), old_cell_type.get().nb_nodes());

        const mesh::sf::SFTag sf_tag =
            mesh::sf::SFTag(old_cell_type_tag.elem_shape(), SFunc::Lagrange,
                            old_cell_type_tag.poly_order(), ModalBasis::Modal);

        // Use Lagrange shape functions for the interpolation
        sf.change_type(old_cell_type_tag, sf_tag);
        math::DenseDMat<Real> const &new_ref_coords_in_cell = new_cell_type[c].get().coordinates();
        sf.get().compute_ref_values(new_ref_coords_in_cell, *interp_mat);

        /*
        std::cout << "Interpolation matrix [" <<
        old_cell_type_tag.as_string()
        << " -> "
                  << new_cell_type_tag.as_string() << "]: " <<
        std::endl; std::cout << *interp_mat << std::endl;
        */

        dof_value_interp_map.insert(std::make_pair(transfer_matrix_key, std::move(interp_mat)));
      }
    }

    // Check if one scratch matrix to store the cell coordinates of the
    // 'old' cell is available
    std::map<Uint, std::unique_ptr<math::DenseDMat<Real>>>::const_iterator cell_dofs_in_iter =
        dof_values_in_map.find(old_cell_type_tag.store_value());

    if (cell_dofs_in_iter == dof_values_in_map.end())
    {
      std::unique_ptr<math::DenseDMat<Real>> old_cell_dofs =
          std::unique_ptr<math::DenseDMat<Real>>(new math::DenseDMat<Real>());

      old_cell_type.change_type(old_cell_type_tag);

      old_cell_dofs->resize(old_cell_type.get().nb_nodes(), m_nb_fields);
      dof_values_in_map.insert(std::pair<Uint, std::unique_ptr<math::DenseDMat<Real>>>(
          old_cell_type_tag.store_value(), std::move(old_cell_dofs)));
    }

    // Check if one scratch matrix to store the cell coordinates of the
    // 'new' cell (with changed polynomial order) is available
    std::map<Uint, std::unique_ptr<math::DenseDMat<Real>>>::const_iterator cell_dofs_out_iter =
        dof_values_out_map.find(new_cell_type_tag.store_value());

    if (cell_dofs_out_iter == dof_values_out_map.end())
    {
      std::unique_ptr<math::DenseDMat<Real>> new_cell_dofs =
          std::unique_ptr<math::DenseDMat<Real>>(new math::DenseDMat<Real>());
      new_cell_dofs->resize(new_cell_type[c].get().nb_nodes(), m_nb_fields);
      dof_values_out_map.insert(std::pair<Uint, std::unique_ptr<math::DenseDMat<Real>>>(
          new_cell_type_tag.store_value(), std::move(new_cell_dofs)));
    }
  }

  // --------------------------------------------------------------------------
  // PASS 2: resize the internal data and recompute values
  // --------------------------------------------------------------------------

  // std::cout << "Number of cells: " << m_cell_type.size() << std::endl;
  // std::cout << "Number of new dofs: " << tot_nb_dofs << std::endl;

  T *new_data = new T[tot_nb_dofs * m_nb_fields];
  for (Uint n = 0; n < tot_nb_dofs * m_nb_fields; ++n)
  {
    new_data[n] = T();
  }

  tot_nb_dofs = 0;

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::PointSetTag old_cell_type_tag =
        cell_dofs.active_cell_std_region_id(mesh::ActiveIdx(c));
    const mesh::PointSetTag new_cell_type_tag = new_cell_type[c].get().pt_set_id();

    const mesh::MeshEntity cell = cell_dofs.active_cell(mesh::ActiveIdx(c));

    // If the 'old' and 'new' cell types differ, we will have to
    // interpolate coordinates
    if (old_cell_type_tag != new_cell_type_tag)
    {
      // Get the transfer matrix between old and new coordinates
      const std::pair<Uint, Uint> transfer_matrix_key =
          std::pair<Uint, Uint>(old_cell_type_tag.store_value(), new_cell_type_tag.store_value());

      interp_mat_map_type::const_iterator interp_mat_iter =
          dof_value_interp_map.find(transfer_matrix_key);

      const math::DenseDMat<Real> &interp_mat = *(interp_mat_iter->second);

      // Get the vector which will hold the values corresponding to the
      // old cell
      std::map<Uint, std::unique_ptr<math::DenseDMat<Real>>>::const_iterator dof_values_in_iter =
          dof_values_in_map.find(old_cell_type_tag.store_value());

      math::DenseDMat<Real> &dofs_in = *(dof_values_in_iter->second);

      // Get the vector which will temporarily hold the new values
      // corresponding to the same cell
      std::map<Uint, std::unique_ptr<math::DenseDMat<Real>>>::const_iterator dof_values_out_iter =
          dof_values_out_map.find(new_cell_type_tag.store_value());

      math::DenseDMat<Real> &dofs_out = *(dof_values_out_iter->second);

      for (Uint n_in = 0; n_in < cell.nb_vert(); ++n_in)
      {
        for (Uint f = 0; f < m_nb_fields; ++f)
        {
          dofs_in(n_in, f) = m_data[cell.vertex(n_in) * m_nb_fields + f];
        }
      }

      dofs_out = interp_mat * dofs_in;

      // Put the coordinates in new geometry
      for (Uint n_out = 0; n_out < dofs_out.rows(); ++n_out)
      {
        for (Uint f = 0; f < m_nb_fields; ++f)
        {
          new_data[tot_nb_dofs * m_nb_fields + f] = dofs_out(n_out, f);
        }
        tot_nb_dofs++;
      }
    }
    else // Otherwise (cell type stays the same) we can just copy the
         // coordinates
    {
      for (Uint n_in = 0; n_in < cell.nb_vert(); ++n_in)
      {
        for (Uint f = 0; f < m_nb_fields; ++f)
        {
          new_data[tot_nb_dofs * m_nb_fields + f] = m_data[cell.vertex(n_in) * m_nb_fields + f];
        }
        tot_nb_dofs++;
      }
    }
  } // Loop over cells which are adapted

  m_nb_entries = tot_nb_dofs;
  delete[] m_data;
  m_data   = new_data;
  new_data = nullptr;

  std::cout << "VectorMeshFunction::p_adapt::finished p-adaptation of function [" << m_name << "]"
            << std::endl;
}

// ----------------------------------------------------------------------------
//             Implementation of VectorMeshFunction operators
// ----------------------------------------------------------------------------

template <typename T>
VectorMeshFunction<T> &VectorMeshFunction<T>::operator=(const VectorMeshFunction<T> &f_rhs)
{
  if (&f_rhs != this)
  {
    m_space_name = f_rhs.m_space_name;
    m_name       = f_rhs.m_name;

    resize(f_rhs.m_nb_fields, f_rhs.m_nb_entries);

    /*
    for (Uint i = 0; i < m_nb_entries * m_nb_fields; ++i)
    {
      m_data[i] = f_rhs.m_data[i];
    }
    */

    detail::VectorMeshFunctionAssign::assign(*this, f_rhs);
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename VF>
inline VectorMeshFunction<T> &VectorMeshFunction<T>::operator=(
    const VectorMeshFunctionBase<VF> &rhs)
{
  const VF &rhs_expr = rhs.wrapped_type();
  for (Uint n = 0; n < m_nb_entries; ++n)
  {
    this->value(n) = rhs_expr.const_value(n);
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename VF>
inline VectorMeshFunction<T> &VectorMeshFunction<T>::operator+=(
    const VectorMeshFunctionBase<VF> &rhs)
{
  const VF &rhs_expr = rhs.wrapped_type();
  for (Uint n = 0; n < m_nb_entries; ++n)
  {
    this->value(n) += rhs_expr.const_value(n);
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename VF>
inline VectorMeshFunction<T> &VectorMeshFunction<T>::operator-=(
    const VectorMeshFunctionBase<VF> &rhs)
{
  const VF &rhs_expr = rhs.wrapped_type();
  for (Uint n = 0; n < m_nb_entries; ++n)
  {
    this->value(n) -= rhs_expr.const_value(n);
  }

  return (*this);
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename SF>
inline VectorMeshFunction<T> &VectorMeshFunction<T>::operator*=(
    const ScalarMeshFunctionBase<SF> &rhs)
{
  const SF &rhs_expr = rhs.wrapped_type();
  for (Uint n = 0; n < m_nb_entries; ++n)
  {
    entry_type node_value = this->value(n);

    for (Uint f = 0; f < m_nb_fields; ++f)
    {
      node_value[f] *= rhs_expr[n];
    }
  }

  return (*this);
}

// ----------------------------------------------------------------------------

// Print the contents of the table:
template <typename T>
std::ostream &operator<<(std::ostream &os, const VectorMeshFunction<T> &function)
{
  for (Uint i = 0; i < function.nb_entries(); ++i)
  {
    os << function.const_value(i) << std::endl;
  }
  return os;
}

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // Namespace pdekit

#endif
