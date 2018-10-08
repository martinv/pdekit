#ifndef PDEKIT_Interpolation_Vector_Mesh_Function_Base_hpp
#define PDEKIT_Interpolation_Vector_Mesh_Function_Base_hpp

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <typename MFT>
class VectorMeshFunctionBase
{
  public:
  /// Typedef for the mesh function type wrapped by this class
  typedef MFT MeshFunctionType;

  inline MeshFunctionType &wrapped_type()
  {
    return static_cast<MeshFunctionType &>(*this);
  }

  inline MeshFunctionType const &wrapped_type() const
  {
    return static_cast<MeshFunctionType const &>(*this);
  }
};

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // Namespace pdekit

#endif
