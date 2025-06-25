#ifndef PDEKIT_Mesh_Containers_View_Traits_hpp
#define PDEKIT_Mesh_Containers_View_Traits_hpp

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Traits for non-constant and constant Views into mesh containers
// ----------------------------------------------------------------------------

struct ViewIsConst
{
};
struct ViewIsNotConst
{
};

} // namespace mesh

} // namespace pdekit

#endif
