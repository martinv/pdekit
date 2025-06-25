#ifndef PDEKIT_Mesh_Mesh_Index_hpp
#define PDEKIT_Mesh_Mesh_Index_hpp

#include "common/PDEKit.hpp"
#include "common/TaggedInt.hpp"

namespace pdekit
{

namespace mesh
{

struct FlatIdxTag
{
};

struct ActiveIdxTag
{
};

typedef common::TaggedInt<Int, FlatIdxTag> FlatIdx;
typedef common::TaggedInt<Int, ActiveIdxTag> ActiveIdx;

} // namespace mesh

} // namespace pdekit

#endif
