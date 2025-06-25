#ifndef PDEKIT_Interpolation_Metric_Flags_hpp
#define PDEKIT_Interpolation_Metric_Flags_hpp

#include "common/TaggedBool.hpp"

namespace pdekit
{

namespace interpolation
{

namespace detail
{

struct RebuildMetricIndexTag
{
};
struct ComputeMetricDerivsTag
{
};

} // namespace detail

// ----------------------------------------------------------------------------

using RebuildMetricIndex  = common::TaggedBool<detail::RebuildMetricIndexTag>;
using ComputeMetricDerivs = common::TaggedBool<detail::ComputeMetricDerivsTag>;

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
