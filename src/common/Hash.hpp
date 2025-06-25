#ifndef PDEKIT_Common_Hash_hpp
#define PDEKIT_Common_Hash_hpp

#include <numeric>

#include "common/PDEKit.hpp"

namespace pdekit
{

// ----------------------------------------------------------------------------

namespace common
{

inline size_t xorrer(size_t const &val_l, size_t const &val_r)
{
  // return val ^ 0x9e3779b9;
  const size_t ret_val = val_l + 0x9e3779b9 + (val_r << 6) + (val_r >> 2);
  return ret_val;
}

// ----------------------------------------------------------------------------

template <typename SeqType>
size_t hash_seq_accum(const SeqType &seq)
{
  return std::accumulate(std::begin(seq), std::end(seq), 0, xorrer);
}

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
