#ifndef PDEKIT_Common_Tagged_Bool_hpp
#define PDEKIT_Common_Tagged_Bool_hpp

// Based on tagged_bool class of Andrzej Krzemienski

namespace pdekit
{

namespace common
{

template <typename Tag>
class TaggedBool
{
  public:
  constexpr explicit TaggedBool(bool v) : m_value{v}
  {
  }

  constexpr explicit TaggedBool(int)    = delete;
  constexpr explicit TaggedBool(double) = delete;
  constexpr explicit TaggedBool(void *) = delete;

  template <typename OtherTag>
  constexpr explicit TaggedBool(TaggedBool<OtherTag> b) : m_value{b.m_value}
  {
  }

  constexpr explicit operator bool() const
  {
    return m_value;
  }
  constexpr TaggedBool operator!() const
  {
    return TaggedBool{!m_value};
  }

  friend constexpr bool operator==(TaggedBool l, TaggedBool r)
  {
    return l.m_value == r.m_value;
  }
  friend constexpr bool operator!=(TaggedBool l, TaggedBool r)
  {
    return l.m_value != r.m_value;
  }

  private:
  bool m_value;

  template <typename /*OtherTag*/>
  friend class TaggedBool;
};

} // namespace common

} // namespace pdekit

#endif
