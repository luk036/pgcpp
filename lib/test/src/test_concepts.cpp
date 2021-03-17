#include <cassert>
#include <concepts>
#include <pgcpp/proj_plane.hpp>

class LA;

class PA
{
  public:
    using dual = LA;
    using value_type = long;

    PA() noexcept = default;
    PA(PA&&) = default;
    PA(PA const&) = delete;
    PA& operator=(PA&&) = default;
    PA& operator=(PA const&) = delete;
    ~PA() = default;

    void operator&() const = delete;
    friend void operator,(PA const&, PA const&) = delete;

    friend bool operator==(PA const&, PA const&) = default;
    LA aux() const noexcept;
    value_type dot(LA const&) const noexcept
    {
        return {};
    }
    value_type operator[](size_t) const noexcept
    {
        return {};
    }
};

class LA
{
  public:
    using dual = PA;
    using value_type = long;

    LA() noexcept = default;
    LA(LA&&) = default;
    LA(LA const&) = delete;
    LA& operator=(LA&&) = default;
    LA& operator=(LA const&) = delete;
    ~LA() = default;

    void operator&() const = delete;
    friend void operator,(LA const&, LA const&) = delete;
    friend bool operator==(LA const&, LA const&) = default;
    PA aux() const noexcept;
    value_type dot(PA const&) const noexcept
    {
        return {};
    }
    value_type operator[](size_t) const noexcept
    {
        return {};
    }
};

// struct RsltP
// {
//     operator PA() const noexcept
//     {
//         return PA{};
//     }
//     RsltP() = default;
//     RsltP(RsltP&&) = delete;
//     RsltP(RsltP const&) = delete;
//     RsltP& operator=(RsltP&&) = default;
//     RsltP& operator=(RsltP const&) = delete;
//     ~RsltP() = default;

//     void operator&() const = delete;
//     LA aux() const { return {}; }
//     friend void operator,(RsltP, RsltP) = delete;
// };


// struct RsltL
// {
//     operator LA() const noexcept
//     {
//         return LA{};
//     }
//     RsltL() = default;
//     RsltL(RsltL&&) = delete;
//     RsltL(RsltL const&) = delete;
//     RsltL& operator=(RsltL&&) = default;
//     RsltL& operator=(RsltL const&) = delete;
//     ~RsltL() = default;

//     void operator&() const = delete;
//     PA aux() const { return {}; }
//     friend void operator,(RsltL, RsltL) = delete;
// };

inline LA operator*(PA const&, PA const&) noexcept
{
    return LA{};
}
inline PA operator*(LA const&, LA const&) noexcept
{
    return PA{};
}
inline LA PA::aux() const noexcept
{
    return LA{};
}
inline PA LA::aux() const noexcept
{
    return PA{};
}
inline PA plucker(const int&, const PA&, const int&, const PA&) noexcept
{
    return PA{};
}
inline LA plucker(const int&, const LA&, const int&, const LA&) noexcept
{
    return LA{};
}
inline bool incident(const PA&, const LA&) noexcept
{
    return true;
}
inline bool incident(const LA&, const PA&) noexcept
{
    return true;
}

using namespace fun;
using PArchetype = PA;
using LArchetype = LA;
static_assert(Projective_plane_coord<PArchetype, LArchetype>);
static_assert(Projective_plane_coord<LArchetype, PArchetype>);

inline void test_concept_usage(PArchetype p, LArchetype l) noexcept
{
    coincident(p * p, p);
    coincident(l * l, l);
    harm_conj(p, p, p);
    harm_conj(l, l, l);
}
