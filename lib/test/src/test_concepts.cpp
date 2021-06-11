#include <cassert>
#include <concepts>
#include <pgcpp/proj_plane.hpp>

class LA;

class PA
{
  public:
    using dual = LA;
    using value_type = long;

    PA() = default;
    PA(PA&&) = default;
    PA(PA const&) = delete;
    auto operator=(PA&&) -> PA& = default;
    auto operator=(PA const&) -> PA& = delete;
    ~PA() = default;

    void operator&() const = delete;
    friend void operator,(PA const&, PA const&) = delete;

    friend auto operator==(PA const&, PA const&) -> bool = default;
    [[nodiscard]] static auto aux() -> LA;
    [[nodiscard]] static auto dot(LA const& /*unused*/) -> value_type
    {
        return {};
    }
    auto operator[](size_t /*unused*/) const -> value_type
    {
        return {};
    }
};

class LA
{
  public:
    using dual = PA;
    using value_type = long;

    LA() = default;
    LA(LA&&) = default;
    LA(LA const&) = delete;
    auto operator=(LA&&) -> LA& = default;
    auto operator=(LA const&) -> LA& = delete;
    ~LA() = default;

    void operator&() const = delete;
    friend void operator,(LA const&, LA const&) = delete;
    friend auto operator==(LA const&, LA const&) -> bool = default;
    [[nodiscard]] static auto aux() -> PA;
    [[nodiscard]] static auto dot(PA const& /*unused*/) -> value_type
    {
        return {};
    }
    auto operator[](size_t /*unused*/) const -> value_type
    {
        return {};
    }
};

// struct RsltP
// {
//     operator PA() const
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
//     operator LA() const
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

inline auto operator*(PA const& /*unused*/, PA const& /*unused*/) -> LA
{
    return LA {};
}
inline auto operator*(LA const& /*unused*/, LA const& /*unused*/) -> PA
{
    return PA {};
}
inline auto PA::aux() -> LA
{
    return LA {};
}
inline auto LA::aux() -> PA
{
    return PA {};
}
inline auto plucker(const int& /*unused*/, const PA& /*unused*/,
    const int& /*unused*/, const PA& /*unused*/) -> PA
{
    return PA {};
}
inline auto plucker(const int& /*unused*/, const LA& /*unused*/,
    const int& /*unused*/, const LA& /*unused*/) -> LA
{
    return LA {};
}
inline auto incident(const PA& /*unused*/, const LA& /*unused*/) -> bool
{
    return true;
}
inline auto incident(const LA& /*unused*/, const PA& /*unused*/) -> bool
{
    return true;
}

using namespace fun;
using PArchetype = PA;
using LArchetype = LA;
static_assert(Projective_plane_coord<PArchetype, LArchetype>);
static_assert(Projective_plane_coord<LArchetype, PArchetype>);

inline void test_concept_usage(PArchetype p, LArchetype l)
{
    coincident(p * p, p);
    coincident(l * l, l);
    harm_conj(p, p, p);
    harm_conj(l, l, l);
}
