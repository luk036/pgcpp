#pragma once

#include "fractions.hpp"
#include "proj_plane_concepts.h"
#include <cassert>
#include <tuple>

/*! @file include/proj_plane.hpp
 *  This is a C++ Library header.
 */

/*!
 @todo: projectivity >=
**/

namespace fun
{

/**
 * @brief Coincident
 * 
 * @tparam[in] L Line
 * @tparam[in] Args points
 * @return true if points are conincident (on a line l)
 * @return false otherwise
 */
template <typename L, typename... Args>
requires (Projective_plane_prim<L, Args> && ...)
constexpr auto coincident(const L& l, const Args&... r) -> bool
{
    return (incident(r, l) && ...);
}

template <typename P>
using Triple = std::tuple<P, P, P>;

/*!
 * @brief
 *
 * @param[in] tri
 * @return auto
 */
template <Projective_plane_prim2 P>
constexpr auto tri_dual(const Triple<P>& tri)

{
    const auto& [a1, a2, a3] = tri;
    assert(!coincident(a2 * a3, a1));
    return std::tuple {a2 * a3, a1 * a3, a1 * a2};
}

/*!
 * @brief
 *
 * @param[in] func
 * @param[in] tri
 * @return auto
 */
template <Projective_plane_prim2 P, typename Fn>
constexpr auto tri_func(Fn&& func, const Triple<P>& tri)

{
    const auto& [a1, a2, a3] = tri;
    return std::tuple{func(a2, a3), func(a1, a3), func(a1, a2)};
}

/*!
 * @brief return whether two triangles are perspective 
 *
 * @param[in] tri1
 * @param[in] tri2
 * @return true
 * @return false
 */
template <Projective_plane_prim2 P>
constexpr auto persp(const Triple<P>& tri1, const Triple<P>& tri2)
-> bool
{
    const auto& [A, B, C] = tri1;
    const auto& [D, E, F] = tri2;
    const auto O = (A * D) * (B * E);
    return incident(O, C * F);
}

/*!
 * @brief
 *
 * @param[in] p
 * @param[in] l
 * @return true
 * @return false
 */
template <typename P, typename L>
requires Projective_plane<P, L>
constexpr auto incident(const P& p, const L& l)
-> bool
{
    return p.dot(l) == Value_type<P>(0);
}

/*!
 * @brief
 *
 * @tparam P
 * @param[in] A
 * @param[in] B
 * @param[in] C
 * @return constexpr P
 */
template <Projective_plane2 P>
constexpr auto harm_conj(const P& A, const P& B, const P& C)
-> P
{
    const auto lC = C * (A * B).aux();
    return plucker(B.dot(lC), A, A.dot(lC), B);
}

/**
 * @brief 
 * 
 * @tparam P 
 * @tparam L 
 */
template <typename P, typename L>
requires Projective_plane<P, L>
class involution
{
    using K = Value_type<P>;

  private:
    const L& _m;
    P _o;
    K _c;

  public:
    /*!
     * @brief Construct a new involution object
     *
     * @param[in] m
     * @param[in] o
     */
    constexpr involution(const L& m, P o) // input mirror and center
        : _m {m}
        , _o {std::move(o)}
        , _c {m.dot(_o)}
    {
    }

    /*!
     * @brief
     *
     * @param[in] p
     * @return P
     */
    constexpr auto operator()(const P& p) const -> P
    {
        return plucker(_c, p, K(-2 * p.dot(_m)), _o);
    }
};

/*!
 * @brief
 *
 * @tparam K
 * @param[in] a
 * @param[in] b
 * @param[in] c
 * @param[in] d
 * @return auto
 */
template <ring K>
constexpr auto ratio_ratio(const K& a, const K& b, const K& c, const K& d)
{
    if constexpr (Integral<K>)
    {
        return Fraction(a, b) / Fraction(c, d);
    }
    else
    {
        return (a * d) / (b * c);
    }
}

/*!
 * @brief Cross Ratio
 *
 * @tparam P
 * @tparam L
 * @param[in] A point \in P
 * @param[in] B point \in P
 * @param[in] l line \in P
 * @param[in] m line \in P
 * @return cross ratio R(A,B;l,m)
 *
 * @todo rewrite by projecting to the y-axis first [:2]
 */
template <typename P, typename L>
requires Projective_plane<P, L>
constexpr auto x_ratio(const P& A, const P& B, const L& l, const L& m)
{
    return ratio_ratio(A.dot(l), A.dot(m), B.dot(l), B.dot(m));
}

/*!
 * @brief
 *
 * @param[in] A
 * @param[in] B
 * @param[in] C
 * @param[in] D
 * @return constexpr auto
 */
template <Projective_plane2 P>
constexpr auto R(const P& A, const P& B, const P& C, const P& D)
{
    const auto O = (C * D).aux();
    return x_ratio(A, B, O * C, O * D);
}

/*!
 * @brief
 *
 * @param[in] A
 * @param[in] B
 * @param[in] C
 * @param[in] D
 * @return constexpr auto
 */
template <Projective_plane_coord2 P>
constexpr auto R0(const P& A, const P& B, const P& C, const P& D)
{
    return ratio_ratio(cross0(A, C), cross0(A, D),
                       cross0(B, C), cross0(B, D)); 
}

/*!
 * @brief
 *
 * @param[in] A
 * @param[in] B
 * @param[in] C
 * @param[in] D
 * @return constexpr auto
 */
template <Projective_plane_coord2 P>
constexpr auto R1(const P& A, const P& B, const P& C, const P& D)
{
    return ratio_ratio(cross1(A, C), cross1(A, D),
                       cross1(B, C), cross1(B, D)); 
}

/*!
 * @brief
 *
 * @param[in] A
 * @param[in] B
 * @param[in] C
 * @param[in] D
 * @return constexpr auto
 */
template <Projective_plane_coord2 P>
constexpr auto R(const P& A, const P& B, const P& C, const P& D)

{
    using K = Value_type<P>;
    if (cross0(A, B) != K(0))
    { // Project points to yz-plane
        return R0(A, B, C, D);
    }
    // Project points to xz-plane
    return R1(A, B, C, D);
}

/*!
 * @brief
 *
 * @param[in] A
 * @param[in] B
 * @param[in] C
 * @param[in] D
 * @return constexpr auto
 */
template <Projective_plane2 P>
constexpr auto is_harmonic(const P& A, const P& B, const P& C, const P& D) -> bool
{
    using K = Value_type<P>;
    return R(A, B, C, D) == K(-1);
}

/*!
 * @brief Check Pappus Theorem
 *
 * @tparam P
 * @tparam L
 * @param[in] co1
 * @param[in] co2
 */
template <Projective_plane_prim2 P>
void check_pappus(const Triple<P>& co1, const Triple<P>& co2)

{
    const auto& [A, B, C] = co1;
    const auto& [D, E, F] = co2;

    const auto G = (A * E) * (B * D);
    const auto H = (A * F) * (C * D);
    const auto I = (B * F) * (C * E);
    assert(coincident(G, H, I));
}

/*!
 * @brief
 *
 * @param[in] tri1
 * @param[in] tri2
 */
template <Projective_plane_prim2 P>
void check_desargue(const Triple<P>& tri1, const Triple<P>& tri2)
{
    const auto trid1 = tri_dual(tri1);
    const auto trid2 = tri_dual(tri2);
    const auto b1 = persp(tri1, tri2);
    const auto b2 = persp(trid1, trid2);
    assert((b1 && b2) || (!b1 && !b2));
}

} // namespace fun
