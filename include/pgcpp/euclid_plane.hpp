#pragma once

#include "fractions.hpp"
#include "pg_common.hpp" // import cross2, dot1
#include "proj_plane.hpp" // import pg_point, involution, tri_func, quad_func, plucker
#include "proj_plane_concepts.h"
#include <type_traits>

namespace fun
{

/*!
 * @brief
 *
 * @param l
 * @return auto
 */
template <Projective_plane_coord2 L> // // and requires p[i]
constexpr auto fB(const L& l) -> typename L::dual
{
    return {l[0], l[1], 0};
}

/*!
 * @brief
 *
 * @param l
 * @param m
 * @return true
 * @return false
 */
template <Projective_plane_coord2 L>
constexpr auto is_perpendicular(const L& l, const L& m) -> bool
{
    return dot1(l, m) == 0;
}

/*!
 * @brief
 *
 * @param l
 * @param m
 * @return true
 * @return false
 */
template <Projective_plane_coord2 L>
constexpr auto is_parallel(const L& l, const L& m) -> bool
{
    return cross2(l, m) == 0;
}

/*!
 * @brief
 *
 * @param a
 * @param l
 * @return L
 */
Projective_plane_coord { P, L }
constexpr auto altitude(const P& a, const L& l) -> L
{
    return a * fB(l);
}

/*!
 * @brief
 *
 * @param tri
 * @return auto
 */
template <Projective_plane_coord2 P>
constexpr auto tri_altitude(const Triple<P>& tri)
{
    const auto& [a1, a2, a3] = tri;

    auto t1 = altitude(a1, a2 * a3);
    auto t2 = altitude(a2, a3 * a1);
    auto t3 = altitude(a3, a1 * a2);
    return std::tuple {std::move(t1), std::move(t2), std::move(t3)};
}

/*!
 * @brief
 *
 * @param tri
 * @return P
 */
template <Projective_plane_coord2 P>
constexpr auto orthocenter(const Triple<P>& tri) -> P
{
    const auto& [a1, a2, a3] = tri;

    const auto t1 = altitude(a1, a2 * a3);
    const auto t2 = altitude(a2, a1 * a3);
    return t1 * t2;
}

/*!
 * @brief
 *
 * @param m
 * @return auto
 */
template <Projective_plane_coord2 L>
constexpr auto reflect(const L& m)
{
    return involution {m, fB(m)};
}

/*!
 * @brief
 *
 * @param a
 * @param b
 * @return P
 */
template <Projective_plane_coord2 P>
constexpr auto midpoint(const P& a, const P& b) -> P
{
    return plucker(b[2], a, a[2], b);
}

/*!
 * @brief
 *
 * @param tri
 * @return auto
 */
template <Projective_plane_coord2 P>
constexpr auto tri_midpoint(const Triple<P>& tri) -> Triple<P>
{
    const auto& [a1, a2, a3] = tri;

    auto m12 = midpoint(a1, a2);
    auto m23 = midpoint(a2, a3);
    auto m13 = midpoint(a1, a3);
    return {std::move(m12), std::move(m23), std::move(m13)};
}

/*!
 * @brief
 *
 * @tparam K
 * @param x1
 * @param z1
 * @param x2
 * @param z2
 * @return auto
 */
template <CommutativeRing K>
constexpr auto quad1(const K& x1, const K& z1, const K& x2, const K& z2)
{
    if constexpr (Integral<K>)
    {
        Fraction<K> res = sq(Fraction<K>(x1, z1) - Fraction<K>(x2, z2));
        return res;
    }
    else
    {
        return sq(x1 / z1 - x2 / z2);
    }
}

/*!
 * @brief
 *
 * @param a1
 * @param a2
 * @return auto
 */
template <Projective_plane_coord2 P>
constexpr auto quadrance(const P& a1, const P& a2)
{
    return quad1(a1[0], a1[2], a2[0], a2[2]) +
        quad1(a1[1], a1[2], a2[1], a2[2]);
}

// Projective_plane2 { L }
// constexpr auto sbase(const L &l1, const L &l2, const Integer &d) {
//     return Fraction(d, omgB(l1, l1)) * Fraction(d, omgB(l2, l2));
// }

/*!
 * @brief
 *
 * @param l1
 * @param l2
 * @param d
 * @return auto
 */
template <Projective_plane_coord2 L>
constexpr auto sbase(const L& l1, const L& l2, auto const& d)
{
    using K = Value_type<L>;
    if constexpr (Integral<K>)
    {
        Fraction<K> res =
            Fraction<K>(d, dot1(l1, l1)) * Fraction<K>(d, dot1(l2, l2));
        return res;
    }
    else
    {
        return (d * d) / (dot1(l1, l1) * dot1(l2, l2));
    }
}

/*!
 * @brief
 *
 * @param l1
 * @param l2
 * @return auto
 */
template <Projective_plane_coord2 L>
constexpr auto spread(const L& l1, const L& l2)
{
    return sbase(l1, l2, cross2(l1, l2));
}

/*!
 * @brief
 *
 * @param triangle
 * @return auto
 */
template <Projective_plane_coord2 P>
constexpr auto tri_quadrance(const Triple<P>& triangle)
{
    const auto& [a1, a2, a3] = triangle;

    auto m1 = quadrance(a2, a3);
    auto m2 = quadrance(a1, a3);
    auto m3 = quadrance(a1, a2);
    return std::tuple {std::move(m1), std::move(m2), std::move(m3)};
}

/*!
 * @brief
 *
 * @param trilateral
 * @return auto
 */
template <Projective_plane_coord2 L>
constexpr auto tri_spread(const Triple<L>& trilateral)
{
    const auto& [a1, a2, a3] = trilateral;

    auto m1 = spread(a2, a3);
    auto m2 = spread(a1, a3);
    auto m3 = spread(a1, a2);
    return std::tuple {std::move(m1), std::move(m2), std::move(m3)};
}

/*!
 * @brief
 *
 * @param l1
 * @param l2
 * @return auto
 */
template <Projective_plane_coord2 L>
constexpr auto cross_s(const L& l1, const L& l2)
{
    return sbase(l1, l2, dot1(l1, l2));
}

/*!
 * @brief
 *
 * @param lda1
 * @param mu1
 * @return P
 */
template <Projective_plane_coord2 P>
constexpr auto uc_point(const Value_type<P>& lda1, const Value_type<P>& mu1)
{
    const auto lda2 = lda1 * lda1;
    const auto mu2 = mu1 * mu1;
    return P {lda2 - mu2, 2 * lda1 * mu1, lda2 + mu2};
}

/*!
 * @brief Archimedes's function
 *
 * @tparam _Q
 * @param a
 * @param b
 * @param c
 * @return auto
 */
template <CommutativeRing _Q>
constexpr auto Ar(const _Q& a, const _Q& b, const _Q& c)
{
    return 4 * a * b - sq(a + b - c);
}

/*!
 * @brief Cyclic quadrilateral quadrea theorem
 *
 * @tparam _Q
 * @param a
 * @param b
 * @param c
 * @param d
 * @return auto
 */
template <CommutativeRing _Q>
constexpr auto cqq(const _Q& a, const _Q& b, const _Q& c, const _Q& d)
{
    const auto t1 = 4 * a * b;
    const auto t2 = 4 * c * d;
    auto m = (t1 + t2) - sq(a + b - c - d);
    auto p = m * m - 4 * t1 * t2;
    return std::tuple {std::move(m), std::move(p)};
}

/*!
 * @brief
 *
 * @tparam _Q
 * @param quad
 * @return auto
 */
constexpr auto Ptolemy(const auto& quad) -> bool
{
    const auto& [Q12, Q23, Q34, Q14, Q13, Q24] = quad;
    return Ar(Q12 * Q34, Q23 * Q14, Q13 * Q24) == 0;
}

#include <cmath>

/*!
 * @brief
 *
 * @param a
 * @param b
 * @return auto
 */
template <Projective_plane_coord2 P>
constexpr auto distance(const P& a, const P& b)
{
    return std::sqrt(double(quadrance(a, b)));
}

/*!
 * @brief
 *
 * @param l
 * @param m
 * @return auto
 */
template <Projective_plane_coord2 L>
constexpr auto angle(const L& l, const L& m)
{
    return std::asin(std::sqrt(double(spread(l, m))));
}

} // namespace fun
