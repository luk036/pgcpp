#ifndef _HOME_UBUNTU_CUBSTORE_PROJ_GEOM_PGCPP_EUCLID_PLANE_HPP
#define _HOME_UBUNTU_CUBSTORE_PROJ_GEOM_PGCPP_EUCLID_PLANE_HPP 1

//#include "proj_plane_concepts.h"
#include "pg_common.hpp"
#include "proj_plane.hpp"
#include "proj_plane_concepts.h"
#include <cassert>
#include <cmath>
#include <type_traits>

namespace fun
{

/**
 * @brief 
 * 
 * @param x 
 * @param y 
 * @return auto 
 */
Projective_plane_coord2 { L } // and requires p[i]
constexpr auto dot1(const L &x, const L &y)
{
    return x[0] * y[0] + x[1] * y[1];
}

/**
 * @brief 
 * 
 * @param l 
 * @return auto 
 */
Projective_plane_coord2 { L } // // and requires p[i]
constexpr auto fB(const L &l)
{
    using P = typename L::dual;
    using K = Value_type<P>;
    return P{l[0], l[1], K(0)};
}

/**
 * @brief 
 * 
 * @param l 
 * @param m 
 * @return true 
 * @return false 
 */
Projective_plane_coord2 { L }
constexpr bool is_perpendicular(const L &l, const L &m)
{
    return dot1(l, m) == 0;
}

/**
 * @brief 
 * 
 * @param l 
 * @param m 
 * @return true 
 * @return false 
 */
Projective_plane_coord2 { L }
constexpr bool is_parallel(const L &l, const L &m)
{
    return l[0] * m[1] == l[1] * m[0];
}

/**
 * @brief 
 * 
 * @param a 
 * @param l 
 * @return L 
 */
Projective_plane_coord { P, L }
constexpr L altitude(const P &a, const L &l) { return a * fB(l); }

/**
 * @brief 
 * 
 * @param tri 
 * @return auto 
 */
Projective_plane_coord2 { P }
constexpr auto tri_altitude(const Triple<P> &tri)
{
    auto const &[l1, l2, l3] = tri_dual(tri);
    auto const &[a1, a2, a3] = tri;
    auto t1 = altitude(a1, l1);
    auto t2 = altitude(a2, l2);
    auto t3 = altitude(a3, l3);
    return std::tuple{std::move(t1), std::move(t2), std::move(t3)};
}

/**
 * @brief 
 * 
 * @param tri 
 * @return P 
 */
Projective_plane_coord2 { P }
constexpr P orthocenter(const Triple<P> &tri)
{
    auto const &[a1, a2, a3] = tri;
    auto t1 = altitude(a1, a2 * a3);
    auto t2 = altitude(a2, a1 * a3);
    return P{t1 * t2};
}

/**
 * @brief 
 * 
 * @param m 
 * @return auto 
 */
Projective_plane_coord2 { L }
constexpr auto reflect(const L &m) { return involution{m, fB(m)}; }

/**
 * @brief 
 * 
 * @param x 
 * @param y 
 * @return auto 
 */
Projective_plane_coord2 { P }
constexpr auto omgB(const P &x, const P &y)
{
    return x[0] * y[0] + x[1] * y[1];
}

/**
 * @brief 
 * 
 * @param x 
 * @param y 
 * @return auto 
 */
Projective_plane_coord2 { P }
constexpr auto det(const P &x, const P &y) { return x[0] * y[1] - x[1] * y[0]; }

/**
 * @brief 
 * 
 * @param a 
 * @param b 
 * @return P 
 */
Projective_plane_coord2 { P }
constexpr P midpoint(const P &a, const P &b)
{
    return plucker(b[2], a, a[2], b);
}

/**
 * @brief 
 * 
 * @param tri 
 * @return auto 
 */
Projective_plane_coord2 { P }
constexpr auto tri_midpoint(const Triple<P> &tri)
{
    auto const &[a1, a2, a3] = tri;
    auto m12 = midpoint(a1, a2);
    auto m23 = midpoint(a2, a3);
    auto m13 = midpoint(a1, a3);
    return Triple<P>{std::move(m12), std::move(m23), std::move(m13)};
}

// Integer { K }
// constexpr auto quad1(const K &x1, const K &z1, const K &x2, const K &z2) {
//     return sq(Fraction(x1, z1) - Fraction(x2, z2));
// }

/**
 * @brief 
 * 
 * @tparam K 
 * @param x1 
 * @param z1 
 * @param x2 
 * @param z2 
 * @return auto 
 */
template <typename K>
constexpr auto quad1(const K &x1, const K &z1,
                     const K &x2, const K &z2)
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

/**
 * @brief 
 * 
 * @param a1 
 * @param a2 
 * @return auto 
 */
Projective_plane_coord2 { P }
constexpr auto quadrance(const P &a1, const P &a2)
{
    return quad1(a1[0], a1[2], a2[0], a2[2]) +
           quad1(a1[1], a1[2], a2[1], a2[2]);
}

// Projective_plane2 { L }
// constexpr auto sbase(const L &l1, const L &l2, const Integer &d) {
//     return Fraction(d, omgB(l1, l1)) * Fraction(d, omgB(l2, l2));
// }

/**
 * @brief 
 * 
 * @param l1 
 * @param l2 
 * @param d 
 * @return auto 
 */
Projective_plane2 { L }
constexpr auto sbase(const L &l1, const L &l2, auto const &d)
{
    using K = Value_type<L>;
    if constexpr (Integral<K>)
    {
        Fraction<K> res =
            Fraction<K>(d, omgB(l1, l1)) * Fraction<K>(d, omgB(l2, l2));
        return res;
    }
    else
    {
        return (d * d) / (omgB(l1, l1) * omgB(l2, l2));
    }
}

/**
 * @brief 
 * 
 * @param l1 
 * @param l2 
 * @return auto 
 */
Projective_plane_coord2 { L }
constexpr auto spread(const L &l1, const L &l2)
{
    return sbase(l1, l2, det(l1, l2));
}

/**
 * @brief 
 * 
 * @param triangle 
 * @return auto 
 */
Projective_plane_coord2 { P }
constexpr auto tri_quadrance(const Triple<P> &triangle)
{
    auto const &[a1, a2, a3] = triangle;
    auto m1 = quadrance(a2, a3);
    auto m2 = quadrance(a1, a3);
    auto m3 = quadrance(a1, a2);
    return std::tuple{std::move(m1), std::move(m2), std::move(m3)};
}

/**
 * @brief 
 * 
 * @param trilateral 
 * @return auto 
 */
Projective_plane_coord2 { L }
constexpr auto tri_spread(const Triple<L> &trilateral)
{
    auto const &[a1, a2, a3] = trilateral;
    auto m1 = spread(a2, a3);
    auto m2 = spread(a1, a3);
    auto m3 = spread(a1, a2);
    return std::tuple{std::move(m1), std::move(m2), std::move(m3)};
}

/**
 * @brief 
 * 
 * @param l1 
 * @param l2 
 * @return auto 
 */
Projective_plane_coord { P, L }
constexpr auto cross_s(const L &l1, const L &l2)
{
    return sbase(l1, l2, omgB(l1, l2));
}

/**
 * @brief 
 * 
 * @param lambda1 
 * @param mu1 
 * @return P 
 */
Projective_plane2 { P }
constexpr P uc_point(const Value_type<P> &lambda1, const Value_type<P> &mu1)
{
    auto lambda2 = lambda1 * lambda1;
    auto mu2 = mu1 * mu1;
    return P(lambda2 - mu2, 2 * lambda1 * mu1, lambda2 + mu2);
}

/**
 * @brief Archimedes's function
 * 
 * @tparam _Q 
 * @param a 
 * @param b 
 * @param c 
 * @return auto 
 */
template <typename _Q>
constexpr auto Ar(const _Q &a, const _Q &b, const _Q &c)
{
    return 4 * a * b - sq(a + b - c);
}

/**
 * @brief Cyclic quadrilateral quadrea theorem
 * 
 * @tparam _Q 
 * @param a 
 * @param b 
 * @param c 
 * @param d 
 * @return auto 
 */
template <typename _Q>
constexpr auto cqq(const _Q &a, const _Q &b,
                   const _Q &c, const _Q &d)
{
    auto t1 = 4 * a * b;
    auto t2 = 4 * c * d;
    auto m = (t1 + t2) - sq(a + b - c - d);
    auto p = m * m - 4 * t1 * t2;
    return std::tuple{m, p};
}

/**
 * @brief 
 * 
 * @tparam _Q 
 * @param quad 
 * @return auto 
 */
template <typename _Q>
constexpr auto Ptolemy(const std::tuple<_Q, _Q, _Q, _Q, _Q, _Q> &quad)
{
    auto const &[Q12, Q23, Q34, Q14, Q13, Q24] = quad;
    return Ar(Q12 * Q34, Q23 * Q14, Q13 * Q24) == 0;
}

/**
 * @brief 
 * 
 * @param a 
 * @param b 
 * @return auto 
 */
Projective_plane_coord2 { P }
constexpr auto distance(const P &a, const P &b)
{
    return std::sqrt(double(quadrance(a, b)));
}

/**
 * @brief 
 * 
 * @param l 
 * @param m 
 * @return auto 
 */
Projective_plane_coord2 { L }
constexpr auto angle(const L &l, const L &m)
{
    return std::asin(std::sqrt(double(spread(l, m))));
}

} // namespace fun

#endif
