#ifndef _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP 1

#include <type_traits> // std::is_base_of_v
#include "pgcpp/pg_common.hpp"
#include "pgcpp/proj_plane.hpp"
#include "pgcpp/proj_plane_concepts.h"

namespace fun
{

/*!
 * @brief
 *
 * @tparam _P
 * @tparam _L
 * @tparam Derived
 */
template<typename _P, typename _L, template<typename P, typename L> class Derived>
requires Projective_plane_prim<_P, _L> // c++20 concept
struct ck
{
    using point_t = _P;
    using line_t  = _L;

    using cDer = const Derived<_P, _L>;
    cDer& self = *static_cast<cDer*>(this);

    constexpr explicit ck()
    {
        static_assert(std::is_base_of_v<ck<_P, _L, Derived>, Derived<_P, _L>>);
    }

    /*!
     * @brief
     *
     * @tparam L
     * @param l
     * @param m
     * @return true
     * @return false
     */
    template<Projective_plane_prim2 L>
    constexpr auto is_perpendicular(const L& l, const L& m) const -> bool
    {
        return incident(m, self.perp(l));
    }

    /*!
     * @brief
     *
     * @param p
     * @param l
     * @return L
     */
    template<typename P, typename L = typename P::dual>
    requires Projective_plane_prim<P, L> // c++20 concept
    constexpr auto altitude(const P& p, const L& l) const -> L
    {
        return p * self.perp(l);
    }

    /*!
     * @brief
     *
     * @param tri
     * @return std::tuple
     */
    constexpr auto tri_altitude(const Triple<Projective_plane_prim2>& tri) const
    {
        auto [l1, l2, l3]   = tri_dual(tri);
        auto&& [a1, a2, a3] = tri;

        auto t1 = this->altitude(a1, l1);
        auto t2 = this->altitude(a2, l2);
        auto t3 = this->altitude(a3, l3);
        return std::tuple{std::move(t1), std::move(t2), std::move(t3)};
    }

    /*!
     * @brief
     *
     * @param tri
     * @return P
     */
    constexpr auto orthocenter(const Triple<Projective_plane_prim2>& tri) const
    {
        auto&& [a1, a2, a3] = tri;

        auto t1 = this->altitude(a1, a2 * a3);
        auto t2 = this->altitude(a2, a1 * a3);
        return t1 * t2;
    }

    auto reflect(const Projective_plane2& m) const { return involution{m, self.perp(m)}; }

    constexpr auto tri_measure(const Triple<Projective_plane2>& tri) const
    {
        auto&& [a1, a2, a3] = tri;

        auto m1 = self.measure(a2, a3);
        auto m2 = self.measure(a1, a3);
        auto m3 = self.measure(a1, a2);
        return std::tuple{std::move(m1), std::move(m2), std::move(m3)};
    }

    constexpr auto quadrance(const _P& p, const _P& q) const { return self.measure(p, q); }

    constexpr auto spread(const _L& l, const _L& m) const { return self.measure(l, m); }

    constexpr auto tri_quadrance(const Triple<_P>& triangle) const
    {
        return this->tri_measure(triangle);
    }

    constexpr auto tri_spread(const Triple<_L>& trilateral) const
    {
        return this->tri_measure(trilateral);
    }
};

/*!
 * @brief
 *
 * @tparam Q_t
 * @tparam S_t
 * @param Q
 * @param S
 * @return true
 * @return false
 */
template<CommutativeRing Q_t>
constexpr bool check_sine_law(const Triple<Q_t>& Q, const Triple<Q_t>& S)
{
    auto&& [q1, q2, q3] = Q;
    auto&& [s1, s2, s3] = S;
    return (s1 * q2 == s2 * q1) && (s2 * q3 == s3 * q2);
}

/*!
 * @brief Elliptic Plane
 *
 * @tparam P
 * @tparam P::dual
 */
template<typename P, typename L = typename P::dual>
requires Projective_plane_prim<P, L> // c++20 concept
struct ellck : ck<P, L, ellck>
{
    constexpr auto perp(const P& v) const -> L { return L(v); }
    constexpr auto perp(const L& v) const -> P { return P(v); }

    template<Projective_plane2 _P>
    constexpr auto measure(const _P& a1, const _P& a2) const
    {
        auto x = x_ratio(a1, a2, this->perp(a2), this->perp(a1));
        return 1 - x;
    }
};

/*!
 * @brief Hyperbolic Plane
 *
 * @tparam P
 * @tparam P::dual
 */
template<typename P, typename L = typename P::dual>
requires Projective_plane_prim<P, L> // c++20 concept
struct hyck : ck<P, L, hyck>
{
    constexpr auto perp(const P& v) const { return L(v[0], v[1], -v[2]); }
    constexpr auto perp(const L& v) const { return P(v[0], v[1], -v[2]); }

    template<Projective_plane2 _P>
    constexpr auto measure(const _P& a1, const _P& a2) const
    {
        auto x = x_ratio(a1, a2, this->perp(a2), this->perp(a1));
        return 1 - x;
    }
};

/*!
 * @brief
 *
 * @tparam K
 * @param Q
 * @return constexpr auto
 */
template<CommutativeRing K>
constexpr auto check_cross_TQF(const Triple<K>& Q)
{
    auto&& [q1, q2, q3] = Q;
    return sq(q1 + q2 + q3) - 2 * (q1 * q1 + q2 * q2 + q3 * q3) - 4 * q1 * q2 * q3;
}

/*!
 * @brief
 *
 * @tparam K
 * @param S
 * @param q3
 * @return constexpr auto
 */
template<CommutativeRing K>
constexpr auto check_cross_law(const Triple<K>& S, const K& q3)
{
    auto&& [s1, s2, s3] = S;
    return sq(s1 * s2 * q3 - (s1 + s2 + s3) + 2) - 4 * (1 - s1) * (1 - s2) * (1 - s3);
}

} // namespace fun

#endif
