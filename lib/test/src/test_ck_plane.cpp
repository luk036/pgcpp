/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "ck_plane.hpp"
#include "persp_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <catch.hpp>
#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>

using namespace fun;

template <typename PG> void chk_int(const PG &myck) {
    using Point = typename PG::point_t;
    using Line = typename PG::line_t;

    Point a1(1, 2, 3);
    Point a2(4, 0, 6);
    Point a3(-7, 1, 2);

    auto triangle = std::tuple{a1, a2, a3};
    auto trilateral = tri_dual(triangle);
    auto [l1, l2, l3] = trilateral;
    CHECK(incident(l1, a2));

    // CHECK(geometry._perp(geometry._perp(a1)) == a1);
    // CHECK(geometry._perp(geometry._perp(l1)) == l1);

    auto [t1, t2, t3] = myck.tri_altitude(triangle);
    CHECK(myck.is_perpendicular(t1, l1));
    CHECK(coincident(t1, t2, t3));

    auto o = myck.orthocenter(triangle);
    CHECK(o == t2 * t3);
    CHECK(a1 == myck.orthocenter(std::tuple{o, a2, a3}));

    auto tau = myck.reflect(l1);
    CHECK(tau(tau(a1)) == a1);

    CHECK(myck.spread(l1, l1) == 0);
    CHECK(myck.quadrance(a1, a1) == 0);

    std::tuple Q = myck.tri_quadrance(triangle);
    std::tuple S = myck.tri_spread(trilateral);
    CHECK(check_sine_law(Q, S));
    CHECK(check_sine_law(S, Q));
}

Projective_plane { P, L }
struct myck : ck<P, L, myck> {
    constexpr L _perp(const P &v) const {
        return L(-2 * v[0], v[1], -2 * v[2]);
    }

    constexpr P _perp(const L &v) const { return P(-v[0], 2 * v[1], -v[2]); }

    Projective_plane2 { _P }
    constexpr auto measure(const _P &a1, const _P &a2) const {
        using K = Value_type<_P>;
        return K(1) - x_ratio(a1, a2, this->_perp(a2), this->_perp(a1));
    }
};

TEST_CASE("CK plane chk_int", "[ck_plane]") {
    chk_int(myck<pg_point<int>, pg_line<int>>());
    chk_int(myck<pg_line<int>, pg_point<int>>());

    chk_int(ellck<pg_point<int>, pg_line<int>>());
    chk_int(ellck<pg_line<int>, pg_point<int>>());

    chk_int(hyck<pg_point<int>, pg_line<int>>());
    chk_int(hyck<pg_line<int>, pg_point<int>>());

    namespace mp = boost::multiprecision;

    auto Ire = pg_point(0, 1, 1);
    auto Iim = pg_point(1, 0, 0);
    auto l_inf = pg_line(0, -1, 1);
    auto P = persp_euclid_plane{Ire, Iim, l_inf};
    chk_int(P);
}
