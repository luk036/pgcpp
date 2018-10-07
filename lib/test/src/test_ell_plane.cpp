/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "ck_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <catch.hpp>
#include <iostream>

using namespace fun;

template <typename PG> void chk_tri(const PG &myck) {
    using Point = typename PG::point_t;
    using Line = typename PG::line_t;

    auto a1 = Point(1, 3, 1);
    auto a2 = Point(4, 2, 1);
    auto a3 = Point(1, 1, -1);

    CHECK(myck._perp(myck._perp(a1)) == a1);

    auto triangle = std::tuple{a1, a2, a3};
    auto trilateral = tri_dual(triangle);
    auto [l1, l2, l3] = trilateral;

    CHECK(myck._perp(myck._perp(l1)) == l1);

    std::tuple Q = myck.tri_quadrance(triangle);
    std::tuple S = myck.tri_spread(trilateral);

    CHECK(check_cross_law(S, std::get<2>(Q)));
    CHECK(check_cross_law(Q, std::get<2>(S)));

    a3 = plucker(2, a1, 3, a2);
    auto collin = std::tuple{a1, a2, a3};
    Q = myck.tri_quadrance(collin);
    CHECK(check_cross_TQF(Q) == 0);
}

TEST_CASE("Elliptic/Hyperbolic plane", "[ell_plane]") {
    chk_tri(ellck<pg_point<int>, pg_line<int>>());
    chk_tri(ellck<pg_line<int>, pg_point<int>>());
    chk_tri(hyck<pg_point<int>, pg_line<int>>());
    chk_tri(hyck<pg_line<int>, pg_point<int>>());
}
