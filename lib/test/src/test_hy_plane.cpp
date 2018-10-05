/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "ck_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <catch.hpp>
#include <iostream>

using namespace fun;

TEST_CASE("Hyperbolic plane", "[hy_plane]") {

    // using P = pg_point<int>;
    // using L = pg_line<int>;
    // using dP = decltype(dualP);
    // using dL = decltype(dualL);

    auto a1 = pg_point(1, 2, 3);
    auto a2 = pg_point(4, -5, 6);
    auto a3 = pg_point(-7, 8, 9);

    auto l1 = join(a2, a3);
    auto l2 = join(a1, a3);
    auto l3 = join(a1, a2);

    using P = decltype(a1);
    using L = decltype(l1);
    auto geometry = hyck<P, L>();

    auto t1 = geometry.altitude(a1, l1);
    auto t2 = geometry.altitude(a2, l2);
    auto t3 = geometry.altitude(a3, l3);
    CHECK(geometry.is_perpendicular(t1, l1));

    auto triangle = std::tuple{a1, a2, a3};
    auto o = geometry.orthocenter(triangle);
    CHECK(o == meet(t2, t3));
    CHECK(a1 == geometry.orthocenter(std::tuple{o, a2, a3}));

    auto tau = geometry.reflect(l1);
    CHECK(tau(tau(a1)) == a1);

    auto q1 = geometry.quadrance(a2, a3);
    auto q2 = geometry.quadrance(a1, a3);
    auto q3 = geometry.quadrance(a1, a2);
    auto s1 = geometry.spread(l2, l3);
    auto s2 = geometry.spread(l1, l3);
    auto s3 = geometry.spread(l1, l2);
    // print(q1/s1, q2/s2, q3/s3)
    CHECK(geometry.spread(l1, l1) == 0);
    CHECK(geometry.quadrance(a1, a1) == 0);
}
