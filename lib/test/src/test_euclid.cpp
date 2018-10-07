/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "euclid_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <catch.hpp>
#include <iostream>

using namespace fun;

TEST_CASE("Euclid plane", "[persp_plane]") {
    auto a1 = pg_point(1, 3, 1);
    auto a2 = pg_point(4, 2, 1);
    auto a3 = pg_point(4, -3, 1);

    auto triangle = std::tuple{a1, a2, a3};
    auto trilateral = tri_dual(triangle);
    auto [l1, l2, l3] = trilateral;

    CHECK(!is_parallel(l1, l2));
    CHECK(!is_parallel(l2, l3));

    auto [t1, t2, t3] = tri_altitude(triangle);
    CHECK(is_perpendicular(t1, l1));
    CHECK(coincident(t1, t2, t3));

    auto o = orthocenter(triangle);
    CHECK(o == t2 * t3);
    CHECK(a1 == orthocenter(std::tuple{o, a2, a3}));

    auto tau = reflect(l1);
    CHECK(tau(tau(a1)) == a1);

    CHECK(spread(l1, l1) == 0);
    CHECK(quadrance(a1, a1) == 0);

    std::tuple Q = tri_quadrance(triangle);
    std::tuple S = tri_spread(trilateral);
    CHECK(check_sine_law(Q, S));
    CHECK(check_sine_law(S, Q));

    auto m12 = midpoint(a1, a2);
    auto m23 = midpoint(a2, a3);
    auto m13 = midpoint(a1, a3);

    t1 = a1 * m23;
    t2 = a2 * m13;
    t3 = a3 * m12;
    CHECK(coincident(t1, t2, t3));

    // auto [q1, q2, q3] = Q;
    // auto tqf = sq(q1 + q2 + q3) - 2*(q1*q1 + q2*q2 + q3*q3);
    // CHECK(tqf == Ar(q1, q2, q3));
}
