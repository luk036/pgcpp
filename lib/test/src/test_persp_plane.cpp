/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "euclid_plane.hpp" // import Ar
#include "persp_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <catch.hpp>
#include <iostream>

using namespace fun;

template <typename PG> void chk_degenerate_int(const PG &myck) {
    using Point = typename PG::point_t;
    using Line = typename PG::line_t;

    auto a1 = Point(-1, 0, 3);
    auto a2 = Point(4, -2, 1);
    auto a3 = Point(3, -1, 1);

    auto triangle = std::tuple{a1, a2, a3};
    auto trilateral = tri_dual(triangle);
    auto [l1, l2, l3] = trilateral;

    CHECK(!myck.is_parallel(l1, l2));
    CHECK(!myck.is_parallel(l2, l3));

    auto m12 = myck.midpoint(a1, a2);
    auto m23 = myck.midpoint(a2, a3);
    auto m13 = myck.midpoint(a1, a3);

    auto t1 = a1 * m23;
    auto t2 = a2 * m13;
    auto t3 = a3 * m12;
    CHECK(coincident(t1, t2, t3));

    auto [q1, q2, q3] = myck.tri_quadrance(triangle);
    auto tqf = sq(q1 + q2 + q3) - (q1 * q1 + q2 * q2 + q3 * q3) * 2;
    CHECK(tqf == Ar(q1, q2, q3));
}

TEST_CASE("Perspective Euclid plane", "[persp_plane]") {
    using boost::multiprecision::cpp_int;

    auto Ire = pg_point<cpp_int>(0, 1, 1);
    auto Iim = pg_point<cpp_int>(1, 0, 0);
    auto l_inf = pg_line<cpp_int>(0, -1, 1);
    auto P = persp_euclid_plane{Ire, Iim, l_inf};
    chk_degenerate_int(P);
}

template <typename PG> void chk_degenerate_float(const PG &myck) {
    using Point = typename PG::point_t;
    using Line = typename PG::line_t;

    auto a1 = Point(-1., 0., 3.);
    auto a2 = Point(4., -2., 1.);
    auto a3 = Point(3., -1., 1.);

    auto triangle = std::tuple{a1, a2, a3};
    auto trilateral = tri_dual(triangle);
    auto [l1, l2, l3] = trilateral;

    CHECK(myck.l_infty().dot(l1 * l2) != Approx(0.0));
    CHECK(myck.l_infty().dot(l2 * l3) != Approx(0.0));

    //auto m12 = myck.midpoint(a1, a2);
    //auto m23 = myck.midpoint(a2, a3);
    //auto m13 = myck.midpoint(a1, a3);
    auto [m12, m23, m13] = myck.tri_midpoint(triangle);
    auto t1 = a1 * m23;
    auto t2 = a2 * m13;
    auto t3 = a3 * m12;
    CHECK(t1.dot(t2 * t3) == Approx(0.0));

    auto [q1, q2, q3] = myck.tri_quadrance(triangle);
    auto tqf = sq(q1 + q2 + q3) - 2 * (q1 * q1 + q2 * q2 + q3 * q3);
    CHECK(Approx(tqf) == Ar(q1, q2, q3));
}

TEST_CASE("Perspective Euclid plane (floating point)", "[persp_plane]") {
    auto Ire = pg_point(0., 1., 1.);
    auto Iim = pg_point(1., 0., 0.);
    auto l_inf = pg_line(0., -1., 1.);
    auto P = persp_euclid_plane{Ire, Iim, l_inf};
    chk_degenerate_float(P);
}
