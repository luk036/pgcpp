/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "ck_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <catch.hpp>
#include <iostream>

using namespace fun;

template <typename T, typename U>
inline bool ApproxEqual(const T &a, const U &b) {
    return a[0] == Approx(b[0])
        && a[1] == Approx(b[1])
        && a[2] == Approx(b[2]); 
}

template <typename PG> void chk_tri_int(const PG &myck) {
    using Point = typename PG::point_t;
    using Line = typename PG::line_t;

    auto a1 = Point(1, 3, 1);
    auto a2 = Point(4, 2, 1);
    auto a3 = Point(1, 1, -1);

    CHECK(myck.perp(myck.perp(a1)) == a1);

    auto triangle = std::tuple{a1, a2, a3};
    auto trilateral = tri_dual(triangle);
    auto [l1, l2, l3] = trilateral;

    CHECK(myck.perp(myck.perp(l1)) == l1);

    std::tuple Q = myck.tri_quadrance(triangle);
    std::tuple S = myck.tri_spread(trilateral);

    CHECK(check_cross_law(S, std::get<2>(Q)) == 0);
    CHECK(check_cross_law(Q, std::get<2>(S)) == 0);

    a3 = plucker(2, a1, 3, a2);
    auto collin = std::tuple{a1, a2, a3};
    Q = myck.tri_quadrance(collin);
    CHECK(check_cross_TQF(Q) == 0);
}

TEST_CASE("Elliptic/Hyperbolic plane", "[ell_plane]") {
    chk_tri_int(ellck<pg_point<int>, pg_line<int>>());
    chk_tri_int(ellck<pg_line<int>, pg_point<int>>());
    chk_tri_int(hyck<pg_point<int>, pg_line<int>>());
    chk_tri_int(hyck<pg_line<int>, pg_point<int>>());
}


template <typename PG> void chk_tri_float(const PG &myck) {
    using Point = typename PG::point_t;
    using Line = typename PG::line_t;

    auto a1 = Point(1., 3., 1.);
    auto a2 = Point(4., 2., 1.);
    auto a3 = Point(1., 1., -1.);

    std::array<double, 3> zero{0, 0, 0};

    CHECK(ApproxEqual(cross(myck.perp(myck.perp(a1)), a1), zero));

    auto triangle = std::tuple{a1, a2, a3};
    auto trilateral = tri_dual(triangle);
    auto [l1, l2, l3] = trilateral;

    CHECK(ApproxEqual(cross(myck.perp(myck.perp(l1)), l1), zero));

    std::tuple Q = myck.tri_quadrance(triangle);
    std::tuple S = myck.tri_spread(trilateral);

    CHECK(check_cross_law(S, std::get<2>(Q)) == Approx(0));
    CHECK(check_cross_law(Q, std::get<2>(S)) == Approx(0));

    a3 = plucker(2, a1, 3, a2);
    auto collin = std::tuple{a1, a2, a3};
    Q = myck.tri_quadrance(collin);
    CHECK(check_cross_TQF(Q) == Approx(0));
}

TEST_CASE("Elliptic/Hyperbolic plane (double)", "[ell_plane]") {
    chk_tri_float(ellck<pg_point<double>, pg_line<double>>());
    chk_tri_float(ellck<pg_line<double>, pg_point<double>>());
    chk_tri_float(hyck<pg_point<double>, pg_line<double>>());
    chk_tri_float(hyck<pg_line<double>, pg_point<double>>());
}
