/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "ck_plane.hpp"
#include "persp_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <catch.hpp>
#include <iostream>
#include "pg_common.hpp"
#include <boost/multiprecision/cpp_int.hpp>

using namespace fun;

template <typename T, typename U>
inline bool ApproxEqual(const T &a, const U &b) {
    return a[0] == Approx(b[0])
        && a[1] == Approx(b[1])
        && a[2] == Approx(b[2]); 
}

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

    // CHECK(geometry.perp(geometry.perp(a1)) == a1);
    // CHECK(geometry.perp(geometry.perp(l1)) == l1);

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
    constexpr L perp(const P &v) const {
        return L(-2 * v[0], v[1], -2 * v[2]);
    }

    constexpr P perp(const L &v) const { return P(-v[0], 2 * v[1], -v[2]); }

    Projective_plane2 { _P }
    constexpr auto measure(const _P &a1, const _P &a2) const {
        using K = Value_type<P>;
        return K(1) - x_ratio(a1, a2, this->perp(a2), this->perp(a1));
    }
};

TEST_CASE("CK plane chk_int", "[ck_plane]") {
    using boost::multiprecision::cpp_int;
    namespace mp = boost::multiprecision;

    chk_int(myck<pg_point<int>, pg_line<int>>());
    chk_int(myck<pg_line<int>, pg_point<int>>());
    chk_int(ellck<pg_point<int>, pg_line<int>>());
    chk_int(ellck<pg_line<int>, pg_point<int>>());
    chk_int(hyck<pg_point<int>, pg_line<int>>());
    chk_int(hyck<pg_line<int>, pg_point<int>>());

    namespace mp = boost::multiprecision;

    auto Ire = pg_point<int>(0, 1, 1);
    auto Iim = pg_point<int>(1, 0, 0);
    auto l_inf= pg_line<int>(0, -1, 1);
    auto P = persp_euclid_plane{Ire, Iim, l_inf};
    chk_int(P);
}



template <typename PG> void chk_float(const PG &myck) {
    using Point = typename PG::point_t;
    using Line = typename PG::line_t;

    Point a1(1., 2., 3.);
    Point a2(4., 0., 6.);
    Point a3(-7., 1., 2.);
    
    auto triangle = std::tuple{a1, a2, a3};
    auto trilateral = tri_dual(triangle);
    auto [l1, l2, l3] = trilateral;
    CHECK(l1.dot(a2) == Approx(0.));

    // CHECK(geometry.perp(geometry.perp(a1)) == a1);
    // CHECK(geometry.perp(geometry.perp(l1)) == l1);

    auto [t1, t2, t3] = myck.tri_altitude(triangle);
    CHECK(l1.dot(myck.perp(t1)) == Approx(0.));
    CHECK(t1.dot(t2 * t3) == Approx(0.));

    std::array<double, 3> zero{0, 0, 0};

    auto o = myck.orthocenter(triangle);
    CHECK(ApproxEqual(cross(o, t2 * t3), zero));
    CHECK(a1 == myck.orthocenter(std::tuple{o, a2, a3}));

    auto tau = myck.reflect(l1);
    CHECK(ApproxEqual(cross(tau(tau(a1)), a1), zero));

    CHECK(myck.measure(l1, l1) == Approx(0.));
    CHECK(myck.measure(a1, a1) == Approx(0.));

    auto [q1, q2, q3] = myck.tri_quadrance(triangle);
    auto [s1, s2, s3] = myck.tri_spread(trilateral);
    // CHECK( q1 * s2 == Approx(q2 * s1) );
    // CHECK( q2 * s3 == Approx(q3 * s2) );
    // CHECK(check_sine_law(Q, S));
    // CHECK(check_sine_law(S, Q));
}

TEST_CASE("CK plane chk_float", "[ck_plane]") {
    chk_float(myck<pg_point<double>, pg_line<double>>());
    chk_float(myck<pg_line<double>, pg_point<double>>());

    chk_float(ellck<pg_point<double>, pg_line<double>>());
    chk_float(ellck<pg_line<double>, pg_point<double>>());

    chk_float(hyck<pg_point<double>, pg_line<double>>());
    chk_float(hyck<pg_line<double>, pg_point<double>>());

    namespace mp = boost::multiprecision;

    auto Ire = pg_point(0., 1., 1.);
    auto Iim = pg_point(1., 0., 0.);
    auto l_inf = pg_line(0., -1., 1.);
    auto P = persp_euclid_plane{Ire, Iim, l_inf};
    chk_float(P);
}
