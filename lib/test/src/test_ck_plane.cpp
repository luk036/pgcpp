/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "ck_plane.hpp"
#include "persp_plane.hpp"
#include "pg_common.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <catch.hpp>
#include <iostream>

using namespace fun;

template <typename T, typename U>
inline bool ApproxEqual(const T &a, const U &b) {
    return a[0] == Approx(b[0]) && a[1] == Approx(b[1]) && a[2] == Approx(b[2]);
}

template <typename PG> void chk_int(const PG &myck) {
    using P = typename PG::point_t;
    using Line = typename PG::line_t;

    // P a1 = {1556, -2535, 3445};
    // P a2 = {4544, 3540, 6754};
    // P a3 = {-7453, 1344, 2534};
    P a1 = {1, -2, 3};
    P a2 = {4, 3, 6};
    P a3 = {-7, 1, 2};

    auto triangle = Triple<P>{std::move(a1), std::move(a2), std::move(a3)};
    auto trilateral = tri_dual(triangle);
    const auto &[l1, l2, l3] = trilateral;
    CHECK(incident(l1, a2));

    // CHECK(geometry.perp(geometry.perp(a1)) == a1);
    // CHECK(geometry.perp(geometry.perp(l1)) == l1);

    auto [t1, t2, t3] = myck.tri_altitude(triangle);
    CHECK(myck.is_perpendicular(t1, l1));
    CHECK(coincident(t1, t2, t3));

    auto o = myck.orthocenter(triangle);
    CHECK(o == t2 * t3);
    CHECK(a1 == myck.orthocenter(
                    std::tuple{std::move(o), std::move(a2), std::move(a3)}));

    auto tau = myck.reflect(l1);
    CHECK(tau(tau(a1)) == a1);

    CHECK(myck.spread(l1, l1) == 0);
    CHECK(myck.quadrance(a1, a1) == 0);

    std::tuple Q = myck.tri_quadrance(triangle);
    std::tuple S = myck.tri_spread(trilateral);
    CHECK(check_sine_law(Q, S));
    CHECK(check_sine_law(S, Q));
}

template <typename P, typename L = typename P::dual>
requires Projective_plane_prim<P, L> // c++20 concept
struct myck : ck<P, L, myck> {
    constexpr L perp(const P &v) const { return L(-2 * v[0], v[1], -2 * v[2]); }

    constexpr P perp(const L &v) const { return P(-v[0], 2 * v[1], -v[2]); }

    Projective_plane2 { _P }
    constexpr auto measure(const _P &a1, const _P &a2) const {
        auto x = x_ratio(a1, a2, this->perp(a2), this->perp(a1));
        using Q_t = decltype(x);
        return Q_t(1) - x;
    }
};

TEST_CASE("CK plane chk_int", "[ck_plane]") {
    // using boost::multiprecision::cpp_int;
    // namespace mp = boost::multiprecision;
    using boost::multiprecision::cpp_int;

    chk_int(myck<pg_point<cpp_int>>());
    chk_int(myck<pg_line<cpp_int>>());
    chk_int(ellck<pg_point<cpp_int>>());
    chk_int(ellck<pg_line<cpp_int>>());
    chk_int(hyck<pg_point<cpp_int>>());
    chk_int(hyck<pg_line<cpp_int>>());

    pg_point<cpp_int> Ire(0, 1, 1);
    pg_point<cpp_int> Iim(1, 0, 0);
    pg_line<cpp_int> l_inf(0, -1, 1);
    auto P = persp_euclid_plane{Ire, Iim, l_inf};
    chk_int(P);
}

template <typename PG> void chk_float(const PG &myck) {
    using Point = typename PG::point_t;
    using Line = typename PG::line_t;

    Point a1 = {1., 2., 3.};
    Point a2(4., 0., 6.);
    Point a3(-7., 1., 2.);

    auto triangle = std::tuple{std::move(a1), std::move(a2), std::move(a3)};
    auto trilateral = tri_dual(triangle);
    const auto &[l1, l2, l3] = trilateral;
    CHECK(l1.dot(a2) == Approx(0.));

    // CHECK(geometry.perp(geometry.perp(a1)) == a1);
    // CHECK(geometry.perp(geometry.perp(l1)) == l1);

    auto [t1, t2, t3] = myck.tri_altitude(triangle);
    CHECK(l1.dot(myck.perp(t1)) == Approx(0.));
    CHECK(t1.dot(t2 * t3) == Approx(0.));

    std::array<double, 3> zero{0, 0, 0};

    auto o = myck.orthocenter(triangle);
    CHECK(ApproxEqual(cross(o, t2 * t3), zero));
    CHECK(a1 == myck.orthocenter(
                    std::tuple{std::move(o), std::move(a2), std::move(a3)}));

    auto tau = myck.reflect(l1);
    CHECK(ApproxEqual(cross(tau(tau(a1)), a1), zero));

    CHECK(myck.measure(l1, l1) == Approx(0.));
    CHECK(myck.measure(a1, a1) == Approx(0.));

    auto [q1, q2, q3] = myck.tri_quadrance(triangle);
    auto [s1, s2, s3] = myck.tri_spread(trilateral);
    double r1 = q1 * s2 - q2 * s1;
    double r2 = q2 * s3 - q3 * s2;
    CHECK(r1 == Approx(0.));
    CHECK(r2 == Approx(0.));
}

TEST_CASE("CK plane chk_float", "[ck_plane]") {
    chk_float(myck<pg_point<double>>());
    chk_float(myck<pg_line<double>>());

    chk_float(ellck<pg_point<float>>());
    chk_float(ellck<pg_line<float>>());

    chk_float(hyck<pg_point<double>>());
    chk_float(hyck<pg_line<double>>());

    // namespace mp = boost::multiprecision;

    auto Ire = pg_point(0., 1., 1.);
    auto Iim = pg_point(1., 0., 0.);
    auto l_inf = pg_line(0., -1., 1.);
    auto P = persp_euclid_plane{Ire, Iim, l_inf};
    chk_float(P);
}
