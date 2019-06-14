/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "pgcpp/ck_plane.hpp"
#include "pgcpp/persp_plane.hpp"
#include "pgcpp/pg_common.hpp"
#include "pgcpp/pg_line.hpp"
#include "pgcpp/pg_point.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <catch.hpp>
#include <iostream>

using namespace fun;

/**
 * @brief 
 * 
 * @param a 
 * @param b 
 * @return true 
 * @return false 
 */
inline auto ApproxEqual(const auto &a, const auto &b) -> bool {
    return a[0] == Approx(b[0]) && a[1] == Approx(b[1]) && a[2] == Approx(b[2]);
}

template <typename PG> void chk_ck(const PG &myck) {
    using P = typename PG::point_t;
    // using Line = typename PG::line_t;
    using K = Value_type<P>;

    // P a1 = {1556, -2535, 3445};
    // P a2 = {4544, 3540, 6754};
    // P a3 = {-7453, 1344, 2534};
    auto a1 = P{1, -2, 3};
    auto a2 = P{4, 0, 6};
    auto a3 = P{-7, 1, 2};
    auto a4 = P{3, 0, 2};
    auto zero = std::array<K, 3>{0, 0, 0};

    auto triangle = std::tuple{std::move(a1), std::move(a2), std::move(a3)};
    auto trilateral = tri_dual(triangle);
    auto &&[l1, l2, l3] = trilateral;
    auto &&[t1, t2, t3] = myck.tri_altitude(triangle);
    auto o = myck.orthocenter(triangle);
    auto tau = myck.reflect(l1);
    auto Q = std::tuple{myck.tri_quadrance(triangle)};
    auto S = std::tuple{myck.tri_spread(trilateral)};

    if constexpr (Integral<K>) {
        CHECK(incident(l1, a2));
        CHECK(myck.is_perpendicular(t1, l1));
        CHECK(coincident(t1, t2, t3));
        CHECK(o == t2 * t3);
        CHECK(a1 == myck.orthocenter(std::tuple{std::move(o), std::move(a2),
                                                std::move(a3)}));
        CHECK(tau(tau(a4)) == a4);
        CHECK(myck.spread(l1, l1) == 0);
        CHECK(myck.quadrance(a1, a1) == 0);
        CHECK(check_sine_law(Q, S));
        CHECK(check_sine_law(S, Q));
    } else {
        CHECK(l1.dot(a2) == Approx(0));
        CHECK(l1.dot(myck.perp(t1)) == Approx(0));
        CHECK(t1.dot(t2 * t3) == Approx(0));
        CHECK(ApproxEqual(cross(o, t2 * t3), zero));
        auto o2 = myck.orthocenter(
            std::tuple{std::move(o), std::move(a2), std::move(a3)});
        CHECK(ApproxEqual(cross(a1, o2), zero));
        CHECK(ApproxEqual(cross(tau(tau(a4)), a4), zero));
        CHECK(myck.measure(l1, l1) == Approx(0));
        CHECK(myck.measure(a1, a1) == Approx(0));
        auto &&[q1, q2, q3] = Q;
        auto &&[s1, s2, s3] = S;
        auto r1 = q1 * s2 - q2 * s1;
        auto r2 = q2 * s3 - q3 * s2;
        CHECK(r1 == Approx(0));
        CHECK(r2 == Approx(0));
    }
}

template <typename P, typename L = typename P::dual>
requires Projective_plane_prim<P, L> // c++20 concept
    struct myck : ck<P, L, myck> {
    constexpr L perp(const P &v) const { return L(-2 * v[0], v[1], -2 * v[2]); }
    constexpr P perp(const L &v) const { return P(-v[0], 2 * v[1], -v[2]); }

    template <Projective_plane2 _P>
    constexpr auto measure(const _P &a1, const _P &a2) const {
        auto x = x_ratio(a1, a2, this->perp(a2), this->perp(a1));
        // using Q_t = decltype(x);
        return 1 - x;
    }
};

TEST_CASE("CK plane chk_ck (int)", "[ck_plane]") {
    // using boost::multiprecision::cpp_int;
    // namespace mp = boost::multiprecision;
    using boost::multiprecision::cpp_int;

    chk_ck(myck<pg_point<cpp_int>>());
    chk_ck(myck<pg_line<cpp_int>>());
    chk_ck(ellck<pg_point<cpp_int>>());
    chk_ck(ellck<pg_line<cpp_int>>());
    chk_ck(hyck<pg_point<cpp_int>>());
    chk_ck(hyck<pg_line<cpp_int>>());

    auto Ire = pg_point<cpp_int>{0, 1, 1};
    auto Iim = pg_point<cpp_int>{1, 0, 0};
    auto l_inf = pg_line<cpp_int>{0, -1, 1};
    auto P =
        persp_euclid_plane{std::move(Ire), std::move(Iim), std::move(l_inf)};
    chk_ck(P);
}

TEST_CASE("CK plane chk_ck (float)", "[ck_plane]") {
    chk_ck(myck<pg_point<double>>());
    chk_ck(myck<pg_line<double>>());
    chk_ck(ellck<pg_point<float>>());
    chk_ck(ellck<pg_line<float>>());
    chk_ck(hyck<pg_point<double>>());
    chk_ck(hyck<pg_line<double>>());

    auto Ire = pg_point{0., 1., 1.};
    auto Iim = pg_point{1., 0., 0.};
    auto l_inf = pg_line{0., -1., 1.};
    auto P =
        persp_euclid_plane{std::move(Ire), std::move(Iim), std::move(l_inf)};
    chk_ck(P);
}
