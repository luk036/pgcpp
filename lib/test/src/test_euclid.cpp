/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "ck_plane.hpp"
#include "euclid_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <catch.hpp>
#include <iostream>

using namespace fun;

TEST_CASE("Euclid plane", "[euclid_plane]")
{
    using boost::multiprecision::cpp_int;
    using K = cpp_int;

    auto a1 = pg_point<cpp_int>{1, 3, 1};
    auto a2 = pg_point<cpp_int>{4, 2, 1};
    auto a3 = pg_point<cpp_int>{4, -3, 1};

    using P = decltype(a1);

    auto triangle = std::tuple{std::move(a1), std::move(a2), std::move(a3)};
    auto trilateral = tri_dual(triangle);
    auto &[l1, l2, l3] = trilateral;

    using P = decltype(a1);
    using L = decltype(l1);

    CHECK(!is_parallel(l1, l2));
    CHECK(!is_parallel(l2, l3));

    auto [t1, t2, t3] = tri_altitude(triangle);
    CHECK(is_perpendicular(t1, l1));
    CHECK(coincident(t1, t2, t3));

    auto o = orthocenter(triangle);
    CHECK(o == t2 * t3);
    CHECK(a1 ==
          orthocenter(std::tuple{std::move(o), std::move(a2), std::move(a3)}));

    auto tau = reflect(l1);
    CHECK(tau(tau(a1)) == a1);

    CHECK(spread(l1, l1) == 0);
    CHECK(quadrance(a1, a1) == 0);

    std::tuple Q = tri_quadrance(triangle);
    std::tuple S = tri_spread(trilateral);
    CHECK(check_sine_law(Q, S));
    CHECK(check_sine_law(S, Q));

    // auto m12 = midpoint(a1, a2);
    // auto m23 = midpoint(a2, a3);
    // auto m13 = midpoint(a1, a3);
    auto [m12, m23, m13] = tri_midpoint(triangle);
    auto tm1 = a1 * m23;
    auto tm2 = a2 * m13;
    auto tm3 = a3 * m12;
    CHECK(coincident(tm1, tm2, tm3));

    auto &[q1, q2, q3] = Q;
    auto &[s1, s2, s3] = S;

    auto tqf = sq(q1 + q2 + q3) - 2 * (q1 * q1 + q2 * q2 + q3 * q3);
    CHECK(tqf == Ar(q1, q2, q3));

    // auto c3 = sq(q1 + q2 - q3) / (4*q1*q2);
    // CHECK( c3 + s3 == 1 ); // get the same

    auto tsf = sq(s1 + s2 + s3) - 2 * (s1 * s1 + s2 * s2 + s3 * s3) - 4 * s1 * s2 * s3;
    CHECK(tsf == 0);

    auto a3p = plucker(3, a1, 4, a2);
    auto q1p = quadrance(a2, a3p);
    auto q2p = quadrance(a1, a3p);
    auto q3p = quadrance(a1, a2);
    auto tqf2 = sq(q1p + q2p + q3p) - 2 * (q1p * q1p + q2p * q2p + q3p * q3p); // get 0
    CHECK(tqf2 == 0);

    auto u1 = uc_point<P>(1, 0);
    auto u2 = uc_point<P>(3, 4);
    auto u3 = uc_point<P>(-1, 2);
    auto u4 = uc_point<P>(0, 1);
    auto q12 = quadrance(u1, u2);
    auto q23 = quadrance(u2, u3);
    auto q34 = quadrance(u3, u4);
    auto q14 = quadrance(u1, u4);
    auto q24 = quadrance(u2, u4);
    auto q13 = quadrance(u1, u3);
    // print(q12, q23, q34, q14, q24, q13)
    auto t = Ar(q12 * q34, q23 * q14, q13 * q24);
    // t = sympy.simplify(t)
    CHECK(t == 0);
    auto okay = Ptolemy(std::tuple{std::move(q12), std::move(q23), std::move(q34), std::move(q14), std::move(q24), std::move(q13)});
    CHECK(okay);
}
