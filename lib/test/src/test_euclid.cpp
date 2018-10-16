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

TEST_CASE("Euclid plane", "[euclid_plane]") {
    using boost::multiprecision::cpp_int;
    using K = cpp_int;

    auto a1 = pg_point<cpp_int>(1, 3, 1);
    auto a2 = pg_point<cpp_int>(4, 2, 1);
    auto a3 = pg_point<cpp_int>(4, -3, 1);

    using P = decltype(a1);

    auto triangle = std::tuple{a1, a2, a3};
    auto trilateral = tri_dual(triangle);
    auto [l1, l2, l3] = trilateral;

    using P = decltype(a1);
    using L = decltype(l1);

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

    auto [q1, q2, q3] = Q;
    auto [s1, s2, s3] = S;

    auto tqf = sq(q1 + q2 + q3) - (q1 * q1 + q2 * q2 + q3 * q3) * 2;
    CHECK(tqf == Ar(q1, q2, q3));

    // auto c3 = sq(q1 + q2 - q3) / (4*q1*q2);
    // CHECK( c3 + s3 == 1 ); // get the same

    auto tsf =
        sq(s1 + s2 + s3) - 2 * (s1 * s1 + s2 * s2 + s3 * s3) - 4 * s1 * s2 * s3;
    CHECK(tsf == 0);

    a3 = plucker(3, a1, 4, a2);
    q1 = quadrance(a2, a3);
    q2 = quadrance(a1, a3);
    q3 = quadrance(a1, a2);
    tqf = sq(q1 + q2 + q3) - 2 * (q1 * q1 + q2 * q2 + q3 * q3); // get 0
    CHECK(tqf == 0);

    a1 = uc_point<P>(1, 0);
    a2 = uc_point<P>(3, 4);
    a3 = uc_point<P>(-1, 2);
    auto a4 = uc_point<P>(0, 1);
    auto q12 = quadrance(a1, a2);
    auto q23 = quadrance(a2, a3);
    auto q34 = quadrance(a3, a4);
    auto q14 = quadrance(a1, a4);
    auto q24 = quadrance(a2, a4);
    auto q13 = quadrance(a1, a3);
    // print(q12, q23, q34, q14, q24, q13)
    auto t = Ar(q12 * q34, q23 * q14, q13 * q24);
    // t = sympy.simplify(t)
    CHECK(t == 0);
    bool b = Ptolemy(std::tuple{q12, q23, q34, q14, q24, q13});
    CHECK(b == true);
}
