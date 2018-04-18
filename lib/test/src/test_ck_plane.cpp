/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include <catch.hpp>
#include "ck_plane.hpp"
#include "pg_point.hpp"
#include "pg_line.hpp"

using namespace fun;

TEST_CASE( "CK plane", "[ck_plane]" ) {
    auto dualP = [](const pg_point<int>& v) {
        return pg_line<int>(-2*v[0], v[1], -2*v[2]);
        };
    auto dualL = [](const pg_line<int>& v) {
        return pg_point<int>(-v[0], 2*v[1], -v[2]);
        };

    // using P = pg_point<int>;
    // using L = pg_line<int>;
    // using dP = decltype(dualP);
    // using dL = decltype(dualL);

    auto a1 = pg_point(1, 2,  3);
    auto a2 = pg_point(4, -5, 6);
    auto a3 = pg_point(-7, 8, 9);
    REQUIRE( dualL(dualP(a1)) == a1 );

    auto l1 = join(a2, a3);
    auto l2 = join(a1, a3);
    auto l3 = join(a1, a2);
    REQUIRE( dualP(dualL(l1)) == l1 );

    auto t1 = altitude(a1, l1, dualL);
    auto t2 = altitude(a2, l2, dualL);
    auto t3 = altitude(a3, l3, dualL);
    REQUIRE( is_perpendicular(t1, l1, dualL) );

    auto o = orthocenter(a1, a2, a3, dualL);
    REQUIRE( o == meet(t2, t3) );
    REQUIRE( a1 == orthocenter(o, a2, a3, dualL) );

    auto tau = line_reflect(l1, dualL);
    REQUIRE( tau(tau(a1)) == a1 );

    auto q1 = quadrance(a2, a3, dualP);
    auto q2 = quadrance(a1, a3, dualP);
    auto q3 = quadrance(a1, a2, dualP);
    auto s1 = spread(l2, l3, dualL);
    auto s2 = spread(l1, l3, dualL);
    auto s3 = spread(l1, l2, dualL);
    // print(q1/s1, q2/s2, q3/s3)
    REQUIRE( spread(l1, l1, dualL) == 0 );
    REQUIRE( quadrance(a1, a1, dualP) == 0 );

}

