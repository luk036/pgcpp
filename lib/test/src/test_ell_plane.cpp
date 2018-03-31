/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include <catch.hpp>
#include "ell_plane.hpp"
#include "pg_point.hpp"
#include "pg_line.hpp"

using namespace fun;
using namespace fun::CK;

TEST_CASE( "Ell plane", "[ell_plane]" ) {
    auto a1 = pg_point(1, 2,  3);
    auto a2 = pg_point(4, -5, 6);
    auto a3 = pg_point(-7, 8, 9);
    REQUIRE( ~(~a1) == a1 );

    auto l1 = join(a2, a3);
    auto l2 = join(a1, a3);
    auto l3 = join(a1, a2);
    REQUIRE( ~(~l1) == l1 );

    auto t1 = CK::altitude(a1, l1);
    auto t2 = CK::altitude(a2, l2);
    auto t3 = CK::altitude(a3, l3);
    REQUIRE( CK::is_perpendicular(t1, l1) );

    auto o = CK::orthocenter(a1, a2, a3);
    REQUIRE( o == meet(t2, t3) );
    REQUIRE( a1 == CK::orthocenter(o, a2, a3) );

    auto tau = CK::line_reflect(l1);
    REQUIRE( tau(tau(a1)) == a1 );

    auto q1 = CK::quadrance(a2, a3);
    auto q2 = CK::quadrance(a1, a3);
    auto q3 = CK::quadrance(a1, a2);
    auto s1 = CK::spread(l2, l3);
    auto s2 = CK::spread(l1, l3);
    auto s3 = CK::spread(l1, l2);
    // print(q1/s1, q2/s2, q3/s3)
    REQUIRE( CK::spread(l1, l1) == 0 );
    REQUIRE( CK::quadrance(a1, a1) == 0 );

}

