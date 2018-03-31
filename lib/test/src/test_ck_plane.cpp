/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include <catch.hpp>
#include "ck_plane.hpp"
#include "pg_point.hpp"
#include "pg_line.hpp"
#include <complex>

using namespace fun;

TEST_CASE( "CK plane", "[ck_plane]" ) {
    auto p = pg_point(1-2j, 3-1j, 2+1j);  // complex number
    auto q = pg_point(-2+1j, 1-3j, -1-1j);
    auto l = p * q;

    REQUIRE( l == q * p );
    REQUIRE( l.incident(p) );
    REQUIRE( l.incident(q) );

    auto pq = plucker(2+1j, p, 3+0j, q);
    REQUIRE( l.incident(pq) );

    auto r = pg_point(2-1j, -2+1j, 1+1j);
    auto s = pg_point(2j, 2-2j, 3+0j);
    auto t = pg_point(2+0j, -2j, 2+0j);

    auto O = meet(join(p, s), join(q, t));
    auto u = plucker(1+0j, O, -1-1j, r);
    check_desargue(p, q, r, s, t, u);
}

