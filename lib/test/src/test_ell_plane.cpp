/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include <catch.hpp>
#include "ell_plane.hpp"
#include "pg_point.hpp"
#include "pg_line.hpp"
#include <iostream>


TEST_CASE( "Ell plane", "[ell_plane]" ) {
    using namespace fun;
    using namespace fun::ELL;

    using P = pg_point<int>;
    using L = pg_line<int>;

    auto a1 = pg_point(1, 2,  3);
    auto a2 = pg_point(4, -5, 6);
    auto a3 = pg_point(-7, 8, 9);
    REQUIRE( dual<L>(dual<P>(a1)) == a1 );

    auto l1 = join(a2, a3);
    auto l2 = join(a1, a3);
    auto l3 = join(a1, a2);
    REQUIRE( dual<P>(dual<L>(l1)) == l1 );

    auto t1 = altitude(a1, l1);
    auto t2 = altitude(a2, l2);
    auto t3 = altitude(a3, l3);
    REQUIRE( is_perpendicular(t1, l1) );

    auto o = orthocenter(a1, a2, a3);
    std::cout << meet(t2, t3) << std::endl;
    std::cout << o << std::endl;

    REQUIRE( o == meet(t2, t3) );
    REQUIRE( a1 == orthocenter(o, a2, a3) );

    auto tau = line_reflect(l1);
    REQUIRE( tau(tau(a1)) == a1 );

    auto q1 = quadrance(a2, a3);
    auto q2 = quadrance(a1, a3);
    auto q3 = quadrance(a1, a2);
    auto s1 = spread(l2, l3);
    auto s2 = spread(l1, l3);
    auto s3 = spread(l1, l2);
    // print(q1/s1, q2/s2, q3/s3)
    REQUIRE( spread(l1, l1) == 0 );
    REQUIRE( quadrance(a1, a1) == 0 );

}

