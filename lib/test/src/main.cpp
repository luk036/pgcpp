/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include "proj_plane.hpp"
#include "pg_point.hpp"
#include "pg_line.hpp"

using namespace fun;

TEST_CASE( "Projective Point", "[proj_plane]" ) {
    auto p = pg_point<int>{{1, 3, 1}};
    auto q = pg_point<int>{{4, 2, 1}};
    // auto a3 = pg_point<int>{{1, 1, -1}};

    REQUIRE( p.incident({p, q}) );
}

TEST_CASE( "Projective Point", "[proj_plane]" ) {
    auto l = pg_line<int>{{1, 3, 1}};
    auto m = pg_line<int>{{4, 2, 1}};

    REQUIRE( l.incident({l, m}) );
}


// int main(int argc, char* argv[]) {
//   //using namespace ModernCppCI;

//   auto result = Catch::Session().run(argc, argv);

//   return (result < 0xff ? result : 0xff);
// }