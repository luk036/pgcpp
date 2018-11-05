/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "fraction.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <catch.hpp>
#include <iostream>

using namespace fun;

TEST_CASE("Fraction", "[Fraction]") {
    using boost::multiprecision::cpp_int;
    static_assert(Integral<cpp_int>);

    cpp_int a = 3, b = 4, c = 5, d = 6;
    cpp_int f = -30, g = 40;
    cpp_int z = 0;
    cpp_int h = -g;

    auto p = Fraction{a, b};

    CHECK( p == Fraction(30, 40) );
    CHECK( p != 0 );
}
