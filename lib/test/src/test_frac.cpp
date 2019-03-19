/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "fractions.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <catch.hpp>
#include <iostream>

using namespace fun;

TEST_CASE("Fraction", "[Fraction]")
{
    using boost::multiprecision::cpp_int;
    static_assert(Integral<cpp_int>);

    cpp_int a = 3, b = 4, c = 5, d = 6;
    cpp_int f = -30, g = 40;
    cpp_int z = 0;
    cpp_int h = -g;

    auto p = Fraction{a, b};
    auto q = 2 * p;

    CHECK(p == Fraction(30, 40));
    CHECK(p != 0);
}

TEST_CASE("Fraction Special Cases", "[Fraction]")
{
    auto p = Fraction{3, 4};
    auto inf = Fraction{1, 0};
    auto nan = Fraction{0, 0};
    auto zero = Fraction{0, 1};

    CHECK(-inf < zero);
    CHECK(zero < inf);
    CHECK(-inf < p);
    CHECK(p < inf);
    CHECK(inf == inf);
    CHECK(-inf < inf);
    CHECK(inf == inf * p);
    CHECK(inf == inf * inf);
    CHECK(inf == p / zero);
    CHECK(inf == inf / zero);
    CHECK(nan == nan);
    CHECK(nan == inf * zero);
    CHECK(nan == -inf * zero);
    CHECK(nan == inf / inf);
    CHECK(nan == nan * zero);
    CHECK(nan == nan * nan);
    CHECK(inf == inf + inf);
    CHECK(nan == inf - inf);
    // CHECK( inf + p == nan ); // ???
    // CHECK( -inf + p == nan ); // ???
}
