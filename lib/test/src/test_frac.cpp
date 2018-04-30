/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "fraction.hpp"
#include <catch.hpp>

using namespace fun;

TEST_CASE("Fraction", "[fraction]") {
  auto p = Fraction(3, 4);
  auto q = Fraction(5, 6);

  REQUIRE(p == Fraction(30, 40));
}
