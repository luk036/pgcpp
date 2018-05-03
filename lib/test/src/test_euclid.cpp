/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "euclid_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <catch.hpp>
#include <iostream>

using namespace fun;

TEST_CASE("Euclid plane", "[persp_plane]") {

  // using P = pg_point<int>;
  // using L = pg_line<int>;
  // using dP = decltype(dualP);
  // using dL = decltype(dualL);

  auto a1 = pg_point(1, 2, 3);
  auto a2 = pg_point(4, -5, 6);
  auto a3 = pg_point(-7, 8, 9);

  auto l1 = join(a2, a3);
  auto l2 = join(a1, a3);
  auto l3 = join(a1, a2);

  using P = decltype(a1);
  using L = decltype(l1);

  auto t1 = altitude(a1, l1);
  auto t2 = altitude(a2, l2);
  auto t3 = altitude(a3, l3);
  REQUIRE(is_perpendicular(t1, l1));

  auto o = orthocenter(a1, a2, a3);
  REQUIRE(o == meet(t2, t3));
  REQUIRE(a1 == orthocenter(o, a2, a3));

  auto tau = reflect(l1);
  REQUIRE(tau(tau(a1)) == a1);

  auto q1 = quadrance(a2, a3);
  auto q2 = quadrance(a1, a3);
  auto q3 = quadrance(a1, a2);
  auto s1 = spread(l2, l3);
  auto s2 = spread(l1, l3);
  auto s3 = spread(l1, l2);
  // print(q1/s1, q2/s2, q3/s3)
  REQUIRE(spread(l1, l1) == 0);
  REQUIRE(quadrance(a1, a1) == 0);
}
