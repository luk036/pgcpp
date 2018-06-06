/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "ck_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <catch.hpp>

using namespace fun;

Projective_plane { P, L }
class myck : public ck<P, L, myck> {
public:
  L _perp(const P &v) const { return L(-2 * v[0], v[1], -2 * v[2]); }
  P _perp(const L &v) const { return P(-v[0], 2 * v[1], -v[2]); }
};

TEST_CASE("CK plane", "[ck_plane]") {
  auto a1 = pg_point(1, 2, 3);
  auto a2 = pg_point(4, -5, 6);
  auto a3 = pg_point(-7, 8, 9);

  auto [l1, l2, l3] = tri(std::tuple{a1, a2, a3});

  using P = decltype(a1);
  using L = decltype(l1);
  auto geometry = myck<P, L>{};

  REQUIRE(geometry._perp(geometry._perp(a1)) == a1);
  REQUIRE(geometry._perp(geometry._perp(l1)) == l1);

  auto [t1, t2, t3] = geometry.tri_altitude(a1, a2, a3);
  CHECK(geometry.is_perpendicular(t1, l1));

  auto o = geometry.orthocenter(a1, a2, a3);
  CHECK(o == t2 * t3);
  CHECK(a1 == geometry.orthocenter(o, a2, a3));

  auto tau = geometry.reflect(l1);
  CHECK(tau(tau(a1)) == a1);

  auto q1 = geometry.quadrance(a2, a3);
  auto q2 = geometry.quadrance(a1, a3);
  auto q3 = geometry.quadrance(a1, a2);
  auto s1 = geometry.spread(l2, l3);
  auto s2 = geometry.spread(l1, l3);
  auto s3 = geometry.spread(l1, l2);
  // print(q1/s1, q2/s2, q3/s3)
  CHECK(geometry.spread(l1, l1) == 0);
  CHECK(geometry.quadrance(a1, a1) == 0);
}
