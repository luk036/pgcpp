/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#define CATCH_CONFIG_MAIN
#include "pg_line.hpp"
#include "pg_point.hpp"
#include "proj_plane.hpp"
#include <catch.hpp>
#include <complex>

using namespace fun;

TEST_CASE("Projective Point", "[proj_plane]") {
  auto p = pg_point(1 - 2j, 3 - 1j, 2 + 1j); // complex number
  auto q = pg_point(-2 + 1j, 1 - 3j, -1 - 1j);
  auto l = p * q;

  REQUIRE(l == q * p);
  REQUIRE(l.incident(p));
  REQUIRE(l.incident(q));

  auto pq = plucker(2 + 1j, p, 3 + 0j, q);
  REQUIRE(l.incident(pq));

  auto r = pg_point(2 - 1j, -2 + 1j, 1 + 1j);
  auto s = pg_point(2j, 2 - 2j, 3 + 0j);
  auto t = pg_point(2 + 0j, -2j, 2 + 0j);

  auto O = meet(join(p, s), join(q, t));
  auto u = plucker(1 + 0j, O, -1 - 1j, r);
  check_desargue(p, q, r, s, t, u);
}

TEST_CASE("Projective Line", "[proj_plane]") {
  auto l = pg_line(1 - 2j, 3 - 1j, 2 + 1j); // complex number
  auto m = pg_line(-2 + 1j, 1 - 3j, -1 - 1j);
  auto A = l * m;
  REQUIRE(A == m * l);
  REQUIRE(A.incident(l));
  REQUIRE(A.incident(m));

  auto lm = plucker(2 + 3j, l, 3 - 4j, m);
  REQUIRE(A.incident(lm));

  // assert coI([l, m, plucker(1, l, 1, m), plucker(1, l, -1, m)])

  auto r = pg_line(2 - 1j, -2 + 1j, 1 + 1j);
  auto s = pg_line(2j, 2 - 2j, 3 + 0j);
  auto t = pg_line(2 + 0j, -2j, 2 + 0j);

  // assert not persp([l, m, l + m], [r, l + r, l])

  auto o = join(meet(l, s), meet(m, t));
  auto u = plucker(1 + 0j, o, -1 - 1j, r);
  check_desargue(l, m, r, s, t, u);
}

// int main(int argc, char* argv[]) {
//   //using namespace ModernCppCI;

//   auto result = Catch::Session().run(argc, argv);

//   return (result < 0xff ? result : 0xff);
// }