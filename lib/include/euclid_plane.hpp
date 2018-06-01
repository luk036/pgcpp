#ifndef _HOME_UBUNTU_CUBSTORE_PROJ_GEOM_PGCPP_EUCLID_PLANE_HPP
#define _HOME_UBUNTU_CUBSTORE_PROJ_GEOM_PGCPP_EUCLID_PLANE_HPP 1

//#include "proj_plane_concepts.h"
#include "ck_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <cassert>
#include <cmath>

namespace fun {

Projective_plane_prim2 { L } // and requires p[i]
constexpr auto dot1(const L &x, const L &y) {
  return x[0] * y[0] + x[1] * y[1];
}

// Projective_plane2 { P }
// constexpr auto cross1(const P &x, const P &y) {
//   return x[0] * y[1] - x[1] * y[0];
// }

Projective_plane2 { L } // // and requires p[i]
constexpr auto fB(const L &l) {
  using P = typename L::dual;
  using K = typename L::value_type;
  return P(l[0], l[1], K(0));
}

Projective_plane2 { L }
constexpr auto is_perpendicular(const L &l, const L &m) {
  using K = typename L::value_type;
  return dot1(l, m) == K(0);
}

Projective_plane2 { L }
constexpr auto is_parallel(const L &l, const L &m) {
  return l[0] * m[1] == l[1] * m[0];
}

Projective_plane { P, L }
constexpr L altitude(const P &a, const L &l) { return a * fB(l); }

Projective_plane2 { P }
constexpr auto tri_altitude(const P &a1, const P &a2, const P &a3) {
  auto [l1, l2, l3] = tri(std::tuple{a1, a2, a3});
  auto t1 = altitude(a1, l1);
  auto t2 = altitude(a2, l2);
  auto t3 = altitude(a3, l3);
  return std::tuple{t1, t2, t3};
}

Projective_plane2 { P }
constexpr P orthocenter(const P &a1, const P &a2, const P &a3) {
  auto t1 = altitude(a1, a2 * a3);
  auto t2 = altitude(a2, a1 * a3);
  return t1 * t2;
}

Projective_plane2 { L }
constexpr auto reflect(const L &m) { return involution(m, fB(m)); }

Projective_plane2 { P }
constexpr auto omgB(const P &x, const P &y) {
  return x[0] * y[0] + x[1] * y[1];
}

Projective_plane2 { P }
constexpr auto det(const P &x, const P &y) { return x[0] * y[1] - x[1] * y[0]; }

Projective_plane2 { P }
constexpr P midpoint(const P &a, const P &b) {
  return plucker(b[2], a, a[2], b);
}

Integer { K }
constexpr auto quad1(const K &x1, const K &z1, const K &x2, const K &z2) {
  return sq(Fraction(x1, z1) - Fraction(x2, z2));
}

template <typename K>
constexpr auto quad1(const K &x1, const K &z1, const K &x2, const K &z2) {
  return sq(x1 / z1 - x2 / z2);
}

Projective_plane2 { P }
constexpr auto quadrance(const P &a1, const P &a2) {
  return quad1(a1[0], a1[2], a2[0], a2[2]) + quad1(a1[1], a1[2], a2[1], a2[2]);
}

Projective_plane2 { L }
constexpr auto sbase(const L &l1, const L &l2, const Integer &d) {
  return Fraction(d, omgB(l1, l1)) * Fraction(d, omgB(l2, l2));
}

Projective_plane2 { L }
constexpr auto sbase(const L &l1, const L &l2, const auto &d) {
  return (d * d) / (omgB(l1, l1) * omgB(l2, l2));
}

Projective_plane2 { L }
constexpr auto spread(const L &l1, const L &l2) {
  return sbase(l1, l2, det(l1, l2));
}

Projective_plane2 { P }
constexpr auto tri_quadrance(const P &a1, const P &a2, const P &a3) {
  return tri_func(quadrance, std::tuple{a1, a2, a3});
}

Projective_plane2 { L }
constexpr auto tri_spread(const L &l1, const L &l2, const L &l3) {
  return tri_func(spread, std::tuple{l1, l2, l3});
}

Projective_plane { P, L }
constexpr auto cross(const L &l1, const L &l2) {
  return sbase(l1, l2, omgB(l1, l2));
}

Projective_plane2 { P }
constexpr P uc_point(const Value_type<P> &lambda1, const Value_type<P> &mu1) {
  auto lambda2 = lambda1 * lambda1;
  auto mu2 = mu1 * mu1;
  return P(lambda2 - mu2, 2 * lambda1 * mu1, lambda2 + mu2);
}

/// Archimedes's function
template <typename _Q>
constexpr auto Ar(const _Q &a, const _Q &b, const _Q &c) {
  return (4 * a * b) - sq(a + b - c);
}

/// Cyclic quadrilateral quadrea theorem
template <typename _Q>
constexpr auto cqq(const _Q &a, const _Q &b, const _Q &c, const _Q &d) {
  auto t1 = 4 * a * b;
  auto t2 = 4 * c * d;
  auto m = (t1 + t2) - sq(a + b - c - d);
  auto p = m * m - 4 * t1 * t2;
  return std::tuple{m, p};
}

template <typename _Q>
constexpr auto Ptolemy(const _Q &Q12, const _Q &Q23, const _Q &Q34,
                       const _Q &Q14, const _Q &Q13, const _Q &Q24) {
  return Ar(Q12 * Q34, Q23 * Q14, Q13 * Q24) == 0;
}

Projective_plane2 { P }
constexpr auto distance(const P &a, const P &b) {
  return std::sqrt(double(quadrance(a, b)));
}

Projective_plane2 { L }
constexpr auto angle(const L &l, const L &m) {
  return std::asin(std::sqrt(double(spread(l, m))));
}

} // namespace fun

#endif
