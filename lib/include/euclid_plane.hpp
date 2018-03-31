#ifndef _HOME_UBUNTU_CUBSTORE_PROJ_GEOM_PGCPP_EUCLID_PLANE_HPP
#define _HOME_UBUNTU_CUBSTORE_PROJ_GEOM_PGCPP_EUCLID_PLANE_HPP 1

//#include "proj_plane_concepts.h"
#include "ck_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <cassert>

namespace fun {

namespace EUG {

template <typename _K> pg_line<_K> operator~(pg_point<_K> p) {
  return pg_line<_K>(0, 0, p[2]);
}

template <typename _K> pg_point<_K> operator~(pg_line<_K> l) {
  return pg_point<_K>(l[0], l[1], 0);
}

template <class L, class P = typename L::dual>
// requires Cayley_Klein_plane<P, L>
bool is_perpendicular(const L &l, const L &m) {
  return m.incident(~l);
}

template <class P, class L>
// requires Cayley_Klein_plane<P, L>
auto altitude(const P &p, const L &l) {
  return L{p * ~l};
}

template <class P, class L = typename P::dual>
// requires Cayley_Klein_plane<P, L>
auto omega(const P &x) {
  return x.dot(~x);
}

template <class P, class L = typename P::dual>
// requires Cayley_Klein_plane<P, L>
auto measure(const P &a1, const P &a2) {
  auto omg = a1.dot(~a2);
  using K = typename P::value_type;
  if constexpr (std::is_integral<K>::value) {
    return 1 - rational<K>(omg, omega(a1)) * rational<K>(omg, omega(a2));
  } else {
    return 1 - (omg * omg) / (omega(a1) * omega(a2));
  }
}

template <class P, class L = typename P::dual>
// requires Cayley_Klein_plane<P, L>
auto quadrance(const P &p, const P &q) {
  return measure(p, q);
}

template <class L, class P = typename L::dual>
// requires Cayley_Klein_plane<P, L>
auto spread(const L &l, const L &m) {
  return measure(l, m);
}

template <typename _K>
auto check_cross_TQF(const _K &q1, const _K &q2, const _K &q3) {
  return sq(q1 + q2 + q3) - 2*(q1*q1 + q2*q2 + q3*q3) - 4*q1*q2*q3;
}

template <typename _K>
bool check_cross_law(const _K &s1, const _K &s2, const _K &s3, const _K &q3) {
  return sq(s1 * s2 * q3 - (s1 + s2 + s3) + 2) ==
         4 * (1 - s1) * (1 - s2) * (1 - s3);
}

} // namespace HG

} // namespace fun

#endif
