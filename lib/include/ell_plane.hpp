#ifndef _HOME_UBUNTU_GITHUB_PGCPP_ELL_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_ELL_PLANE_HPP 1

//#include "proj_plane_concepts.h"
#include "ck_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <cassert>

namespace fun {

namespace ELL {

template <class P, class L = typename P::dual>
  requires Projective_plane<P, L>
auto dual = [](const P& v) { return L(v); };

template <class L, class P = typename L::dual>
  requires Projective_plane<P, L>
bool is_perpendicular(const L & l, const L & m) {
  return fun::is_perpendicular(l, m, dual<L>);
}

template <class P, class L = typename P::dual>
  requires Projective_plane<P, L>
auto altitude(const P & p, const L & l) {
  return fun::altitude(p, l, dual<L>);
}

template <class P, class L = typename P::dual>
  requires Projective_plane<P, L>
auto orthocenter(const P & a1, const P & a2, const P & a3) {
  return fun::orthocenter(a1, a2, a3, dual<L>);
}

template <class L, class P = typename L::dual>
  requires Projective_plane<P, L>
auto line_reflect(const L& m) {
  return fun::line_reflect(m, dual<L>);
}

template <class P, class L = typename P::dual>
  requires Projective_plane<P, L>
auto measure(const P &a1, const P &a2) {
  return fun::measure(a1, a2, dual<P>);
}

template <class P, class L = typename P::dual>
  requires Projective_plane<P, L>
auto quadrance(const P &p, const P &q) {
  return fun::quadrance(p, q, dual<P>);
}

template <class L, class P = typename L::dual>
  requires Projective_plane<P, L>
auto spread(const L &l, const L &m) {
  return fun::spread(l, m, dual<L>);
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

} // namespace ELL

} // namespace fun

#endif
