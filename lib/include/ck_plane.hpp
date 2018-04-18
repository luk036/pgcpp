#ifndef _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP 1

#include "proj_plane.hpp"
#include "proj_plane_concepts.h"
#include <cassert>

namespace fun {

template <typename dL, class L, class P = typename L::dual>
  requires Projective_plane<P, L>
bool is_perpendicular(const L & l, const L & m, dL dualL) {
  return m.incident(dualL(l));
}

template <typename dL, class P, class L = typename P::dual>
  requires Projective_plane<P, L>
auto altitude(const P & p, const L & l, dL dualL) {
  return p * dualL(l);
}

template <typename dL, class P, class L = typename P::dual>
  requires Projective_plane<P, L>
auto orthocenter(const P & a1, const P & a2, const P & a3, dL dualL) {
  auto t1 = altitude<dL, P, L>(a1, a2 * a3, dualL);
  auto t2 = altitude<dL, P, L>(a2, a1 * a3, dualL);
  return t1 * t2;
}

template <typename dL, class L, class P = typename L::dual>
  requires Projective_plane<P, L>
auto line_reflect(const L& m, dL dualL) {
  return involution(m, dualL(m));
}

// Projective_plane2{P}
// auto omega(const P & x) {
//   return x.dot(~x);
// }

template <typename dP, class P, class L = typename P::dual>
  requires Projective_plane<P, L>
auto measure(const P &a1, const P &a2, dP dualP) {
  return 1 - x_ratio(a1, a2, dualP(a2), dualP(a1));
}

template <typename dP, class P, class L = typename P::dual>
  requires Projective_plane<P, L>
auto quadrance(const P &p, const P &q, dP dualP) {
  return measure<dP, P, L>(p, q, dualP);
}

template < typename dL, class L, class P = typename L::dual>
  requires Projective_plane<P, L>
auto spread(const L &l, const L &m, dL dualL) {
  return measure<dL, L, P>(l, m, dualL);
}

template <typename _K>
bool check_sine_law(const _K &s1, const _K &q1, const _K &s2, const _K &q2) {
  return s1*q2 == s2*q1;
}


} // namespace fun

#endif
