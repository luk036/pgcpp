#ifndef _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP 1

#include "proj_plane.hpp"
#include "proj_plane_concepts.h"
#include <cassert>

namespace fun {

namespace CK {

Cayley_Klein_plane2{L}
bool is_perpendicular(const L & l, const L & m) {
  return m.incident(~l);
}

Cayley_Klein_plane{P, L}
auto altitude(const P & p, const L & l) {
  return p * ~l;
}

//template <class P, class L = typename P::dual>
//requires Cayley_Klein_plane<P, L>
Cayley_Klein_plane2{P}
auto orthocenter(const P & a1, const P & a2, const P & a3) {
  auto t1 = altitude(a1, a2 * a3);
  auto t2 = altitude(a2, a1 * a3);
  return t1 * t2;
}

Cayley_Klein_plane2{L}
auto line_reflect(const L& m) {
  return involution(m, ~m);
}

// Cayley_Klein_plane2{P}
// auto omega(const P & x) {
//   return x.dot(~x);
// }

template <class P, class L = typename P::dual>
requires Cayley_Klein_plane<P, L>
auto measure(const P &a1, const P &a2) {
  return 1 - x_ratio(a1, a2, ~a2, ~a1);  
}

Cayley_Klein_plane2{P}
auto quadrance(const P &p, const P &q) {
  return measure(p, q);
}

Cayley_Klein_plane2{L}
auto spread(const L &l, const L &m) {
  return measure(l, m);
}

template <typename _K>
bool check_sine_law(const _K &s1, const _K &q1, const _K &s2, const _K &q2) {
  return s1*q2 == s2*q1;
}


}  // namespace CK

} // namespace fun

#endif
