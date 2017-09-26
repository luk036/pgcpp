#ifndef _HOME_UBUNTU_CUBSTORE_PROJ_GEOM_PGCPP_CK_PLANE_HPP
#define _HOME_UBUNTU_CUBSTORE_PROJ_GEOM_PGCPP_CK_PLANE_HPP 1

#include "proj_plane_concepts.h"
#include "proj_plane.hpp"
#include <cassert>

namespace fun {

Cayley_Klein_plane2{L}
bool is_perpendicular(const L &l, const L &m) {
  return m.incident(~l);
}

Cayley_Klein_plane{P, L}
auto altitude(const P &p, const L &l) {
  return L{p, ~l};
}

//template <class P, class L = typename P::dual>
//requires Cayley_Klein_plane<P, L>
Cayley_Klein_plane{P, L}
auto orthocenter(const P &a1, const P &a2, const P &a3) {
  auto t1 = altitude(a1, L{a2, a3});
  auto t2 = altitude(a2, L{a1, a3});
  return P{t1, t2};
}

template <typename _K>
bool check_sine_law(const _K &s1, const _K &q1, const _K &s2, const _K &q2) {
  return s1*q2 == s2*q1;
}

namespace CK {

Cayley_Klein_plane2{P}
auto omega(const P &x) {
  return x.dot(~x);
}

/*
template <class P, class L = typename P::dual>
requires Cayley_Klein_plane<P, L>
auto measure(const P &a1, const P &a2) {
  auto omg = omega(L{a1, a2});
  using K = typename P::value_type;
  if constexpr (std::is_integral<K>::value) {
    return rational<K>(omg, omega(a1) * omega(a2));
  } 
    return omg / (omega(a1) * omega(a2));
  
}

template <class P, class L>
requires Cayley_Klein_plane<P, L>
auto quadrance(const P &p, const P &q) {
  return measure(p, q);
}

template <class P, class L>
requires Cayley_Klein_plane<P, L>
auto spread(const L &l, const L &m) {
  return measure(l, m);
}
*/

}

} // namespace fun

#endif
