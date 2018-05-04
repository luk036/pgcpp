#ifndef _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP 1

#include "proj_plane.hpp"
#include "proj_plane_concepts.h"
#include <cassert>

namespace fun {

Projective_plane { _P, _L }
class ck {
public:
  virtual _L _perp(const _P &) const = 0; // virtual

  virtual _P _perp(const _L &) const = 0; // virtual

  Projective_plane2 { L }
  bool is_perpendicular(const L &l, const L &m) const {
    return m.incident(_perp(l));
  }

  Projective_plane { P, L }
  auto altitude(const P &p, const L &l) const { return p * _perp(l); }

  Projective_plane2 { P }
  auto orthocenter(const P &a1, const P &a2, const P &a3) const {
    auto t1 = altitude(a1, a2 * a3);
    auto t2 = altitude(a2, a1 * a3);
    return t1 * t2;
  }

  Projective_plane2 { L }
  auto reflect(const L &m) const { return involution(m, _perp(m)); }

  // Projective_plane2{P}
  // auto omega(const P & x) {
  //   return x.dot(~x);
  // }

  Projective_plane2 { P }
  auto measure(const P &a1, const P &a2) const {
    using K = typename P::value_type;
    return K(1) - x_ratio(a1, a2, _perp(a2), _perp(a1));
  }

  Projective_plane2 { P }
  auto tri_measure(const Triple<P> &triangle) {
    return tri_func(this->measure, triangle);
  }

  auto quadrance(const _P &p, const _P &q) const { return measure(p, q); }

  auto spread(const _L &l, const _L &m) const { return measure(l, m); }

  auto tri_quadrance(const _P &a1, const _P &a2, const _P &a3) {
    return this->tri_measure(std::tuple{a1, a2, a3});
  }

  auto tri_spread(const _L &l1, const _L &l2, const _L &l3) {
    return this->tri_measure(std::tuple{l1, l2, l3});
  }
};

template <typename _K>
bool check_sine_law(const _K &s1, const _K &q1, const _K &s2, const _K &q2) {
  return s1 * q2 == s2 * q1;
}

Projective_plane { P, L }
class ellck : public ck<P, L> {
public:
  virtual L _perp(const P &v) const final { return L(v); }
  virtual P _perp(const L &v) const final { return P(v); }
};

Projective_plane { P, L }
class hyck : public ck<P, L> {
public:
  virtual L _perp(const P &v) const final { return L(v[0], v[1], -v[2]); }
  virtual P _perp(const L &v) const final { return P(v[0], v[1], -v[2]); }
};

template <typename _K>
auto check_cross_TQF(const _K &q1, const _K &q2, const _K &q3) {
  return sq(q1 + q2 + q3) - 2 * (q1 * q1 + q2 * q2 + q3 * q3) -
         4 * q1 * q2 * q3;
}

template <typename _K>
bool check_cross_law(const _K &s1, const _K &s2, const _K &s3, const _K &q3) {
  return sq(s1 * s2 * q3 - (s1 + s2 + s3) + 2) ==
         4 * (1 - s1) * (1 - s2) * (1 - s3);
}

} // namespace fun

#endif
