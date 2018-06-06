#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PERSP_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PERSP_PLANE_HPP 1

#include "ck_plane.hpp"
#include "pg_common.hpp"
#include <cassert>

namespace fun {

Projective_plane { P, L }
class persp_euclid_plane : public ck<P, L, persp_euclid_plane> {
  using K = typename P::value_type;

private:
  P _Ire;
  P _Iim;
  L _l_infty;

public:
  persp_euclid_plane(const P &Ire, const P &Iim, const L &l_infty)
      : _Ire{Ire}, _Iim{Iim}, _l_infty{l_infty} {}

  L _perp(const P &x) const { return _l_infty; }

  P _perp(const L &x) const {
    return plucker(x.dot(_Ire), _Ire, x.dot(_Iim), _Iim);
  }

  bool is_parallel(const L &l, const L &m) const {
    return incident(_l_infty, l * m);
  }

  P midpoint(const P &a, const P &b) const {
    return plucker(b.dot(_l_infty), a, a.dot(_l_infty), b);
  }

  auto omega(const P &x) const { return sq(x.dot(_l_infty)); }

  auto omega(const L &x) const { return sq(x.dot(_Ire)) + sq(x.dot(_Iim)); }

  Projective_plane2 { _P }
  auto measure(const _P &a1, const _P &a2) const {
    auto omg = this->omega(a1 * a2);
    auto den = this->omega(a1) * this->omega(a2);
    if constexpr (std::is_integral<K>::value) {
      return Fraction<K>(omg, den);
    } else {
      return omg / den;
    }
  }

  auto quadrance(const P &a1, const P &a2) const {
    return this->measure(a1, a2);
  }

  auto spread(const L &l1, const L &l2) const { return this->measure(l1, l2); }

  auto tri_quadrance(const P &a1, const P &a2, const P &a3) {
    return tri_func(this->quadrance, std::tuple{a1, a2, a3});
  }

  auto tri_spread(const L &l1, const L &l2, const L &l3) {
    return tri_func(this->spread, std::tuple{l1, l2, l3});
  }

  auto cross(const L &l1, const L &l2) const {
    return 1 - this->spread(l1, l2); // ???
  }
};

template <typename _Q> auto Ar(const _Q &a, const _Q &b, const _Q &c) {
  return (4 * a * b) - sq(a + b - c);
}

} // namespace fun

#endif