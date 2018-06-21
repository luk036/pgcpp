#ifndef _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP 1

#include "proj_plane.hpp"
#include "proj_plane_concepts.h"
#include <cassert>
#include <type_traits>

namespace fun {

template <typename _P, typename _L,
          template <typename P, typename L> class Derived>
    requires Projective_plane_prim<_P, _L>
struct ck {
    using cDer = const Derived<_P, _L>;
    cDer &self = *static_cast<cDer *>(this);

    explicit ck() {
        static_assert(
            std::is_base_of<ck<_P, _L, Derived>, Derived<_P, _L>>::value);
    }

    Projective_plane_prim2 { L }
    bool is_perpendicular(const L &l, const L &m) const {
        return incident(m, self._perp(l));
    }

    Projective_plane_prim{P, L}
    L altitude(const P &p, const L &l) const {
        return p * self._perp(l);
    }

    Projective_plane_prim2 { P }
    constexpr auto tri_altitude(const P &a1, const P &a2, const P &a3) {
        auto [l1, l2, l3] = tri(std::tuple{a1, a2, a3});
        auto t1 = altitude(a1, l1);
        auto t2 = altitude(a2, l2);
        auto t3 = altitude(a3, l3);
        return std::tuple{t1, t2, t3};
    }

    Projective_plane_prim2 { P }
    auto orthocenter(const P &a1, const P &a2, const P &a3) const {
        auto t1 = altitude(a1, a2 * a3);
        auto t2 = altitude(a2, a1 * a3);
        return t1 * t2;
    }

    Projective_plane2 { L }
    auto reflect(const L &m) const {
        static_assert(
            std::is_base_of<ck<_P, _L, Derived>, Derived<_P, _L>>::value);
        return involution(m, self._perp(m));
    }

    // Projective_plane2{P}
    // auto omega(const P & x) {
    //   return x.dot(~x);
    // }

    Projective_plane2 { P }
    auto measure(const P &a1, const P &a2) const {
        using K = Value_type<P>;
        return K(1) - x_ratio(a1, a2, self._perp(a2), self._perp(a1));
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

template <typename _Q>
bool check_sine_law(const _Q &s1, const _Q &q1, const _Q &s2, const _Q &q2) {
    return s1 * q2 == s2 * q1;
}

Projective_plane_prim { P, L } // and requires vector computations
struct ellck : ck<P, L, ellck> {
    L _perp(const P &v) const { return L(v); }
    P _perp(const L &v) const { return P(v); }
};

Projective_plane_prim { P, L } // and requires vector computations
struct hyck : ck<P, L, hyck> {
    L _perp(const P &v) const { return L(v[0], v[1], -v[2]); }
    P _perp(const L &v) const { return P(v[0], v[1], -v[2]); }
};

template <typename _Q>
auto check_cross_TQF(const _Q &q1, const _Q &q2, const _Q &q3) {
    return sq(q1 + q2 + q3) - 2 * (q1 * q1 + q2 * q2 + q3 * q3) -
           4 * q1 * q2 * q3;
}

template <typename _Q>
bool check_cross_law(const _Q &s1, const _Q &s2, const _Q &s3, const _Q &q3) {
    return sq(s1 * s2 * q3 - (s1 + s2 + s3) + 2) ==
           4 * (1 - s1) * (1 - s2) * (1 - s3);
}

} // namespace fun

#endif
