#ifndef _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_CK_PLANE_HPP 1

#include "proj_plane.hpp"
#include "proj_plane_concepts.h"
#include <cassert>
#include <type_traits>

namespace fun {

template <typename _P, typename _L,
          template <typename P, typename L> class Derived>
requires Projective_plane_prim<_P, _L> // c++20 concept
    struct ck {
    using point_t = _P;
    using line_t = _L;

    using cDer = const Derived<_P, _L>;
    cDer &self = *static_cast<cDer *>(this);

    constexpr explicit ck() {
        static_assert(
            std::is_base_of<ck<_P, _L, Derived>, Derived<_P, _L>>::value);
    }

    Projective_plane_prim2 { L }
    constexpr bool is_perpendicular(const L &l, const L &m) const {
        return incident(m, self.perp(l));
    }

    Projective_plane_prim { P, L }
    constexpr L altitude(const P &p, const L &l) const {
        return p * self.perp(l);
    }

    Projective_plane_prim2 { P }
    constexpr auto tri_altitude(const Triple<P> &tri) const {
        auto &&[l1, l2, l3] = tri_dual(tri);
        auto &&[a1, a2, a3] = tri;
        using L = typename P::dual;
        L &&t1 = altitude(a1, l1);
        L &&t2 = altitude(a2, l2);
        L &&t3 = altitude(a3, l3);
        return std::tuple{t1, t2, t3};
    }

    Projective_plane_prim2 { P }
    constexpr P orthocenter(const Triple<P> &tri) const {
        auto &&[a1, a2, a3] = tri;
        using L = typename P::dual;
        L &&t1 = altitude(a1, a2 * a3);
        L &&t2 = altitude(a2, a1 * a3);
        return t1 * t2;
    }

    Projective_plane2 { L }
    auto reflect(const L &m) const {
        static_assert(
            std::is_base_of<ck<_P, _L, Derived>, Derived<_P, _L>>::value);
        return involution(m, self.perp(m));
    }

    Projective_plane2 { P }
    constexpr auto tri_measure(const Triple<P> &tri) const {
        auto &&[a1, a2, a3] = tri;
        using ret_t = decltype(self.measure(a1, a2));
        ret_t &&m1 = self.measure(a2, a3);
        ret_t &&m2 = self.measure(a1, a3);
        ret_t &&m3 = self.measure(a1, a2);
        return Triple<ret_t>{m1, m2, m3};
    }

    constexpr auto quadrance(const _P &p, const _P &q) const {
        return self.measure(p, q);
    }

    constexpr auto spread(const _L &l, const _L &m) const {
        return self.measure(l, m);
    }

    constexpr auto tri_quadrance(const Triple<_P> &triangle) const {
        return this->tri_measure(triangle);
    }

    constexpr auto tri_spread(const Triple<_L> &trilateral) const {
        return this->tri_measure(trilateral);
    }
};

template <typename Q_t, typename S_t>
constexpr bool check_sine_law(const Q_t &Q, const S_t &S) {
    auto &&[q1, q2, q3] = Q;
    auto &&[s1, s2, s3] = S;
    if (s1 * q2 != s2 * q1) {
        return false;
    }
    if (s2 * q3 != s3 * q2) {
        return false;
    }
    return true;
}

Projective_plane_prim { P, L } // and requires vector computations
struct ellck : ck<P, L, ellck> {
    constexpr L perp(const P &v) const { return L(v); }

    constexpr P perp(const L &v) const { return P(v); }

    Projective_plane2 { _P }
    constexpr auto measure(const _P &a1, const _P &a2) const {
        using K = Value_type<_P>;
        return K(1) - x_ratio(a1, a2, this->perp(a2), this->perp(a1));
    }
};

Projective_plane_prim { P, L } // and requires vector computations
struct hyck : ck<P, L, hyck> {
    constexpr L perp(const P &v) const { return L(v[0], v[1], -v[2]); }

    constexpr P perp(const L &v) const { return P(v[0], v[1], -v[2]); }

    Projective_plane2 { _P }
    constexpr auto measure(const _P &a1, const _P &a2) const {
        using K = Value_type<_P>;
        return K(1) - x_ratio(a1, a2, this->perp(a2), this->perp(a1));
    }
};


template <typename Q_t> constexpr auto check_cross_TQF(const Q_t &Q) {
    auto &&[q1, q2, q3] = Q;
    auto &&sum = q1 + q2 + q3;
    return sum*sum - 2 * (q1 * q1 + q2 * q2 + q3 * q3) -
           4 * q1 * q2 * q3;
}

template <typename S_t, typename K>
constexpr auto check_cross_law(const S_t &S, const K &q3) {
    auto &&[s1, s2, s3] = S;
    auto &&temp = s1 * s2 * q3 - (s1 + s2 + s3) + 2;
    return temp * temp - 4 * (1 - s1) * (1 - s2) * (1 - s3);
}

} // namespace fun

#endif
