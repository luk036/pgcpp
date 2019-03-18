#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PERSP_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PERSP_PLANE_HPP 1

#include "ck_plane.hpp"
#include "pg_common.hpp"
#include <cassert>

namespace fun {

/**
 * @brief
 *
 * @tparam P
 * @tparam P::dual
 */
template <typename P, typename L = typename P::dual>
requires Projective_plane_prim<P, L> // c++20 concept
class persp_euclid_plane : public ck<P, L, persp_euclid_plane> {
    using K = Value_type<P>;

  private:
    P _Ire;
    P _Iim;
    L _l_infty;

  public:
    /**
     * @brief Construct a new persp euclid plane object
     *
     * @param Ire
     * @param Iim
     * @param l_infty
     */
    constexpr persp_euclid_plane(P &&Ire, P &&Iim, L &&l_infty)
        : _Ire{std::move(Ire)}, //
          _Iim{std::move(Iim)}, //
          _l_infty{std::move(l_infty)} {}

    /**
     * @brief Construct a new persp euclid plane object
     *
     * @param Ire
     * @param Iim
     * @param l_infty
     */
    constexpr persp_euclid_plane(const P &Ire, const P &Iim, const L &l_infty)
        : _Ire{Ire}, _Iim{Iim}, _l_infty{l_infty} {}

    /**
     * @brief
     *
     * @param x
     * @return const L&
     */
    constexpr const L &perp(const P &x) const { return _l_infty; }

    /**
     * @brief
     *
     * @return const L&
     */
    constexpr const L &l_infty() const { return _l_infty; }

    /**
     * @brief
     *
     * @param x
     * @return P
     */
    constexpr P perp(const L &x) const {
        return plucker(x.dot(_Ire), _Ire, x.dot(_Iim), _Iim);
    }

    /**
     * @brief
     *
     * @param l
     * @param m
     * @return true
     * @return false
     */
    constexpr bool is_parallel(const L &l, const L &m) const {
        return incident(_l_infty, l * m);
    }

    /**
     * @brief
     *
     * @param a
     * @param b
     * @return P
     */
    constexpr P midpoint(const P &a, const P &b) const {
        return plucker(b.dot(_l_infty), a, a.dot(_l_infty), b);
    }

    /**
     * @brief
     *
     * @param tri
     * @return auto
     */
    constexpr auto tri_midpoint(const Triple<P> &tri) const {
        auto const &[a1, a2, a3] = tri;
        auto m12 = this->midpoint(a1, a2);
        auto m23 = this->midpoint(a2, a3);
        auto m13 = this->midpoint(a1, a3);
        return Triple<P>{std::move(m12), std::move(m23), std::move(m13)};
    }

    /**
     * @brief
     *
     * @param x
     * @return K
     */
    constexpr K omega(const P &x) const { return sq(x.dot(_l_infty)); }

    /**
     * @brief
     *
     * @param x
     * @return K
     */
    constexpr K omega(const L &x) const {
        return sq(x.dot(_Ire)) + sq(x.dot(_Iim));
    }

    /**
     * @brief
     *
     * @param a1
     * @param a2
     * @return auto
     */
    Projective_plane2 { _P }
    constexpr auto measure(const _P &a1, const _P &a2) const {
        auto omg = K(this->omega(a1 * a2));
        auto den = K(this->omega(a1) * this->omega(a2));
        if constexpr (Integral<K>) {
            return Fraction<K>(omg, den);
        } else {
            return omg / den;
        }
    }

    /**
     * @brief
     *
     * @param l1
     * @param l2
     * @return auto
     */
    constexpr auto cross_s(const L &l1, const L &l2) const {
        return 1 - this->spread(l1, l2); // ???
    }
};

// Projective_plane{P, L} persp_euclid_plane(const P &, const P &, const L &)
//     ->persp_euclid_plane<P, L>;

} // namespace fun

#endif