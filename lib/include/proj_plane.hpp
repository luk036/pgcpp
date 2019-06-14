#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PROJ_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PROJ_PLANE_HPP 1

#include "fractions.hpp"
#include "proj_plane_concepts.h"
#include <cassert>
#include <tuple>

/** @file include/proj_plane.hpp
 *  This is a C++ Library header.
 */

/**
 @todo: projectivity >=
**/

namespace fun {

/**
 * @brief Coincident
 *
 * @tparam P Point
 * @tparam L Line
 * @param p p
 * @param q q
 * @param r r
 * @return true if three points are conincident
 * @return false otherwise
 */
template <Projective_plane_prim2 P>
constexpr bool coincident(const P &p, const P &q, const P &r) {
    return incident(r, p * q);
}

/**
 * @brief Coincident
 *
 * @tparam P Point
 * @tparam L Line
 * @param l line
 * @param seq Sequence of points
 * @return true if all points are incident with l
 * @return false otherwise
 */
Projective_plane_prim { P, L }
constexpr bool coincident(const L &l, const Sequence &seq) {
    for (auto const &p : seq) {
        if (!incident(l, p))
            return false;
    }
    return true;
}

template <typename P> using Triple = std::tuple<P, P, P>;

/**
 * @brief
 *
 * @param tri
 * @return auto
 */
template <Projective_plane_prim2 P>
constexpr auto tri_dual(const Triple<P> &tri) {
    auto &&[a1, a2, a3] = tri;
    return std::tuple{a2 * a3, a1 * a3, a1 * a2};
}

/**
 * @brief
 *
 * @param func
 * @param tri
 * @return auto
 */
template <Projective_plane_prim2 P>
constexpr auto tri_func(auto func, const Triple<P> &tri) {
    auto &&[a1, a2, a3] = tri;
    auto m1 = func(a2, a3);
    auto m2 = func(a1, a3);
    auto m3 = func(a1, a2);
    return std::tuple{std::move(m1), std::move(m2), std::move(m3)};
}

/**
 * @brief
 *
 * @param tri1
 * @param tri2
 * @return true
 * @return false
 */
template <Projective_plane_prim2 P>
constexpr bool persp(const Triple<P> &tri1, const Triple<P> &tri2) {
    auto &&[A, B, C] = tri1;
    auto &&[D, E, F] = tri2;
    auto O = (A * D) * (B * E);
    return incident(O, C * F);
}

/**
 * @brief
 *
 * @param p
 * @param l
 * @return true
 * @return false
 */
Projective_plane { P, L }
constexpr bool incident(const P &p, const L &l) { return p.dot(l) == 0; }

/**
 * @brief 
 * 
 * @tparam P 
 * @param A 
 * @param B 
 * @param C 
 * @return constexpr P 
 */
template <Projective_plane2 P>
constexpr auto harm_conj(const P &A, const P &B, const P &C) -> P {
    auto lC = C * (A * B).aux();
    return plucker(B.dot(lC), A, A.dot(lC), B);
}

Projective_plane { P, L }
class involution {
    using K = Value_type<P>;

  private:
    const L &_m;
    P _o;
    K _c;

  public:
    /**
     * @brief Construct a new involution object
     *
     * @param m
     * @param o
     */
    constexpr involution(const L &m, P &&o)
        : // input mirror and center
          _m{m}, _o{std::move(o)}, _c{m.dot(_o)} {}

    /**
     * @brief Construct a new involution object
     *
     * @param m
     * @param o
     */
    constexpr involution(const L &m, const P &o)
        : // input mirror and center
          _m{m}, _o{o}, _c{m.dot(o)} {}

    /**
     * @brief
     *
     * @param p
     * @return P
     */
    constexpr auto operator()(const P &p) const -> P {
        return plucker(_c, p, K(-2 * p.dot(_m)), _o);
    }
};

/**
 * @brief
 *
 * @tparam K
 * @param a
 * @param b
 * @param c
 * @param d
 * @return auto
 */
template <CommutativeRing K>
constexpr auto ratio_ratio(const K &a, const K &b, const K &c, const K &d) {
    if constexpr (Integral<K>) {
        return Fraction(a, b) / Fraction(c, d);
    } else {
        return (a * d) / (b * c);
    }
}

/**
 * @brief Cross Ratio
 *
 * @tparam P
 * @tparam L
 * @param A point \in P
 * @param B point \in P
 * @param l line \in P
 * @param m line \in P
 * @return cross ratio R(A,B;l,m)
 *
 * @todo rewrite by projecting to the y-axis first [:2]
 */
Projective_plane { P, L }
constexpr auto x_ratio(const P &A, const P &B, const L &l, const L &m) {
    auto dAl = A.dot(l);
    auto dAm = A.dot(m);
    auto dBl = B.dot(l);
    auto dBm = B.dot(m);
    return ratio_ratio(dAl, dAm, dBl, dBm);
}

/**
 * @brief
 *
 * @param A
 * @param B
 * @param C
 * @param D
 * @return constexpr auto
 */
template <Projective_plane2 P>
constexpr auto R(const P &A, const P &B, const P &C, const P &D) {
    auto O = (C * D).aux();
    return x_ratio(A, B, O * C, O * D);
}

/**
 * @brief
 *
 * @param A
 * @param B
 * @param C
 * @param D
 * @return constexpr auto
 */
template <Projective_plane_coord2 P>
constexpr auto R0(const P &A, const P &B, const P &C, const P &D) {
    auto ac = cross0(A, C);
    auto ad = cross0(A, D);
    auto bc = cross0(B, C);
    auto bd = cross0(B, D);
    return ratio_ratio(ac, ad, bc, bd);
}

/**
 * @brief
 *
 * @param A
 * @param B
 * @param C
 * @param D
 * @return constexpr auto
 */
template <Projective_plane_coord2 P>
constexpr auto R1(const P &A, const P &B, const P &C, const P &D) {
    auto ac = cross1(A, C);
    auto ad = cross1(A, D);
    auto bc = cross1(B, C);
    auto bd = cross1(B, D);
    return ratio_ratio(ac, ad, bc, bd);
}

/**
 * @brief
 *
 * @param A
 * @param B
 * @param C
 * @param D
 * @return constexpr auto
 */
template <Projective_plane_coord2 P>
constexpr auto R(const P &A, const P &B, const P &C, const P &D) {
    if (cross0(A, B) != 0) { // Project points to yz-plane
        return R0(A, B, C, D);
    }
    // Project points to xz-plane
    return R1(A, B, C, D);
}

/**
 * @brief
 *
 * @param A
 * @param B
 * @param C
 * @param D
 * @return constexpr auto
 */
template <Projective_plane2 P>
constexpr auto is_harmonic(const P &A, const P &B, const P &C, const P &D) {
    return R(A, B, C, D) == -1;
}

/**
 * @brief Check Pappus Theorem
 *
 * @tparam P
 * @tparam L
 * @param co1
 * @param co2
 */
template <Projective_plane_prim2 P>
void check_pappus(const Triple<P> &co1, const Triple<P> &co2) {
    auto &&[A, B, C] = co1;
    auto &&[D, E, F] = co2;
    auto G = (A * E) * (B * D);
    auto H = (A * F) * (C * D);
    auto I = (B * F) * (C * E);
    assert(coincident(G, H, I));
}

/**
 * @brief
 *
 * @param tri1
 * @param tri2
 */
template <Projective_plane_prim2 P>
void check_desargue(const Triple<P> &tri1, const Triple<P> &tri2) {
    auto trid1 = tri_dual(tri1);
    auto trid2 = tri_dual(tri2);
    auto b1 = persp(tri1, tri2);
    auto b2 = persp(trid1, trid2);
    assert((b1 && b2) || (!b1 && !b2));
}

} // namespace fun

#endif
