#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PROJ_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PROJ_PLANE_HPP 1

#include "fraction.hpp"
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
Projective_plane_prim2 { P }
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
    for (P &&p : seq) {
        if (!incident(l, p))
            return false;
    }
    return true;
}

template <typename P> using Triple = std::tuple<P, P, P>;

Projective_plane_prim2 { P }
constexpr auto tri_dual(const Triple<P> &tri) {
    auto &&[a1, a2, a3] = tri;
    using Line = typename P::dual;
    Line &&l1 = a2 * a3;
    Line &&l2 = a1 * a3;
    Line &&l3 = a1 * a2;
    return Triple<Line>{l1, l2, l3};
}

Projective_plane_prim2 { P }
constexpr auto tri_func(auto func, const Triple<P> &tri) {
    auto &&[a1, a2, a3] = tri;
    using ret_t = decltype(func(a1, a2));
    ret_t &&m1 = func(a2, a3);
    ret_t &&m2 = func(a1, a3);
    ret_t &&m3 = func(a1, a2);
    return Triple<ret_t>{m1, m2, m3};
}

Projective_plane_prim2 { Point }
constexpr bool persp(const Triple<Point> &tri1, const Triple<Point> &tri2) {
    auto &&[A, B, C] = tri1;
    auto &&[D, E, F] = tri2;
    Point &&O = (A * D) * (B * E);
    return incident(O, C * F);
}

Projective_plane { P, L }
constexpr bool incident(const P &p, const L &l) {
    using K = typename P::value_type;
    return p.dot(l) == K(0);
}

Projective_plane2 { P }
constexpr P harm_conj(const P &A, const P &B, const P &C) {
    using Line = typename P::dual;
    Line &&lC = C * (A * B).aux();
    return plucker(B.dot(lC), A, A.dot(lC), B);
}

Projective_plane { P, L }
class involution {
    using K = typename P::value_type;

  private:
    L _m;
    P _o;
    K _c;

  public:
    constexpr involution(const L &m, const P &o)
        : // input mirror and center
          _m{m}, _o{o}, _c{m.dot(o)} {}

    constexpr P operator()(const P &p) const {
        return plucker(_c, p, -K(2) * p.dot(_m), _o);
    }
};

/**
 * @brief Basic ratio of ratio
 *
 * @tparam K
 * @param a
 * @param b
 * @param c
 * @param d
 * @return (a/b) / (c/d)
 */
// Integer { K }
// constexpr Fraction<K> ratio_ratio(const K &a, const K &b, const K &c,
//                                   const K &d) {
//     return Fraction<K>(a, b) / Fraction<K>(c, d);
// }

template <typename K>
constexpr auto ratio_ratio(const K &a, const K &b, const K &c, const K &d) {
    if constexpr (std::is_integral<K>::value) {
        return Fraction<K>(a, b) / Fraction<K>(c, d);
    }
    else {
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
    using ret_t = decltype(A.dot(l));
    ret_t &&dAl = A.dot(l);
    ret_t &&dAm = A.dot(m);
    ret_t &&dBl = B.dot(l);
    ret_t &&dBm = B.dot(m);
    return ratio_ratio(dAl, dAm, dBl, dBm);
}

/* @todo

def R(A, B, C, D):
    # not sure???
    if A[1]*B[2] != B[1]*A[2]:
        # Project points to yz-plane
        a, b, c, d = A[1:], B[1:], C[1:], D[1:]
    else:
        # Project points to xz-plane
        a, b, c, d = A[(0, 2)], B[(0, 2)], C[(0, 2)], D[(0, 2)]
    return R1(a, b, c, d)

def isharmonic(A, B, C, D):
    return R(A, B, C, D) == -1

*/

/**
 * @brief Check Pappus Theorem
 *
 * @tparam P
 * @tparam L
 * @param co1
 * @param co2
 */
Projective_plane_prim2 { Point }
void check_pappus(const Triple<Point> &co1, const Triple<Point> &co2) {
    auto &&[A, B, C] = co1;
    auto &&[D, E, F] = co2;
    Point &&G = (A * E) * (B * D);
    Point &&H = (A * F) * (C * D);
    Point &&I = (B * F) * (C * E);
    assert(coincident(G, H, I));
}

// template <class P, class L = typename P::dual>
// requires Projective_plane<P, L>
Projective_plane_prim2 { P }
void check_desargue(const Triple<P> &tri1, const Triple<P> &tri2) {
    auto &&trid1 = tri_dual(tri1);
    auto &&trid2 = tri_dual(tri2);
    bool b1 = persp(tri1, tri2);
    bool b2 = persp(trid1, trid2);
    assert((b1 && b2) || (!b1 && !b2));
}

} // namespace fun

#endif
