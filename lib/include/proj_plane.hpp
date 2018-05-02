#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PROJ_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PROJ_PLANE_HPP 1

#include "proj_plane_concepts.h"
// #include <boost/rational.hpp>
#include "fraction.hpp"
#include <tuple>
#include <cassert>

/** @file include/proj_plane.hpp
 *  This is a C++ Library header.
 */

/**
 @todo: projectivity >=
**/

//using boost::rational;

namespace fun {

// Projective_plane{P, L}


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
template <class P, class L>
  requires Projective_plane<P, L>
constexpr bool coincident(const P & p, const P & q, const P & r) {
  return r.incident(p * q);
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
template <class P, class L>
  requires Projective_plane<P, L>
constexpr bool coincident(const L & l, const Sequence & seq) {
  for (const P &p : seq) {
    if (!l.incident(p))
      return false;
  }
  return true;
}

Projective_plane{P, L}
class involution {
  using K = typename P::value_type;

private:
  L _m;
  P _o;
  K _c;

public:
  involution(const L &m, const P &o): _m{m}, _o{o}, _c{m.dot(o)} {}

  auto operator()(const P &p) const {
    return plucker(_c, p, -K(2)*p.dot(_m), _o);
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
template <typename K>
auto ratio_ratio(const K &a, const K &b, const K &c, const K &d) {
  if constexpr (std::is_integral<K>::value) {
    return Fraction<K>(a, b) / Fraction<K>(c, d);
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
template <class P, class L>
  requires Projective_plane<P, L>
auto x_ratio(const P & A, const P & B, const L & l, const L & m) {
  auto dAl = A.dot(l);
  auto dAm = A.dot(m);
  auto dBl = B.dot(l);
  auto dBm = B.dot(m);
  return ratio_ratio(dAl, dAm, dBl, dBm);
}

template <typename P>
using Triple = std::tuple<P, P, P>;


/**
 * @brief Check Pappus Theorem
 *
 * @tparam P
 * @tparam L
 * @param co1
 * @param co2
 */
template <class P, class L = typename P::dual>
  requires Projective_plane<P, L>
void check_pappus(const Triple<P>& co1, const Triple<P>& co2) {
  auto [A, B, C] = co1;
  auto [D, E, F] = co2;
  // auto G = P{L{A*E} * L{B*D}};
  // auto H = P{L{A*F} * L{C*D}};
  // auto I = P{L{B*F} * L{C*E}};
  auto G = (A*E) * (B*D);
  auto H = (A*F) * (C*D);
  auto I = (B*F) * (C*E);
  assert (coincident(G, H, I));
}

//Projective_plane{P, L}
/**
 * @brief
 *
 * @tparam P
 * @tparam L
 * @param A
 * @param B
 * @param C
 * @param D
 * @param E
 * @param F
 * @return true
 * @return false
 */
template <class P, class L = typename P::dual>
  requires Projective_plane<P, L>
bool persp(P & A, P & B, P & C, P & D, P & E, P & F) {
  auto O = (A*D) * (B*E);
  return O.incident(C*F);
}

//template <class P, class L = typename P::dual>
//requires Projective_plane<P, L>
Projective_plane2{P}
void check_desargue(P & A, P & B, P & C, P & D, P & E, P & F) {
  auto a = B * C;
  auto b = A * C;
  auto c = B * A;
  auto d = E * F;
  auto e = D * F;
  auto f = E * D;

  auto b1 = persp(A, B, C, D, E, F);
  auto b2 = persp(a, b, c, d, e, f);
  assert((b1 && b2) || (!b1 && !b2));
}

} // namespace fun

#endif
