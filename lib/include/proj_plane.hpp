#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PROJ_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PROJ_PLANE_HPP 1

#include "proj_plane_concepts.h"
#include <boost/rational.hpp>
#include <cassert>

/** @file include/proj_plane.hpp
 *  This is a C++ Library header.
 */

/**
 @todo: projectivity >=
**/

using boost::rational;

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
  return r.incident(L{p, q});
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

//template <class P, class L = typename P::dual>
//requires Projective_plane<P, L>

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
auto x_ratio(P & A, P & B, L & l, L & m) {
  auto dAl = A.dot(l);
  auto dAm = A.dot(m);
  auto dBl = B.dot(l);
  auto dBm = B.dot(m);
  using K = typename P::value_type;
  if constexpr (std::is_integral<K>::value) {
    return rational<K>(dAl, dAm) / rational<K>(dBl, dBm);
  } else {
    return dAl * dBm / (dAm * dBl);
  }
}

/**
 * @brief Check Pappus Theorem
 * 
 * @tparam P 
 * @tparam L 
 * @param A 
 * @param B 
 * @param C 
 * @param D 
 * @param E 
 * @param F 
 */
template <class P, class L = typename P::dual>
  requires Projective_plane<P, L>
void check_pappus(P & A, P & B, P & C, P & D, P & E, P & F) {
  auto G = P{L{A, E}, L{B, D}};
  auto H = P{L{A, F}, L{C, D}};
  auto I = P{L{B, F}, L{C, E}};
  assert(G.incident(L{H, I}) && "Pappus Theorem failed");
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
  auto O = P{L{A, D}, L{B, E}};
  return O.incident(L{C, F});
}

//template <class P, class L = typename P::dual>
//requires Projective_plane<P, L>
Projective_plane{P, L}
void check_desargue(P & A, P & B, P & C, P & D, P & E, P & F) {
  auto a = L{B, C};
  auto b = L{A, C};
  auto c = L{B, A};
  auto d = L{E, F};
  auto e = L{D, F};
  auto f = L{E, D};
  auto b1 = persp(A, B, C, D, E, F);
  auto b2 = persp(a, b, c, d, e, f);
  assert((b1 && b2) || (!b1 && !b2));
}

} // namespace fun

#endif
