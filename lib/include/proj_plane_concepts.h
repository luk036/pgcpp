/** @file include/proj_plane_concepts.h
 *  This is a C++ Library header.
 */


#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PROJ_PLANE_CONCEPTS_H
#define _HOME_UBUNTU_GITHUB_PGCPP_PROJ_PLANE_CONCEPTS_H 1

#include "common_concepts.h"

/**
 @todo: projectivity >=
**/

namespace fun {

/**
 * @brief Projective plane Concept (half)
 * 
 * @tparam P Point
 * @tparam L Line
 */
template <class P, class L>
concept bool Projective_plane_h = 
  Equality_comparable<P> && requires(P p, P q, L l) {
  { P(p) } -> P; // copyable
  { p.incident(l) } -> bool; // incidence
  { L(p, q) } -> L; // join or meet
  { p.aux() } -> L; // line not incident with p
  { p.dot(l) } -> typename P::value_type; 
};

/**
 * @brief Projective plane Concept (full)
 * 
 * @tparam P Point
 * @tparam L Line 
 */
template <class P, class L = typename P::dual>
concept bool Projective_plane =
  Projective_plane_h<P, L> && Projective_plane_h<L, P>;

/*
axiom(P p, P q, P r, L l) {
  l == L{p, q} => I(p, l) && I(q, l);
}
*/


/**
 * @brief Shorthand Notation of Projective_plane
 * 
 * @tparam P Point
 */
template <class P>
concept bool Projective_plane2 = Projective_plane<P>; // Make the compiler happy

/**
 * @brief Cayley Klein plane Concept (half)
 * 
 * @tparam P Point
 * @tparam L Line
 */
template <class P, class L>
concept bool Cayley_Klein_plane_h =
  Projective_plane_h<P, L> && requires(P p, P q, L l) {
  { ~p } -> L; // line not incident with p
};

/**
 * @brief Cayley Klein plane Concept (full)
 * 
 * @tparam P Point
 * @tparam L Line 
 */
template <class P, class L = typename P::dual>
concept bool Cayley_Klein_plane =
  Cayley_Klein_plane_h<P, L> && Cayley_Klein_plane_h<L, P>;

/* (non-degenerate) 
axiom(P p) {
  ~(~p) == p
}
*/

/**
 * @brief Shorthand Notation of Cayley_Klein_plane
 * 
 * @tparam P Point
 */
template <class P>
concept bool Cayley_Klein_plane2 = Cayley_Klein_plane<P>; // Make the compiler happy


} // namespace fun

#endif
