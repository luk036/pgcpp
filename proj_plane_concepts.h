#ifndef _PROJ_PLANE_CONCEPTS_H
#define _PROJ_PLANE_CONCEPTS_H 1

#include "common_concepts.h"

/**
 @todo: projectivity >=
**/

namespace fun {

template <class P, class L>
concept bool Projective_plane_h = 
  Equality_comparable<P> && requires(P p, P q, L l) {
  { P(p) } -> P; // copyable
  { p.incident(l) } -> bool; // incidence
  { L(p, q) } -> L; // join or meet
  { p.aux() } -> L; // line not incident with p
  { p.dot(l) } -> typename P::value_type; 
};

/*
axiom(P p, P q, P r, L l) {
  l == {p, q} => I(p, l) && I(q, l);
}
*/

template <class P, class L = typename P::dual>
concept bool Projective_plane =
  Projective_plane_h<P, L> && Projective_plane_h<L, P>;

template <class P>
concept bool Projective_plane2 = Projective_plane<P>; // Make the compiler happy

///  Return true if @a p, @a q, @a r are collinear
Projective_plane2{P}
constexpr bool coincident(const P &p, const P &q, const P &r) {
  using L = typename P::dual;
  return r.incident(L{p, q});
}

template <class P, class L>
concept bool Cayley_Klein_plane_h =
  Projective_plane_h<P, L> && requires(P p, P q, L l) {
  { ~p } -> L; // line not incident with p
};

/*
axiom(P p, P q, P r, L l) {
  l == {p, q} => I(p, l) && I(q, l);
}
*/

template <class P, class L = typename P::dual>
concept bool Cayley_Klein_plane =
  Cayley_Klein_plane_h<P, L> && Cayley_Klein_plane_h<L, P>;

template <class P>
concept bool Cayley_Klein_plane2 = Cayley_Klein_plane<P>; // Make the compiler happy


} // namespace fun

#endif
