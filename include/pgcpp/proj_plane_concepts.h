/*! @file include/proj_plane_concepts.h
 *  This is a C++ Library header.
 */

#pragma once

#include "common_concepts.h"

/*!
 @todo: projectivity >=
**/

namespace fun
{

/*!
 * @brief Projective plane Concept (half)
 *
 * @tparam P Point
 * @tparam L Line
 */
template <class P, class L>
concept bool Projective_plane_prim_h = Equality_comparable<P>&& requires(
    P& p, P& q, L& l)
{
    // { P(p) } -> P; // copyable
    {
        incident(p, l)
    }
    ->bool; // incidence
    {
        p* q
    }
    ->L; // join or meet
    {
        p.aux()
    }
    ->L; // line not incident with p
};

/*!
 * @brief Projective plane Concept (full)
 *
 * @tparam P Point
 * @tparam L Line
 */
template <class P, class L = typename P::dual>
concept bool Projective_plane_prim =
    Projective_plane_prim_h<P, L>&& Projective_plane_prim_h<L, P>;


/*!
 * @brief Shorthand Notation of Projective_plane
 *
 * @tparam P Point
 */
template <class P>
concept bool Projective_plane_prim2 =
    Projective_plane_prim<std::remove_reference_t<P>>; // Make the compiler
                                                       // happy

/*!
 * @brief Projective plane Concept (half)
 *
 * @tparam P Point
 * @tparam L Line
 */
template <class P, class L>
concept bool Projective_plane_h = Equality_comparable<P>&& requires(
    P& p, P& q, L& l, Value_type<P> a)
{
    typename Value_type<P>;
    // { P(p) } -> P; // copyable
    // { incident(p, l) } -> bool; // incidence
    {
        p* q
    }
    ->L; // join or meet
    {
        p.aux()
    }
    ->L; // line not incident with p
    {
        p.dot(l)
    }
    ->Value_type<P>; // for measurement
    {
        plucker(a, p, a, q)
    }
    ->P; // vector computation
};

/*!
 * @brief Projective plane Concept (full)
 *
 * @tparam P Point
 * @tparam L Line
 */
template <class P, class L = typename P::dual>
concept bool Projective_plane =
    Projective_plane_h<P, L>&& Projective_plane_h<L, P>;

/*
axiom(P p, P q, P r, L l) {
  l == L{p, q} => I(p, l) and I(q, l);
}
*/


/*!
 * @brief Shorthand Notation of Projective_plane
 *
 * @tparam P Point
 */
template <class P>
concept bool Projective_plane2 =
    Projective_plane<std::remove_reference_t<P>>; // Make the compiler happy

/*!
 * @brief Projective plane Concept (half)
 *
 * @tparam P Point
 * @tparam L Line
 */
template <class P, class L>
concept bool Projective_plane_coord_h = Projective_plane_h<P, L>&& requires(
    P& p, size_t idx)
{
    {
        p[idx]
    }
    ->Value_type<P>; // for coordinate acess
};

/*!
 * @brief Projective plane Concept (full)
 *
 * @tparam P Point
 * @tparam L Line
 */
template <class P, class L = typename P::dual>
concept bool Projective_plane_coord =
    Projective_plane_coord_h<P, L>&& Projective_plane_coord_h<L, P>;

/*!
 * @brief Shorthand Notation of Projective_plane
 *
 * @tparam P Point
 */
template <class P>
concept bool Projective_plane_coord2 =
    Projective_plane_coord<std::remove_reference_t<P>>; // Make the compiler
                                                        // happy

} // namespace fun
