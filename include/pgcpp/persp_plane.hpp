#pragma once

#include "ck_plane.hpp"
#include "fractions.hpp"
// #include "pg_common.hpp"
#include "proj_plane.hpp" // import pg_point, involution, tri_func,

namespace fun
{

/*!
 * @brief
 *
 * @tparam P
 * @tparam P::dual
 */
template <typename P, typename L = typename P::dual>
requires Projective_plane_prim<P, L> // c++20 concept
class persp_euclid_plane : public ck<P, L, persp_euclid_plane>
{
    using K = Value_type<P>;

  private:
    P _Ire;
    P _Iim;
    L _l_infty;

  public:
    /*!
     * @brief Construct a new persp euclid plane object
     *
     * @param Ire
     * @param Iim
     * @param l_infty
     */
    constexpr persp_euclid_plane(P Ire, P Iim, L l_infty)
        : _Ire {std::move(Ire)}
        , _Iim {std::move(Iim)}
        , _l_infty {std::move(l_infty)}
    {
    }

    // /*!
    //  * @brief Construct a new persp euclid plane object
    //  *
    //  * @param Ire
    //  * @param Iim
    //  * @param l_infty
    //  */
    // constexpr persp_euclid_plane(const P &Ire, const P &Iim, const L
    // &l_infty)
    //     : _Ire{Ire}, _Iim{Iim}, _l_infty{l_infty} {}

    // /*!
    //  * @brief
    //  *
    //  * @param x
    //  * @return const L&
    //  */
    // constexpr const L &perp(const P &x) const { return _l_infty; }

    /*!
     * @brief
     *
     * @return const L&
     */
    constexpr const L& l_infty() const
    {
        return this->_l_infty;
    }

    /*!
     * @brief
     *
     * @param x
     * @return P
     */
    constexpr auto perp(const L& v) const -> P
    {
        const auto alpha = v.dot(this->_Ire);
        const auto beta = v.dot(this->_Iim);
        return plucker(alpha, this->_Ire, beta, this->_Iim);
    }

    /*!
     * @brief
     *
     * @param l
     * @param m
     * @return true
     * @return false
     */
    constexpr auto is_parallel(const L& l, const L& m) const -> bool
    {
        return incident(this->_l_infty, l * m);
    }

    /*!
     * @brief
     *
     * @param a
     * @param b
     * @return P
     */
    constexpr auto midpoint(const P& a, const P& b) const -> P
    {
        const auto alpha = a.dot(this->_l_infty);
        const auto beta = b.dot(this->_l_infty);
        return plucker(alpha, a, beta, b);
    }

    /*!
     * @brief
     *
     * @param tri
     * @return auto
     */
    constexpr auto tri_midpoint(const Triple<P>& tri) const
    {
        const auto& [a1, a2, a3] = tri;

        auto m12 = this->midpoint(a1, a2);
        auto m23 = this->midpoint(a2, a3);
        auto m13 = this->midpoint(a1, a3);
        return Triple<P> {std::move(m12), std::move(m23), std::move(m13)};
    }

    /*!
     * @brief
     *
     * @param x
     * @return K
     */
    constexpr auto omega(const P& x) const -> K
    {
        return sq(x.dot(this->_l_infty));
    }

    /*!
     * @brief
     *
     * @param x
     * @return K
     */
    constexpr auto omega(const L& x) const -> K
    {
        return sq(x.dot(this->_Ire)) + sq(x.dot(this->_Iim));
    }

    /*!
     * @brief
     *
     * @param a1
     * @param a2
     * @return auto
     */
    template <Projective_plane2 _P>
    constexpr auto measure(const _P& a1, const _P& a2) const
    {
        const auto omg = K(this->omega(a1 * a2));
        const auto den = K(this->omega(a1) * this->omega(a2));
        if constexpr (Integral<K>)
        {
            return Fraction<K>(omg, den);
        }
        else
        {
            return omg / den;
        }
    }

    // /*!
    //  * @brief
    //  *
    //  * @param l1
    //  * @param l2
    //  * @return auto
    //  */
    // constexpr auto cross_s(const L &l1, const L &l2) const {
    //     return 1 - this->spread(l1, l2); // ???
    // }
};

} // namespace fun
