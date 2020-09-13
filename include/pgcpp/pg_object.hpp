// The template and inlines for the -*- C++ -*- pg object classes.
// Initially implemented by Wai-Shing Luk <luk036@gmail.com>
//

/*! @file include/pg_object.hpp
 *  This is a C++ Library header.
 */

#pragma once

#include "pg_common.hpp"

namespace fun
{

/*!
 * @brief Projective object
 *
 * @tparam _K Type of object elements
 * @tparam _dual
 */
template <CommutativeRing _K, typename _dual>
class pg_object : public std::array<_K, 3>
{
    /// Value typedef.
    using _Base = std::array<_K, 3>;
    using _Self = pg_object<_K, _dual>;

  public:
    using value_type = _K;
    using dual = _dual;

    // pg_object(_Self &&) = default;

    /*!
     * @brief Construct a new pg object object
     *
     * @param[in] a array of coordinates
     */
    constexpr explicit pg_object(const _Base& a)
        : _Base {a}
    {
    }

    /*!
     * @brief Construct a new pg object
     *
     */
    pg_object(const _Self&) = delete;

    /*!
     * @brief
     *
     * @return _Self&
     */
    _Self& operator=(const _Self&) = delete;

    /*!
     * @brief Construct a new pg object
     *
     */
    pg_object(_Self&&) = default;

    /*!
     * @brief
     *
     * @return _Self&
     */
    _Self& operator=(_Self&&) = default;

    // Operators:

    /*!
     * @brief Equal to
     *
     * @param[in] rhs
     * @return true if this object is equivalent to the rhs
     * @return false otherwise
     */
    friend constexpr bool operator==(const _Self& lhs, const _Self& rhs)
    {
        if (&lhs == &rhs)
        {
            return true;
        }
        return cross(lhs, rhs) == _Base {_K(0), _K(0), _K(0)};
    }

    /*!
     * @brief Not equal to
     *
     * @param[in] rhs
     * @return true if this object is not equivalent to the rhs
     * @return false otherwise
     */
    friend constexpr bool operator!=(const _Self& lhs, const _Self& rhs)
    {
        return !(lhs == rhs);
    }

    /*!
     * @brief Equal to
     *
     * @param[in] rhs
     * @return true if this object is equivalent to the rhs
     * @return false otherwise
     */
    constexpr bool is_NaN() const
    {
        const _Base& base = *this;
        return base == _Base {_K(0), _K(0), _K(0)};
    }


    /*!
     * @brief the dot product
     *
     * @param[in] l
     * @return _K
     */
    constexpr _K dot(const dual& l) const
    {
        return fun::dot_c(*this, l);
    }

    /*!
     * @brief Generate a new line not incident with p
     *
     * @return dual
     */
    constexpr dual aux() const
    {
        return dual(*this);
    }

    /*!
     * @brief Join or meet
     *
     * @param[in] rhs
     * @return true if this point is equivalent to the rhs
     * @return false otherwise
     */
    friend constexpr dual operator*(const _Self& lhs, const _Self& rhs)
    {
        return dual(cross(lhs, rhs));
    }
};

/*!
 * @brief
 *
 * @tparam P
 * @tparam Value_type<P>
 * @param[in] lda1
 * @param[in] p
 * @param[in] mu1
 * @param[in] q
 * @return P
 */
template <typename P, CommutativeRing _K = Value_type<P>>
constexpr P plucker(const _K& lda1, const P& p, const _K& mu1, const P& q)
{
    return P {plucker_c(lda1, p, mu1, q)};
}

/*!
 * @brief
 *
 * @tparam _K
 * @tparam _dual
 * @tparam _Stream
 * @param[in] os
 * @param[in] p
 * @return _Stream&
 */
template <CommutativeRing _K, typename _dual, class _Stream>
_Stream& operator<<(_Stream& os, const pg_object<_K, _dual>& p)
{
    os << '(' << p[0] << ':' << p[1] << ':' << p[2] << ')';
    return os;
}

} // namespace fun
