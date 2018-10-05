// The template and inlines for the -*- C++ -*- 3d object classes.
// Initially implemented by Wai-Shing Luk <luk036@gmail.com>
//

/** @file include/pg_object.hpp
 *  This is a C++ Library header.
 */

#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PG_OBJECT_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PG_OBJECT_HPP 1

#include "pg_common.hpp"

namespace fun {

/**
 * @brief Projective object
 *
 * @tparam _K Type of object elements
 * @tparam _dual
 */
template <typename _K, typename _dual>
class pg_object : public std::array<_K, 3> {
    /// Value typedef.
    using _Base = std::array<_K, 3>;
    using _Self = pg_object<_K, _dual>;

  public:
    using value_type = _K;
    using dual = _dual;

    /**
     * @brief Construct a new pg object object
     *
     * @param a array of coordinates
     */
    constexpr explicit pg_object(const _Base &a) : _Base{a} {}

    // Operators:

    /**
     * @brief Equal to
     *
     * @param rhs
     * @return true if this object is equivalent to the rhs
     * @return false otherwise
     */
    constexpr bool operator==(const _Self &rhs) const {
        if (this == &rhs) {
            return true;
        }
        return cross(*this, rhs) == _Base{_K(0), _K(0), _K(0)};
    }

    /**
     * @brief Not equal to
     *
     * @param rhs
     * @return true if this object is not equivalent to the rhs
     * @return false otherwise
     */
    constexpr bool operator!=(const _Self &rhs) const {
        return !(*this == rhs);
    }

    /// Return the dot product
    constexpr auto dot(const dual &l) const { return fun::dot_c(*this, l); }

    /**
     * @brief Join or meet
     *
     * @param rhs
     * @return true if this point is equivalent to the rhs
     * @return false otherwise
     */
    constexpr dual operator*(const _Self &rhs) const {
        return dual(cross(*this, rhs));
    }

    ///  Return new line not incident with p
    constexpr dual aux() const { return dual(*this); }
};

template <typename P, typename _K = typename P::value_type>
constexpr P plucker(const _K &lambda1, const P &p, const _K &mu1, const P &q) {
    return P(plucker_c(lambda1, p, mu1, q));
}

///  Insertion operator for object values.
template <typename _K, typename _dual, class _Stream>
_Stream &operator<<(_Stream &os, const pg_object<_K, _dual> &p) {
    os << '(' << p[0] << ':' << p[1] << ':' << p[2] << ')';
    return os;
}

} // namespace fun

#endif
