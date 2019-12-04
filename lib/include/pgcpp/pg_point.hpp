// The template and inlines for the -*- C++ -*- pg point classes.
// Initially implemented by Wai-Shing Luk <luk036@gmail.com>
//

/*! @file include/pg_point.hpp
 *  This is a C++ Library header.
 */

#pragma once

#include "pg_common.hpp"
#include "pg_object.hpp"

namespace fun
{

// Forward declarations.
template <CommutativeRing _K>
class pg_line;

template <CommutativeRing _K>
struct pg_point : pg_object<_K, pg_line<_K>>
{
    /// Value typedef.
    using _Base = pg_object<_K, pg_line<_K>>;
    using _Base2 = std::array<_K, 3>;
    using value_type = _K;

    /*!
     * @brief Construct a new pg point object
     *
     */
    pg_point(const pg_point<_K>&) = delete;

    /*!
     * @brief Construct a new pg point object
     *
     */
    pg_point(pg_point<_K>&&) = default;

    /*!
     * @brief
     *
     * @return pg_point<_K>&
     */
    pg_point<_K>& operator=(const pg_point<_K>&) = delete;

    /*!
     * @brief
     *
     * @return pg_point<_K>&
     */
    pg_point<_K>& operator=(pg_point<_K>&&) = default;

    /*!
     * @brief Construct a new pg object object
     *
     * @param a array of coordinates
     */
    constexpr explicit pg_point(const _Base2& a)
        : _Base {a}
    {
    }

    /*!
     * @brief Construct a new pg_object object
     *
     * @param x
     * @param y
     * @param z
     */
    constexpr pg_point(const _K& x, const _K& y, const _K& z)
        : _Base {_Base2 {x, y, z}}
    {
    }
};

/*!
 * @brief Return join of two points.
 *
 * @param p
 * @param q
 * @return pg_line<_K>
 */
template <CommutativeRing _K>
constexpr pg_line<_K> join(const pg_point<_K>& p, const pg_point<_K>& q)
{
    return p * q;
}

// template deduction guides (C++17)
// CommutativeRing{_K} pg_point(const std::array<_K, 3> &)->pg_point<_K>;

// CommutativeRing{_K} pg_point(const _K &, const _K &, const _K
// &)->pg_point<_K>;

} // namespace fun
