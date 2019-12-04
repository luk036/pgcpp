// The template and inpoints for the -*- C++ -*- 3d line classes.
// Initially implemented by Wai-Shing Luk <luk036@gmail.com>
//

/*! @file include/pg_line.hpp
 *  This is a C++ Library header.
 */

#pragma once

#include "pg_common.hpp"
#include "pg_object.hpp"

namespace fun
{

// Forward declarations.

template <CommutativeRing _K>
class pg_point;

/*!
 * @brief Projective line: two dimensional subspace of K^3
 *
 * @tparam  _K  Type of line elements
 */
template <CommutativeRing _K>
struct pg_line : pg_object<_K, pg_point<_K>>
{
    /// Value typedef.
    using _Base = pg_object<_K, pg_point<_K>>;
    using _Base2 = std::array<_K, 3>;
    using value_type = _K;

    pg_line(const pg_line<_K>&) = delete;
    pg_line(pg_line<_K>&&) = default;

    /*!
     * @brief Construct a new pg object object
     *
     * @param a array of coordinates
     */
    constexpr explicit pg_line(const _Base2& a)
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
    constexpr pg_line(const _K& x, const _K& y, const _K& z)
        : _Base {_Base2 {x, y, z}}
    {
    }
};

/// Return meet of two lines.
template <CommutativeRing _K>
constexpr pg_point<_K> meet(const pg_line<_K>& l, const pg_line<_K>& m)
{
    return l * m;
}

// template deduction guides (C++17)
// template <CommutativeRing _K>
// pg_line(const std::array<_K, 3>)->pg_line<_K>;

// CommutativeRing{_K} pg_line(const _K &, const _K &, const _K &)->pg_line<_K>;

} // namespace fun
