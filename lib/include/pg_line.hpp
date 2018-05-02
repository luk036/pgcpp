// The template and inpoints for the -*- C++ -*- 3d line classes.
// Initially implemented by Wai-Shing Luk <luk036@gmail.com>
//

/** @file include/pg_line.hpp
 *  This is a C++ Library header.
 */

#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PG_LINE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PG_LINE_HPP 1

#include "pg_common.hpp"
#include "pg_object.hpp"

namespace fun {

// Forward declarations.

template <typename _K> class pg_point;

/**
 * @brief Projective line: two dimensional subspace of K^3
 *
 * @tparam  _K  Type of line elements
 */
template <typename _K> class pg_line : public pg_object<_K, pg_point<_K>> {
  /// Value typedef.
  using _Base = pg_object<_K, pg_point<_K>>;
  using _Base2 = std::array<_K, 3>;

public:
  using value_type = _K;

  /**
   * @brief Construct a new pg object object
   *
   * @param a array of coordinates
   */
  constexpr explicit pg_line(const _Base2 &a) : _Base{a} {}

  /**
   * @brief Construct a new pg_object object
   *
   * @param x
   * @param y
   * @param z
   */
  constexpr pg_line(const _K &x, const _K &y, const _K &z) : _Base{x, y, z} {}
};

/// Return meet of two lines.
template <typename _K> auto meet(const pg_line<_K> &l, const pg_line<_K> &m) {
  return l * m;
}

// template deduction guides (C++17)
template <typename _K> pg_line(const std::array<_K, 3>)->pg_line<_K>;

template <typename _K> pg_line(const _K &, const _K &, const _K &)->pg_line<_K>;

} // namespace fun

#endif
