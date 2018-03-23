// The template and inlines for the -*- C++ -*- 3d point classes.
// Initially implemented by Wai-Shing Luk <luk036@gmail.com>
//

/** @file include/pg_point.hpp
 *  This is a C++ Library header.
 */

#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PG_POINT_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PG_POINT_HPP 1

#include "pg_common.hpp"
#include "pg_line.hpp"

namespace fun {

// Forward declarations.

template <typename _K> class pg_line;

/**
 * @brief Projective point: one dimensional subspace of K^3
 * 
 * @tparam _K Type of point elements
 */
template <typename _K> class pg_point : public std::array<_K, 3> {
  /// Value typedef.
  typedef std::array<_K, 3> _Base;

public:
  using value_type = _K;
  using dual = pg_line<_K>;

  /// Construct by meet of two lines @a l and @a m. (p. 53)
  constexpr explicit pg_point(const _Base &a) : _Base{a} {}

  /**
   * @brief Construct by data
   * 
   */

  /**
   * @brief Construct a new pg_point object
   * 
   * @param x 
   * @param y 
   * @param z 
   */
  constexpr explicit pg_point(const _K &x, const _K &y, const _K &z) 
      : _Base{x, y, z} {}

  /**
   * @brief Construct a new pg_point object by meet of two lines (p. 53)
   * 
   * @param l 
   * @param m 
   */
  constexpr pg_point(const pg_line<_K> &l, const pg_line<_K> &m)
      : _Base{cross(l, m)} {}

  // Operators:

  /**
   * @brief Equal to
   * 
   * @param rhs
   * @return true if this point is equivalent to the rhs
   * @return false otherwise
   */
  constexpr bool operator==(const pg_point<_K> &rhs) const {
    return cross(*this, rhs) == _Base({0, 0, 0});
  }

  /// Return the dot product
  constexpr auto dot(const pg_line<_K> &l) const {
    return fun::dot_c(*this, l);
  }

  /// Return true if a line @a l incident with point @a p
  constexpr bool incident(const pg_line<_K> &l) const {
    return this->dot(l) == 0;
  }

  ///  Return new line not incident with p
  constexpr pg_line<_K> aux() { return pg_line<_K>{*this}; }
};

template <typename _K>
auto plucker(const _K &l, const pg_point<_K> &p1, const _K &m,
             const pg_point<_K> &p2) {
  return pg_point<_K>{plucker_c(l, p1, m, p2)};
}


///  Insertion operator for point values.
template <typename _K, class _Stream>
_Stream &operator<<(_Stream &os, const pg_point<_K> &p) {
  os << '(' << p[0] << ':' << p[1] << ':' << p[2] << ')';
  return os;
}

// template deduction guides (C++17)
template <typename _K>
pg_point(const std::array<_K, 3> &) -> pg_point<_K>;

template <typename _K>
pg_point(const _K &, const _K &, const _K &) -> pg_point<_K>; 

} // namespace fun

#endif
