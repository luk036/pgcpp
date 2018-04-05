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
  using _Base = std::array<_K, 3>;
  using _Self = pg_point<_K>;

public:
  using value_type = _K;
  using dual = pg_line<_K>;
  
  /**
   * @brief Construct a new pg point object
   * 
   * @param a array of coordinates
   */
  constexpr explicit pg_point(const _Base &a) : _Base{a} {}

  /**
   * @brief Construct a new pg_point object
   * 
   * @param x 
   * @param y 
   * @param z 
   */
  constexpr pg_point(const _K &x, const _K &y, const _K &z) 
      : _Base{x, y, z} {}

  /**
   * @brief Construct a new pg_point object by meet of two lines (p. 53)
   * 
   * @param l 
   * @param m 
   */
  // constexpr pg_point(const pg_line<_K> &l, const pg_line<_K> &m)
  //    : _Base{cross(l, m)} {}

  // Operators:

  /**
   * @brief Equal to
   * 
   * @param rhs
   * @return true if this point is equivalent to the rhs
   * @return false otherwise
   */
  constexpr bool operator==(const _Self &rhs) const {
    return cross(*this, rhs) == _Base({0, 0, 0});
  }

  /**
   * @brief Not Equal to
   * 
   * @param rhs
   * @return true if this point is not equivalent to the rhs
   * @return false otherwise
   */
  constexpr bool operator!=(const _Self &rhs) const {
    return ! (*this == rhs);
  }

  /// Return the dot product
  constexpr auto dot(const dual &l) const {
    return fun::dot_c(*this, l);
  }

  /// Return true if a line @a l incident with point @a p
  constexpr bool incident(const dual &l) const {
    return this->dot(l) == 0;
  }

  /**
   * @brief Equal to
   * 
   * @param rhs
   * @return true if this point is equivalent to the rhs
   * @return false otherwise
   */
  constexpr dual operator*(const _Self &rhs) const {
    return dual{cross(*this, rhs)};
  }

  ///  Return new line not incident with p
  //constexpr dual aux() { return dual{*this}; }
};


/// Return join of two points.
template <typename _K>
auto join(const pg_point<_K> &p, const pg_point<_K> &q) {
    return p * q;
}

template <typename _K>
auto plucker(const _K &lambda1, const pg_point<_K> &p, 
             const _K &mu1, const pg_point<_K> &q) {
  return pg_point<_K>{plucker_c(lambda1, p, mu1, q)};
}

///  Insertion operator for point values.
template <typename _K, class _Stream>
_Stream &operator<<(_Stream &os, const pg_point<_K> &p) {
  os << '[' << p[0] << ',' << p[1] << ',' << p[2] << ']';
  return os;
}

// template deduction guides (C++17)
template <typename _K>
pg_point(const std::array<_K,3> ) -> pg_point<_K>;


template <typename _K>
pg_point(const _K &, const _K &, const _K &) -> pg_point<_K>; 

} // namespace fun

#endif
