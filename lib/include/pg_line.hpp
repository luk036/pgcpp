// The template and inpoints for the -*- C++ -*- 3d line classes.
// Initially implemented by Wai-Shing Luk <luk036@gmail.com>
//

/** @file include/pg_line.hpp
 *  This is a C++ Library header.
 */

#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PG_LINE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PG_LINE_HPP 1

#include "pg_common.hpp"
#include "pg_point.hpp"

namespace fun {

// Forward declarations.

template <typename _K> class pg_point;

/**
 * @brief Projective line: two dimensional subspace of K^3
 * 
 * @tparam  _K  Type of line elements 
 */
template <typename _K> class pg_line : public std::array<_K, 3> {
  /// Value typedef.
  using _Base = std::array<_K, 3>;
  using _Self = pg_line<_K>;

public:
  using value_type = _K;
  using dual = pg_point<_K>;

  /**
   * @brief Construct by an array
   * 
   */
  constexpr explicit pg_line(const _Base &a) : _Base{a} {}

  /**
   * @brief Construct by data
   * 
   */
  // constexpr explicit pg_line(const _K &a, const _K &b, const _K &c) 
  //    : _Base{a, b, c} {}

  /// Construct by meet of two points @a p and @a q. (p. 53)
  // constexpr pg_line(const dual &p, const dual &q)
  //     : _Base{cross(p, q)} {}

  // Operators:

  /// Return true if @a p is equivalent to @a q (in projective sense).
  constexpr bool operator==(const _Self &other) const {
    return cross(*this, other) == _Base({0, 0, 0});
  }

  /// Return true if @a p is equivalent to @a q (in projective sense).
  constexpr bool operator!=(const _Self &other) const {
    return !(*this == other);
  }

  /// Return the dot product
  constexpr auto dot(const dual &l) const {
    return dot_c(*this, l);
  }

  /// Return true if a point @a l incident with line @a p
  constexpr bool incident(const dual &p) const {
    return this->dot(p) == 0;
  }

  /// Return meet of two line.
  constexpr dual operator*(const _Self &other) const {
    return dual{cross(*this, other)};
  }

  ///  Return new point not incident with p
  // constexpr dual aux() { return dual{*this}; }
};


/// Return meet of two lines.
template <typename _K>
auto meet(pg_line<_K> &l, pg_line<_K> &m) {
    return l * m;
}


template <typename _K>
auto plucker(const _K &lambda1, const pg_line<_K> &l, 
             const _K &mu1, const pg_line<_K> &m) {
  return pg_line<_K>{plucker_c(lambda1, l, mu1, m)};
}


///  Insertion operator for line values.
template <typename _K, class _Stream>
_Stream &operator<<(_Stream &os, const pg_line<_K> &l) {
  os << '<' << l[0] << ':' << l[1] << ':' << l[2] << '>';
  return os;
}


// template deduction guides (C++17)
template <typename _K>
pg_line(const std::array<_K, 3>) -> pg_line<_K>;


//template <typename _K>
//pg_line(const _K &, const _K &, const _K &) -> pg_line<_K>; 

} // namespace fun

#endif
