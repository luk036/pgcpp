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

/**
 * @defgroup 3d_line 3d line in projective geometry
 * @ingroup geometry
 *
 * Classes and functions for 3d line.
 * @{
 */

// Forward declarations.

template <typename _K> class pg_point;

/**
 * @brief Projective line: one dimensional subspace of K^3
 * 
 * @tparam  _K  Type of line elements 
 */
template <typename _K> class pg_line : public std::array<_K, 3> {
  /// Value typedef.
  typedef std::array<_K, 3> _Base;

public:
  using value_type = _K;
  using dual = pg_point<_K>;

  // Lets the compiler synthesize the copy constructor
  // pg_line (const pg_line<_K>&);
  /// Return the base class.
  //constexpr _Base base() { return *this; }

  // Lets the compiler synthesize the copy constructor
  // pg_point (const pg_point<_K>&);
  /// Return the base class.
  //constexpr const _Base &base() const { return *this; }

  /**
   * @brief Construct by an array
   * 
   */
  constexpr explicit pg_line(const _Base &a) : _Base{a} {}

  /**
   * @brief Construct by data
   * 
   */
  constexpr explicit pg_line(const _K &a, const _K &b, const _K &c) 
      : _Base{a, b, c} {}

  /// Construct by meet of two points @a p and @a q. (p. 53)
  constexpr pg_line(const pg_point<_K> &p, const pg_point<_K> &q)
      : _Base{cross(p, q)} {}

  // Operators:

  /// Return true if @a p is equivalent to @a q (in projective sense).
  constexpr bool operator==(const pg_line<_K> &other) const {
    return cross(*this, other) == _Base({0, 0, 0});
  }

  /// Return the dot product
  constexpr auto dot(const pg_point<_K> &l) const {
    return dot_c(*this, l);
  }

  /// Return true if a point @a l incident with line @a p
  constexpr bool incident(const pg_point<_K> &p) const {
    return dot_c(*this, p) == 0;
  }

  ///  Return new point not incident with p
  constexpr pg_point<_K> aux() { return pg_point<_K>{*this}; }
};

/**
 * @brief Pl\u{"}cker
 * 
 * @tparam _K coordinate type
 * @param a 
 * @param l1 
 * @param b 
 * @param l2 
 * @return auto 
 */
template <typename _K>
auto plucker(const _K &a, const pg_line<_K> &l1, const _K &b,
             const pg_line<_K> &l2) {
  return pg_line<_K>{plucker_c(a, l1, b, l2)};
}

///  Insertion operator for line values.
template <typename _K, class _Stream>
_Stream &operator<<(_Stream &os, const pg_line<_K> &l) {
  os << '<' << l[0] << ':' << l[1] << ':' << l[2] << '>';
  return os;
}

// template deduction guides (C++17)
template <typename _K>
pg_line(const std::array<_K, 3> &) -> pg_line<_K>;

template <typename _K>
pg_line(const _K &, const _K &, const _K &) -> pg_line<_K>; 

} // namespace fun

#endif
