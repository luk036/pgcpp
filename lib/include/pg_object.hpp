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
 * @brief Projective object: one dimensional subspace of K^3
 *
 * @tparam _K Type of object elements
 */
template <typename _K, typename dual>
class pg_object : public std::array<_K, 3> {
  /// Value typedef.
  typedef std::array<_K, 3> _Base;

public:
  using value_type = _K;

  /**
   * @brief Construct a new pg object object
   *
   * @param a array of coordinates
   */
  explicit pg_object(const _Base &a) : _Base{a} {}

  /**
   * @brief Construct a new pg_object object
   *
   * @param x
   * @param y
   * @param z
   */
  constexpr explicit pg_object(const _K &x, const _K &y, const _K &z)
      : _Base{x, y, z} {}

  /**
   * @brief Construct a new pg_object object by meet of two lines (p. 53)
   *
   * @param l
   * @param m
   */
  constexpr pg_object(const dual &l, const dual &m)
      : _Base{cross(l, m)} {}

  // Operators:

  /**
   * @brief Equal to
   *
   * @param rhs
   * @return true if this object is equivalent to the rhs
   * @return false otherwise
   */
  constexpr bool operator==(const pg_object<_K> &rhs) const {
    return cross(*this, rhs) == _Base({0, 0, 0});
  }

  /// Return the dot product
  constexpr auto dot(const pg_object<_K> &l) const {
    return fun::dot_c(*this, l);
  }

  /// Return true if a line @a l incident with object @a p
  constexpr bool incident(const dual &l) const {
    return this->dot(l) == 0;
  }

  ///  Return new line not incident with p
  constexpr dual aux() { return dual{*this}; }
};

template <typename _K>
auto plucker(const _K &l, const pg_object<_K> &p1, const _K &m,
             const pg_object<_K> &p2) {
  return pg_object<_K>{plucker_c(l, p1, m, p2)};
}


///  Insertion operator for object values.
template <typename _K, class _Stream>
_Stream &operator<<(_Stream &os, const pg_object<_K> &p) {
  os << '(' << p[0] << ':' << p[1] << ':' << p[2] << ')';
  return os;
}

// template deduction guides (C++17)
template <typename _K>
pg_object(const std::array<_K, 3> &) -> pg_object<_K>;

template <typename _K>
pg_object(const _K &, const _K &, const _K &) -> pg_object<_K>;

} // namespace fun

#endif
