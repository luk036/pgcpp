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
template <typename _K, typename pg_dual>
class pg_object : public std::array<_K, 3> {
  /// Value typedef.
  using _Base = std::array<_K, 3>;
  using _Self = pg_object<_K, pg_dual>;

public:
  using value_type = _K;
  using dual = pg_dual;

  /**
   * @brief Construct a new pg object object
   *
   * @param a array of coordinates
   */
  constexpr explicit pg_object(const _Base &a) : _Base{a} {}

  /**
   * @brief Construct a new pg_object object
   *
   * @param x
   * @param y
   * @param z
   */
  constexpr pg_object(const _K &x, const _K &y, const _K &z)
      : _Base{x, y, z} {}

  /**
   * @brief Construct a new pg_object object by meet of two lines (p. 53)
   *
   * @param l
   * @param m
   */
  // constexpr pg_object(const dual &l, const dual &m)
  //     : _Base{cross(l, m)} {}

  // Operators:

  /**
   * @brief Equal to
   *
   * @param rhs
   * @return true if this object is equivalent to the rhs
   * @return false otherwise
   */
  constexpr bool operator==(const _Self &rhs) const {
    if (this == &rhs) return true;
    return cross(*this, rhs) == _Base({0, 0, 0});
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
  constexpr auto dot(const dual &l) const {
    return fun::dot_c(*this, l);
  }

  /// Return true if a line @a l incident with object @a p
  constexpr bool incident(const dual &l) const {
    return this->dot(l) == 0;
  }

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
  constexpr dual aux() { return dual(*this); }
};


// template <typename _K, typename pg_dual>
// auto plucker(const _K &l, const pg_object<_K, pg_dual> &p1, const _K &m,
//              const pg_object<_K, pg_dual> &p2) {
//   return pg_object<_K, pg_dual>{plucker_c(l, p1, m, p2)};
// }
template <typename P, typename _K = typename P::value_type>
auto plucker(const _K &lambda1, const P &p, const _K &mu1, const P &q) {
  return P(plucker_c(lambda1, p, mu1, q));
}


///  Insertion operator for object values.
template <typename _K, typename pg_dual, class _Stream>
_Stream &operator<<(_Stream &os, const pg_object<_K, pg_dual> &p) {
  os << '(' << p[0] << ':' << p[1] << ':' << p[2] << ')';
  return os;
}


} // namespace fun

#endif
