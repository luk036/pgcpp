// The template and inlines for the -*- C++ -*- 3d point classes.
// Initially implemented by Wai-Shing Luk <luk036@gmail.com>
//

/** @file include/pg_point.hpp
 *  This is a C++ Library header.
 */

#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PG_POINT_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PG_POINT_HPP 1

#include "pg_common.hpp"
#include "pg_object.hpp"


namespace fun {

// Forward declarations.

template <typename _K> class pg_line;

template <typename _K>
class pg_point : public pg_object< _K, pg_line<_K> > {
  /// Value typedef.
  using _Base = pg_object< _K, pg_line<_K> >;
  using _Base2 = std::array<_K, 3>;
  using _Self = pg_point<_K>;

public:
  using value_type = _K;

  /**
   * @brief Construct a new pg object object
   *
   * @param a array of coordinates
   */
  constexpr explicit pg_point(const _Base2 &a) : _Base{a} {}

  /**
   * @brief Construct a new pg_object object
   *
   * @param x
   * @param y
   * @param z
   */
  constexpr pg_point(const _K &x, const _K &y, const _K &z)
      : _Base{x, y, z} {}
};

/// Return join of two points.
template <typename _K>
auto join(const pg_point<_K> &p, const pg_point<_K> &q) {
    return p * q;
}

// template <typename _K>
// auto plucker(const _K &lambda1, const pg_point<_K> &p,
//              const _K &mu1, const pg_point<_K> &q) {
//   return pg_point<_K>(plucker_c(lambda1, p, mu1, q));
// }

// ///  Insertion operator for point values.
// template <typename _K, class _Stream>
// _Stream &operator<<(_Stream &os, const pg_point<_K> &p) {
//   os << '[' << p[0] << ',' << p[1] << ',' << p[2] << ']';
//   return os;
// }

// template deduction guides (C++17)
template <typename _K>
pg_point(const std::array<_K,3> &) -> pg_point<_K>;

template <typename _K>
pg_point(const _K &, const _K &, const _K &) -> pg_point<_K>;

} // namespace fun

#endif
