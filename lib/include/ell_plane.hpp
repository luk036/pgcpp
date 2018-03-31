#ifndef _HOME_UBUNTU_GITHUB_PGCPP_ELL_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_ELL_PLANE_HPP 1

//#include "proj_plane_concepts.h"
#include "ck_plane.hpp"
#include "pg_line.hpp"
#include "pg_point.hpp"
#include <cassert>

namespace fun {

namespace CK {

template <typename _K> pg_line<_K> operator~(const pg_point<_K>& p) {
  return pg_line<_K>(p);
}

template <typename _K> pg_point<_K> operator~(const pg_line<_K>& l) {
  return pg_point<_K>(l);
}

namespace ELL {

template <typename _K>
auto check_cross_TQF(const _K &q1, const _K &q2, const _K &q3) {
  return sq(q1 + q2 + q3) - 2*(q1*q1 + q2*q2 + q3*q3) - 4*q1*q2*q3;
}

template <typename _K>
bool check_cross_law(const _K &s1, const _K &s2, const _K &s3, const _K &q3) {
  return sq(s1 * s2 * q3 - (s1 + s2 + s3) + 2) ==
         4 * (1 - s1) * (1 - s2) * (1 - s3);
}

} // namespace ELL

} // namespace CK

} // namespace fun

#endif
