#ifndef _HOME_UBUNTU_CUBSTORE_PROJ_GEOM_PGCPP_PG_COMMON_HPP
#define _HOME_UBUNTU_CUBSTORE_PROJ_GEOM_PGCPP_PG_COMMON_HPP 1

#include <array>
#include <tuple>

namespace fun {

template <typename _K>
auto cross0 (const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return v[1] * w[2] - w[1] * v[2];
}

template <typename _K>
auto cross1(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return v[0] * w[2] - w[0] * v[2];
}

template <typename _K>
auto cross2(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return v[0] * w[1] - w[0] * v[1];
}

template <typename _K>
auto cross(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return std::array<_K, 3>({cross0(v,w), -cross1(v,w), cross2(v,w)});
}

template <typename _K>
auto dot_c(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  auto [x1, y1, z1] = v;
  auto [x2, y2, z2] = w;
  return x1*x2 + y1*y2 + z1*z2;
}

template <typename _K>
auto plucker_c(const _K &ld, const std::array<_K, 3> &v, const _K &mu,
               const std::array<_K, 3> &w) {
  auto [x1, y1, z1] = v;
  auto [x2, y2, z2] = w;
  //_K x1 = v[0], y1 = v[1], z1 = v[2];
  //_K x2 = w[0], y2 = w[1], z2 = w[2];
  return std::array<_K, 3>(
      {ld * x1 + mu * x2, ld * y1 + mu * y2, ld * z1 + mu * z2});
}

template <typename _K>
auto dot1(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return v[0] * w[0] + v[1] * w[1];
}

template <typename _K>
auto dot2(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return v[0] * w[0] + v[2] * w[2];
}


template <typename T>
constexpr auto sq(T&& a) { return a*a; }

} // namespace fun 

#endif