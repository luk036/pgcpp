#ifndef _HOME_UBUNTU_GITHUB_PGCPP_PG_COMMON_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_PG_COMMON_HPP 1

#include <array>
#include <tuple>

namespace fun {

/**
 * @brief 1st term of Cross product
 * 
 * @tparam _K 
 * @param v 
 * @param w 
 * @return 1st term of Cross product
 */
template <typename _K>
auto cross0 (const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return v[1] * w[2] - w[1] * v[2];
}

/**
 * @brief 2nd term of Cross product
 * 
 * @tparam _K 
 * @param v 
 * @param w 
 * @return 2nd term of Cross product
 */
template <typename _K>
auto cross1(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return v[0] * w[2] - w[0] * v[2];
}

/**
 * @brief 3rd term of Cross product
 * 
 * @tparam _K 
 * @param v 
 * @param w 
 * @return 3rd term of Cross product 
 */
template <typename _K>
auto cross2(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return v[0] * w[1] - w[0] * v[1];
}

/**
 * @brief Cross product
 * 
 * @tparam _K 
 * @param v 
 * @param w 
 * @return Cross product 
 */
template <typename _K>
auto cross(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return std::array<_K, 3>({cross0(v,w), -cross1(v,w), cross2(v,w)});
}

/**
 * @brief Dot product
 * 
 * @tparam _K 
 * @param v 
 * @param w 
 * @return auto 
 */
template <typename _K>
auto dot_c(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  auto [x1, y1, z1] = v;
  auto [x2, y2, z2] = w;
  return x1*x2 + y1*y2 + z1*z2;
}

/**
 * @brief generic Plucker function
 * 
 * @tparam _K data type
 * @param ld lamda
 * @param v 
 * @param mu 
 * @param w 
 * @return lamda*v + mu*w
 */
template <typename _K>
auto plucker_c(const _K &ld, const std::array<_K, 3> &v1, const _K &mu,
               const std::array<_K, 3> &v2) {
  auto [x1, y1, z1] = v1;
  auto [x2, y2, z2] = v2;
  return std::array<_K, 3>(
      {ld * x1 + mu * x2, ld * y1 + mu * y2, ld * z1 + mu * z2});
}

/**
 * @brief dot product of the (0,1)-component of two vectors 
 * 
 * @tparam _K 
 * @param v 
 * @param w 
 * @return auto 
 */
template <typename _K>
auto dot1(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return v[0] * w[0] + v[1] * w[1];
}

/**
 * @brief dot product of the (0,2)-component of two vectors 
 * 
 * @tparam _K 
 * @param v 
 * @param w 
 * @return auto 
 */
template <typename _K>
auto dot2(const std::array<_K, 3> &v, const std::array<_K, 3> &w) {
  return v[0] * w[0] + v[2] * w[2];
}

/**
 * @brief Square function
 * 
 * @tparam T data type
 * @param a input value
 * @return a^2
 */
template <typename T>
constexpr auto sq(T&& a) { return a*a; }

} // namespace fun 

#endif