#ifndef _HOME_UBUNTU_GITHUB_PGCPP_FRACTION_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_FRACTION_HPP 1

#include <cassert>
#include <type_traits>
#include <cmath>

namespace fun {

template <typename Z>
  requires std::is_integral<Z>::value
inline constexpr Z gcd(const Z &a, const Z &b) noexcept {
  return b == Z(0) ? std::abs(a) : gcd(b, a % b);
}

template <typename Z>
//  requires std::is_integral<Z>::value
inline constexpr Z lcm(const Z &a, const Z &b) noexcept {
  return a / gcd(a, b) * b;
}

template <typename Z> class Fraction {
  using _Self = Fraction<Z>;

public:
  Z _numerator;
  Z _denominator;

public:
  constexpr Fraction(const Z &numerator, const Z &denominator) {
    auto common = gcd(numerator, denominator);
    _numerator = numerator / common;
    _denominator = denominator / common;
  }

  constexpr auto operator+(const _Self &frac) const {
    auto common = lcm(_denominator, frac._denominator);
    auto n = common / _denominator * _numerator +
        common / frac._denominator * frac._numerator;
    return _Self(n, common);
  }

  constexpr auto operator-(const _Self &frac) const {
    return *this + (-frac);
  }

  constexpr auto operator-() const {
    return _Self(-_numerator, _denominator);
  }

  constexpr auto abs() const {
    return _Self(std::abs(_numerator), std::abs(_denominator));
  }

  constexpr auto operator*(const _Self &frac) const {
    return _Self(_numerator * frac._numerator,
                    _denominator * frac._denominator);
  }

  constexpr auto operator/(const _Self &frac) const {
    return *this * frac.reciprocal();
  }

  constexpr auto reciprocal() const {
    return _Self(_denominator, _numerator);
  }

  /**
   * @brief Three way comparison
   * 
   * @param frac 
   * @return constexpr auto 
   */
  constexpr auto cmp(const _Self &frac) const {
    return _numerator * frac._denominator - _denominator * frac._numerator;
  }

  constexpr bool operator==(const _Self &frac) const {
    return this->cmp(frac) == 0;
  }

  constexpr bool operator<(const _Self &frac) const {
    return this->cmp(frac) < 0;
  }

  constexpr auto cmp(const Z &c) const {
    return _numerator - _denominator * c;
  }

  constexpr bool operator==(const Z &c) const {
    return this->cmp(c) == 0;
  }

  constexpr bool operator<(const Z &c) const {
    return this->cmp(c) < 0;
  }
};

template <typename Z>
  requires std::is_integral<Z>::value
constexpr auto operator-(const Z &c, const Fraction<Z> &frac) {
  return Fraction<Z>(frac._denominator*c - frac._numerator, frac._denominator);
}

// template <typename _Stream, typename Z>
//   requires std::is_integral<Z>::value
// _Stream operator<<(_Stream& os, const Fraction<Z> &frac) {
//   os << frac._numerator << "/" << frac._denominator;
//   return os;
// }

// For template deduction
template <typename Z> Fraction(const Z &, const Z &) -> Fraction<Z>;

} // namespace fun

#endif