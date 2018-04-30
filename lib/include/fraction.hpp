#ifndef _HOME_UBUNTU_GITHUB_PGCPP_FRACTION_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_FRACTION_HPP 1

#include <cassert>
#include <type_traits>

namespace fun {

template <typename Z>
  requires std::is_integral<Z>::value
inline constexpr Z gcd(const Z &a, const Z &b) noexcept {
  return b == Z(0) ? abs(a) : gcd(b, a % b);
}

template <typename Z>
//  requires std::is_integral<Z>::value
inline constexpr Z lcm(const Z &a, const Z &b) noexcept {
  return a / gcd(a, b) * b;
}

template <typename Z> class Fraction {
  using _Self = Fraction<Z>;

private:
  Z _numerator;
  Z _denominator;

public:
  constexpr Fraction(const Z &numerator, const Z &denominator) {
    auto common = gcd(numerator, denominator);
    _numerator = numerator / common;
    _denominator = denominator / common;
  }

  constexpr auto operator+(const _Self &frac) const {
    auto common = lcm(_denominator, frac.denominator);
    auto n = common / _denominator * _numerator +
        common / frac.denominator * frac.numerator;
    return _Self(n, common);
  }

  constexpr auto operator-(const _Self &frac) const {
    return *this + (-frac);
  }

  constexpr auto operator-() const {
    return _Self(-_numerator, _denominator);
  }

  constexpr auto abs() const {
    return _Self(abs(_numerator), abs(_denominator));
  }

  constexpr auto operator*(const _Self &frac) const {
    return _Self(_numerator * frac.numerator,
                    _denominator * frac.denominator);
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


  // constexpr auto operator float() {
  //     return float(_numerator / _denominator);
  // }

  // constexpr auto operator int() {
  //     return int(_numerator / _denominator);
  // }
};

// For template deduction
template <typename Z> Fraction(const Z &, const Z &) -> Fraction<Z>;

} // namespace fun

#endif