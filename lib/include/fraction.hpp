#ifndef _HOME_UBUNTU_GITHUB_PGCPP_FRACTION_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_FRACTION_HPP 1

#include <cassert>
#include <cmath>
#include <type_traits>
#include <numeric>

template <typename Z> bool concept Integer = requires {
    std::is_integral<Z>::value;
};

namespace fun {

// Integer { Z }
// inline constexpr Z std::gcd(const Z &a, const Z &b) noexcept {
//     return b == Z(0) ? std::abs(a) : std::gcd(b, a % b);
// }

// Integer { Z }
// inline constexpr Z std::lcm(const Z &a, const Z &b) noexcept {
//     return a / std::gcd(a, b) * b;
// }

Integer { Z }
class Fraction {
    using _Self = Fraction<Z>;

  private:
    Z _numerator;
    Z _denominator;

  public:
    constexpr Fraction(const Z &numerator, const Z &denominator) {
        auto common = std::gcd(numerator, denominator);
        _numerator = numerator / common;
        _denominator = denominator / common;
    }

    constexpr const Z &numerator() const { return _numerator; }

    constexpr const Z &denominator() const { return _denominator; }

    constexpr auto operator+(const _Self &frac) const {
        auto common = std::lcm(_denominator, frac._denominator);
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
        return _numerator * frac._denominator - 
               _denominator * frac._numerator;
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

    constexpr bool operator==(const Z &c) const { return this->cmp(c) == 0; }

    constexpr bool operator<(const Z &c) const { return this->cmp(c) < 0; }

    operator double() { return double(_numerator) / _denominator; }
};

Integer { Z }
constexpr auto operator-(const Z &c, const Fraction<Z> &frac) {
    return Fraction<Z>(frac.denominator() * c - frac.numerator(),
                       frac.denominator());
}

template <typename _Stream, typename Z>
_Stream &operator<<(_Stream &os, const Fraction<Z> &frac) {
    os << frac.numerator() << "/" << frac.denominator();
    return os;
}

// For template deduction
Integer { Z } Fraction(const Z &, const Z &)->Fraction<Z>;

} // namespace fun

#endif