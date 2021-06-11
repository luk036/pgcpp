// -*- coding: utf-16 -*-
#pragma once

/*! @file include/fractions.hpp
 *  This is a C++ Library header.
 */

#include <boost/operators.hpp>
// #include <cmath>
#include <numeric>
#include <type_traits>
#include "common_concepts.h"

namespace fun
{

template <typename T>
requires std::is_unsigned_v<T>
constexpr auto abs(const T& a) -> T
{
    return a;
}

template <typename T>
requires (!std::is_unsigned_v<T> && ordered_ring<T>)
constexpr auto abs(const T& a) -> T
{
    return (a < T(0))? -a : a;
}

template <typename T>
requires Integral<T> 
constexpr auto gcd_recur(const T& a, const T& b) -> T
{
    // if (a == T(0)) return abs(b);
    if (b == T(0)) return abs(a);
    return gcd_recur(b, a % b);
}

template <typename T>
requires Integral<T> 
constexpr auto gcd(const T& a, const T& b) -> T
{
    if (a == T(0)) return abs(b);
    // if (b == T(0)) return abs(a);
    return gcd_recur(a, b);
}

/*!
 * @brief Least common multiple
 *
 * @tparam _Mn
 * @param[in] __m
 * @param[in] __n
 * @return _Mn
 */
template <Integral _Mn>
inline constexpr auto lcm(_Mn __m, _Mn __n) -> _Mn
{
    if (__m == 0 || __n == 0) return 0;
    return (abs(__m) / gcd(__m, __n)) * abs(__n);
}


template <Integral Z>
struct Fraction : boost::totally_ordered<Fraction<Z>,
                      boost::totally_ordered2<Fraction<Z>, Z,
                          boost::multipliable2<Fraction<Z>, Z,
                              boost::dividable2<Fraction<Z>, Z>>>>
{
    Z _num;
    Z _den;

    /*!
     * @brief Construct a new Fraction object
     *
     * @param[in] num
     * @param[in] den
     */
    constexpr Fraction(Z&& num, Z&& den)
        : _num {std::move(num)}
        , _den {std::move(den)}
    {
        this->normalize();
    }

    /*!
     * @brief Construct a new Fraction object
     *
     * @param[in] num
     * @param[in] den
     */
    constexpr Fraction(const Z& num, const Z& den)
        : _num {num}
        , _den {den}
    {
        this->normalize();
    }

    constexpr void normalize()
    {
        Z common = gcd(this->_num, this->_den);
        if (common == Z(1) || common == Z(0))
        {
            return;
        }
        if (this->_den < Z(0))
        {
            common = -common;
        }
        this->_num /= common;
        this->_den /= common;
    }

    /*!
     * @brief Construct a new Fraction object
     *
     * @param[in] num
     */
    constexpr explicit Fraction(Z&& num)
        : _num {std::move(num)}
        , _den(Z(1))
    {
    }

    /*!
     * @brief Construct a new Fraction object
     *
     * @param[in] num
     */
    constexpr explicit Fraction(const Z& num)
        : _num {num}
        , _den(1)
    {
    }

    /*!
     * @brief Construct a new Fraction object
     *
     * @param[in] num
     */
    constexpr Fraction()
        : _num (0)
        , _den (1)
    {
    }

    /*!
     * @brief
     *
     * @return const Z&
     */
    constexpr auto num() const noexcept -> const Z&
    {
        return _num;
    }

    /*!
     * @brief
     *
     * @return const Z&
     */
    constexpr auto den() const noexcept -> const Z&
    {
        return _den;
    }

    /*!
     * @brief
     *
     * @return Fraction
     */
    constexpr auto abs() const -> Fraction
    {
        return Fraction(abs(this->_num), abs(this->_den));
    }

    /*!
     * @brief
     *
     */
    constexpr void reciprocal() noexcept(std::is_nothrow_swappable_v<Z>)
    {
        std::swap(this->_num, this->_den);
    }

    /*!
     * @brief
     *
     * @return Fraction
     */
    constexpr auto operator-() const -> Fraction
    {
        auto res = Fraction(*this);
        res._num = -res._num;
        return res;
    }

    /*!
     * @brief
     *
     * @param[in] frac
     * @return Fraction
     */
    constexpr auto operator+(const Fraction& frac) const -> Fraction
    {
        if (this->_den == frac._den)
        {
            return Fraction(this->_num + frac._num, this->_den);
        }
        auto d = this->_den * frac._den;
        auto n = frac._den * this->_num + this->_den * frac._num;
        return Fraction(n, d);
    }

    /*!
     * @brief
     *
     * @param[in] frac
     * @return Fraction
     */
    constexpr auto operator-(const Fraction& frac) const -> Fraction
    {
        return *this + (-frac);
    }

    /*!
     * @brief
     *
     * @param[in] frac
     * @return Fraction
     */
    constexpr auto operator*(const Fraction& frac) const -> Fraction
    {
        auto n = this->_num * frac._num;
        auto d = this->_den * frac._den;
        return Fraction(std::move(n), std::move(d));
    }

    /*!
     * @brief
     *
     * @param[in] frac
     * @return Fraction
     */
    constexpr auto operator/(Fraction frac) const -> Fraction
    {
        frac.reciprocal();
        return *this * frac;
    }

    /*!
     * @brief
     *
     * @param[in] i
     * @return Fraction
     */
    constexpr auto operator+(const Z& i) const -> Fraction
    {
        auto n = this->_num + this->_den * i;
        return Fraction(std::move(n), this->_den);
    }

    /*!
     * @brief
     *
     * @param[in] i
     * @return Fraction
     */
    constexpr auto operator-(const Z& i) const -> Fraction
    {
        return *this + (-i);
    }

    // /*!
    //  * @brief
    //  *
    //  * @param[in] i
    //  * @return Fraction
    //  */
    // constexpr Fraction operator*(const Z& i) const noexcept
    // {
    //     auto n = _num * i;
    //     return Fraction(n, _den);
    // }

    // /*!
    //  * @brief
    //  *
    //  * @param[in] i
    //  * @return Fraction
    //  */
    // constexpr Fraction operator/(const Z& i) const noexcept
    // {
    //     auto d = _den * i;
    //     return Fraction(_num, d);
    // }

    /*!
     * @brief
     *
     * @param[in] frac
     * @return Fraction
     */
    constexpr auto operator+=(const Fraction& frac) -> Fraction&
    {
        return *this = *this + frac;
    }

    /*!
     * @brief
     *
     * @param[in] frac
     * @return Fraction
     */
    constexpr auto operator-=(const Fraction& frac) -> Fraction&
    {
        return *this = *this - frac;
    }

    /*!
     * @brief
     *
     * @param[in] frac
     * @return Fraction
     */
    constexpr auto operator*=(const Fraction& frac) -> Fraction&
    {
        return *this = *this * frac;
    }

    /*!
     * @brief
     *
     * @param[in] frac
     * @return Fraction
     */
    constexpr auto operator/=(const Fraction& frac) -> Fraction&
    {
        return *this = *this / frac;
    }

    /*!
     * @brief
     *
     * @param[in] i
     * @return Fraction
     */
    constexpr auto operator+=(const Z& i) -> Fraction&
    {
        return *this = *this + i;
    }

    /*!
     * @brief
     *
     * @param[in] i
     * @return Fraction
     */
    constexpr auto operator-=(const Z& i) -> Fraction&
    {
        return *this = *this - i;
    }

    /*!
     * @brief
     *
     * @param[in] i
     * @return Fraction
     */
    constexpr auto operator*=(const Z& i) -> Fraction&
    {
        const auto common = gcd(i, this->_den);
        if (common == Z(1))
        {
            this->_num *= i;
        }
        // else if (common == Z(0)) [[unlikely]] // both i and den are zero
        // {
        //     this->_num = Z(0);
        // }
        else
        {
            this->_num *= (i / common);
            this->_den /= common;
        }
        return *this;
    }

    /*!
     * @brief
     *
     * @param[in] i
     * @return Fraction
     */
    constexpr auto operator/=(const Z& i) -> Fraction&
    {
        const auto common = gcd(this->_num, i);
        if (common == Z(1))
        {
            this->_den *= i;
        }
        // else if (common == Z(0)) [[unlikely]] // both i and num are zero
        // {
        //     this->_den = Z(0);
        // }
        else
        {
            this->_den *= (i / common);
            this->_num /= common;
        }
        return *this;
    }

    /*!
     * @brief Three way comparison
     *
     * @param[in] frac
     * @return auto
     */
    template <typename U>
    constexpr auto cmp(const Fraction<U>& frac) const -> Fraction&
    {
        if (this->_den == frac._den)
        {
            return (this->_num - frac._num) * this->_den;
        }
        return this->_num * frac._den - this->_den * frac._num;
    }

    template <typename U>
    constexpr auto operator==(const Fraction<U>& rhs) const -> bool
    {
        if (this->_den == rhs._den)
        {
            return this->_num == rhs._num;
        }

        return this->_num * rhs._den == this->_den * rhs._num;
    }

    template <typename U>
    constexpr auto operator<(const Fraction<U>& rhs) const -> bool
    {
        if (this->_den == rhs._den)
        {
            return this->_num < rhs._num;
        }

        return this->_num * rhs._den < this->_den * rhs._num;
    }

    /**
     * @brief
     *
     */
    constexpr auto operator==(const Z& rhs) const -> bool
    {
        return this->_den == Z(1) && this->_num == rhs;
    }

    /**
     * @brief
     *
     */
    constexpr auto operator<(const Z& rhs) const -> bool
    {
        return this->_num < this->_den * rhs;
    }

    /**
     * @brief
     *
     */
    constexpr auto operator>(const Z& rhs) const -> bool
    {
        return this->_num > this->_den * rhs;
    }

    // /*!
    //  * @brief
    //  *
    //  * @return double
    //  */
    // constexpr explicit operator double())
    // {
    //     return double(_num) / _den;
    // }

    // /**
    //  * @brief
    //  *
    //  */
    // friend constexpr bool operator<(const Z& lhs, const Fraction<Z>& rhs))
    // {
    //     return lhs * rhs.den() < rhs.num();
    // }
};


/*!
 * @brief
 *
 * @param[in] c
 * @param[in] frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr auto operator+(const Z& c, const Fraction<Z>& frac) -> Fraction<Z>
{
    return frac + c;
}

/*!
 * @brief
 *
 * @param[in] c
 * @param[in] frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr auto operator-(const Z& c, const Fraction<Z>& frac) -> Fraction<Z>
{
    return c + (-frac);
}

// /*!
//  * @brief
//  *
//  * @param[in] c
//  * @param[in] frac
//  * @return Fraction<Z>
//  */
// template <typename Z>
// constexpr Fraction<Z> operator*(const Z& c, const Fraction<Z>& frac)
// {
//     return frac * c;
// }

/*!
 * @brief
 *
 * @param[in] c
 * @param[in] frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr auto operator+(int&& c, const Fraction<Z>& frac) -> Fraction<Z>
{
    return frac + c;
}

/*!
 * @brief
 *
 * @param[in] c
 * @param[in] frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr auto operator-(int&& c, const Fraction<Z>& frac) -> Fraction<Z>
{
    return (-frac) + c;
}

/*!
 * @brief
 *
 * @param[in] c
 * @param[in] frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr auto operator*(int&& c, const Fraction<Z>& frac) -> Fraction<Z>
{
    return frac * c;
}

/*!
 * @brief
 *
 * @tparam _Stream
 * @tparam Z
 * @param[in] os
 * @param[in] frac
 * @return _Stream&
 */
template <typename _Stream, typename Z>
auto operator<<(_Stream& os, const Fraction<Z>& frac) -> _Stream&
{
    os << frac.num() << "/" << frac.den();
    return os;
}

// For template deduction
// Integral{Z} Fraction(const Z &, const Z &) noexcept -> Fraction<Z>;

} // namespace fun
