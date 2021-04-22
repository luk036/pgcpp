#pragma once

#include <utility>
#include <cmath>
#include <numeric>
#include <type_traits>
// #include <concepts>
// #include <ranges>

namespace fun
{

template <typename T>
using Value_type = typename T::value_type;

template <typename T>
using Element_type = typename std::decay<decltype(back(std::declval<T>()))>::type;

/**
 * @brief Sequence
 * 
 * @tparam T 
 */
template <typename T>
concept Sequence = requires(T t, Element_type<T> x)
{
    { t.size() }  -> std::convertible_to<std::size_t>;
    { t.empty() } -> std::convertible_to<bool>;
    { t.back() }  -> std::same_as<Element_type<T> >;
    { t.push_back(x) };
};


template <typename K>
concept ring = std::equality_comparable<K> && requires(K a, K b)
{
    { a +  b } -> std::convertible_to<K>;
    { a -  b } -> std::convertible_to<K>;
    { a *  b } -> std::convertible_to<K>;
    { a += b } -> std::same_as<K&>;
    { a -= b } -> std::same_as<K&>;
    { a *= b } -> std::same_as<K&>;
    { -a }     -> std::convertible_to<K>;
    { K(a) }   -> std::convertible_to<K>;
    { K(0) }   -> std::convertible_to<K>;
};

template <typename K>
concept ordered_ring = ring<K> && std::totally_ordered<K>;


template <typename Z>
concept Integral = ordered_ring<Z> && requires(Z a, Z b)
{
    { a %  b } -> std::convertible_to<Z>;
    { a /  b } -> std::convertible_to<Z>;
    { a %= b } -> std::same_as<Z&>;
    { a /= b } -> std::same_as<Z&>;
};

} // namespace fun
