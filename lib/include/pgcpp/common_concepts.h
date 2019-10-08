#pragma once

#include <utility>

namespace fun
{

template <typename T>
using Value_type = typename T::value_type;

template <typename T>
using Element_type = decltype(back(std::declval<T>()));

template <typename T>
using Iterator_type = decltype(begin(std::declval<T>()));

/**
 * @brief Input iterator concept
 * 
 * @tparam I 
 */
template <typename I>
concept bool Input_iter = requires(I iter)
{
    ++iter;
};

/**
 * @brief Sequence
 * 
 * @tparam T 
 */
template <typename T>
concept bool Sequence = requires(T t, Element_type<T> x)
{
    { t.size() } -> int;
    { t.empty() } -> bool;
    { t.back() } -> Element_type<T>;
    { t.push_back(x) }
};

template <typename T>
concept bool Range = requires(T range)
{
    typename Iterator_type<T>;
    { std::begin(range) } -> Iterator_type<T>;
    { std::end(range) } -> Iterator_type<T>;
};

template <typename P, typename... Args>
concept bool Predicate()
{
    return requires(P pred, Args... args)
    {
        { pred(args...) } ->bool;
    };
}

template <typename T>
concept bool Equality_comparable = requires(T a, T b)
{
    { a == b } -> bool;
    { a != b } -> bool;
};
// ( T, == ) must be reflective, symmetric, and transitive.

template <Equality_comparable K>
concept bool CommutativeRing = requires(K a, K b)
{
    { a + b } -> K;
    { a - b } -> K;
    { a * b } -> K;
    { -a } -> K;
    { K(a) } -> K;
    { K(0) } -> K;
};

template <CommutativeRing Z>
concept bool Integral = requires(Z a, Z b)
{
    { a % b } -> Z;
    { a / b } -> Z;
};

} // namespace fun
