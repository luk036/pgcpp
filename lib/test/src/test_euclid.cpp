/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "pgcpp/ck_plane.hpp"
#include "pgcpp/euclid_plane.hpp"
#include "pgcpp/pg_line.hpp"
#include "pgcpp/pg_point.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <doctest.h>
#include <iostream>

using namespace fun;

static auto Zero = doctest::Approx(0).epsilon(0.01);

/*!
 * @brief
 *
 * @param a
 * @return true
 * @return false
 */
inline auto ApproxZero(const auto& a) -> bool
{
    return a[0] == Zero and a[1] == Zero and a[2] == Zero;
}

/*!
 * @brief
 *
 * @tparam T
 * @param triangle
 */
template <Projective_plane_prim2 P>
void chk_euclid(const Triple<P>& triangle)
{
    auto trilateral = tri_dual(triangle);

    const auto& [a1, a2, a3] = triangle;
    const auto& [l1, l2, l3] = trilateral;

    // using P = decltype(a1);
    using L = decltype(l1);
    using K = Value_type<P>;
    static_assert(Projective_plane_prim2<L>);
    static_assert(CommutativeRing<K>);

    // auto zero = std::array<K, 3> {0, 0, 0};

    auto [t1, t2, t3] = tri_altitude(triangle);

    auto t4 = harm_conj(t1, t2, t3);
    auto o = orthocenter(triangle);
    auto tau = reflect(l1);
    auto Q = tri_quadrance(triangle);
    auto S = tri_spread(trilateral);

    auto [m12, m23, m13] = tri_midpoint(triangle);

    auto mt1 = a1 * m23;
    auto mt2 = a2 * m13;
    auto mt3 = a3 * m12;

    const auto& [q1, q2, q3] = Q;
    const auto& [s1, s2, s3] = S;

    auto tqf = sq(q1 + q2 + q3) - 2 * (q1 * q1 + q2 * q2 + q3 * q3);
    auto tsf =
        sq(s1 + s2 + s3) - 2 * (s1 * s1 + s2 * s2 + s3 * s3) - 4 * s1 * s2 * s3;
    auto c3 = sq(q1 + q2 - q3) / (4 * q1 * q2);

    auto a3p = plucker(3, a1, 4, a2);
    auto q1p = quadrance(a2, a3p);
    auto q2p = quadrance(a1, a3p);
    auto q3p = quadrance(a1, a2);
    auto tqf2 = Ar(q1p, q2p, q3p); // get 0

    if constexpr (Integral<K>)
    {
        CHECK(not is_parallel(l1, l2));
        CHECK(not is_parallel(l2, l3));
        CHECK(is_perpendicular(t1, l1));
        CHECK(spread(t1, l1) == 1);
        CHECK(coincident(t1, t2, t3));
        CHECK(R(t1, t2, t3, t4) == -1);
        CHECK(o == t2 * t3);
        CHECK(tau(tau(a1)) == a1);
        CHECK(spread(l1, l1) == 0);
        CHECK(quadrance(a1, a1) == 0);
        CHECK(check_sine_law(Q, S));
        CHECK(check_sine_law(S, Q));
        CHECK(coincident(mt1, mt2, mt3));
        // CHECK(cross_s(l1, l2) == c3);
        CHECK((c3 + s3) == 1);
        CHECK(tqf == Ar(q1, q2, q3));
        CHECK(tsf == 0);
        CHECK(tqf2 == 0);
        // auto o2 = orthocenter(
        //               std::tuple {std::move(o), std::move(a2),
        //               std::move(a3)});
        // CHECK(a1 == o2);
    }
    else
    {
        CHECK(cross2(l1, l2) != Zero);
        CHECK(cross2(l2, l3) != Zero);
        CHECK(dot1(t1, l1) == Zero);
        CHECK(spread(t1, l1) - 1 == Zero);
        CHECK(t1.dot(t2 * t3) == Zero);
        CHECK(R(t1, t2, t3, t4) + 1 == Zero);
        CHECK(ApproxZero(cross(meet(t2, t3), o)));
        CHECK(ApproxZero(cross(tau(tau(a1)), a1)));
        CHECK(mt1.dot(mt2 * mt3) == Zero);
        CHECK(spread(l1, l1) == Zero);
        CHECK(quadrance(a1, a1) == Zero);
        CHECK(angle(l1, l1) == Zero);
        CHECK(distance(a1, a1) == Zero);
        // CHECK(cross_s(l1, l2) == doctest::Approx(c3).epsilon(0.01));
        CHECK((c3 + s3) - 1 == Zero);
        CHECK(tqf - Ar(q1, q2, q3) == Zero);
        CHECK(tsf == Zero);
        CHECK(tqf2 == Zero);
        // CHECK(ApproxEqual(a1, orthocenter(std::tuple{std::move(o), // not
        // quite accurate
        //                     std::move(a2), std::move(a3)})));
    }
}

template <typename T>
void chk_cyclic(const T& quadangle)
{
    auto& [u1, u2, u3, u4] = quadangle;

    auto q12 = quadrance(u1, u2);
    auto q23 = quadrance(u2, u3);
    auto q34 = quadrance(u3, u4);
    auto q14 = quadrance(u1, u4);
    auto q24 = quadrance(u2, u4);
    auto q13 = quadrance(u1, u3);

    using P = decltype(u1);
    using K = Value_type<P>;

    if constexpr (Integral<K>)
    {
        auto okay = Ptolemy(std::tuple {std::move(q12), std::move(q23),
            std::move(q34), std::move(q14), std::move(q24), std::move(q13)});
        CHECK(okay);
    }
    else
    {
        auto t = Ar(q12 * q34, q23 * q14, q13 * q24);
        CHECK(t == Zero);
    }
}

TEST_CASE("Euclid plane (cpp_int)")
{
    using boost::multiprecision::cpp_int;

    auto a1 = pg_point<cpp_int> {1, 3, 1};
    auto a2 = pg_point<cpp_int> {4, 2, 1};
    auto a3 = pg_point<cpp_int> {4, -3, 1};

    auto triangle = std::tuple {std::move(a1), std::move(a2), std::move(a3)};
    chk_euclid(triangle);
}

TEST_CASE("Euclid plane (floating point)")
{
    auto a1 = pg_point {1., 3., 1.};
    auto a2 = pg_point {4., 2., 1.};
    auto a3 = pg_point {4., -3., 1.};

    auto triangle = std::tuple {std::move(a1), std::move(a2), std::move(a3)};
    chk_euclid(triangle);
}

TEST_CASE("Euclid Cyclic Points (cpp_int)")
{
    using boost::multiprecision::cpp_int;
    using P = pg_point<cpp_int>;

    auto u1 = uc_point<P>(1, 0);
    auto u2 = uc_point<P>(3, 4);
    auto u3 = uc_point<P>(-1, 2);
    auto u4 = uc_point<P>(0, 1);

    auto quadangle =
        std::tuple {std::move(u1), std::move(u2), std::move(u3), std::move(u4)};
    chk_cyclic(quadangle);
}

TEST_CASE("Euclid Cyclic Points (double)")
{
    using P = pg_point<double>;

    auto u1 = uc_point<P>(1, 0);
    auto u2 = uc_point<P>(3, 4);
    auto u3 = uc_point<P>(-1, 2);
    auto u4 = uc_point<P>(0, 1);

    auto quadangle =
        std::tuple {std::move(u1), std::move(u2), std::move(u3), std::move(u4)};
    chk_cyclic(quadangle);
}
