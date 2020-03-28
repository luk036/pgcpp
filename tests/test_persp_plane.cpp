/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "pgcpp/euclid_plane.hpp" // import Ar
#include "pgcpp/persp_plane.hpp"
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
 * @tparam PG
 * @param myck
 */
template <typename PG>
void chk_degenerate(const PG& myck)
{
    using Point = typename PG::point_t;
    // using Line = typename PG::line_t;
    using K = Value_type<Point>;

    auto a1 = Point {-1, 0, 3};
    auto a2 = Point {4, -2, 1};
    auto a3 = Point {3, -1, 1};

    auto triangle = std::tuple {std::move(a1), std::move(a2), std::move(a3)};
    const auto trilateral = tri_dual(triangle);

    // const auto& [a1, a2, a3] = triangle;
    const auto& [l1, l2, l3] = trilateral;

    const auto m12 = myck.midpoint(a1, a2);
    const auto m23 = myck.midpoint(a2, a3);
    const auto m13 = myck.midpoint(a1, a3);
    const auto t1 = a1 * m23;
    const auto t2 = a2 * m13;
    const auto t3 = a3 * m12;

    const auto [q1, q2, q3] = myck.tri_quadrance(triangle);
    const auto [s1, s2, s3] = myck.tri_spread(trilateral);

    const auto tqf = sq(q1 + q2 + q3) - 2 * (q1 * q1 + q2 * q2 + q3 * q3);
    const auto tsf =
        sq(s1 + s2 + s3) - 2 * (s1 * s1 + s2 * s2 + s3 * s3) - 4 * s1 * s2 * s3;
    auto a4 = plucker(3, a1, 4, a2);
    const auto tri2 = std::tuple {std::move(a1), std::move(a2), std::move(a4)};

    const auto [qq1, qq2, qq3] = myck.tri_quadrance(tri2);

    const auto tqf2 = Ar(qq1, qq2, qq3); // get 0

    if constexpr (Integral<K>)
    {
        CHECK(not myck.is_parallel(l1, l2));
        CHECK(not myck.is_parallel(l2, l3));
        CHECK(coincident(t1, t2, t3));
        CHECK(tqf == Ar(q1, q2, q3));
        CHECK(tsf == 0);
        CHECK(tqf2 == 0);
    }
    else
    {
        CHECK(myck.l_infty().dot(l1 * l2) != Zero);
        CHECK(myck.l_infty().dot(l2 * l3) != Zero);
        CHECK(t1.dot(t2 * t3) == Zero);
        CHECK(tqf - Ar(q1, q2, q3) == Zero);
        CHECK(tsf == Zero);
        CHECK(tqf2 == Zero);
    }
}

TEST_CASE("Perspective Euclid plane (cpp_int)")
{
    using boost::multiprecision::cpp_int;

    auto Ire = pg_point<cpp_int>(0, 1, 1);
    auto Iim = pg_point<cpp_int>(1, 0, 0);
    auto l_inf = pg_line<cpp_int>(0, -1, 1);

    const auto P =
        persp_euclid_plane {std::move(Ire), std::move(Iim), std::move(l_inf)};
    chk_degenerate(P);
}

TEST_CASE("Perspective Euclid plane (floating point)")
{
    auto Ire = pg_point {0., 1., 1.};
    auto Iim = pg_point {1., 0., 0.};
    auto l_inf = pg_line {0., -1., 1.};

    const auto P =
        persp_euclid_plane {std::move(Ire), std::move(Iim), std::move(l_inf)};
    chk_degenerate(P);
}
