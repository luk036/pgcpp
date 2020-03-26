/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "pgcpp/ck_plane.hpp"
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
 * @tparam PG
 * @param myck
 */
template <typename PG>
void chk_tri(const PG& myck)
{
    using Point = typename PG::point_t;
    using K = Value_type<Point>;

    auto a1 = Point {1, 3, 1};
    auto a2 = Point {4, 2, 1};
    auto a3 = Point {1, 1, -1};

    // auto zero = std::array<K, 3> {0, 0, 0};

    const auto triangle =
        std::tuple {std::move(a1), std::move(a2), std::move(a3)};
    const auto trilateral = tri_dual(triangle);

    const auto& [l1, l2, l3] = trilateral;

    const auto Q = std::tuple {myck.tri_quadrance(triangle)};
    const auto S = std::tuple {myck.tri_spread(trilateral)};

    auto a4 = plucker(2, a1, 3, a2);
    const auto collin =
        std::tuple {std::move(a1), std::move(a2), std::move(a4)};
    const auto Q2 = myck.tri_quadrance(collin);

    if constexpr (Integral<K>)
    {
        CHECK(myck.perp(myck.perp(a1)) == a1);
        CHECK(myck.perp(myck.perp(l1)) == l1);
        CHECK(myck.perp(myck.perp(l2)) == l2);
        CHECK(myck.perp(myck.perp(l3)) == l3);
        CHECK(check_cross_law(S, std::get<2>(Q)) == 0);
        CHECK(check_cross_law(Q, std::get<2>(S)) == 0);
        CHECK(check_cross_TQF(Q2) == 0);
    }
    else
    {
        CHECK(ApproxZero(cross(myck.perp(myck.perp(a1)), a1)));
        CHECK(ApproxZero(cross(myck.perp(myck.perp(l1)), l1)));
        CHECK(ApproxZero(cross(myck.perp(myck.perp(l2)), l2)));
        CHECK(ApproxZero(cross(myck.perp(myck.perp(l3)), l3)));
        CHECK(check_cross_law(S, std::get<2>(Q)) == Zero);
        CHECK(check_cross_law(Q, std::get<2>(S)) == Zero);
        CHECK(check_cross_TQF(Q2) == Zero);
    }
}

TEST_CASE("Elliptic/Hyperbolic plane")
{
    using boost::multiprecision::cpp_int;

    chk_tri(ellck<pg_point<cpp_int>>());
    chk_tri(ellck<pg_line<cpp_int>>());
    chk_tri(hyck<pg_point<cpp_int>>());
    chk_tri(hyck<pg_line<cpp_int>>());
}

TEST_CASE("Elliptic/Hyperbolic plane (double)")
{
    chk_tri(ellck<pg_point<double>>());
    chk_tri(ellck<pg_line<double>>());
    chk_tri(hyck<pg_point<double>>());
    chk_tri(hyck<pg_line<double>>());
}
