/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch.hpp>
#include <iostream>
#include "pgcpp/ck_plane.hpp"
#include "pgcpp/pg_line.hpp"
#include "pgcpp/pg_point.hpp"

using namespace fun;

static auto Zero = Approx(0).margin(0.01);

/*!
 * @brief
 *
 * @param a
 * @param b
 * @return true
 * @return false
 */
inline auto ApproxEqual(const auto& a, const auto& b) -> bool
{
    return a[0] - b[0] == Zero && a[1] - b[1] == Zero && a[2] - b[2] == Zero;
}

/*!
 * @brief
 *
 * @tparam PG
 * @param myck
 */
template<typename PG>
void chk_tri(const PG& myck)
{
    using Point = typename PG::point_t;
    // using Line = typename PG::line_t;
    using K = Value_type<Point>;

    auto a1 = Point{1, 3, 1};
    auto a2 = Point{4, 2, 1};
    auto a3 = Point{1, 1, -1};

    auto zero = std::array<K, 3>{0, 0, 0};

    auto triangle   = std::tuple{std::move(a1), std::move(a2), std::move(a3)};
    auto trilateral = tri_dual(triangle);

    auto&& [l1, l2, l3] = trilateral;

    auto Q      = std::tuple{myck.tri_quadrance(triangle)};
    auto S      = std::tuple{myck.tri_spread(trilateral)};
    auto a4     = plucker(2, a1, 3, a2);
    auto collin = std::tuple{std::move(a1), std::move(a2), std::move(a4)};
    auto Q2     = myck.tri_quadrance(collin);

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
        CHECK(ApproxEqual(cross(myck.perp(myck.perp(a1)), a1), zero));
        CHECK(ApproxEqual(cross(myck.perp(myck.perp(l1)), l1), zero));
        CHECK(ApproxEqual(cross(myck.perp(myck.perp(l2)), l2), zero));
        CHECK(ApproxEqual(cross(myck.perp(myck.perp(l3)), l3), zero));
        CHECK(check_cross_law(S, std::get<2>(Q)) == Zero);
        CHECK(check_cross_law(Q, std::get<2>(S)) == Zero);
        CHECK(check_cross_TQF(Q2) == Zero);
    }
}

TEST_CASE("Elliptic/Hyperbolic plane", "[ell_plane]")
{
    using boost::multiprecision::cpp_int;

    chk_tri(ellck<pg_point<cpp_int>>());
    chk_tri(ellck<pg_line<cpp_int>>());
    chk_tri(hyck<pg_point<cpp_int>>());
    chk_tri(hyck<pg_line<cpp_int>>());
}

TEST_CASE("Elliptic/Hyperbolic plane (double)", "[ell_plane]")
{
    chk_tri(ellck<pg_point<double>>());
    chk_tri(ellck<pg_line<double>>());
    chk_tri(hyck<pg_point<double>>());
    chk_tri(hyck<pg_line<double>>());
}
