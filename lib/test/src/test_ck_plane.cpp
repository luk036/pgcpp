/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "pgcpp/ck_plane.hpp"
#include "pgcpp/persp_plane.hpp"
#include "pgcpp/pg_common.hpp"
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

template <typename PG>
void chk_ck(const PG& myck)
{
    using P = typename PG::point_t;
    using K = Value_type<P>;

    auto a1 = P {1, -2, 3};
    auto a2 = P {4, 0, 6};
    auto a3 = P {-7, 1, 2};
    const auto triangle =
        std::tuple {std::move(a1), std::move(a2), std::move(a3)};
    const auto trilateral = tri_dual(triangle);
    const auto& [l1, l2, l3] = trilateral;
    const auto [t1, t2, t3] = myck.tri_altitude(triangle);

    auto o = myck.orthocenter(triangle);
    const auto tau = myck.reflect(l1);
    const auto Q = std::tuple {myck.tri_quadrance(triangle)};
    const auto S = std::tuple {myck.tri_spread(trilateral)};

    const auto a4 = P {3, 0, 2};

    if constexpr (Integral<K>)
    {
        CHECK(incident(l1, a2));
        CHECK(myck.is_perpendicular(t1, l1));
        CHECK(coincident(t1, t2, t3));
        CHECK(o == t2 * t3);
        CHECK(a1 ==
            myck.orthocenter(
                std::tuple {std::move(o), std::move(a2), std::move(a3)}));
        CHECK(tau(tau(a4)) == a4);
        CHECK(myck.spread(l2, l2) == 0);
        CHECK(myck.spread(l3, l3) == 0);
        CHECK(myck.quadrance(a1, a1) == 0);
        CHECK(check_sine_law(Q, S));
        CHECK(check_sine_law(S, Q));
    }
    else
    {
        CHECK(l1.dot(a2) == Zero);
        CHECK(l1.dot(myck.perp(t1)) == Zero);
        CHECK(t1.dot(t2 * t3) == Zero);
        CHECK(ApproxZero(cross(o, t2 * t3)));
        const auto o2 = myck.orthocenter(
            std::tuple {std::move(o), std::move(a2), std::move(a3)});
        CHECK(ApproxZero(cross(a1, o2)));
        CHECK(ApproxZero(cross(tau(tau(a4)), a4)));
        CHECK(myck.measure(l2, l2) == Zero);
        CHECK(myck.measure(l3, l3) == Zero);
        CHECK(myck.measure(a1, a1) == Zero);
        const auto& [q1, q2, q3] = Q;
        const auto& [s1, s2, s3] = S;

        const auto r1 = q1 * s2 - q2 * s1;
        const auto r2 = q2 * s3 - q3 * s2;
        CHECK(r1 == Zero);
        CHECK(r2 == Zero);
    }
}

template <typename P, typename L = typename P::dual>
requires Projective_plane_prim<P, L> // c++20 concept
    struct myck : ck<P, L, myck>
{
    constexpr L perp(const P& v) const
    {
        return L(-2 * v[0], v[1], -2 * v[2]);
    }

    constexpr P perp(const L& v) const
    {
        return P(-v[0], 2 * v[1], -v[2]);
    }

    template <Projective_plane2 _P>
    constexpr auto measure(const _P& a1, const _P& a2) const
    {
        auto x = x_ratio(a1, a2, this->perp(a2), this->perp(a1));
        // using Q_t = decltype(x);
        return 1 - x;
    }
};

TEST_CASE("CK plane chk_ck (int)")
{
    // using boost::multiprecision::cpp_int;
    // namespace mp = boost::multiprecision;
    using boost::multiprecision::cpp_int;

    chk_ck(myck<pg_point<cpp_int>>());
    chk_ck(myck<pg_line<cpp_int>>());
    chk_ck(ellck<pg_point<cpp_int>>());
    chk_ck(ellck<pg_line<cpp_int>>());
    chk_ck(hyck<pg_point<cpp_int>>());
    chk_ck(hyck<pg_line<cpp_int>>());

    auto Ire = pg_point<cpp_int> {0, 1, 1};
    auto Iim = pg_point<cpp_int> {1, 0, 0};
    auto l_inf = pg_line<cpp_int> {0, -1, 1};

    auto P =
        persp_euclid_plane {std::move(Ire), std::move(Iim), std::move(l_inf)};
    chk_ck(P);
}

TEST_CASE("CK plane chk_ck (float)")
{
    chk_ck(myck<pg_point<double>>());
    chk_ck(myck<pg_line<double>>());
    chk_ck(ellck<pg_point<float>>());
    chk_ck(ellck<pg_line<float>>());
    chk_ck(hyck<pg_point<double>>());
    chk_ck(hyck<pg_line<double>>());

    auto Ire = pg_point {0., 1., 1.};
    auto Iim = pg_point {1., 0., 0.};
    auto l_inf = pg_line {0., -1., 1.};

    auto P =
        persp_euclid_plane {std::move(Ire), std::move(Iim), std::move(l_inf)};
    chk_ck(P);
}
