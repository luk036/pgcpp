/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "pgcpp/pg_line.hpp"
#include "pgcpp/pg_point.hpp"
#include "pgcpp/proj_plane.hpp"
#include <catch.hpp>
#include <complex>

using namespace fun;

TEST_CASE("Projective Point", "[proj_plane]") {
    auto p = pg_point{1 - 2j, 3 - 1j, 2 + 1j}; // complex number
    auto q = pg_point{-2 + 1j, 1 - 3j, -1 - 1j};
    auto l = p * q;

    CHECK(l == q * p);
    CHECK(!(l == q));
    CHECK(incident(l, p));
    CHECK(incident(l, q));

    auto pq = plucker(2 + 1j, p, 3 + 0j, q);
    CHECK(incident(l, pq));

    auto h = harm_conj(p, q, pq);
    CHECK(is_harmonic(p, q, pq, h));

    auto r = pg_point{2 - 1j, -2 + 1j, 1 + 1j};
    auto s = pg_point{2j, 2 - 2j, 3 + 0j};
    auto t = pg_point{2 + 0j, -2j, 2 + 0j};

    auto O = meet(join(p, s), join(q, t));
    auto u = plucker(1 + 0j, O, -1 - 1j, r);
    // check_desargue(std::tuple{p, q, r}, std::tuple{s, t, u});
}

TEST_CASE("Projective Line", "[proj_plane]") {
    auto l = pg_line{1 - 2j, 3 - 1j, 2 + 1j}; // complex number
    auto m = pg_line{-2 + 1j, 1 - 3j, -1 - 1j};
    auto A = l * m;
    CHECK(A == m * l);
    CHECK(incident(A, l));
    CHECK(incident(A, m));

    auto lm = plucker(2 + 3j, l, 3 - 4j, m);
    CHECK(incident(A, lm));

    auto r = pg_line{2 - 1j, -2 + 1j, 1 + 1j};
    auto s = pg_line{2j, 2 - 2j, 3 + 0j};
    auto t = pg_line{2 + 0j, -2j, 2 + 0j};

    // assert not persp([l, m, l + m], [r, l + r, l])

    auto o = join(meet(l, s), meet(m, t));
    auto u = plucker(1 + 0j, o, -1 - 1j, r);
    // check_desargue(std::tuple{l, m, r}, std::tuple{s, t, u});
}

TEST_CASE("Projective Point (Special case)", "[proj_plane]") {
    auto p = pg_point{1, 3, 2};    
    auto l = pg_line{-2, 3, 1};    
    // auto l_inf = pg_line{0, 0, 1}; 
    auto l_nan = pg_line{0, 0, 0}; 
    auto p_nan = pg_point{0, 0, 0};

    CHECK(l_nan == l_nan);
    CHECK(l_nan == p * p); // join two equal points
    CHECK(p_nan == l * l);
    CHECK(l_nan == p_nan * p);
    CHECK(p_nan == l_nan * l);
    CHECK(incident(l_nan, p));
    CHECK(incident(p_nan, l));
    CHECK(incident(p_nan, l_nan));
}
