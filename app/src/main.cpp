/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include "pg_line.hpp"
#include "pg_point.hpp"
#include "proj_plane.hpp"
#include "proj_plane_concepts.h"
#include <iostream>

/**
 * @brief Main program
 *
 * @param argc
 * @param argv
 * @return int
 */
int main(int argc, char *argv[]) {
  using namespace fun;

  auto a1 = pg_point{1., 3., 1.};
  auto a2 = pg_point{4., 2., 1.};
  auto a3 = pg_point{1., 1., -1.};

  auto l1 = pg_line(a2 * a3);
  auto l2 = pg_line(a1 * a3);
  auto l3 = pg_line(a1 * a2);

  return 0;
}
