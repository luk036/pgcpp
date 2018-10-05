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

  pg_point a1{1., 3., 1.};
  pg_point a2{4., 2., 1.};
  pg_point a3{1., 1., -1.};

  pg_line l1{a2 * a3};
  pg_line l2{a1 * a3};
  pg_line l3{a1 * a2};

  return 0;
}
