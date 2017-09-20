#include "proj_plane_concepts.h"
#include "hyper_plane.hpp"
#include <iostream>

using namespace fun;
using namespace std;
using namespace fun::HG;

int main()
{
  auto a1 = pg_point<int>({1, 3, 1});
  auto a2 = pg_point<int>({4, 2, 1});
  auto a3 = pg_point<int>({1, 1, -1});

  auto l1 = pg_line<int>(a2, a3);
  auto l2 = pg_line<int>(a1, a3);
  auto l3 = pg_line<int>(a1, a2);
  auto s1 = spread(l2, l3);
  auto s2 = spread(l1, l3);
  auto s3 = spread(l1, l2);

  auto q1 = quadrance(a2, a3);
  auto q2 = quadrance(a1, a3);
  auto q3 = quadrance(a1, a2);

  cout << s1 << ',' << s2 << ',' << s3 << '\n';
  cout << q1 << ',' << q2 << ',' << q3 << '\n';

  auto t12 = q1 * s2 - q2 * s1;
  // t12 = sympy.simplify(t12)

  cout << t12 << '\n';

  cout << check_cross_law(s1, s2, s3, q3) << '\n';
  cout << check_cross_law(q1, q2, q3, s3) << '\n';
}