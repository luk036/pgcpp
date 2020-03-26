#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include <array>
#include <gsl/span>
#include <vector>

TEST_CASE("Test span")
{
    auto k = std::array {3, 4, 5};
    auto s = gsl::span<int> {k};

    std::cout << s[1] << '\n';
}
