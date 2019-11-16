#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do
                          // this in one cpp file
#include "catch2/catch.hpp"
#include <array>
#include <gsl/span>
#include <vector>

TEST_CASE("Test span", "[testmain]")
{
    auto k = std::array {3, 4, 5};
    auto s = gsl::span<int> {k};

    std::cout << s[1] << '\n';
}
