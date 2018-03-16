/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#define CATCH_CONFIG_RUNNER
#include <catch.hpp>


int main(int argc, char* argv[]) {
  //using namespace ModernCppCI;

  auto result = Catch::Session().run(argc, argv);

  return (result < 0xff ? result : 0xff);
}