cmake_minimum_required (VERSION 2.8.11)
project (pgcpp)
set (CMAKE_BUILD_TYPE Release )
add_definitions ( -Wall -fconcepts -std=c++1z )
add_executable (Main Main.cpp)
target_link_libraries (Main)


