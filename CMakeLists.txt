cmake_minimum_required (VERSION 3.10)
project (ALPHA_STRONG)

find_package(Boost 1.66.0)

include_directories(${BOOST_INCLUDE_DIRS})

add_executable(ALPHA_STRONG main.cpp)
